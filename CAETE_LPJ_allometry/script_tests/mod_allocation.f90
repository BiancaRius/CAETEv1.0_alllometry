module types
    implicit none
 
    ! FOR THE GNU FORTRAN COMPILER
    integer,parameter,public :: l_1 = 2  ! standart Logical type
    integer,parameter,public :: i_2 = 2  ! 16 bits integer
    integer,parameter,public :: i_4 = 4  ! 32 bits integer
    integer,parameter,public :: r_4 = 4  ! 32 bits float
    integer,parameter,public :: r_8 = 8  ! 64 bits float
 
 end module types
 
 ! program allocation !module to test allocation logic module of LPJmL-Fire
 subroutine allocation(gc_area,lm_test, cw_test, rm_test, dwood_test, sla_test, nind_test,&
     &bminc_test, height_test, cl_inc,cw_inc,ch_inc, cs_inc, cr_inc) !, cr_inc,  cs_inc, ch_inc, ctotal_inc)    
     use types
     use iso_fortran_env, only : output_unit
 
     !VARIABLES [INPUT] - vem do módulo de estabelecimento/FP
     real(r_8), intent(in) :: gc_area
     real, intent(in) :: lm_test, cw_test, rm_test, dwood_test
     real, intent(in) :: sla_test, nind_test, bminc_test, height_test
 
     !VARIABLES [OUTPUT] - o incremento de cada compartimento de carbono neste ano
     real(r_8), intent(out) :: cl_inc !leaf increment (gC)
     real(r_8), intent(out) :: cw_inc !wood increment (gC)
     ! real(r_8), intent(out) :: cr_inc !root increment (gC)
     ! real(r_8), intent(out) :: cs_inc !sapwood increment (gC)
     ! real(r_8), intent(out) :: ch_inc !heartwood increment (gC)
     ! real(r_8), intent(out) :: ctotal_inc !total carbon increment (gC)
 
     integer(i_4),parameter :: npls = 3000
     integer(i_4),parameter :: nseg = 20 ! segment number to bisection loop
     integer(i_4), parameter :: time = 1000
     real(r_8),parameter :: pi   =  3.14159265
     real(r_8),parameter :: xacc =  0.1     !x-axis precision threshold for the allocation solution
     real(r_8),parameter :: yacc =  1.e-10  !y-axis precision threshold for the allocation solution
     integer(i_4),parameter :: ntl=365
     integer, parameter :: stdout = output_unit
 
 
     integer(i_4) :: pls
     integer(i_4) :: xmin, xmax
     real(r_8) :: allom1 = 100
     real(r_8) :: allom2 = 40.0
     real(r_8) :: allom3 = 0.5
     real(r_8) :: latosa = 6000.0
     real(r_8) :: reinickerp = 1.6
     real(r_8) :: ltor = 0.77302587552347657 !leaf:root from Philip
     real(r_8) :: grid_area = 1000 !m2
     real(r_8)  :: dwood   !in g/cm3 but must be in g/m2
     real(r_8)  :: sla !m2/g
     real(r_8)  :: nind !m2
     real(r_8)  :: bminc !kgC/m2 - !total biomass increment this year for area
     
 
     real(r_8)  :: lm_ind !in kgC/m2 but use in gC/ind [transformation below]
     real(r_8)  :: sm_ind !in kgC/m2 but use in gC/ind [transformation below]
     real(r_8)  :: hm_ind !in kgC/m2 but use in gC/ind [transformation below]
     real(r_8)  :: rm_ind !in kgC/m2 but use in gC/ind [transformation below]
     real(r_8)  :: cw_ind !to calculate sap and heartwood (in kgC/m2 but use in gC/ind)
     real(r_8)  :: cwood
     real(r_8)  :: litter_ag_fast
     real(r_8)  :: litter_ag_slow
     real(r_8)  :: litter_bg
     real(r_8)  :: lai_ind
     real(r_8)  :: crown_area
     real(r_8)  :: crown_area_ind
     real(r_8)  :: diameter
     real(r_4)  :: height
     real(r_8)  :: fpc_grid
     real(r_8)  :: fpc_ind
     real(r_8)  :: fpc_inc
     real(r_8)  :: fpc_grid_old
 
 
     !Local Variables
     real(r_8)  :: bminc_ind !gC/ind - individual total biomass increment this year 
     real(r_8)  :: lm2rm          !ratio of leafmass to fine rootmass
     real(r_8)  :: lminc_ind      !individual leafmass increment this year
     real(r_8)  :: rminc_ind      !individual fineroot mass increment this year
     real(r_8)  :: lminc_ind_min  !min leafmass increment to maintain current sapwood
     real(r_8)  :: rminc_ind_min  !min rootmass increment to support new leafmass
     real(r_8)  :: sap_xsa        !cross sectional area of sapwood  
     real(r_8)  :: sminc_ind      !individual sapmass increment this year
     real(r_8) :: fpc_inc_tree     !this years total FPC increment for tree PFTs
     real(r_8) :: fpc_tree_total   !total grid FPC for tree PFTs      
     real(r_8) :: excess                           !total tree FPC or grass cover to be reduced
     real(r_8) :: fpc_tree_max
     real(r_8)  :: nind_kill        !reduction in individual density to reduce tree FPC to permitted maximum (indiv/m2)
     real(r_8)  :: rm_kill          !reduction in grass PFT root mass to reduce grass cover to permitted maximum (gC)  
     real(r_8)  :: lm_kill          !reduction in grass PFT leaf mass to reduce grass cover to permitted maximum (gC)
     real(r_8)  :: lm_old
 
     real(r_8)  :: lm !leaf mass
     real(r_8)  :: sm !sapwood mass
     real(r_8)  :: hm !heartwood mass
     real(r_8)  :: rm !root mass
 
     real(r_8)  :: x1             !working vars in bisection
     real(r_8)  :: x2
     real(r_8)  :: rtbis
     real(r_8)  :: dx
     real(r_8)  :: xmid
     real(r_8)  :: root1, root2, root3
     real(r_8)  :: sign
     real(r_8) :: wooddens = 2.e5
     logical  :: normal
 
     real(r_8)  :: fx1
     real(r_8)  :: fmid
 
 
     real(r_8)  :: lm1     !allometric leafmass requirement (leafmass req'd to keep sapwood alive; gC ind-1)
 
     integer :: i
 
     ! print*, 'cl2 inside alloc', lm_test, cw_test, cr_test, dwood_test,&
      ! &sla_test, nind_test, bminc_test
     !Arrays with values to some variables (generic values)
     dwood = dwood_test*1000000 !(/0.74,0.73,0.59,0.52,0.41,0.44,0.86,0.42,0.64,0.69,0.92,&
     ! &0.60,0.36,0.99,0.59,0.52,0.41,0.44,0.86,0.42/) !atenção para a unidade
     bminc = bminc_test !(/2.15,2.,2.18,2.6,2.5,1.8,2.3,2.,1.8,2.84,2.25,3.,2.2,1.7,&
     !&1.18,2.6,3.5,2.8,3.3,2./)
     sla = sla_test! (/0.002,0.018,0.009,0.023,0.013,0.039,0.040,0.0028,0.0025,&
     ! &0.027,0.032,0.007,0.013,0.025,0.002,0.008,0.004,0.016,0.023,0.015/)
     nind = nind_test !(/1.,2.,8.,6.,5.,9.,3.,4.,7.,1.,2.,8.,5.,3.,6.,4.,5.,8.,9.,3./)
     height = height_test !(/5.,9.,15.,10.9,11.5,18.9,12.6,2.5,14.9,22.5,28.7,23.6,&
     ! &28.8,19.6,13.3,27.6,29.5,21.6,30.,2./)
     lm_ind = lm_test!(/2.15,2.,1.18,1.6,1.5,1.8,0.3,2.,0.8,.84,0.25,1.,0.2,1.7,&
     ! &1.18,1.6,1.5,1.8,0.3,2./)
     cw_ind = cw_test !(/7.,12.,7.2,8.3,8.8,9.7,7.5,11.5,10.,8.6,7.3,10.3,6.8,9.9,&
     ! &5.3,9.2,15.,12.6,10.7,11.4/)
     rm_ind = rm_test !(/0.63,0.8,0.9,0.5,1.3,0.9,0.4,1.0,0.56,0.87,0.33,0.97,0.31,&
     ! &0.55,0.2,0.8,0.4,0.66,0.23,1.5/)
 
     !-----------------------------------------------------------------
    ! print*,'inside', height
 
    
 
         lm  = (lm_ind)! /nind )*1.D3 !já vem como gC/ ind
         sm  = (cw_ind *0.05) !/nind )*1.D3
         hm  = (cw_ind *0.95)!/nind )*1.D3
         rm  = rm_ind !/nind )*1.D3
         bminc_ind  = bminc ! bminc ! /nind )*1.D3
         
         !checking entering values
         ! print*, 'LM=', lm ,'SM',sm ,'HM', hm ,'RM', rm, 'BMINC', bminc_ind 
 
         ! ====== TREE ALLOCATION ======
 
         lm1  = (latosa*sm /(dwood)*height *sla )  !allometric leaf mass requirement *****ATENÇÃO*****
         !print*, 'LM1', lm1, latosa, sm, dwood, height, sla, lm 
 
         ! lm1  = 1000.0  !valor arbitrario colocado para rever a unidade do dwood
 
         lminc_ind_min  = lm1  - lm   !eqn (27)
         ! print*, 'LM MIN', lminc_ind_min , pls
 
         ! lminc_ind_min  = 0.6
     
         !calculate minimum root production to support this leaf mass (i.e. lm_ind + lminc_ind_min)
         !May be negative following a reduction in soil water limitation (increase in lm2rm) relative to last year.
 
         rminc_ind_min  = lm1  / ltor - rm       !eqn (30)
         ! print*, 'RM MIN', rminc_ind_min 
 
         rminc_ind_min  = (latosa*sm /(dwood)*height *sla *ltor) - rm       !eqn (30)
         ! print*, 'RM MIN teste', rminc_ind_min 
 
 
         if (rminc_ind_min  .gt. 0. .and. lminc_ind_min  .gt. 0. .and. &
             &(rminc_ind_min  + lminc_ind_min ) .le. bminc_ind ) then
 
             !Normal allocation (positive increment to all living C compartments)
            !  print*, 'normal'
             normal = .true.
 
             !Calculation of leaf mass increment (lminc_ind) that satisfies Eqn (22)
             !Since this is normal allocation, we set the lower bound for the leafmass allocation (x1)
             !to its allometric minimum, because it should be able to be fulfilled, i.e.:
 
             !Start to find root procedure (relate to bisection method)
 
             x1  = lminc_ind_min 
             x2  = (bminc_ind  - (lm  / ltor - rm )) / (1. + 1. / ltor)
             
             dx  = x2  - x1 
 
             if (dx  < 0.01) then
 
                 !there seems to be rare cases where lminc_ind_min (x1) is almost equal to x2. In this case,
                 !assume that the leafmass increment is equal to the midpoint between the values and skip 
                 !the root finding procedure
 
                 lminc_ind  = x1  + 0.5 * dx 
 
             else
                 !Find a root for non-negative lminc_ind, rminc_ind and sminc_ind using Bisection Method (Press et al 1986, p 346)
                 !There should be exactly one solution (no proof presented, but Steve has managed one).
                     
                 dx  = dx /nseg
 
                 !! ===== FIND ROOT FUNCTION ===== [**must be a function**]
 
                 pi4 = pi/4
                 a1 = 2./allom3
                 a2 = 1. + a1
                 a3 = allom2**a1
 
 
                 root1  = a3*((sm +bminc_ind -x1 -((lm +x1 )/ltor)+&
                         &rm +hm )/dwood )/pi4-((sm +bminc_ind -x1 -&
                         &((lm +x1 )/ltor)+rm )/((lm +x1 )*sla *&
                         &(dwood )/latosa))**a2
 
                 ! ======================================================
 
                 !evaluate f(x1) = LHS of eqn (22) at x1
 
                 fx1  = root1 
 
                 !Find approximate location of leftmost root on the interval (x1,x2).
                 !Subdivide (x1,x2) into nseg equal segments seeking change in sign of f(xmid) relative to f(x1).
 
                 fmid  = fx1 
                 xmid  = x1 
 
                 i = 1
 
                 do           
 
                     xmid  = xmid  + dx 
 
                     ! root2  = a3*((sm +bminc_ind -xmid -((lm +xmid )/ltor)+&
                     ! &rm +hm )/wooddens)/pi4-((sm +bminc_ind -xmid -&
                     ! &((lm +xmid )/ltor)+rm )/((lm +xmid )*sla *&
                     ! &(wooddens)/latosa))**a2
 
                     root2  = a3*((sm +bminc_ind -xmid -((lm +xmid )/ltor)+&
                     &rm +hm )/dwood)/pi4-((sm +bminc_ind -xmid -&
                     &((lm +xmid )/ltor)+rm )/((lm +xmid )*sla *&
                     &(dwood)/latosa))**a2
 
                     fmid  = root2 
 
                     if ((fmid *fx1 ) .le. 0. .or. xmid  .ge. x2 ) exit  !sign has changed or we are over the upper bound
 
                     if (i > 20) write(stdout,*)'first alloc loop flag',i,pls,fmid *fx1 ,&
                          &xmid ,x1 ,x2 ,dx ,bminc_ind 
 
                     if (i > 50) stop 'Too many iterations allocmod'
 
                     i = i + 1
 
                 enddo
 
                 !the interval that brackets zero in f(x) becomes the new bounds for the root search
 
                 x1  = xmid  - dx 
                 x2  = xmid 
 
                 !Apply bisection method to find root on the new interval (x1,x2)
 
                 fx1  = root1 
 
                 if (fx1  .ge. 0.) then
                     sign  = -1.
                 else
                     sign  =  1.
                 end if
 
                 rtbis  = x1 
                 dx     = x2  - x1 
 
                 !Bisection loop: search iterates on value of xmid until xmid lies within xacc of the root,
                 !i.e. until |xmid-x| < xacc where f(x) = 0. the final value of xmid with be the leafmass increment
 
                 i = 1
 
                 do               
 
                     dx    = 0.5 * dx 
                     xmid  = rtbis  + dx 
 
                     !calculate fmid = f(xmid) [eqn (22)]
 
                     ! root3  = a3*((sm +bminc_ind -xmid -((lm +xmid )/ltor)+&
                     ! &rm +hm )/wooddens)/pi4-((sm +bminc_ind -xmid -&
                     ! &((lm +xmid )/ltor)+rm )/((lm +xmid )*sla *&
                     ! &wooddens/latosa))**a2
 
                     root3  = a3*((sm +bminc_ind -xmid -((lm +xmid )/ltor)+&
                     &rm +hm )/dwood)/pi4-((sm +bminc_ind -xmid -&
                     &((lm +xmid )/ltor)+rm )/((lm +xmid )*sla *&
                     &dwood/latosa))**a2
 
                     fmid  = root3 
 
                     if (fmid  * sign  .le. 0.) rtbis  = xmid 
 
                     if (dx  < xacc .or. abs(fmid ) <= yacc) exit
 
                     if (i > 20) write(stdout,*)'second alloc loop flag',i,pls,dx ,abs(fmid )
                     if (i > 50) stop 'Too many iterations allocmod'
 
                     i = i + 1
                 enddo
               
 
                 !Now rtbis contains numerical solution for lminc_ind given eqn (22)
 
                 lminc_ind  = rtbis
                !  print*, 'lminc', lminc_ind 
 
             endif
 
             !Calculate increments in other compartments using allometry relationships
 
             rminc_ind  = (lm  + lminc_ind ) / ltor - rm        !eqn (9)
 
             sminc_ind  = bminc_ind  - rminc_ind  - lminc_ind   !eqn (1)
 
             ! print*, 'LEAF_INC (gC/ind)', (lminc_ind /1.D3), 'ROOT_INC (gC/ind)', (rminc_ind /1.D3),&
             ! & 'SAP_INC(gC/ind)', (sminc_ind /1.D3), pls, 'NORMAL'
 
         else 
 
             !Abnormal allocation: reduction in some C compartment(s) to satisfy allometry
             
             normal = .false.
 
             !Attempt to distribute this year's production among leaves and roots only
 
             lminc_ind  = (bminc_ind -lm /ltor+rm )/(1.+1./ltor)  !eqn (33)
             !print*,'anormal', lminc_ind
 
             if (lminc_ind  > 0.) then
 
                 !Positive allocation to leafmass
 
                 rminc_ind  = bminc_ind  - lminc_ind   !eqn (31)
                 
                 !Add killed roots (if any) to below-ground litter
 
                 if (rminc_ind  < 0.) then
 
                     lminc_ind  = bminc_ind 
                     rminc_ind  = (lm  + lminc_ind ) / ltor - rm 
 
                     litter_bg  = litter_bg  + abs(rminc_ind ) * nind 
 
                 end if
                 
                 i = 1
 
             else
 
                 !Negative allocation to leaf mass
 
                 rminc_ind  = bminc_ind 
                 lminc_ind  = (rm  + rminc_ind ) * ltor - lm   !from eqn (9)
 
                 !Add killed leaves to litter
 
                 litter_ag_fast  = litter_ag_fast  + abs(lminc_ind ) * nind 
                 
                 i = 2
 
             endif
 
             !Calculate sminc_ind (must be negative)
       
             sminc_ind  = (lm  + lminc_ind ) * sla  /&
             & latosa * 2.e5 * height  - sm   !eqn (35)
 
             !Convert killed sapwood to heartwood
 
             hm  = hm  + abs(sminc_ind )
 
             ! print*, pls, 'ANNORMAL'
 
 
         endif !normal/abnormal allocation
 
         !Increment on C compartments - OUTPUT FINAL (gC)
         cl_inc = lminc_ind
         cr_inc = rminc_ind
         cs_inc = sminc_ind
         ch_inc = hm 
         cw_inc = cs_inc + ch_inc
         ctotal_inc = cl_inc + cr_inc + cs_inc + ch_inc 
 
         ! lm_ind  = ((lm  + lminc_ind )*nind )/1.D3
         ! rm_ind  = ((rm  + rminc_ind )*nind )/1.D3 
         ! sm_ind  = ((sm  + sminc_ind )*nind )/1.D3
         ! hm_ind  = (hm *nind )/1.D3
         ! cwood  = sm_ind +hm_ind 
         ! print*, 'leaf carbon',lm_ind 
         ! print*, 'LMINC', lminc_ind , pls
         !print*, 'LM', lm_ind , 'RM', rm_ind , 'SM', sm_ind , 'HM', hm_ind , 'CWOOD', cwood , pls
 
         !ALLOMETRY EQUATIONS
 
         !DIAMETER (m)
         diameter  = (4*(cwood *1.0D3)/(dwood *1D7)*pi*allom2)&
         &**(1/(2+allom3))
         
         !CROWN AREA (m2)
         crown_area  = allom1*(diameter **1.6)
         crown_area_ind  = (crown_area /nind )
 
         !LAI (m2/m2)
         lai_ind =(((lm_ind /nind )*sla )/crown_area_ind )
         
         !ALTURA (m)
         ! height  = allom2*(diameter **allom3)
 
         
         fpc_ind  = 1. - exp(-0.5 * 8.)
         fpc_grid  = crown_area  * nind  * fpc_ind 
         fpc_grid_old  = fpc_grid 
         fpc_inc  = max(fpc_grid  - fpc_grid_old ,0.)
 
         !SELF-THINNING LOGIC 
         fpc_tree_max = grid_area*0.95
 
         ! fpc_inc_tree    = sum(fpc_inc(:))
         fpc_inc_tree    = fpc_inc            
         ! fpc_tree_total  = sum(fpc_grid, mask = present .and. tree)
 
 
 
         ! print*, 'FPC_IND', fpc_ind , 'FPC GRID', fpc_grid 
    
 end subroutine allocation
 
 ! ! program allocation !module to test allocation logic module of LPJmL-Fire
 ! subroutine allocation(gc_area,lm_test, cw_test, rm_test, dwood_test, sla_test, nind_test,&
 !     &bminc_test, height_test)    
 !     use types
 !     use iso_fortran_env, only : output_unit
 
 !     !VARIABLES [INPUT] - Determinadas arbitrariamente
 !     real(r_8), intent(in) :: gc_area
 !     real, intent(in) :: lm_test, cw_test, rm_test, dwood_test
 !     real, intent(in) :: sla_test, nind_test, bminc_test, height_test
 
 
 !     integer(i_4),parameter :: npls = 3000
 !     integer(i_4),parameter :: nseg = 20 ! segment number to bisection loop
 !     integer(i_4), parameter :: time = 1000
 !     real(r_8),parameter :: pi   =  3.14159265
 !     real(r_8),parameter :: xacc =  0.1     !x-axis precision threshold for the allocation solution
 !     real(r_8),parameter :: yacc =  1.e-10  !y-axis precision threshold for the allocation solution
 !     integer(i_4),parameter :: ntl=365
 !     integer, parameter :: stdout = output_unit
 
 
 !     integer(i_4) :: pls
 !     integer(i_4) :: xmin, xmax
 !     real(r_8) :: allom1 = 100
 !     real(r_8) :: allom2 = 40.0
 !     real(r_8) :: allom3 = 0.5
 !     real(r_8) :: latosa = 6000.0
 !     real(r_8) :: reinickerp = 1.6
 !     real(r_8) :: ltor = 0.77302587552347657 !leaf:root from Philip
 !     real(r_8) :: grid_area = 1000 !m2
 !     real(r_8), dimension(npls) :: dwood   !in g/cm3 but must be in g/m2
 !     real(r_8), dimension(npls) :: sla !m2/g
 !     real(r_8), dimension(npls) :: nind !m2
 !     real(r_8), dimension(npls) :: bminc !kgC/m2 - !total biomass increment this year for area
     
 
 !     real(r_8), dimension(npls) :: lm_ind !in kgC/m2 but use in gC/ind [transformation below]
 !     real(r_8), dimension(npls) :: sm_ind !in kgC/m2 but use in gC/ind [transformation below]
 !     real(r_8), dimension(npls) :: hm_ind !in kgC/m2 but use in gC/ind [transformation below]
 !     real(r_8), dimension(npls) :: rm_ind !in kgC/m2 but use in gC/ind [transformation below]
 !     real(r_8), dimension(npls) :: cw_ind !to calculate sap and heartwood (in kgC/m2 but use in gC/ind)
 !     real(r_8), dimension(npls) :: cwood
 !     real(r_8), dimension(npls) :: litter_ag_fast
 !     real(r_8), dimension(npls) :: litter_ag_slow
 !     real(r_8), dimension(npls) :: litter_bg
 !     real(r_8), dimension(npls) :: lai_ind
 !     real(r_8), dimension(npls) :: crown_area
 !     real(r_8), dimension(npls) :: crown_area_ind
 !     real(r_8), dimension(npls) :: diameter
 !     real(r_4), dimension(npls) :: height
 !     real(r_8), dimension(npls) :: fpc_grid
 !     real(r_8), dimension(npls) :: fpc_ind
 !     real(r_8), dimension(npls) :: fpc_inc
 !     real(r_8), dimension(npls) :: fpc_grid_old
 
 
 !     !Local Variables
 !     real(r_8),dimension(npls) :: bminc_ind !gC/ind - individual total biomass increment this year 
 !     real(r_8),dimension(npls) :: lm2rm          !ratio of leafmass to fine rootmass
 !     real(r_8),dimension(npls) :: lminc_ind      !individual leafmass increment this year
 !     real(r_8),dimension(npls) :: rminc_ind      !individual fineroot mass increment this year
 !     real(r_8),dimension(npls) :: lminc_ind_min  !min leafmass increment to maintain current sapwood
 !     real(r_8),dimension(npls) :: rminc_ind_min  !min rootmass increment to support new leafmass
 !     real(r_8),dimension(npls) :: sap_xsa        !cross sectional area of sapwood  
 !     real(r_8),dimension(npls) :: sminc_ind      !individual sapmass increment this year
 !     real(r_8) :: fpc_inc_tree     !this years total FPC increment for tree PFTs
 !     real(r_8) :: fpc_tree_total   !total grid FPC for tree PFTs      
 !     real(r_8) :: excess                           !total tree FPC or grass cover to be reduced
 !     real(r_8) :: fpc_tree_max
 !     real(r_8),dimension(npls) :: nind_kill        !reduction in individual density to reduce tree FPC to permitted maximum (indiv/m2)
 !     real(r_8),dimension(npls) :: rm_kill          !reduction in grass PFT root mass to reduce grass cover to permitted maximum (gC)  
 !     real(r_8),dimension(npls) :: lm_kill          !reduction in grass PFT leaf mass to reduce grass cover to permitted maximum (gC)
 !     real(r_8),dimension(npls) :: lm_old
 
 !     real(r_8),dimension(npls) :: lm !leaf mass
 !     real(r_8),dimension(npls) :: sm !sapwood mass
 !     real(r_8),dimension(npls) :: hm !heartwood mass
 !     real(r_8),dimension(npls) :: rm !root mass
 
 !     real(r_8),dimension(npls) :: x1             !working vars in bisection
 !     real(r_8),dimension(npls) :: x2
 !     real(r_8),dimension(npls) :: rtbis
 !     real(r_8),dimension(npls) :: dx
 !     real(r_8),dimension(npls) :: xmid
 !     real(r_8),dimension(npls) :: root1, root2, root3
 !     real(r_8),dimension(npls) :: sign
 !     real(r_8) :: wooddens = 2.e5
 !     logical  :: normal
 
 !     real(r_8),dimension(npls) :: fx1
 !     real(r_8),dimension(npls) :: fmid
 
 
 !     real(r_8),dimension(npls) :: lm1     !allometric leafmass requirement (leafmass req'd to keep sapwood alive; gC ind-1)
 
 !     integer :: i
 
 !     ! print*, 'cl2 inside alloc', lm_test, cw_test, cr_test, dwood_test,&
 !         ! &sla_test, nind_test, bminc_test
 !     !Arrays with values to some variables (generic values)
 !     dwood = dwood_test !(/0.74,0.73,0.59,0.52,0.41,0.44,0.86,0.42,0.64,0.69,0.92,&
 !     ! &0.60,0.36,0.99,0.59,0.52,0.41,0.44,0.86,0.42/) !atenção para a unidade
 !     bminc = bminc_test !(/2.15,2.,2.18,2.6,2.5,1.8,2.3,2.,1.8,2.84,2.25,3.,2.2,1.7,&
 !     !&1.18,2.6,3.5,2.8,3.3,2./)
 !     sla = sla_test! (/0.002,0.018,0.009,0.023,0.013,0.039,0.040,0.0028,0.0025,&
 !     ! &0.027,0.032,0.007,0.013,0.025,0.002,0.008,0.004,0.016,0.023,0.015/)
 !     nind = nind_test !(/1.,2.,8.,6.,5.,9.,3.,4.,7.,1.,2.,8.,5.,3.,6.,4.,5.,8.,9.,3./)
 !     height= height_test !(/5.,9.,15.,10.9,11.5,18.9,12.6,2.5,14.9,22.5,28.7,23.6,&
 !     ! &28.8,19.6,13.3,27.6,29.5,21.6,30.,2./)
 !     lm_ind= lm_test!(/2.15,2.,1.18,1.6,1.5,1.8,0.3,2.,0.8,.84,0.25,1.,0.2,1.7,&
 !     ! &1.18,1.6,1.5,1.8,0.3,2./)
 !     cw_ind = cw_test !(/7.,12.,7.2,8.3,8.8,9.7,7.5,11.5,10.,8.6,7.3,10.3,6.8,9.9,&
 !     ! &5.3,9.2,15.,12.6,10.7,11.4/)
 !     rm_ind= rm_test !(/0.63,0.8,0.9,0.5,1.3,0.9,0.4,1.0,0.56,0.87,0.33,0.97,0.31,&
 !     ! &0.55,0.2,0.8,0.4,0.66,0.23,1.5/)
 
 !     !-----------------------------------------------------------------
 
 
 !     do pls = 1,npls
 
 !         lm(pls) = (lm_ind(pls)/nind(pls))*1.D3 !PRECISA COLOCAR VALORES INICIAIS
 !         sm(pls) = ((cw_ind(pls)*0.05)/nind(pls))*1.D3
 !         hm(pls) = ((cw_ind(pls)*0.95)/nind(pls))*1.D3
 !         rm(pls) = (rm_ind(pls)/nind(pls))*1.D3
 
 !         ! print*, 'LM=', lm(pls),'SM',sm(pls),'HM', hm(pls),'RM', rm(pls)
         
 !         bminc_ind(pls) = (bminc(pls)/nind(pls))*1.D3
 
 !         ! ====== TREE ALLOCATION ======
 
 !         lm1(pls) = (latosa*sm(pls)/(dwood(pls)*1000)*height(pls)*sla(pls))  !allometric leaf mass requirement *****ATENÇÃO*****
 !         ! print*, 'LM1', lm1(pls)
 
 !         ! lm1(pls) = 1000.0  !valor arbitrario colocado para rever a unidade do dwood
 
 !         lminc_ind_min(pls) = lm(pls) - lm1(pls)  !eqn (27)
 !         ! print*, 'LM MIN', lminc_ind_min(pls), pls
 
 !         ! lminc_ind_min(pls) = 0.6
     
 !         !calculate minimum root production to support this leaf mass (i.e. lm_ind + lminc_ind_min)
 !         !May be negative following a reduction in soil water limitation (increase in lm2rm) relative to last year.
 
 !         rminc_ind_min(pls) = lm1(pls) / ltor - rm(pls)      !eqn (30)
 !         ! print*, 'RM MIN', rminc_ind_min(pls)
 
 !         rminc_ind_min(pls) = (latosa*sm(pls)/(dwood(pls)*1000)*height(pls)*sla(pls)*ltor) - rm(pls)      !eqn (30)
 !         ! print*, 'RM MIN teste', rminc_ind_min(pls)
 
 
 !         if (rminc_ind_min(pls) .gt. 0. .and. lminc_ind_min(pls) .gt. 0. .and. &
 !             &(rminc_ind_min(pls) + lminc_ind_min(pls)) .le. bminc_ind(pls)) then
 
 !             !Normal allocation (positive increment to all living C compartments)
 !             print*, 'normal'
 !             normal = .true.
 
 !             !Calculation of leaf mass increment (lminc_ind) that satisfies Eqn (22)
 !             !Since this is normal allocation, we set the lower bound for the leafmass allocation (x1)
 !             !to its allometric minimum, because it should be able to be fulfilled, i.e.:
 
 !             !Start to find root procedure (relate to bisection method)
 
 !             x1(pls) = lminc_ind_min(pls)
 !             x2(pls) = (bminc_ind(pls) - (lm(pls) / ltor - rm(pls))) / (1. + 1. / ltor)
             
 !             dx(pls) = x2(pls) - x1(pls)
 
 !             if (dx(pls) < 0.01) then
 
 !                 !there seems to be rare cases where lminc_ind_min (x1) is almost equal to x2. In this case,
 !                 !assume that the leafmass increment is equal to the midpoint between the values and skip 
 !                 !the root finding procedure
 
 !                 lminc_ind(pls) = x1(pls) + 0.5 * dx(pls)
 
 !             else
 !                 !Find a root for non-negative lminc_ind, rminc_ind and sminc_ind using Bisection Method (Press et al 1986, p 346)
 !                 !There should be exactly one solution (no proof presented, but Steve has managed one).
                     
 !                 dx(pls) = dx(pls)/nseg
 
 !                 !! ===== FIND ROOT FUNCTION ===== [**must be a function**]
 
 !                 pi4 = pi/4
 !                 a1 = 2./allom3
 !                 a2 = 1. + a1
 !                 a3 = allom2**a1
 
 
 !                 root1(pls) = a3*((sm(pls)+bminc_ind(pls)-x1(pls)-((lm(pls)+x1(pls))/ltor)+&
 !                         &rm(pls)+hm(pls))/dwood(pls))/pi4-((sm(pls)+bminc_ind(pls)-x1(pls)-&
 !                         &((lm(pls)+x1(pls))/ltor)+rm(pls))/((lm(pls)+x1(pls))*sla(pls)*&
 !                         &(dwood(pls))/latosa))**a2
 
 !                 ! ======================================================
 
 !                 !evaluate f(x1) = LHS of eqn (22) at x1
 
 !                 fx1(pls) = root1(pls)
 
 !                 !Find approximate location of leftmost root on the interval (x1,x2).
 !                 !Subdivide (x1,x2) into nseg equal segments seeking change in sign of f(xmid) relative to f(x1).
 
 !                 fmid(pls) = fx1(pls)
 !                 xmid(pls) = x1(pls)
 
 !                 i = 1
 
 !                 do
 
 !                     xmid(pls) = xmid(pls) + dx(pls)
 
 !                     root2(pls) = a3*((sm(pls)+bminc_ind(pls)-xmid(pls)-((lm(pls)+xmid(pls))/ltor)+&
 !                     &rm(pls)+hm(pls))/wooddens)/pi4-((sm(pls)+bminc_ind(pls)-xmid(pls)-&
 !                     &((lm(pls)+xmid(pls))/ltor)+rm(pls))/((lm(pls)+xmid(pls))*sla(pls)*&
 !                     &(wooddens)/latosa))**a2
 
 !                     fmid(pls) = root2(pls)
 
 !                     if ((fmid(pls)*fx1(pls)) .le. 0. .or. xmid(pls) .ge. x2(pls)) exit  !sign has changed or we are over the upper bound
 
 !                     if (i > 20) write(stdout,*)'first alloc loop flag',i,pls,fmid(pls)*fx1(pls),&
 !                          &xmid(pls),x1(pls),x2(pls),dx(pls),bminc_ind(pls)
 
 !                     if (i > 50) stop 'Too many iterations allocmod'
 
 !                     i = i + 1
 
 !                 end do
 
 !                 !the interval that brackets zero in f(x) becomes the new bounds for the root search
 
 !                 x1(pls) = xmid(pls) - dx(pls)
 !                 x2(pls) = xmid(pls)
 
 !                 !Apply bisection method to find root on the new interval (x1,x2)
 
 !                 fx1(pls) = root1(pls)
 
 !                 if (fx1(pls) .ge. 0.) then
 !                     sign(pls) = -1.
 !                 else
 !                     sign(pls) =  1.
 !                 end if
 
 !                 rtbis(pls) = x1(pls)
 !                 dx(pls)    = x2(pls) - x1(pls)
 
 !                 !Bisection loop: search iterates on value of xmid until xmid lies within xacc of the root,
 !                 !i.e. until |xmid-x| < xacc where f(x) = 0. the final value of xmid with be the leafmass increment
 
 !                 i = 1
 
 !                 do 
 
 !                     dx(pls)   = 0.5 * dx(pls)
 !                     xmid(pls) = rtbis(pls) + dx(pls)
 
 !                     !calculate fmid = f(xmid) [eqn (22)]
 
 !                     root3(pls) = a3*((sm(pls)+bminc_ind(pls)-xmid(pls)-((lm(pls)+xmid(pls))/ltor)+&
 !                     &rm(pls)+hm(pls))/wooddens)/pi4-((sm(pls)+bminc_ind(pls)-xmid(pls)-&
 !                     &((lm(pls)+xmid(pls))/ltor)+rm(pls))/((lm(pls)+xmid(pls))*sla(pls)*&
 !                     &wooddens/latosa))**a2
 
 !                     fmid(pls) = root3(pls)
 
 !                     if (fmid(pls) * sign(pls) .le. 0.) rtbis(pls) = xmid(pls)
 
 !                     if (dx(pls) < xacc .or. abs(fmid(pls)) <= yacc) exit
 
 !                     if (i > 20) write(stdout,*)'second alloc loop flag',i,pls,dx(pls),abs(fmid(pls))
 !                     if (i > 50) stop 'Too many iterations allocmod'
 
 !                     i = i + 1
 
 !                 end do
 
 !                 !Now rtbis contains numerical solution for lminc_ind given eqn (22)
 
 !                 lminc_ind(pls) = rtbis(pls)
 
 !             endif
 
 !             !Calculate increments in other compartments using allometry relationships
 
 !             rminc_ind(pls) = (lm(pls) + lminc_ind(pls)) / ltor - rm(pls)       !eqn (9)
 
 !             sminc_ind(pls) = bminc_ind(pls) - rminc_ind(pls) - lminc_ind(pls)  !eqn (1)
 
 !             ! print*, 'LEAF_INC (gC/ind)', (lminc_ind(pls)/1.D3), 'ROOT_INC (gC/ind)', (rminc_ind(pls)/1.D3),&
 !             ! & 'SAP_INC(gC/ind)', (sminc_ind(pls)/1.D3), pls, 'NORMAL'
 
 !         else 
 
 !             !Abnormal allocation: reduction in some C compartment(s) to satisfy allometry
             
 !             normal = .false.
 
 !             !Attempt to distribute this year's production among leaves and roots only
 
 !             lminc_ind(pls) = (bminc_ind(pls)-lm(pls)/ltor+rm(pls))/(1.+1./ltor)  !eqn (33)
 
 
 !             if (lminc_ind(pls) > 0.) then
 
 !                 !Positive allocation to leafmass
 
 !                 rminc_ind(pls) = bminc_ind(pls) - lminc_ind(pls)  !eqn (31)
                 
 !                 !Add killed roots (if any) to below-ground litter
 
 !                 if (rminc_ind(pls) < 0.) then
 
 !                     lminc_ind(pls) = bminc_ind(pls)
 !                     rminc_ind(pls) = (lm(pls) + lminc_ind(pls)) / ltor - rm(pls)
 
 !                     litter_bg(pls) = litter_bg(pls) + abs(rminc_ind(pls)) * nind(pls)
 
 !                 end if
                 
 !                 i = 1
 
 !             else
 
 !                 !Negative allocation to leaf mass
 
 !                 rminc_ind(pls) = bminc_ind(pls)
 !                 lminc_ind(pls) = (rm(pls) + rminc_ind(pls)) * ltor - lm(pls)  !from eqn (9)
 
 !                 !Add killed leaves to litter
 
 !                 litter_ag_fast(pls) = litter_ag_fast(pls) + abs(lminc_ind(pls)) * nind(pls)
                 
 !                 i = 2
 
 !             endif
 
 !             !Calculate sminc_ind (must be negative)
       
 !             sminc_ind(pls) = (lm(pls) + lminc_ind(pls)) * sla(pls) /&
 !             & latosa * 2.e5 * height(pls) - sm(pls)  !eqn (35)
 
 !             !Convert killed sapwood to heartwood
 
 !             hm(pls) = hm(pls) + abs(sminc_ind(pls))
 
 !             print*, pls, 'ANNORMAL'
 
 
 !         endif !normal/abnormal allocation
 
 !         !Increment C compartments - OUTPUT FINAL (kgC/m²)
 
 !         lm_ind(pls) = ((lm(pls) + lminc_ind(pls))*nind(pls))/1.D3
 !         rm_ind(pls) = ((rm(pls) + rminc_ind(pls))*nind(pls))/1.D3 
 !         sm_ind(pls) = ((sm(pls) + sminc_ind(pls))*nind(pls))/1.D3
 !         hm_ind(pls) = (hm(pls)*nind(pls))/1.D3
 !         cwood(pls) = sm_ind(pls)+hm_ind(pls)
 !         ! print*, 'leaf carbon',lm_ind(pls)
 !         ! print*, 'LMINC', lminc_ind(pls), pls
 !         !print*, 'LM', lm_ind(pls), 'RM', rm_ind(pls), 'SM', sm_ind(pls), 'HM', hm_ind(pls), 'CWOOD', cwood(pls), pls
 
 !         !ALLOMETRY EQUATIONS
 
 !         !DIAMETER (m)
 !         diameter(pls) = (4*(cwood(pls)*1.0D3)/(dwood(pls)*1D7)*pi*allom2)&
 !         &**(1/(2+allom3))
         
 !         !CROWN AREA (m2)
 !         crown_area(pls) = allom1*(diameter(pls)**1.6)
 !         crown_area_ind(pls) = (crown_area(pls)/nind(pls))
 
 !         !LAI (m2/m2)
 !         lai_ind(pls)=(((lm_ind(pls)/nind(pls))*sla(pls))/crown_area_ind(pls))
         
 !         !ALTURA (m)
 !         ! height(pls) = allom2*(diameter(pls)**allom3)
 
         
 !         fpc_ind(pls) = 1. - exp(-0.5 * 8.)
 !         fpc_grid(pls) = crown_area(pls) * nind(pls) * fpc_ind(pls)
 !         fpc_grid_old(pls) = fpc_grid(pls)
 !         fpc_inc(pls) = max(fpc_grid(pls) - fpc_grid_old(pls),0.)
 
 !         !SELF-THINNING LOGIC 
 !         fpc_tree_max = grid_area*0.95
 
 !         fpc_inc_tree    = sum(fpc_inc(:))
 !         ! fpc_tree_total  = sum(fpc_grid, mask = present .and. tree)
 
 
 
 !         ! print*, 'FPC_IND', fpc_ind(pls), 'FPC GRID', fpc_grid(pls)
 !     enddo
 ! end subroutine allocation

