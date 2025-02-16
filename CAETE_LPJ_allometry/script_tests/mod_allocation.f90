 
 ! program allocation !module to test allocation logic module of LPJmL-Fire
 subroutine allocation(lm_1, wm_1,hm_1, sm_1, rm_1, dwood_1, sla_1, nind_1,&
    &bminc_1, height_1, cl_inc,cw_inc,ch_inc, cs_inc, cr_inc, ctotal_inc,&
    &lm_2, ch_2, cs_2, rm_2, cw_2) 

    
    use iso_fortran_env, only : output_unit
    use types


    !VARIABLES [INPUT] - vem do módulo de estabelecimento/FPC
    real(r_8), intent(in) :: lm_1, wm_1, rm_1, hm_1, sm_1, height_1 !valor atualizado do carbono do ano passado.
    real(r_8), intent(in) :: sla_1, nind_1, dwood_1 !variable traits
    real(r_8), intent(in) :: bminc_1 !valor aleatório usado tbm no estabelecimento efbmin FPC posterior deve ser a NPP acumulada do ano
    

    !VARIABLES [OUTPUT] - o incremento de cada compartimento de carbono neste ano
    real(r_8), intent(out) :: cl_inc !leaf increment (gC)
    real(r_8), intent(out) :: cw_inc !wood increment (gC)
    real(r_8), intent(out) :: cr_inc !root increment (gC)
    real(r_8), intent(out) :: cs_inc !sapwood increment (gC)
    real(r_8), intent(out) :: ch_inc !heartwood increment (gC)
    real(r_8), intent(out) :: ctotal_inc !total carbon increment (gC)
    real(r_8), intent(out) :: lm_2, ch_2, cs_2, rm_2, cw_2

    integer(i_4),parameter :: nseg = 20 ! segment number to bisection loop
    integer(i_4), parameter :: time = 1000
    real(r_8),parameter :: pi   =  3.14159265
    real(r_8),parameter :: xacc =  0.1     !x-axis precision threshold for the allocation solution
    real(r_8),parameter :: yacc =  1.e-10  !y-axis precision threshold for the allocation solution
    integer(i_4),parameter :: ntl=365
    integer(i_4), parameter :: stdout = output_unit


    integer(i_4) :: pls
    integer(i_4) :: xmin, xmax
    real(r_8) :: allom1 = 100
    real(r_8) :: allom2 = 40.0
    real(r_8) :: allom3 = 0.5
    real(r_8) :: latosa = 8000.0
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
    real(r_8)  :: height
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
    logical  :: normal

    real(r_8)  :: fx1
    real(r_8)  :: fmid


    real(r_8)  :: lm1     !allometric leafmass requirement (leafmass req'd to keep sapwood alive; gC ind-1)

    integer(i_4) :: i

   !print*, 'inside alloc', lm_1, wm_1,hm_1, sm_1 
   
   !initializing variables
    cl_inc = 0.
    cw_inc = 0. !wood increment (gC)
    cr_inc = 0. !root increment (gC)
    cs_inc = 0.!sapwood increment (gC)
    ch_inc = 0.!heartwood increment (gC)
    ctotal_inc = 0.!total carbon increment (gC)
    lm_2 = 0.
    ch_2 = 0.
    cs_2 = 0.
    rm_2 = 0.
    cw_2 = 0.

    dwood = dwood_1*1000000 !*1e6 transforms dwood to gC/m3

    bminc = bminc_1 

    sla = sla_1

    nind = nind_1 

    height = height_1 
 
    lm_ind = lm_1

    cw_ind = wm_1

    hm_ind = hm_1

    sm_ind = sm_1

    rm_ind = rm_1 

    !-----------------------------------------------------------------
   ! print*,'inside', height

   

        lm  = (lm_ind)
        sm  = sm_ind !(cw_ind *0.1) !provisorio
        hm  = hm_ind !(cw_ind *0.9) !provisorio
        rm  = rm_ind
        bminc_ind  = bminc 
     
        ! ====== TREE ALLOCATION ======

        lm1 = (latosa*sm)/(dwood*height*sla) !allometric leaf mass requirement *****ATENÇÃO*****
       !  print*, 'LM1', lm1, 'lm', lm, sm, dwood, height, sla, lm, latosa
       !  print*, 'latosa*sm', latosa*sm, 'denominador', dwood*height*sla, (latosa*sm)/(dwood*height*sla)
       !  print*, 'lm1', lm1
        
        lminc_ind_min  = lm1 - lm  !eqn (27)
       !  if(lminc_ind_min.ge.0.)then
       !     print*, 'LM MIN', lminc_ind_min, lm, lm1
       !  endif
       !   print*, 'LM MIN', lminc_ind_min


    
        !calculate minimum root production to support this leaf mass (i.e. lm_ind + lminc_ind_min)
        !May be negative following a reduction in soil water limitation (increase in lm2rm) relative to last year.

        rminc_ind_min  = lm1  / ltor - rm       !eqn (30)
       !  print*, 'RM MIN', rminc_ind_min 

       ! rminc_ind_min  = (latosa*sm /(dwood)*height *sla *ltor) - rm       !eqn (30)
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
           !  print*, 'x1', x2
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

                   !  if (i > 20) write(stdout,*)'first alloc loop flag',i,pls,fmid *fx1 ,&
                   !       &xmid ,x1 ,x2 ,dx ,bminc_ind 

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

                   !  if (i > 20) write(stdout,*)'second alloc loop flag',i,pls,dx ,abs(fmid )
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
            !print*, 'anormal'
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
               !  print*, 'lminc', lminc_ind
                !Add killed leaves to litter

                litter_ag_fast  = litter_ag_fast  + abs(lminc_ind ) * nind 
                
                i = 2

            endif

            !Calculate sminc_ind (may be negative)
      
            sminc_ind  = (lm  + lminc_ind ) * sla  /&
            & latosa * dwood * height  - sm   !eqn (35)

            !Convert killed sapwood to heartwood

            hm  = hm  + abs(sminc_ind )

           !  print*, pls, 'ANNORMAL'


        endif !normal/abnormal allocation

        !Increment on C compartments 
        cl_inc = lminc_ind
        !print*, 'clinc', cl_inc !GOTOLITTER???
        cr_inc = rminc_ind
        !print*, 'cRinc', cr_inc
        cs_inc = sminc_ind
        !print*, 'csinc', cs_inc
        ch_inc = hm
        !print*, 'chinc', hm
        cw_inc = cs_inc + ch_inc

        
       !  print*, 'ctotal_inc', ctotal_inc
        if(cl_inc.le.0.)then
           cl_inc = 0. !this leaf goes to litter
        endif

        if(cs_inc.le.0.) then
           cs_inc = 0. !this sap goes to heart (I think it already done in the code)
        endif
        ctotal_inc = cl_inc + cr_inc + cs_inc + ch_inc

        !Compartments plus increment
        
        lm_2 = lm  + cl_inc

        ch_2 = hm + ch_inc

        cs_2 = sm + cs_inc

        rm_2 = rm + cr_inc

        cw_2 = cs_2 + ch_2
        
       !  print*, 'after alloc with inc'
       !  print*, 'l', lm_2
       !  print*, 'h', ch_2
       !  print*, 's', cs_2
       !  print*, 'r', rm_2
       !  if(cl_inc.lt.0) then
       !     print*, cl_inc, lm_2, lm
       !  endif

       !  if(ch_inc.lt.0) thenhm
       !     print*, ch_inc, ch_2, hm
       !  endif

       !  if(cr_inc.lt.0) then
       !     print*, cr_inc, rm_2, rm
       !  endif

       ! if(cs_inc.lt.0) then
       !     print*, cs_inc, cs_2, sm, cl_inc
       !  endif

      
        

end subroutine allocation    
        