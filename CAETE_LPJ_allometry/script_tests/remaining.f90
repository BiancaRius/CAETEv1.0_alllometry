program self_thinning

    ! ================= VARIABLES TO USE DECLARATION ===================== !
    integer :: j,k
    integer, parameter:: npls = 20
    real, dimension(npls), allocatable :: lai (:) !Leaf Area Index (m2/m2)
    real, dimension(npls), allocatable :: diam (:) !Tree diameter in m. (Smith et al., 2001 - Supplementary)
    real, dimension(npls), allocatable :: crown_area (:) !Tree crown area (m2) (Sitch et al., 2003)
    real, dimension(npls) :: dens_1 != 2  !densidade de indivíduos em cada PLS (0.5 é um número genérico)
    real, dimension(npls) :: dens_2 != 2  !densidade de indivíduos em cada PLS (0.5 é um número genérico)

   
    
    real, dimension(npls) :: est_pls !establishment for a specific PLS
    real, allocatable :: FPC_ind (:) !Foliage projective cover for each average individual of a PLS (Stich et al., 2003)
    real, allocatable :: FPC_pls_1 (:) !Total Foliage projective cover of a PLS (Stich et al., 2003)
    real, allocatable :: FPC_pls_2 (:) !Total Foliage projective cover of a PLS (Stich et al., 2003)
    real, allocatable :: FPC_pls (:) !Total Foliage projective cover of a PLS (Stich et al., 2003)
    real, allocatable :: FPC_grid_t1 (:) !Fractional projective cover for each PLS in time 1 (Sitch et al., 2003)
    real, allocatable :: FPC_grid_t2 (:) !Fractional projective cover for each PLS in time 2 (Sitch et al., 2003)
    real, allocatable :: fpc_dec (:) !'decrease FPC' in other words: the excedend of FPC in eac year [LPJ-GUESS - Phillip's video]
    real, allocatable :: fpc_dec_prop (:) !proportion of the decrease of a PLS
    real, allocatable :: mort (:) !equivalente ao 'mort_shade' no LPJ-GUESS [Phillip's video]
    real, allocatable :: mort_greff (:) ! motallity from growth efficiency (Sitch et al 2003)
    real, allocatable :: greff (:) ! motallity from growth efficiency (Sitch et al 2003)
    real, allocatable :: remaining (:) !taxa de redução
    real, allocatable :: FPC_inc (:)
    real, allocatable :: FPC_inc_cont (:)
    real, allocatable :: carbon_increment (:) ! used to calculate mort greff (Sitch et al 2003)
    
    real :: FPC_total_1 = 0.0 !sum of FPC_grid 
    real :: FPC_total_2 = 0.0 !sum of FPC_grid in
    
    real :: FPC_total_accu_1 = 0.0
    real :: FPC_total_accu_2 = 0.0

    real :: gc_area = 15. !grid cell size - 15 m2 FOR TESTING PURPOSE (the real value will be 1ha or 10000 m2)
    
    real :: fpc_max_tree !95% of grid-cell (in m2)
    real :: exc_area
    

    !Parameters and constants
    real :: k_allom1 = 100. !allometric constant (Table 3; Sitch et al., 2003)
    real :: krp = 1.6 !allometric constant (Table 3; Sitch et al., 2003)
    real :: ltor = 0.77302587552347657
    real :: k_est = 0.06 !establishment constant !Smith et al 2001 - Table A1
    real :: leaf_allocation = 0.35 !% of NPP allocaed to leaves
    real :: wood_allocation = 0.3  !% of NPP allocated to wood
    real :: root_allocation = 0.35 !% of NPP allocated to roots
    real :: k_mort1 = 0.01 !mortality parameter from Sitch et al 2003
    real :: k_mort2 = 0.3

    !Variables to allocation prototype
    real, dimension(npls) :: npp1 !KgC/ano
    real, dimension(npls) :: npp_inc  !incremento anual de C para cada PLS
    real, dimension(npls) :: annual_npp !quantidade de NPP com os incrementos.
    real, dimension(npls) :: cl2 !carbon on leaves after allocation
    real, dimension(npls) :: cw2 !carbon on wood after allocation
    real, dimension(npls) :: cr2 !carbon on wood after allocation

    ! Variables with generic values for testing the logic code
    real, dimension(npls) :: dwood !wood density (g/cm-3) *Fearnside, 1997 - aleatory choices
    real, dimension(npls) :: cw1 !KgC/m2 (Cheart + Csap)
    real, dimension(npls) :: cl1 !KgC/m2 
    real, dimension(npls) :: cr1 !KgC/m2
    real, dimension(npls) :: spec_leaf !m2/gC
    real, dimension(npls) :: leaf_inc !kgC/ ind
    real, dimension(npls) :: wood_inc !kgC/ ind
    real, dimension(npls) :: root_inc !kgC/ ind
    real, dimension(npls) :: diameter !
    

    
    !ARRAYS WITH GENERIC VALUES OF SOME VARIABLES
    dwood=(/0.24,0.53,0.39,0.32,0.31,0.44,0.66,0.42,0.74,0.39,0.82,0.40,0.26,0.79,0.39,0.52,0.41,0.44,0.86,0.42/) !atenção para a unidade
    cl1=(/.7,1.,0.3,1.6,1.10,1.8,0.3,0.2,0.8,0.84,0.25,1.,0.2,1.7,0.4,.6,.5,.8,0.3,1.8/)
    spec_leaf=(/0.0153,0.0101,0.0107,0.0112,0.012,0.0141,0.0137,0.0115,0.0122,0.010,0.012,0.011,&
    &0.013,0.014,0.0112,0.012,0.0141,0.0137,0.0115,0.0122/)
    cw1=(/30.,22.,34.,28.3,20.2,19.7,27.5,19.5,20.,28.6,24.3,19.3,26.8,22.,18.3,22.,15.,22.6,10.7,21.4/)
    cr1=(/0.63,0.8,0.9,1.5,1.3,0.9,0.4,1.0,0.56,0.87,0.33,0.97,0.31,0.55,0.2,0.8,0.4,0.66,0.23,1.5/)
    npp1 = (/0.5,0.8,1.5,1.2,1.9,1.3,1.7,0.8,0.6,2.0,0.7,1.1,1.9,1.85,1.96,1.77,1.33,1.54,1.62,0.55/)
    diameter = (/0.16,0.45,0.17,0.25,0.34,0.4,0.23,0.49,0.37,0.5,0.53,0.12,0.75,0.22,0.63,0.31,0.41,0.63,0.52,0.15/)
    dens_1 = (/1.,2.,5.,3.3, 1.3, 7., 2.8, 3.,4.5,1.7,3.6,9.,4.,2.45,5.27,4.6,8.2,9.29,3.,4.8/)

  

    allocate (FPC_ind(1:npls))
    allocate (FPC_pls_1(1:npls))
    allocate (FPC_pls_2(1:npls))
    allocate (FPC_pls(1:npls))
    allocate (FPC_grid_t1(1:npls))
    allocate (FPC_grid_t2(1:npls))
    allocate (fpc_dec(1:npls))
    allocate (fpc_dec_prop(1:npls))
    allocate (FPC_inc(1:npls))
    allocate (FPC_inc_cont(1:npls))
    allocate (remaining(1:npls))
    allocate (mort(1:npls))
    allocate (mort_greff(1:npls))
    allocate (greff(1:npls))
    allocate (diam(1:npls))
    allocate (crown_area(1:npls))
    allocate (lai(1:npls))
    allocate (carbon_increment(1:npls))
   

    ! ================= END VARIABLES DECLARATION ===================== !

    ! ==================== ALLOMETRY EQUATIONS =========================!
    !        Increment of carbon on tissues per individual 

   

    ! do k = 1,4 !looping for years
    !     print*,'k', k
    !     do j = 1, npls !looping for pls'
    !         cl2(j) = cl1(j)*3
    !         print*, 'cl2', cl2(j)
    !         print*, 'dwood', dwood(j)
            

    !     enddo
    !     cl1 = cl2  !updating variable for the next year
    !     print*, 'cl1_atualizado', cl1
    ! enddo


!!!------------------------------------------------------
        !Transforms from kgC to gC (as in LPJ)

    do j = 1, npls ! print*, j
        ! print*, 'FPC_pls2', FPC_pls_2(j), j, 'dens', dens_1(j)

        cl1(j) = cl1(j)*1000.

        cw1(j) = cw1(j)*1000.

        cr1(j) = cr1(j)*1000.

        npp1(j) = npp1(j)*1000.

       

!----------------------------------------------------------
        !Define a general value for NPP increment per year 
        !Transforms from kgC to gC (as in LPJ)

        npp_inc(j) = 0.1*1000.

!Annual NPP available to allocation (??????? é essa NPP ou a NPP inc?)
        
        annual_npp(j) = ((npp1(j)/dens_1(j)) + npp_inc(j))

        ! print*, 'annual npp', annual_npp(j)/1000.

         !-------------------------------------------------------------------------------
         ! !Increments to each compartments per individual. Here, the NPP proportions allocated
         ! to each compartment is being used for testing purpose. The actual values will be calculated
         ! in allocation routine.

        leaf_inc(j) = leaf_allocation * annual_npp(j)

        root_inc(j) = root_allocation * annual_npp(j) 

        wood_inc(j) = wood_allocation * annual_npp(j)  

        carbon_increment(j) = leaf_inc(j) + root_inc(j) + wood_inc(j)
        print*, '1st year', carbon_increment(j)/1000.

!!---------------------------------------------------
        
!!!-------------------------------------------------------
!----------------------------------------------------------
        !Define a general value for FPC in order to initialize and 
        !use to calculate the FPC increments

        FPC_pls_1(j) = .1

        FPC_total_1 = FPC_total_1 + FPC_pls_1(j)
        !print*, 'FPC_total_1', FPC_total_1
       
        if (j.eq.npls) then
            FPC_total_accu_1 = FPC_total_1
            ! print*, 'FPC_total_accu_1', FPC_total_accu_1
        endif

    enddo

   
    do k = 1, 10

        print*, '**********************************************************'
        print*, '                                                           '
        print*, '                                                           '
        print*, '                                                           '
        print*, '                                                           '
        print*, 'year',k
       

        
        
        FPC_ind = 0.
        FPC_pls_2 = 0.
        FPC_total_2 = 0.
        FPC_total_accu_2 = 0.
        fpc_max_tree = 0.
        exc_area = 0.
        FPC_inc = 0.
        FPC_inc_cont = 0.
        FPC_dec = 0.
        fpc_dec_prop = 0.
        dens_2 = 0.
        cl2 = 0.
        cw2 = 0.
        cr2 = 0.
        lai = 0.
        crown_area = 0.
        diam = 0.
        annual_npp = 0.
        leaf_inc = 0.
        wood_inc = 0.
        root_inc = 0.
      
        

        do j = 1, npls
        
        !--------------------------------------------------------------------------
        !transforming the carbon content from gC/m2 to gc/average individual 
        !(the carbon divided by dens gives the individual carbon, as in LPJ)

            ! print*, '1st cl', cl1(j)/1000.
            if(cl1(j).eq.0) then
                cl2(j) = 0.

                cw2(j) =0.
    
                cr2(j) = 0.                             
                              
                diam(j) = 0.
            
                crown_area(j) = 0.
            
                lai(j) = 0.
           
                FPC_ind(j) = 0.
                
                FPC_pls_2(j) = 0.

                dens_1(j) = 0.

                dens_2(j) = 0.
               
            else
                cl2(j) = (cl1(j)/dens_1(j)) 

                cw2(j) = (cw1(j)/dens_1(j)) 

                cr2(j) = (cr1(j)/dens_1(j)) 

        
                  !----------------------------------------------------------------------------
                 !Structuring PLSs [diameter, crown area and leaf area index]

                diam(j) = ((4*(cw2(j)))/((dwood(j)*1000000.)*3.14*36))**(1/(2+0.22)) !nessa equação dwood deve estar em *g/m3*
            
                crown_area(j) = k_allom1*(diam(j)**krp)
            
                lai(j) = (cl2(j)*spec_leaf(j))/crown_area(j) 

                
                !------------------------------------------------------------------------------
                !---------------------------------------------------------------------------
                !Calculatin Foliage Projective Cover of average individual(FPC_ind), of the PLS(FPC_pls)
                ! and of the grid cell (FPC_total)

                FPC_ind(j) = (1-(exp(-0.5*lai(j))))
                ! print*, "FPC_ind", FPC_ind(j)
                
            
                FPC_pls_2(j) = crown_area(j) * dens_1(j) * FPC_ind(j) 

                
            endif

            ! print*, 'FPC_pls_2', FPC_pls_2(j),j

            FPC_total_2 = FPC_total_2 + (FPC_pls_2(j)) !accumulate the values in the variable FPC_total.
                                                        !the actual value will only be obtained when j = npls

            if (j.eq.npls) then   !take the value accumulated until the last pls
              
                FPC_total_accu_2 = FPC_total_2
                ! print*, 'FPC_total_accu_2', FPC_total_accu_2, j

            endif  
        
        enddo

        !--------------------------------------------------------------------------- 
        !Verifying if FPCs together occupy more than 95% of the plot area
               
        fpc_max_tree = gc_area*0.95 !utilizaremos 1 ha !! 5% é destinado ao novo estabelecimento
                
        !print*, 'fpc_max_tree', fpc_max_tree
        
        ! print*, 'FPC_total_accu_2', FPC_total_accu_2, 'FPC_pls_2', FPC_pls_2
        if (FPC_total_accu_2 .gt. fpc_max_tree) then
                    
            print*, 'ultrapassou'
            ! Excedent area
        
            exc_area = FPC_total_accu_2 - fpc_max_tree

            ! print*, exc_are

            !------------------------------------------------------------

       
            do j = 1, npls
                ! print*, 'FPC_pls2', FPC_pls_2(j)
                if(FPC_pls_2(j).eq.0.)then
                    FPC_inc(j) = 0.
                    FPC_inc_cont(j) = 0.
                    fpc_dec(j) = 0.                   
                    fpc_dec_prop(j) = 0.               
                    greff(j) = 0.
                    mort(j) = 0.
                    ! print*, 'dead PLS', j
                 
                else
                    ! print*, j
                    !Calculating the increment of PLS from a year to the next
            
                

                    FPC_inc(j) = FPC_pls_2(j) - FPC_pls_1(j)
               
                    !Calculating the relative contribution to total increment considering all PLSs

                    FPC_inc_cont(j) = (FPC_inc(j)/(FPC_total_accu_2-FPC_total_accu_1))
                
               
                    !    Calculating the percentage of FPC reduction of each PLS in relation to the area excedent

                    fpc_dec(j) = (exc_area)*(FPC_inc_cont(j))

       
                    !calculating shade mortality
                    !!!ATENTION: include the other mortality sources

               
                    fpc_dec_prop(j) = (((FPC_pls_2(j) - fpc_dec(j))/FPC_pls_2(j)))

                    
                    
                
                    greff(j) = carbon_increment(j)/(cl2(j)*spec_leaf(j))

                    mort_greff(j) = k_mort1/(1+(k_mort2*greff(j)))
                    ! print*, 'mort_greff', mort_greff(j), j

                    mort(j) = 1.0 - (fpc_dec_prop(j)+mort_greff(j))
                    print*, 'mort', mort(j)
                   
                endif    
              

                if (mort(j).lt.0.)then !maximum mortality in this case
                    
                    mort(j) = 1.
                                  
                endif        
 

                

            enddo

        else
            print*, 'n ultrapassou'
            !if the occupation is smaller than the stand area the mortality is defined only by
            !the growth efficiency and the loss of carbon through turnover
            do j=1, npls
                greff(j) = carbon_increment(j)/(cl2(j)*spec_leaf(j))

                mort_greff(j) = k_mort1/(1+(k_mort2*greff(j)))
                
                mort(j) = 1.0 - mort_greff(j)
                ! print*, 'mort_greff', mort_greff(j), j
                ! print*, 'mort', mort(j)
                !fpc_dec(j) = 0.
            enddo
            
            !print*, 'NAO ULTRAPASSOU ==', 'este FPC_total_accu_2 tem que ser igual ao valor anterior==', FPC_total_accu_2

            !exc_area = 0.


        endif
        
        
        
        do j=1,npls
            ! print*, 'mort',mort(j)   
            remaining(j) = 1.0 - mort(j)
            print*, remaining(j)
            if (remaining(j) .le. 0.) then
                ! print*, 'PLS dead===============================================================',j
                ! goto 10 
                dens_2(j) = 0.
                cl2(j) = 0.
                cw2(j) = 0.
                cr2(j) = 0.
                FPC_pls_2(j) = 0.
                FPC_total_accu_2 = 0. 
                npp_inc(j) = 0
                annual_npp(j) = 0.
                leaf_inc(j) = 0.
                root_inc(j) = 0.
                wood_inc(j) = 0.
                cl1(j) = 0.
                cw1(j) = 0.
                cr1(j) = 0.
            endif

            ! print*, remaining(j) , j

            dens_2(j) = dens_1(j) * remaining(j)

            ! print*, 'dens_2 (pos remaining)', dens_2(j), 'dens_1', dens_1(j), remaining(j), mort(j)
            ! print*, '                            '
            ! print*, '                            '
            ! print*, '                            '

            cl2(j) = cl2(j) * remaining(j)
            print*, 'cl2', cl2(j)/1000.

            cw2(j) = cw2(j) * remaining(j)

            cr2(j) = cr2(j) * remaining(j)

            !!carbon to litter
            ! carbon_litter_leaf(j) = cl2(j) - ((cl2(j))*remaining(j))
            ! print*,  'litter', carbon_litter_leaf(j), cl2(j), remaining(j)
            !add litter of previous time step
            !the litter when pls dies is equal to cl, cw, cr
            
            ! print*, 'cl2 remaining', cl2(j)/1000., 'remaining', remaining(j)
            ! print*, j
     
            ! print*, 'cw2', cw2(j)/1000., 'cw1', cw1(j)/1000.
         
            ! print*, 'cr2', cr2(j)/1000., 'cr1', cr1(j)/1000.
            ! print*, j
        enddo
       

        ! !---------------------------------------------------------------------------
       
        !----------------------------------------------------------------------------
        !updating the variables for the next year

        cl1 = cl2
        
        ! print*, 'cl1 atualizado p/ a aloca', cl1/1000.

        cw1 = cw2

        ! print*, 'cw1 atualizado', cw1/1000.
        
        cr1 = cr2

        ! print*, 'cr1 atualizado', cr1/1000.

        !FPCs

        FPC_pls_1 = FPC_pls_2

        ! print*, 'FPC_atualizado', FPC_pls_1


        FPC_total_accu_1 = FPC_total_accu_2
        

        dens_1 = dens_2

    

        
        !-----------------------------------------------------------------------------
        !!!---------------------------------------------------------------------------
        !!!Fictitious allocation process in order to test the logic developed
         !------------------------------------------------------------------------------
        !NPP increment (NPPt-NPPt-1); for testing purpose a general value was defined
        !Transforms NPP increment from m2 to NPP increment for each averge individual
        
        do j=1,npls

           
            if(dens_1(j).le.0.) then
                npp_inc(j) = 0.
            
                annual_npp(j) = 0.

                leaf_inc(j) = 0.

                root_inc(j) = 0.

                wood_inc(j) = 0.

                cl1(j) = 0.

                cw1(j) = 0.

                cr1(j) = 0.

                ! delta_carbon_pls(j) = 0.

                ! print*, 'delta', delta_carbon_pls(j), j
            
            else
            
            
            
                npp_inc(j) = npp_inc(j)/dens_1(j)

           

            !-------------------------------------------------------------------------------
            !Annual NPP available to allocation (??????? é essa NPP ou a NPP inc?)
        
                annual_npp(j) = ((npp1(j)/dens_1(j)) + npp_inc(j))

            ! print*, 'annual npp', annual_npp(j)/1000.

             !-------------------------------------------------------------------------------
             ! !Increments to each compartments per individual. Here, the NPP proportions allocated
             ! to each compartment is being used for testing purpose. The actual values will be calculated
             ! in allocation routine.

                leaf_inc(j) = leaf_allocation * annual_npp(j)

                root_inc(j) = root_allocation * annual_npp(j) 

                wood_inc(j) = wood_allocation * annual_npp(j)  

                carbon_increment(j) = leaf_inc(j) + root_inc(j) + wood_inc(j)
                ! print*, 'final', carbon_increment(j)/1000.

                cl1(j) = cl1(j) + leaf_inc(j)

                cw1(j) = cw1(j) + wood_inc(j)

                cr1(j) = cr1(j) + root_inc(j)
            endif

            ! print*, 'cl1 com incremento após aloca', cl1(j)/1000., j, npp_inc(j), dens_1(j)
            !print*, j

            !print*, 'densidade p/ ano seguinte =======', dens_1(j)

            cl1(j) = cl1(j) * dens_1(j)
            cw1(j) = cw1(j) * dens_1(j)
            cr1(j) = cr1(j) * dens_1(j)

            npp_inc(j) = npp_inc(j) * dens_1(j)

            carbon_increment(j) = carbon_increment(j)

            ! delta_carbon_pls(j) = delta_carbon_pls(j)

            ! print*, '==============================='
            ! print*,'cl final', cl1(j)/1000., j

        enddo

        ! FPC_total_2 = 0.
        !FPC_total_accu_2 = 0.
        !dens_1 = dens_2

    enddo

end program self_thinning 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
!   do k = 1, 3 !Loop dos anos
!         if (k .eq. 1) then !*usando o time step t1 - INITIAL ALLOMETRY*

!             !Carbon on tissues (wood and leaf) per average-individual (this considers the individual density [dens])
!             cw1 = (cw1/dens)*1000. !*1000 transforma de kgC para gC - carbono
            
!             cl1 = (cl1/dens)*1000.  !*1000 transforma de kgC para gC - carbono
           
!             cr1 = (cr1/dens)*1000. !*1000 transforma de kgC para gC - carbono

!             !PLS structure [diam, crown area and leaf area index]
!             diam = ((4*(cw1))/((dwood*1000000.)*3.14*36))**(1/(2+0.22)) !nessa equação dwood deve estar em *g/m3*
!             crown_area = k_allom1*(diam**krp)
!             lai = (cl1*spec_leaf)/crown_area 
            
!             !print*, crown_area, lai
!             !Net Primary Productive logic (may be per average-individual). NPP ano atual - NPP ano anterior
!             ! #1 NPP increment to each average-individual. Cada individuo médio de cada PLS vai ter um incremento
!             npp_inc = 0.15/dens !this increment is fixed in 0.15kgC/yr (valor generalizado)

!             ! #2 Annual NPP disponible to alloc 
!             annual_npp = ((npp1/dens) + npp_inc)*1000 !*1000 transforma de kgC para gC

!             !==================================================
!             !INCREMENTS TO LEAF AND WOOD TISSUES PER INDIVIDUO

!             leaf_inc = 0.35*npp_inc ![35% of NPP] - esse valor de 0.35 é calculado pela alocação
!             root_inc = 0.35*npp_inc ![35% of NPP - Functional balance]
!             wood_inc = 0.3*npp_inc  ![30% of NPP]

!             !==================================================
!             !CARBON TISSUES (atualização dos tecidos)

!             cl2 = cl1 + leaf_inc !cl1 e leaf_inc já está dividido peela densidade
!             cw2 = cw1 + wood_inc !cw1 e wood_inc já está dividido peela densidade
!             cr2 = cr1 + root_inc

!             !==================================================
!             !Foliage Projective Cover (FPC_ind) & Fractional Projective Cover (FPC_grid) - calculado após
!             ! a alocação e a estruturação do PLS
            
!             FPC_ind = (1-exp(-0.5*lai)) !FPC ind médio [m2]
!             FPC_grid_t1 = crown_area*dens*FPC_ind !FPC of pls [occupation on grid cell considers all average-individual of PLS; m2]
!             FPC_total_t1 = sum(FPC_grid_t1)


!         else !*usando time step 2*

!             !PLS structure [diam, crown area and leaf area index]
!             diam = ((4*(cw2))/((dwood*1000000.)*3.14*36))**(1/(2+0.22)) !nessa equação dwood deve estar em *g/m3*
!             crown_area = k_allom1*(diam**krp)
!             lai = (cl2*spec_leaf)/crown_area 
!             !print*, 'diametro_old', diam

!             !==================================================
!             !Foliage Projective Cover (FPC_ind) & Fractional Projective Cover (FPC_grid)               
            
!             FPC_ind = (1-exp(-0.5*lai)) !FPC ind médio [m2]
!             FPC_grid_t2 = crown_area*dens*FPC_ind !FPC of PLS [occupation on grid cell considers all average-individual of PLS; m2]
!             FPC_total_t2 = sum(FPC_grid_t2)
!             ! print*, 'inc summing t1', FPC_total_t2
            

!         endif

!         ! print*, k, 'diametro_old', diam, 'cl2_old', cl2, 'cw2_old', cw2, 'dens_ind_old', dens

!         fpc_max_tree = gc_area*0.95 !utilizaremos 1 ha !! 5% é destinado ao novo estabelecimento

!         if (FPC_total_t2 .gt. fpc_max_tree) then
!             !Self-Thinning imposed when total tree cover above 'FPC MAX TREE' [max. of occupation]
!             !partitioned among tree PLS in proportion to this year's FPC increment

!             !Area excedent
!             exc_area = FPC_total_t2 - fpc_max_tree
!             ! print*, 'exc_area', exc_area

!             !Contribuição excedente de cada PLS
!             FPC_inc = FPC_grid_t2 - FPC_grid_t1 !!incremento de cada PLS ! delta entre tempo 1 e 2
!             ! print*, 'soma dos incrementos', sum(FPC_inc), 'diferença t2t1', FPC_total_t2-FPC_total_t1

!             FPC_inc_cont = (FPC_inc/(FPC_total_t2-FPC_total_t1))  !contribuição relativa de cada PLS para o incremento total
!             ! print*, 'FPC_INC_pls cont====', FPC_inc_cont

!             fpc_dec = (exc_area)*(FPC_inc_cont) !porcentagem em relação ao exc de area que cada pls tem que reduzir o fpc
!             mort = 1.0 - (((FPC_grid_t2-fpc_dec)/FPC_grid_t2)) 
            
!             !incluir outros fatores da mortalidade
            
!             do j=1,npls
!                 if (mort(j).gt.1.)then
!                     mort(j) = 1.
!                 endif
!             enddo
!             ! print*, 'mort total', mort, 'mort', ((FPC_grid_t2-fpc_dec)/FPC_grid_t2), 'fpd dec',fpc_dec, sum(FPC_dec)
!         else
!             mort = 0.0
!         endif
!         remaining = 1.0-mort
!         ! print*, 'remainig', remaining
!         dens = dens*remaining 
!         cl2 = cl2*remaining 
!         cw2 = cw2*remaining 
!         cr2 = cr2*remainig
        
!         carbon_pls = cl2 + cw2 + cr2

!         carbon_grid_cell = sum(carbon_pls)

!         ! print*, 'carbon_pls', carbon_pls, 'total carbon', carbon_grid_cell/1000.

      


!         !print*,'dens update', dens, 'cl2 updt', cl2, 'cw2 updt', cw2
        
!         do j = 1, npls
!             if (remaining(j) .le. 0.) then
                   
!                 FPC_grid_t2(j) = 0.0
!                     ! dens(j) = 0.0
!                     ! diam(j) = 0.0
!                     ! crown_area(j) = 0.0
!                     ! lai(j) = 0.0
!                     ! cl2(j) = 0.0
!                     ! cw2 (j) = 0.0
!                     ! FPC_ind(j) = 0.0
!                     ! remaining(j) = 0.0 
! ! 
!             else
! !          recalculate  PLS structure [diam, crown area and leaf area index]
!                 diam(j) = ((4*(cw2(j)))/((dwood(j)*1000000.)*3.14*36))**(1/(2+0.22)) !nessa equação dwood deve estar em *g/m3*
!                 crown_area(j) = k_allom1*(diam(j)**krp)
!                 lai(j) = (cl2(j)*spec_leaf(j))/crown_area(j) 
! ! 
!                     ! ==================================================
!                     ! Foliage Projective Cover (FPC_ind) & Fractional Projective Cover (FPC_grid)               
!                     ! 
!                 FPC_ind(j) = (1-exp(-0.5*lai(j))) !FPC ind médio [m2]
!                 FPC_grid_t2(j) = crown_area(j)*dens(j)*FPC_ind(j) !FPC of PLS [occupation on grid cell considers all average-individual of PLS; m2]
!                 FPC_grid_total_t2 = sum(FPC_grid_t2)
! ! 
!             endif
!         enddo
!         ! print*, k, 'SOMA_TOTAL',sum(FPC_grid_t2),'FPC pls', FPC_grid_t2, '95% area', fpc_max_tree, 'diametro', diam
!         ! print*, 'exc_area', exc_area
!         FPC_perc = FPC_grid_total_t2/gc_area
        
        
        
!        ! gen_est = k_est*(1-exp(-5*(1-FPC_perc)))*(1-FPC_perc) !Eqn 17 in Smith et al 2001 (se ocupação > 90%)
        
!         gen_est = k_est*(1-FPC_perc)


        
!         !print*, 'general establish', gen_est, FPC_grid_total_t2, 'FPC_perc', FPC_perc

!         est_pls = gen_est*(carbon_pls/carbon_grid_cell)*FPC_grid_t2*(1-FPC_perc) !saplings/yr
        
!         !após o estabelecimento, os saplings devem ser adicionados à densidade de individuos
!         !N_new = N_old + est_pls
!         !update os pools de C

!         ! if (k.ne.1) then
!         ! !   print*, 'EST_PLS', est_pls, sum(est_pls), '+ DENS', dens+est_pls
!         ! endif

!         !update dens
        
!         dens_new = dens + est_pls

!         !updating carbon pools according to the establishment of saplings (Eq. 20 in Smith 2001) -
!         !ps: como nosso modelo ainda não tem crescimento nós optamos por utilizar o valor de carbono do indivíduo médio
!         !do PLS para calcular quanto de carbono é incrementado para este compartimento
        
!         !!!!ATENÇÃO - DÚVIDA: isso mantém o balanço de carbono??????


!     enddo
!     !incremento da densidade populacional (Eqn 19 in SMith 2001)
!     !update of biomass compartments (Eqn 20 in Smith 2001)
!     !update of biomass compartments with establishment
!     !veja no vídeo do Philip no minuto 45:06
!     !inserir quantidade de carbono que vai pra liteira

