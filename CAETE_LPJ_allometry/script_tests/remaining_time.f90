program self_thinning

    ! ================= VARIABLES TO USE DECLARATION ===================== !
    integer :: j,k
    integer, parameter :: npls = 20 !40 !20
    integer, parameter :: time = 100
    real, dimension(npls,time) :: lai !Leaf Area Index (m2/m2)
    real, dimension(npls,time) :: diam !Tree diameter in m. (Smith et al., 2001 - Supplementary)
    real, dimension(npls,time) :: crown_area !Tree crown area (m2) (Sitch et al., 2003)

   
    
    ! real, dimension(npls) :: est_pls !establishment for a specific PLS
    real, dimension(npls,time) :: FPC_ind !Foliage projective cover for each average individual of a PLS (Stich et al., 2003)
    real, dimension(npls,time) :: FPC_pls_1  !Total Foliage projective cover of a PLS (Stich et al., 2003)
    real, dimension(npls,time) :: FPC_pls_2  !Total Foliage projective cover of a PLS (Stich et al., 2003)
    real, dimension(npls) :: fpc_dec  !'decrease FPC' in other words: the excedend of FPC in eac year [LPJ-GUESS - Phillip's video]
    real, dimension(npls) :: fpc_dec_prop  !proportion of the decrease of a PLS
    real, dimension(npls) :: mort  !equivalente ao 'mort_shade' no LPJ-GUESS [Phillip's video]
    real, dimension(npls) :: mort_greff  ! motallity from growth efficiency (Sitch et al 2003)
    real, dimension(npls) :: greff  ! motallity from growth efficiency (Sitch et al 2003)
    real, dimension(npls) :: remaining  !taxa de redução
    real, dimension(npls) :: FPC_inc 
    real, dimension(npls) :: FPC_inc_cont 
    real, dimension(npls) :: carbon_increment  ! used to calculate mort greff (Sitch et al 2003)

    
    real, dimension(time) :: FPC_total_1 = 0.0 !sum of FPC_grid 
    real, dimension(time) :: FPC_total_2 = 0.0 !sum of FPC_grid in
    
    real, dimension(time):: FPC_total_accu_1 = 0.0
    real, dimension(time) :: FPC_total_accu_2 = 0.0

    real :: gc_area = 1000. !grid cell size - 15 m2 FOR TESTING PURPOSE (the real value will be 1ha or 10000 m2)
    
    real :: fpc_max_tree !95% of grid-cell (in m2)
    real :: exc_area
    

    !Parameters and constants
    real :: k_allom1 = 100. !allometric constant (Table 3; Sitch et al., 2003)
    real :: krp = 1.6 !allometric constant (Table 3; Sitch et al., 2003)
    real :: ltor = 0.77302587552347657
    real :: k_est = 0.06 !establishment constant !Smith et al 2001 - Table A1
    real :: leaf_allocation = 0.4 !% of NPP allocaed to leaves
    real :: wood_allocation = 0.3  !% of NPP allocated to wood
    real :: root_allocation = 0.3 !% of NPP allocated to roots
    real :: k_mort1 = 0.01 !mortality parameter from Sitch et al 2003
    real :: k_mort2 = 0.3
    real :: res_time_leaf = 4 !general residence time value for testing purpose
    real :: res_time_root = 4
    real :: res_time_wood = 60 !ATENÇÃO! ESSE NUMERO PRECISA SER REVISADO POIS EM SITCH ET AL 2003 APENAS O SAPWOOD É PERDIDO POR TURNOVER

    !Variables to allocation prototype
    real, dimension(npls,time) :: npp_inc  !incremento anual de C para cada PLS
    real, dimension(npls,time) :: annual_npp !quantidade de NPP com os incrementos.
    real, dimension(npls,time) :: cl2 !carbon on leaves after allocation
    real, dimension(npls,time) :: cw2 !carbon on wood after allocation
    real, dimension(npls,time) :: cr2 !carbon on wood after allocation
    real, dimension(npls,time) :: dens2

    real, dimension(npls,time) :: cleaf_avg_ind
    real, dimension(npls,time) :: cwood_avg_ind
    real, dimension(npls,time) :: croot_avg_ind

    real, dimension(npls,time) :: cw1 !KgC/m2 (Cheart + Csap)
    real, dimension(npls,time) :: cl1 !KgC/m2 
    real, dimension(npls,time) :: cr1 !KgC/m2
    real, dimension(npls,time) :: dens1

    ! Variables with generic values for testing the logic code
    real, dimension(npls,time) :: dwood !wood density (g/cm-3) *Fearnside, 1997 - aleatory choices
    real, dimension(npls,time) :: spec_leaf !m2/gC
    real, dimension(npls) :: leaf_inc !kgC/ ind
    real, dimension(npls) :: wood_inc !kgC/ ind
    real, dimension(npls) :: root_inc !kgC/ ind

    !variables with initial values
    real, dimension(npls,time) :: cl1_initial
    real, dimension(npls,time) :: cw1_initial
    real, dimension(npls,time) :: cr1_initial
    real, dimension(npls,time) :: npp1_initial
    real, dimension(npls) :: FPC_pls_initial
    real, dimension(npls,time) :: dens1_initial
    

    !auxiliary variables for outputs
    real, dimension (npls,time) :: cl1_aux
    real, dimension (npls,time) :: cw1_aux
    real, dimension (npls,time) :: cr1_aux
    real, dimension (npls,time) :: FPC_pls_1_aux
    real, dimension (npls,time) :: dens1_aux


    !creating random numbers for npp increment
    
    real:: x(npls,time)
    
            
   
    
    ! dens_1 = (/1.,2.,5.,3.3, 1.3, 7., 2.8, 3.,4.5,1.7,3.6,9.,4.,2.45,5.27,4.6,8.2,9.29,3.,4.8/)!,&
    !&1.,2.,5.,3.3, 1.3, 7., 2.8, 3.,4.5,1.7,3.6,9.,4.,2.45,5.27,4.6,8.2,9.29,3.,4.8/)

  

  
   

    ! ================= END VARIABLES DECLARATION ===================== !

    ! ==================== ALLOMETRY EQUATIONS =========================!
    !        Increment of carbon on tissues per individual 


    
!!!------------------------------------------------------

! !creating value for initial density


    xmin = 1.7
    xmax = 9.5
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
            ! print*,x(j,k),j,k
        enddo
    enddo

    do j = 1, npls      

        dens1_initial(j,:) =x(j,:)
        ! print*, ' dens initial', dens1_initial(j,:)

    enddo
!________________________________________________________________
!!!------------------------------------------------------

! !creating value for initial cleaf


    xmin = 0.2
    xmax = 1.8
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
            ! print*,x(j,k),j,k
        enddo
    enddo

!________________________________________________________________
!!!!!Transforms from kgC to gC (as in LPJ)

    do j = 1, npls      

        cl1_initial(j,:) =x(j,:)*1000.
        ! print*, 'initial', cl1_initial(j,:)/1000.

    enddo

!_______________________________________________
!!    creating value for initial cwood
    xmin = 10.
    xmax = 35.
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
            ! print*, 'x cwood',x(j,k),j,k
        enddo
    enddo

    do j = 1, npls      

        cw1_initial(j,:) =x(j,:)*1000.
        ! print*, ' cw initial', cw1_initial(j,:)/1000.

    enddo

!____________________________________________________

!!    creating value for initial croot
    xmin = 0.2
    xmax = 1.5
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
            ! print*, 'x croot',x(j,k),j,k
        enddo
    enddo

    do j = 1, npls      

        cr1_initial(j,:) =x(j,:)*1000.
        ! print*, ' cr initial', cr1_initial(j,:)/1000.

    enddo

!____________________________________________________
   !!    creating value for initial npp
    xmin = 0.5
    xmax = 2.
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
            ! print*, 'x npp',x(j,k),j,k
        enddo
    enddo

    do j = 1, npls      

        npp1_initial(j,:) =x(j,:)*1000.
        ! print*, ' npp initial', npp1_initial(j,:)/1000.

    enddo

!____________________________________________________
   !!    creating value for dwood
    xmin = 0.24
    xmax = 0.86
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
            ! print*, 'x dwood',x(j,k),j,k
        enddo
    enddo

    do j = 1, npls      

        dwood(j,:) =x(j,:)
       

    enddo

!___________________________________________________________
 
!____________________________________________________
   !!    creating value for spec_leaf
    xmin = 0.0101
    xmax = 0.0153
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
            ! print*, 'x dwood',x(j,k),j,k
        enddo
    enddo

    do j = 1, npls      

        spec_leaf(j,:) =x(j,:)
       

    enddo

!___________________________________________________________
   
   
  
  
   
  

!______________________________________________
!!!!creating value for initial npp_inc
!__________________________________________

    xmin = 0.1
    xmax =  2.5
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
            ! print*,x(j,k),j,k
        enddo
    enddo  



!Annual NPP available to allocation (??????? é essa NPP ou a NPP inc?)
    do j = 1, npls
        npp_inc(j,:) = x(j,:)*1000
        ! print*,'___________________________________'
        ! print*, 'npp_inc', npp_inc(j,:)/1000.

        annual_npp(j,:) = ((npp1_initial(j,:)/dens1_initial(j,:)) + npp_inc(j,:))

        ! print*, 'annual npp', annual_npp(j,:)/1000.

         !-------------------------------------------------------------------------------
         ! !Increments to each compartments per individual. Here, the NPP proportions allocated
         ! to each compartment is being used for testing purpose. The actual values will be calculated
         ! in allocation routine.

        leaf_inc(j) = leaf_allocation * annual_npp(j,k)

        root_inc(j) = root_allocation * annual_npp(j,k) 

        wood_inc(j) = wood_allocation * annual_npp(j,k)  

        carbon_increment(j) = leaf_inc(j) + root_inc(j) + wood_inc(j)
        ! print*, '1st year', carbon_increment(j)/1000.,j

!!---------------------------------------------------
        
!!!-------------------------------------------------------
!----------------------------------------------------------
        !Define a general value for FPC in order to initialize and 
        !use to calculate the FPC increments

        FPC_pls_initial(j) = .1

        FPC_total_1(:) = FPC_total_1(:) + FPC_pls_initial(j)
        !print*, 'FPC_total_1', FPC_total_1
       
        if (j.eq.npls) then
            FPC_total_accu_1(:) = FPC_total_1(:)
            print*, 'FPC_total_accu_1', FPC_total_accu_1(k)
        endif

    enddo

   
    do k = 1, time
                                                   
        print*, 'year',k
       

        
        
        FPC_ind(:,k) = 0.
        FPC_pls_2(:,k) = 0.
        FPC_total_2(k) = 0.
        FPC_total_accu_2(k) = 0.
        fpc_max_tree = 0.
        exc_area = 0.
        FPC_inc = 0.
        FPC_inc_cont = 0.
        FPC_dec = 0.
        fpc_dec_prop = 0.

        cleaf_avg_ind(:,k) = 0.
        cwood_avg_ind(:,k) = 0.
        croot_avg_ind(:,k) = 0.


        cl2(:,k) = 0.
        cw2(:,k) = 0.
        cr2(:,k) = 0.
        lai(:,k) = 0.
        dens2(:,k) = 0.
        crown_area(:,k) = 0.
        diam(:,k) = 0.
        annual_npp(:,k) = 0.
        leaf_inc = 0.
        wood_inc = 0.
        root_inc = 0.
      
        if(k.eq.1)then
            cl1(:,k) = cl1_initial(:,k)
            cw1(:,k) = cw1_initial(:,k)
            cr1(:,k) = cr1_initial(:,k)
            FPC_pls_1(:,k) = FPC_pls_initial(:)
            dens1(:,k) = dens1_initial(:,k)
            ! print*,'cl previous yr', cl1(:,k)/1000.
        else

        cl1(:,k) = cl1_aux(:,k-1)
        cw1(:,k) = cw1_aux(:,k-1)
        cr1(:,k) = cr1_aux(:,k-1)
        dens1(:,k) = dens1_aux(:,k-1)
        FPC_pls_1(:,k) = FPC_pls_1_aux(:, k-1)
        ! print*, 'cl1 nxt year', cl1(1,k)/1000.
        
        endif

        do j = 1, npls
        
        !--------------------------------------------------------------------------
        !transforming the carbon content from gC/m2 to gc/average individual 
        !(the carbon divided by dens gives the individual carbon, as in LPJ)
            

            ! print*, '1st cl', cl1(j)/1000.
            if(cl1(j,k).eq.0) then
                ! print*, 'cl1 eq 0'
                cl2(j,k) = 0.

                cw2(j,k) =0.
    
                cr2(j,k) = 0.                             
                              
                diam(j,k) = 0.
            
                crown_area(j,k) = 0.
            
                lai(j,k) = 0.
           
                FPC_ind(j,k) = 0.
                
                FPC_pls_2(j,k) = 0.

                dens1(j,k) = 0.

                dens2(j,k) = 0.
               
            else
                ! print*, 'cl1 ne 0'
                cl2(j,k) = (cl1(j,k)/dens1(j,k)) 

                cw2(j,k) = (cw1(j,k)/dens1(j,k)) 

                cr2(j,k) = (cr1(j,k)/dens1(j,k)) 

        
                  !----------------------------------------------------------------------------
                 !Structuring PLSs [diameter, crown area and leaf area index]

                diam(j,k) = ((4*(cw2(j,k)))/((dwood(j,k)*1000000.)*3.14*36))**(1/(2+0.22)) !nessa equação dwood deve estar em *g/m3*
            
                crown_area(j,k) = k_allom1*(diam(j,k)**krp)
            
                lai(j,k) = (cl2(j,k)*spec_leaf(j,k))/crown_area(j,k) 

                
                !------------------------------------------------------------------------------
                !---------------------------------------------------------------------------
                !Calculatin Foliage Projective Cover of average individual(FPC_ind), of the PLS(FPC_pls)
                ! and of the grid cell (FPC_total)

                FPC_ind(j,k) = (1-(exp(-0.5*lai(j,k))))
                ! print*, "FPC_ind", FPC_ind(j)
                
            
                FPC_pls_2(j,k) = crown_area(j,k) * dens1(j,k) * FPC_ind(j,k) 

                
            endif

            
                ! print*,cl2(j), j
            
            
            ! print*, 'FPC_pls_2', FPC_pls_2(j),j

            FPC_total_2(k) = FPC_total_2(k) + (FPC_pls_2(j,k)) !accumulate the values in the variable FPC_total.
                                                        !the actual value will only be obtained when j = npls

            if (j.eq.npls) then   !take the value accumulated until the last pls
              
                FPC_total_accu_2(k) = FPC_total_2(k)
                ! print*, 'FPC_total_accu_2', FPC_total_accu_2, j

            endif  
        
        enddo

        !--------------------------------------------------------------------------- 
        !Verifying if FPCs together occupy more than 95% of the plot area
               
        fpc_max_tree = gc_area*0.95 !utilizaremos 1 ha !! 5% é destinado ao novo estabelecimento
                
        print*, 'fpc_max_tree', fpc_max_tree
        
        ! print*, 'FPC_total_accu_2', FPC_total_accu_2, 'FPC_pls_2', FPC_pls_2
        if (FPC_total_accu_2(k) .gt. fpc_max_tree) then
                    
            print*, 'ultrapassou'
            ! Excedent area
        
            exc_area = FPC_total_accu_2(k) - fpc_max_tree

            ! print*, exc_are

            !------------------------------------------------------------

       
            do j = 1, npls
                ! print*, 'FPC_pls2', FPC_pls_2(j)
                if(FPC_pls_2(j,k).eq.0.)then
                    FPC_inc(j) = 0.
                    FPC_inc_cont(j) = 0.
                    fpc_dec(j) = 0.                   
                    fpc_dec_prop(j) = 0.               
                    greff(j) = 0.
                    mort(j) = 1.
                    mort_greff(j) = 0.
                    ! print*, 'dead PLS', j
                 
                else
                    ! print*, j
                    !Calculating the increment of PLS from a year to the next
            
                

                    FPC_inc(j) = FPC_pls_2(j,k) - FPC_pls_1(j,k)
                    
                    if(FPC_inc(j).lt.0..or.FPC_total_accu_1(k).gt.FPC_total_accu_2(k))then
                        FPC_inc(j) = 0.
                        FPC_inc_cont(j) = 0.
                        fpc_dec(j) = 0.                   
                        fpc_dec_prop(j) = 0.               
                        greff(j) = 0.
                        mort(j) = 1.
                        mort_greff(j) = 0.
                        ! print*, 'dead PLSSSSSSSSSS', j
                    
                    else
                        !Calculating the relative contribution to total increment considering all PLSs

                        FPC_inc_cont(j) = (FPC_inc(j)/(FPC_total_accu_2(k)-FPC_total_accu_1(k)))
                        ! print*, 'inc_cont', FPC_inc_cont(j), j
                        ! print*,''
                        ! print*, 'FPC inc',FPC_inc(j),j
                        ! print*,''
                        ! print*, 'fpc total accu 2', FPC_total_accu_2
                        ! print*,''
                        ! print*, 'fpc total accu 1', FPC_total_accu_1
               
                        !    Calculating the percentage of FPC reduction of each PLS in relation to the area excedent

                        fpc_dec(j) = (exc_area)*(FPC_inc_cont(j))

       
                    
                        !!!ATENTION: include the other mortality sources

               
                        fpc_dec_prop(j) = (((FPC_pls_2(j,k) - fpc_dec(j))/FPC_pls_2(j,k))) !calculating shade mortality

                                      
                        greff(j) = carbon_increment(j)/(cl2(j,k)*spec_leaf(j,k)) !growth efficiency

                        mort_greff(j) = k_mort1/(1+(k_mort2*greff(j))) !mortality by gowth efficiency

                        ! print*, 'mort_greff', mort_greff(j), j

                        mort(j) = (fpc_dec_prop(j)+mort_greff(j)) !sum of all mortality
                        ! print*, 'mort', mort(j), 'fpc_decprop', fpc_dec_prop(j),'fpc_dec',fpc_dec(j),j
                        ! print*, 'fpc inc',FPC_inc(j),'FPCpls2', FPC_pls_2(j),'mort_greff', mort_greff(j), j
                   
                    endif

                endif    
              

                if (mort(j).lt.0.)then !maximum mortality in this case
                    
                    mort(j) = 1.
                                  
                endif        
 

                

            enddo

        else
            print*, 'n ultrapassou', FPC_total_accu_2(k)
            !if the occupation is smaller than the stand area the mortality is defined only by
            !the growth efficiency and the loss of carbon through turnover
            do j=1, npls
                if(cl2(j,k).eq.0.) then
                    greff(j) = 0.
                    mort_greff(j) = 0.
                    mort(j) = 1.
                    ! print*,'cl2 eq 0'

                else    
                    greff(j) = carbon_increment(j)/(cl2(j,k)*spec_leaf(j,k))

                    mort_greff(j) = k_mort1/(1+(k_mort2*greff(j)))
                
                    mort(j) = mort_greff(j)
                    ! print*, 'greff', greff(j), carbon_increment(j)/1000., cl2(j)/1000., spec_leaf(j)
                    ! print*, 'mort_greff', mort_greff(j), j
                    ! print*, 'mort', mort(j)
                    !fpc_dec(j) = 0.
                endif
            enddo
            
            !print*, 'NAO ULTRAPASSOU ==', 'este FPC_total_accu_2 tem que ser igual ao valor anterior==', FPC_total_accu_2

            !exc_area = 0.


        endif
        
        
        
        do j=1,npls

            if(mort(j).gt.1.) then !maximum mortality is equal to 1
                mort(j) = 1.
            else
                mort(j) = mort(j)
            endif

            ! print*, 'mort',mort(j)   
            remaining(j) = 1.0 - mort(j)
           
            ! print*, 'remaining', remaining(j), 'mort', mort(j), j
           
            if (remaining(j) .le. 0.) then
                ! print*, 'PLS dead===============================================================',j
                ! goto 10 
                dens2(j,k) = 0.
                cl2(j,k) = 0.
                cw2(j,k) = 0.
                cr2(j,k) = 0.
                FPC_pls_2(j,k) = 0.
                FPC_total_accu_2(k) = 0. 
                npp_inc(j,k) = 0
                annual_npp(j,k) = 0.
                leaf_inc(j) = 0.
                root_inc(j) = 0.
                wood_inc(j) = 0.
                cl1(j,k) = 0.
                cw1(j,k) = 0.
                cr1(j,k) = 0.
            endif

            ! print*, remaining(j) , j

            dens2(j,k) = dens1(j,k) * remaining(j)

            ! print*, 'dens_2 (pos remaining)', dens_2(j), 'dens_1', dens_1(j), remaining(j), mort(j)
            ! print*, '                            '
            ! print*, '                            '
            ! print*, '                            '

            cl2(j,k) = cl2(j,k) * remaining(j)         

            cw2(j,k) = cw2(j,k) * remaining(j)

            cr2(j,k) = cr2(j,k) * remaining(j)



            !! Loss of carbon through residence time
           

            cl2(j,k) = cl2(j,k) - (cl2(j,k)/res_time_leaf)
            ! print*, 'cl2 after restime', cl2(j)/1000., (cl1(j)/res_time)/1000.
            ! print*, ''

            cw2(j,k) = cw2(j,k) - (cw2(j,k)/res_time_wood)

            cr2(j,k) = cr2(j,k) - (cr2(j,k)/res_time_root)


            !!carbon to litter
            ! carbon_litter_leaf(j) = cl2(j) - ((cl2(j))*remaining(j))
            ! print*,  'litter', carbon_litter_leaf(j), cl2(j), remaining(j)
            !add litter of previous time step
            !the litter when pls dies is equal to cl, cw, cr
        
        enddo


         !----------------------------------------------------------------------------
        !updating the variables for the next year

        cl1_aux(:,k) = cl2(:,k)
        
        ! print*, 'cl1_aux', cl1_aux(:,k)/1000.

        cw1_aux(:,k) = cw2(:,k)

        ! print*, 'cw1 atualizado', cw1/1000.
        
        cr1_aux(:,k) = cr2(:,k)

        ! print*, 'cr1 atualizado', cr1/1000.


        dens1_aux(:,k) = dens2(:,k)


        !FPCs

        FPC_pls_1_aux(:,k) = FPC_pls_2(:,k)

        ! print*, 'FPC_atualizado', FPC_pls_1


        FPC_total_accu_1(k)= FPC_total_accu_2(k)
        ! print*, 'FPC_atualizado', FPC_pls_1

       
        !-----------------------------------------------------------------------------
        !!!---------------------------------------------------------------------------
        !!!Fictitious allocation process in order to test the logic developed
         !------------------------------------------------------------------------------
        !NPP increment (NPPt-NPPt-1); for testing purpose a general value was defined
        !Transforms NPP increment from m2 to NPP increment for each averge individual
        
        do j=1,npls

           
            if(dens1_aux(j,k).le.0.) then
                npp_inc(j,k) = 0.
            
                annual_npp(j,k) = 0.

                leaf_inc(j) = 0.

                root_inc(j) = 0.

                wood_inc(j) = 0.

                cl1_aux(j,k) = 0.

                cw1_aux(j,k) = 0.

                cr1_aux(j,k) = 0.

                ! delta_carbon_pls(j) = 0.

                ! print*, 'delta', delta_carbon_pls(j), j
            
            else
            
            
            
                npp_inc(j,k) = npp_inc(j,k)/dens1_aux(j,k)

           

            !-------------------------------------------------------------------------------
            !Annual NPP available to allocation (??????? é essa NPP ou a NPP inc?)
        
                annual_npp(j,k) = ((npp1_initial(j,k)/dens1_aux(j,k)) + npp_inc(j,k))

            ! print*, 'annual npp', annual_npp(j)/1000.

             !-------------------------------------------------------------------------------
             ! !Increments to each compartments per individual. Here, the NPP proportions allocated
             ! to each compartment is being used for testing purpose. The actual values will be calculated
             ! in allocation routine.

                leaf_inc(j) = leaf_allocation * annual_npp(j,k)

                root_inc(j) = root_allocation * annual_npp(j,k) 

                wood_inc(j) = wood_allocation * annual_npp(j,k)  

                carbon_increment(j) = leaf_inc(j) + root_inc(j) + wood_inc(j)
                ! print*, 'final', carbon_increment(j)/1000.

                cl1_aux(j,k) = cl1_aux(j,k) + leaf_inc(j)
                
                cw1_aux(j,k) = cw1_aux(j,k) + wood_inc(j)
                
                cr1_aux(j,k) = cr1_aux(j,k) + root_inc(j)

                !saving value for avg individual for outputs
                cleaf_avg_ind(j,k) = cl1_aux(j,k)
                cwood_avg_ind(j,k) = cw1_aux(j,k)
                croot_avg_ind(j,k) = cr1_aux(j,k)
                
            endif

            ! print*, 'cl1 com incremento após aloca', cl1_aux(j,k)/1000., j
            print*, 'cl avg ind', cleaf_avg_ind(j,k)/1000.
            print*, 'cw avg ind', cwood_avg_ind(j,k)/1000.
            print*, 'cr avg ind', croot_avg_ind(j,k)/1000.



            ! print*, ''
            ! print*, 'cw1 com incremento após aloca', cw1(j)/1000., j
            ! print*, ''
            ! print*, 'cr1 com incremento após aloca', cr1(j)/1000., j

            !print*, 'densidade p/ ano seguinte =======', dens_1(j)

            cl1_aux(j,k) = cl1_aux(j,k) * dens1_aux(j,k)
            print*, 'cl * dens', cl1_aux(j,k)/1000
            cw1_aux(j,k) = cw1_aux(j,k) * dens1_aux(j,k)
            print*, 'cw * dens', cw1_aux(j,k)/1000
            cr1_aux(j,k) = cr1_aux(j,k) * dens1_aux(j,k)
            print*, 'cr * dens', cr1_aux(j,k)/1000

            npp_inc(j,k) = npp_inc(j,k) * dens1_aux(j,k)

            carbon_increment(j) = carbon_increment(j)

            ! delta_carbon_pls(j) = delta_carbon_pls(j)

            ! print*, '==============================='
            ! print*,'cl final', cl1(j)/1000., j

        enddo


    enddo

    open(unit=1,file='carbon_pools_time.csv',status='unknown')
    do k=1, time
        do j = 1,npls

       
            write(1,*) cl1_aux(j,k)/1000.,',',cw1_aux(j,k)/1000.,',', cr1_aux(j,k)/1000.,',','pls',j,',',k !newline
        enddo
    enddo    

    close(1)

    open(unit=2,file='FPC_total_time.csv',status='unknown')
    do k=1, time
        

       
            write(2,*) FPC_total_accu_2(k),',',k !newline
        
    enddo    

    close(2)


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

