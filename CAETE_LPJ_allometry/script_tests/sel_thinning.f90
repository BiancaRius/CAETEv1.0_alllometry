!Este módulo calcula as variáveis para indivíduos médios de cada PLS
!Designa-se também a calcular o FPC (foliar projective cover) através
!da estruturação alométrica dos PLSs. Além disso, calcula a mortalidade por espaço,
!ou seja, considerando uma área disponível de 1ha, os PLSs lenhosos só podem ocupar
!95% desta área. Caso ultrapassem esse montante, haverá penalização (redução na densidade 
!de indivíduos) que impactará os compartimentos de carbono. Caso não ultrapassem, será 
!possível que novos indivíduos se estabeleçam. Tais indivíduos são chamados de "sapling",
!apresentam alometria própria, mas não há cohort, de modo que ao incorporá-los na população do 
!PLSs, deve haver um ajustamento no balanço de carbono, o "shrink", que representa uma reestruturação
!dos indivíduos médios.
!Este código é baseado principalmente no modelo LPJ (Sitch et al 2003 e Smith et al 2001)
program self_thinning

    
    ! use csv_file
    use establish
    use types
    ! ================= VARIABLES TO USE DECLARATION ===================== !
    
    integer(i_4) :: j !dimensão de PLSs
    integer(i_4) :: k !dimensão de tempo
    ! logical :: grass !gramíneas não estão sendo incorporadas no momento
    
    integer(i_4) :: file_unit
    integer(i_4), parameter :: npls = 3000
    integer(i_4), parameter :: time = 500

    ! integer, parameter :: grassess = 0.1*npls !gramíneas não estão sendo incorporadas no momento
    real(r_8), dimension(npls,time) :: lai !Leaf Area Index (m2/m2)
    real(r_8), dimension(npls,time) :: diam !Tree diameter in m. (Smith et al., 2001 - Supplementary)
    real(r_8), dimension(npls,time) :: crown_area !Tree crown area (m2) (Sitch et al., 2003)
    real(r_8), dimension(npls,time) :: height !Height (m) (Sitch et al., 2003)

    
    ! real(r_8), dimension(npls) :: est_pls !establishment for a specific PLS
    real(r_8), dimension(npls,time) :: FPC_ind    !Foliage projective cover for each average individual of a PLS (Stich et al., 2003)
    real(r_8), dimension(npls,time) :: FPC_pls_1  !Total Foliage projective cover of a PLS (Stich et al., 2003)
    real(r_8), dimension(npls,time) :: FPC_pls_2  !Total Foliage projective cover of a PLS (Stich et al., 2003)
    real(r_8), dimension(npls,time) :: FPC_dec    !Quantidade de área que deve ser reduzida de um PLS em cada ano 
    real(r_8), dimension(npls,time) :: nind_kill_FPC !Número de indivíduos que devem morrer para que o PLS cumpra sua redução
    real(r_8), dimension(npls,time) :: nind_kill_greff !número de indivíduos a serem mortos por ineficiência de crescimento
    real(r_8), dimension(npls,time) :: nind_kill_total !Número total de indivíduos que devem morrer (somando todas as fontes de mortalidade)
    real(r_8), dimension(npls,time) :: mort  !% a ser descontada dos compartimentos de carbono por conta da mortalidade
    real(r_8), dimension(npls,time) :: mort_greff  ! mortality from growth efficiency takes into account mort by wd (Sakschewski et al 2015)
    real(r_8), dimension(npls,time) :: mort_wd  ! fator de mortalidade que entra na mort_greff considerando o wood density (Sakschewski et al 2015)
    real(r_8), dimension(npls,time) :: greff  ! cálculo de eficiência de crescimento (Sitch et al 2003)
    real(r_8), dimension(npls,time) :: remaining  !1 - mort (% multiplicada pelos estados dos PLSs para saber quanto do que existe deve permanecer)  real(r_8), dimension(npls,time) :: mort_wd
    real(r_8), dimension(npls,time) :: FPC_inc !incremento em área para um PLS entre um ano e outro
   
    real(r_8), dimension(npls,time) :: carbon_increment_initial  ! used to calculate mort greff (Sitch et al 2003)
    real(r_8), dimension(npls,time) :: carbon_increment ! used to calculate mort greff (Sitch et al 2003)
    ! real(r_8), dimension (npls) :: FPC_inc_grid
    
    real(r_8), dimension(time) :: FPC_total_initial = 0.0 !sum of FPC_grid
    real(r_8), dimension(time) :: FPC_total_accu_initial = 0.0 !sum of FPC_grid  

    real(r_8), dimension(time) :: FPC_total_2 = 0.0 !sum of FPC_grid in
    real(r_8) :: dead_pls, alive_pls, count_pls
    
    real(r_8), dimension(time):: FPC_total_accu_1 = 0.0
    real(r_8), dimension(time) :: FPC_total_accu_2 = 0.0

    real(r_8) :: gc_area = 10000!grid cell size - 15 m2 FOR TESTING PURPOSE (the real(r_8) value will be 1ha or 10000 m2)

    real(r_8), dimension(time) :: gc_available

    real(r_8) :: fpc_max_tree !95% of grid-cell (in m2)
    real(r_8), dimension(time) :: exc_area   
    

    !Parameters and constants
    real(r_8) :: k_allom1 = 100. !allometric constant (Table 3; Sitch et al., 2003)
    real(r_8) :: k_allom2 = 40.0
    real(r_8) :: k_allom3 = 0.85
    real(r_8) :: krp = 1.6 !allometric constant (Table 3; Sitch et al., 2003)
    real(r_8) :: ltor = 0.77302587552347657
    real(r_8) :: k_est = 0.06 !establishment constant !Smith et al 2001 - Table A1
    real(r_8) :: leaf_allocation = 0.4 !% of NPP allocaed to leaves
    real(r_8) :: wood_allocation = 0.3  !% of NPP allocated to wood
    real(r_8) :: root_allocation = 0.3 !% of NPP allocated to roots
    real(r_8) :: k_mort1 = 0.01 !mortality parameter from Sitch et al 2003
    real(r_8) :: k_mort2 = 0.5
    real(r_8) :: res_time_leaf = 2 !general residence time value for testing purpose
    real(r_8) :: res_time_root = 2
    real(r_8) :: res_time_sap = 2
    real(r_8) :: res_time_wood = 40 !ATENÇÃO! ESSE NUMERO PRECISA SER REVISADO POIS EM SITCH ET AL 2003 APENAS O SAPWOOD É PERDIDO POR TURNOVER
    real(r_8) :: crown_area_max = 30 !m2 !number from lplmfire code (establishment.f90)
    real(r_8) :: pi = 3.1415

    !Variables to allocation prototype
    real(r_8), dimension(npls,time) :: npp_inc, npp_inc2  !incremento anual de C para cada PLS
    real(r_8), dimension(npls,time) :: npp_inc_init  !incremento anual de C para cada PLS

    real(r_8), dimension(npls,time) :: annual_npp !quantidade de NPP com os incrementos.
    real(r_8), dimension(npls,time) :: cl2 !carbon on leaves after allocation
    real(r_8), dimension(npls,time) :: cw2 !carbon on wood after allocation
    real(r_8), dimension(npls,time) :: cs2 !carbon on sapwood after allocation
    real(r_8), dimension(npls,time) :: ch2 !carbon on heartwood after allocation
    real(r_8), dimension(npls,time) :: cr2 !carbon on wood after allocation
    real(r_8), dimension(npls,time) :: dens2

    real(r_8), dimension(npls,time) :: cleaf_avg_ind
    real(r_8), dimension(npls,time) :: csap_avg_ind
    real(r_8), dimension(npls,time) :: cheart_avg_ind
    real(r_8), dimension(npls,time) :: cwood_avg_ind
    real(r_8), dimension(npls,time) :: croot_avg_ind

    real(r_8), dimension(npls,time) :: cw1 !KgC/m2 (Cheart + Csap)
    real(r_8), dimension(npls,time) :: cs1 !KgC/m2 (Cheart + Csap)
    real(r_8), dimension(npls,time) :: ch1 !KgC/m2 (Cheart + Csap)
    real(r_8), dimension(npls,time) :: cl1 !KgC/m2 
    real(r_8), dimension(npls,time) :: cr1 !KgC/m2
    real(r_8), dimension(npls,time) :: dens1

    ! Variables with generic values for testing the logic code
    real(r_8), dimension(npls,time) :: dwood !wood density (g/cm-3) *Fearnside, 1997 - aleatory choices
    real(r_8), dimension(npls,time) :: spec_leaf !m2/gC
    real(r_8), dimension(npls) :: leaf_inc !kgC/ ind
    real(r_8), dimension(npls) :: wood_inc !kgC/ ind
    real(r_8), dimension(npls) :: root_inc !kgC/ ind

    !variables with initial values
    real(r_8), dimension(npls,time) :: cl1_initial
    real(r_8), dimension(npls,time) :: cw1_initial
    real(r_8), dimension(npls,time) :: cs1_initial
    real(r_8), dimension(npls,time) :: ch1_initial
    real(r_8), dimension(npls,time) :: cr1_initial
    real(r_8), dimension(npls,time) :: npp1_initial
    real(r_8), dimension(npls,time) :: FPC_pls_initial
    real(r_8), dimension(npls,time) :: dens1_initial
    

    !auxiliary variables for outputs
    real(r_8), dimension (npls,time) :: cl1_aux
    real(r_8), dimension (npls,time) :: cw1_aux
    real(r_8), dimension (npls,time) :: ch1_aux
    real(r_8), dimension (npls,time) :: cs1_aux
    real(r_8), dimension (npls,time) :: cr1_aux
    real(r_8), dimension (npls,time) :: FPC_pls_1_aux
    real(r_8), dimension (npls,time) :: dens1_aux
    real(r_8), dimension (time) :: FPC_total_accu_1_aux


    !creating random numbers for npp increment
    
    real(r_8):: x(npls,time)

    !variables for module of establishment

    real(r_8), dimension (time) :: est
    real(r_8), dimension (npls,time) :: est_pls
    real(r_8), dimension (npls,time) :: cleaf_sapl
    real(r_8), dimension (npls,time) :: csap_sapl
    real(r_8), dimension (npls,time) :: cheart_sapl
    real(r_8), dimension (npls,time) :: croot_sapl
    real(r_8), dimension (npls,time) :: dens_est
    real(r_8), dimension (npls,time) :: cleaf_new
    real(r_8), dimension (npls,time) :: cwood_new
    real(r_8), dimension (npls,time) :: csap_new
    real(r_8), dimension (npls,time) :: cheart_new
    real(r_8), dimension (npls,time) :: croot_new

     !variables for module allocation
    real(r_8), dimension (npls,time) :: cl_inc !leaf increment from allocation (gC, average_in)
    real(r_8), dimension (npls,time) :: cr_inc !root increment from allocation (gC, average_in)
    real(r_8), dimension (npls,time) :: cw_inc !wood increment from allocation (gC, average_in)
    real(r_8), dimension (npls,time) :: ch_inc !heart increment from allocation (gC, average_in)
    real(r_8), dimension (npls,time) :: cs_inc !sap increment from allocation (gC, average_in)
    real(r_8), dimension (npls,time) :: ctotal_inc !total increment from allocation (gC, average_in)
    real(r_8), dimension (npls,time) :: cl2_aux = 0. !leaf carbon after allocation (gC, avera_in)
    real(r_8), dimension (npls,time) :: cw2_aux = 0. !wood carbon after allocation (gC, avera_in)
    real(r_8), dimension (npls,time) :: ch2_aux = 0. !heartwood carbon after allocation (gC, avera_in)
    real(r_8), dimension (npls,time) :: cs2_aux = 0.!sapwood carbon after allocation (gC, avera_in)
    real(r_8), dimension (npls,time) :: cr2_aux = 0. !root carbon after allocation (gC, avera_in)
   

   
  
    

    ! ================= END VARIABLES DECLARATION ===================== !

    ! ==================== ALLOMETRY EQUATIONS =========================!
    !        Increment of carbon on tissues per individual 

!!!------------------------------------------------------

! !creating value for initial density
!! this value is used to calculate the inc in tissues for the 1st day


    xmin = 0.05
    xmax = 5.
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
        enddo
    enddo

    do j = 1, npls      
        dens1_initial(j,:) =x(j,:)
    enddo
!________________________________________________________________
!!!------------------------------------------------------

! !creating value for initial cleaf


    xmin = 0.2
    xmax = 2.5
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
        enddo
    enddo

!________________________________________________________________
!!!!!Transforms from kgC to gC (as in LPJ)

    do j = 1, npls      

        cl1_initial(j,:) =x(j,:)*1000.
       
    enddo

!_______________________________________________
!!    creating value for initial cheartwood
    xmin = 5.
    xmax = 80.
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)

        enddo
    enddo

!!!!!Transforms from kgC to gC (as in LPJ)
    do j = 1, npls      

        ch1_initial(j,:) =x(j,:)*1000.
        
    enddo

!____________________________________________________
!____________________________________________________
!    creating value for initial csapwood
    xmin = 0.2
    xmax = 2.5
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)

        enddo
    enddo
!!!!!Transforms from kgC to gC (as in LPJ)
    do j = 1, npls      

        cs1_initial(j,:) =x(j,:)*1000.
        
    enddo


!Total carbon in wood tissues is the sum of csap and cheart
    do j = 1, npls      

        cw1_initial(j,:) = cs1_initial(j,:) + ch1_initial(j,:)
        
    enddo
!____________________________________________________


!!    creating value for initial croot
    xmin = 0.2
    xmax = 2.5
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
            
        enddo
    enddo
!!!!!Transforms from kgC to gC (as in LPJ)
    do j = 1, npls      

        cr1_initial(j,:) =x(j,:)*1000.
        
    enddo

!____________________________________________________
   !!    creating value for initial npp
    xmin = 0.5
    xmax = 3.5
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
           
        enddo
    enddo
!!!!!Transforms from kgC to gC (as in LPJ)
    do j = 1, npls      

        npp1_initial(j,:) =x(j,:)*1000.
        
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
            
        enddo
    enddo

    do j = 1, npls      

        spec_leaf(j,:) =x(j,:)
       

    enddo

!___________________________________________________________
 
 

!______________________________________________
!!!!creating value for initial npp_inc_init
!__________________________________________
!____________________________
    !!calculates increment npp to be allocated
    xmin = 0.5
    xmax = 3.5
     
    x(:,:) = 0.
    call random_number(x)

    do k = 1, time
        do j = 1, npls
            x(j,k) = xmin + (xmax-xmin)*x(j,k)
            
        enddo
    enddo  
!!!!!Transforms from kgC to gC (as in LPJ)
    do j = 1, npls
        npp_inc_init(j,:) = (x(j,:)*1000)
    enddo
        
        ! annual_npp(j,:) = (npp1_initial(j,:)/dens1_initial(j,:)) + (npp_inc_init(j,:)/dens1_initial(j,:))

       

!-------------------------------------------------------------------------------
    ! !Increments to each compartments per individual. Here, the NPP proportions allocated
    ! to each compartment is being used for the first day. For the rest of the days it is
    ! calculated from allocation routine

    do j=1, npls
        leaf_inc(j) = (leaf_allocation * npp_inc_init(j,k))/dens1_initial(j,k)

        root_inc(j) = (root_allocation * npp_inc_init(j,k))/dens1_initial(j,k)
 
    !here wood_inc do not consider sap and heart separetely    
        wood_inc(j) = (wood_allocation * npp_inc_init(j,k))/dens1_initial(j,k)
  

        carbon_increment_initial(j,k) = leaf_inc(j) + root_inc(j) + wood_inc(j)
    enddo   

!!---------------------------------------------------
        
!!!-------------------------------------------------------
!----------------------------------------------------------
        !Define a general value for FPC in order to initialize and 
        !use to calculate the FPC increments
    do k = 1, time
        do j=1,npls

            FPC_pls_initial(j,k) = 1.

            FPC_total_initial(k) = FPC_total_initial(k) + FPC_pls_initial(j,k)

        
            if (j.eq.npls) then
                FPC_total_accu_initial(k) = FPC_total_initial(k)
                ! print*, 'FPC_total_accu_initial', FPC_total_accu_initial(k),k
            endif

        enddo

        ! if(k.eq.1)then
        !     print*,'FPC_total_accu_initial', FPC_total_accu_initial(k),k
        ! endif

    enddo

    dead_pls = 0.
    alive_pls = 0. 
    count_pls = 0.
    do k = 1, time
                                                  
        ! print*, 'year',k
       
        
        
        
        FPC_ind(:,k) = 0.
        FPC_pls_2(:,k) = 0.
        FPC_total_2(k) = 0.
        FPC_total_accu_2(k) = 0.
        fpc_max_tree = 0.
        exc_area(k) = 0.
        FPC_inc(:,k) = 0.
        FPC_dec (:,k) = 0.
        nind_kill_FPC (:,k) = 0.
        nind_kill_greff(:, k) = 0.

        cleaf_avg_ind(:,k) = 0.
        csap_avg_ind(:,k) = 0.
        cheart_avg_ind(:,k) = 0.
        cwood_avg_ind(:,k) = 0.
        croot_avg_ind(:,k) = 0.


        cl2(:,k) = 0.
        cw2(:,k) = 0.
        cs2(:,k) = 0.
        ch2(:,k) = 0.
        cr2(:,k) = 0.
        lai(:,k) = 0.
        dens2(:,k) = 0.
        crown_area(:,k) = 0.
        diam(:,k) = 0.
        annual_npp(:,k) = 0.
        leaf_inc = 0.
        wood_inc = 0.
        root_inc = 0.
        mort(:,k) = 0.
        ! nind_kill_prop(:,k) = 0.
        npp_inc2(:,k) = 0.
        cl2_aux(:,k) = 0.
        cw2_aux(:,k) = 0.
        cs2_aux(:,k) = 0.
        ch2_aux(:,k) = 0.
        cr2_aux(:,k) = 0.
      
        if(k.eq.1)then
            cl1(:,k) = cl1_initial(:,k)
            cw1(:,k) = cw1_initial(:,k)
            cs1(:,k) = cs1_initial(:,k)
            ch1(:,k) = ch1_initial(:,k)
            cr1(:,k) = cr1_initial(:,k)
            FPC_pls_1(:,k) = FPC_pls_initial(:,k)
            dens1(:,k) = dens1_initial(:,k)
            FPC_total_accu_1(k) = FPC_total_accu_initial(k)
            ! print*, 'fpc total accu 1, entrando no modelo', FPC_total_accu_1(k)
            npp_inc(:,k) = npp_inc_init(:,k)
            ! print*, 'npp inc init', npp_inc_init(:,k)
            carbon_increment(:,k) = carbon_increment_initial(:,k)
           
        else

            cl1(:,k) = cl1_aux(:,k-1)
            cw1(:,k) = cw1_aux(:,k-1)
            cs1(:,k) = cs1_aux(:,k-1) !- provisorio
            ch1(:,k) = ch1_aux(:,k-1) !- provisorio
            cr1(:,k) = cr1_aux(:,k-1)
            dens1(:,k) = dens1_aux(:,k-1)
            FPC_pls_1(:,k) = FPC_pls_1_aux(:, k-1)
            FPC_total_accu_1(k) = FPC_total_accu_1_aux(k-1)
            !!!atenção
            carbon_increment(:,k) = carbon_increment(:,k-1)

        
            npp_inc(:,k)=0. !reinitializing for a new sampling
        
            xmin = 0.1
            xmax = 3.5
     
            x(:,:) = 0.
            call random_number(x)

           
            do j = 1, npls
                x(j,k) = xmin + (xmax-xmin)*x(j,k)
                
            enddo
            npp_inc(:,k) = x(:,k)*1000
        
        endif

           

        do j = 1, npls
                
                !--------------------------------------------------------------------------
        !transforming the carbon content from gC/m2 to gc/average individual 
        !(the carbon divided by dens gives the individual carbon, as in LPJ)
            
            !if(cl1(j,k).le.0) then
                !print*, 'cl1 0',cl1(j,k),dens1(j,k), cr1(j,k), diam(j,k), height(j,k)
            if(dens1(j,k).le.1.e-1) then !densidade mínima de indivíduos (from LPJmfire code)
                !print*, cl1(j,k), cr1(j,k), cw1(j,k), FPC_pls_2(j,k)
            ! if(dens1(j,k).lt.0.001) then !densidade mínima de indivíduos (from LPJmfire code)    
                                    !quase uma mortalidade por tamanho máximo
                cl2(j,k) = 0.

                cw2(j,k) =0.

                ch2(j,k) =0.

                cs2(j,k) =0.
    
                cr2(j,k) = 0.                             
                              
                diam(j,k) = 0.
            
                crown_area(j,k) = 0.
            
                lai(j,k) = 0.

                height(j,k) = 0.
           
                FPC_ind(j,k) = 0.
                
                FPC_pls_2(j,k) = 0.

                dens1(j,k) = 0.

                dens2(j,k) = 0.

                npp_inc2(j,k) = 0.
                
               
            else
               
                cl2(j,k) = (cl1(j,k)/dens1(j,k)) 

                ! cw2(j,k) = (cw1(j,k)/dens1(j,k))
                
                ch2(j,k) = (ch1(j,k)/dens1(j,k))

                cs2(j,k) = (cs1(j,k)/dens1(j,k))

                cw2(j,k) = ch2(j,k) + cs2(j,k)

                cr2(j,k) = (cr1(j,k)/dens1(j,k)) 

                npp_inc2(j,k) = (npp_inc(j,k)/dens1(j,k)) 
        
                  !----------------------------------------------------------------------------
                 !Structuring PLSs [diameter, crown area and leaf area index]
                diam(j,k) = (4*(cs2(j,k)+ch2(j,k)) / (dwood(j,k)*1000000.) / pi / k_allom2)**(1./(2. + k_allom3)) !Eqn 9
                !print*, 'diam lpjmfire', diam(j,k)*100
                

                crown_area(j,k) = min(crown_area_max,k_allom1 * diam(j,k)**krp)
                !crown_area(j,k)=k_allom1 * diam(j,k)**krp
                !print*,'ca lpj', crown_area(j,k)

                lai(j,k) = (cl2(j,k)*spec_leaf(j,k))/crown_area(j,k)
                !print*, 'lai', lai(j,k)
                
                height(j,k) = k_allom2 *diam(j,k)**k_allom3

                ! if(height(j,k).gt.30.) then
                
                !     print*, height(j,k), lai(j,k), crown_area(j,k), diam(j,k)*100
                ! endif
                !------------------------------------------------------------------------------
                !---------------------------------------------------------------------------
                !Calculatin Foliage Projective Cover of average individual(FPC_ind), of the PLS(FPC_pls)
                ! and of the grid cell (FPC_total)

                FPC_ind(j,k) = (1-(exp(-0.5*lai(j,k))))
                !print*, 'fpcind', FPC_ind(j,k)
                
                
            
                FPC_pls_2(j,k) = (crown_area(j,k) * dens1(j,k) * FPC_ind(j,k)) 
               
            endif
            
            
            
            FPC_total_2(k) = FPC_total_2(k) + (FPC_pls_2(j,k)) !accumulate the values in the variable FPC_total.
                                                        !the actual value will only be obtained when j = npls
            
            if (j.eq.npls) then   !take the value accumulated until the last pls
              
                FPC_total_accu_2(k) = FPC_total_2(k)

                ! if(height(j,k).le.0.)then
                !     print*, height(j,k), FPC_pls_2(j,k), FPC_inc(j,k), dens1(j,k), diam(j,k), cl2(j,k)
    
                ! endif
               
            endif
            ! if(k.eq.1)then
            !     print*, 'fpc_total_accu2', FPC_total_accu_2(k)
            ! endif
            
        enddo

        !--------------------------------------------------------------------------- 
        
        

        ! FPC_inc_grid(k) = FPC_total_accu_2(k) - FPC_total_accu_1(k)
        
        ! if(FPC_inc_grid(k).le.0.)FPC_inc_grid(k)=0!print*, 'inc le ', FPC_inc_grid(k), k,FPC_total_accu_2(k), FPC_total_accu_1(k) 

        dead_pls = 0.
        do j=1, npls

            
            if (FPC_pls_2(j,k).le.0..or.dens1(j,k).lt.1.e-1) then !REVER ESSA MORTALIDADE POR DENSIDADE
                dead_pls = dead_pls +1
                
            endif
                    
            
            
            if (j.eq.npls) then
                dead_pls = dead_pls
                ! print*,'d lt 0', dead_pls
                alive_pls = npls - dead_pls
               
                
                ! PRINT*,'dead2',dead_pls,'alive', alive_pls
            endif

            ! FPC_inc(j,k) = FPC_pls_2(j,k) - FPC_pls_1(j,k)
                      

            if(FPC_pls_2(j,k).le.0..or.dens1(j,k).lt.1.e-1)then
                FPC_inc(j,k) = 0.
                FPC_dec(j,k) = 0.                   
                nind_kill_FPC(j,k) = 0.               
                greff(j,k) = 0.
                mort(j,k) = 1.
                mort_greff(j,k) = 0.
                mort_wd(j,k) = 0.
                dens2(j,k) = 0.
                nind_kill_greff(j,k) = 0
                nind_kill_total(j,k) = 0
            endif
            ! if (FPC_inc(j,k).lt.0.)then
            !     print*, 'fpc negativo'
            ! endif
            
        enddo

        !Verifying if FPCs together occupy more than 95% of the plot area

        fpc_max_tree = gc_area*0.95 !utilizaremos 1 ha !! 5% é destinado ao novo estabelecimento
        !---------------------------------------------------------------------------------------
                
        if (FPC_total_accu_2(k) .gt. fpc_max_tree) then
            print*, 'ULTRAPASSSSSOOOOUUUUUUUUUUUUUUUUUUUU', FPC_total_accu_2(k),k, FPC_total_accu_2(k)-fpc_max_tree
            
            ! print*, FPC_inc_grid(k), 'inc grid ultrapassou', k
       
           
            est_pls(:,k) = 0.0 !if the total FPC (considering all PLS) is grater than fpc_max_tree there is no new establishment
           
            ! Excedent area
           
            exc_area(k) = FPC_total_accu_2(k) - fpc_max_tree !quanto foi o excedente em área

            !------------------------------------------------------------
            !LPJmlFire scheme
            do j = 1, npls
                if(FPC_pls_2(j,k).gt.0) then !se a ocupação é maior q zero = PLS vivo.
                    FPC_dec(j,k) = min(FPC_pls_2(j,k), exc_area(k)*(FPC_pls_2(j,k)/FPC_total_accu_2(k)))
                    !quanto q a ocupação de um pls vai ter que reduzir. Essa redução é proporcional ao FPC (ocupação) do PLS em relação ao FPC total (soma de todos os FPC da célula-grade)
                    ! print*, 'dec new', FPC_dec(j,k), FPC_pls_2(j,k), j
                else
                    FPC_dec(j,k) = 0
                    ! print*,'0'
                endif

                if(FPC_dec(j,k).le.0.) then
                    nind_kill_FPC(j,k) = 0.
                    ! print*, '000000000000'
                    
                else
                    nind_kill_FPC(j,k) = (dens1(j,k) * FPC_dec(j,k))/FPC_pls_2(j,k) !NIND_KILL.
                    ! print*, 'prop', FPC_dec(j,k)/FPC_pls_2(j,k)
                    !numero de ind. que vão morrer (ind/m2) devido ocupação maior que 95%
                    ! print*, nind_kill_FPC(j,k), dens1(j,k), 'FPCs', ((FPC_dec(j,k))/(FPC_pls_2(j,k)))
                    ! nind_kill_prop(j,k) = 1 - (dens1(j,k)-FPC_dec_prop(j,k))/dens1(j,k)
                    ! print*, 'new kill', FPC_dec_prop(j,k), dens1(j,k), FPC_dec(j,k),j,nind_kill_prop(j,k)
                endif
    
                ! print*, 'dec', FPC_dec_prop(j,k)
            enddo            

            do j = 1, npls
               
                if(FPC_pls_2(j,k).le.0..or.dens1(j,k).lt.1e-1)then
                    ! print*, 'MORTALITYY'
                    FPC_dec(j,k) = 0.                   
                    nind_kill_FPC(j,k) = 0.               
                    greff(j,k) = 0.
                    mort(j,k) = 1.
                    mort_greff(j,k) = 0.
                    mort_wd(j,k) = 0.
                    dens1(j,k) = 0.
                    ! nind_kill_prop(j,k) = 0.        
                else
                

                    greff(j,k) = carbon_increment(j,k)/(cl2(j,k)*spec_leaf(j,k)) !growth efficiency in m2/gC

                    mort_wd(j,k) = exp(-2.66+(0.255/dwood(j,k))) !
                    ! print*, 'mort dwood', mort_wd(j,k)
                    
                    ! mort_greff(j,k) = k_mort1/(1+(k_mort2*greff(j,k))) !mortality by gowth efficiency
                        ! print*, 'mort_greff', mort_greff(j), j

                    mort_greff(j,k) = mort_wd(j,k)/(1+k_mort2*greff(j,k))
                    ! print*, 'mort greff dwood', mort_greff(j,k)

                    
                    ! mort(j,k) = max((FPC_dec_prop(j,k) + mort_greff(j,k),1)) !sum of all mortality
                    ! print*, 'mort_total', mort(j,k)

                    nind_kill_greff(j,k) = (dens1(j,k) * mort_greff(j,k)) !valores absolutos de ind.
                    ! print*, 'NIND_KILL', nind_kill_greff(j,k)

                    !soma nind_kill
                    nind_kill_total(j,k) = nind_kill_FPC(j,k) + nind_kill_greff(j,k) !em ind/m2
                    ! print*, 'NIND_KILL TOTAL', nind_kill_total(j,k), 'KILL_FPC',&
                    ! nind_kill_FPC(j,k), 'KILL GREFF', nind_kill_greff(j,k), 'DENS',dens1(j,k)

                    !mort(j,k) = ((dens1(j,k)-nind_kill_total(j,k))/dens1(j,k)) !quanto vai morrer em relação a densidade atual
                    mort(j,k) = nind_kill_total(j,k)/dens1(j,k)
                    ! print*,'mort fpc', mort(j,k),nind_kill_total(j,k),dens1(j,k)
                    ! print*, 'quem ficou', mort(j,k)
                    ! mort(j,k) = 1 - mort(j,k)
                    ! print*, 'quem morreu de vdd', mort(j,k)
                    
                    cleaf_new(j,k) = cl2(j,k)
                   
                    cwood_new(j,k) = cw2(j,k)

                    cheart_new (j,k) = ch2(j,k)

                    csap_new (j,k) = cs2(j,k)

                    croot_new(j,k) = cr2(j,k)

                    dens_est(j,k) = dens1(j,k)


                endif    
              

                if (mort(j,k).lt.0.)then !maximum mortality in this case
                    
                    mort(j,k) = 1.
                                  
                endif        
               

            enddo

        else !total FPC of all PLS is smaller than fpc_max_tree
            
            gc_available(k) = fpc_max_tree - FPC_total_accu_2(k)
            ! print*, 'gc aavailable', gc_available

            ! print*, FPC_inc_grid(k), 'inc grid ', k, FPC_total_accu_2(k)
            !if the occupation is smaller than the stand area the mortality is defined only by
            !the growth efficiency and the loss of carbon through turnover
            count_pls = 0.
            ! if(FPC_total_accu_2(k).lt.5000.) then
            !     print*,FPC_total_accu_2(k), k
            ! endif 
            print*, 'n ultrapassou', FPC_total_accu_2(k), k
            do j=1, npls
                
                
                FPC_inc(j,k) = FPC_pls_2(j,k) - FPC_pls_1(j,k)

                ! if (height(j,k).le.0.) then
                !     print*, height(j,k), cl2(j,k), cw2(j,k), cr2(j,k), FPC_pls_2(j,k), diam(j,k), FPC_inc(j,k)
                ! endif
                ! print*, 'alive_pls', alive_pls
                call establishment(j,gc_available(k),alive_pls, FPC_total_accu_2(k),gc_area, est(k),est_pls(j,k),&
            &       FPC_pls_2(j,k))
                ! pint*,'establishment', FPC_total_accu_2(k), est(k),j,k, est_pls(j,k)
                call sapling_allometry(alive_pls,cleaf_sapl(j,k),csap_sapl(j,k),cheart_sapl(j,k),croot_sapl(j,k))
                
                call shrink(spec_leaf(j,k),dwood(j,k),cl2(j,k),ch2(j,k),cs2(j,k),cw2(j,k),cr2(j,k),est_pls(j,k),dens1(j,k),&
            &      cleaf_sapl(j,k),csap_sapl(j,k),cheart_sapl(j,k),croot_sapl(j,k),&
            &      dens_est(j,k),cleaf_new(j,k),cwood_new(j,k),cheart_new(j,k),&
            &      csap_new(j,k),croot_new(j,k))
                ! if(cleaf_new(j,k).gt.0.) then
                !     print*, 'after shrink',cleaf_new(j,k)/1000,cwood_new(j,k)/1000,croot_new(j,k)/1000, height(j,k)
                ! endif

                !PRINT*, 'as', cl2(j,k)/1000, cw2(j,k)/1000, cr2(j,k)/1000, dens_est(j,k),height(j,k)
                ! endif            
            
                ! cl2(j,k) = cleaf_new(j,k)
                ! cw2(j,k) = cwood_new(j,k)
                ! ch2(j,k) = cheart_new(j,k)
                ! cs2(j,k) = csap_new(j,k)
                ! cr2(j,k) = croot_new(j,k)
                ! dens1(j,k) = dens_est(j,k)
                
                if(FPC_pls_2(j,k).le.0..or.dens_est(j,k).lt.(1.e-1)) then
                    ! print*, 'LT 0  ', FPC_pls_2(j,k),j
                    cleaf_new(j,k) = 0.
                    cwood_new(j,k) = 0.
                    cheart_new(j,k) = 0.
                    csap_new(j,k) = 0.
                    croot_new(j,k) = 0.
                    dens_est(j,k) = 0.
                    count_pls = count_pls+1
                endif
                
                if(j.eq.npls)then
                    count_pls = count_pls
                    ! print*,'count_pls', count_pls
                endif
            
                if(cleaf_new(j,k).le.0..or.cl2(j,k).le.0.) then
                    
                    greff(j,k) = 0.
                    mort_greff(j,k) = 0.
                    mort(j,k) = 1.
                    
                else    
                    
                    greff(j,k) = carbon_increment(j,k)/(cleaf_new(j,k)*spec_leaf(j,k))

                    mort_wd(j,k) = exp(-2.66+(0.255/dwood(j,k))) !

                    mort_greff(j,k) = mort_wd(j,k)/(1+k_mort2*greff(j,k))

                    ! mort_greff(j,k) = k_mort1/(1+(k_mort2*greff(j,k)))

                    nind_kill_greff(j,k) = (dens_est(j,k) * mort_greff(j,k))

                    nind_kill_total(j,k) = nind_kill_greff(j,k)
                
                    mort(j,k) = nind_kill_total(j,k)/dens_est(j,k)
                    ! mort(j,k) = mort_greff(j,k)
                    ! mort(j,k) = (1-(dens1(j,k)-nind_kill_total(j,k))/dens1(j,k)) 
                   !mort(j,k) = 1 - mort(j,k) !adicionando mortalidade pra quando nãoultrapassa
                    ! print*, 'greff', greff(j), carbon_increment(j)/1000., cl2(j)/1000., spec_leaf(j)
                    !print*, 'mort_greff', mort_greff(j), j
                    ! print*, 'mort else', mort(j,k),nind_kill_total(j,k)
                    !fpc_dec(j) = 0.

                    ! print*, nind_kill_FPC(j,k), dens1(j,k), 'FPCs', ((FPC_dec(j,k))/(FPC_pls_2(j,k)))
                endif


            enddo
            
        endif
       
        
        
        
        do j=1,npls
            
            if(mort(j,k).ge.1.) then !maximum mortality is equal to 1
                mort(j,k) = 1.
                
            else
                mort(j,k) = mort(j,k)
            endif

            ! print*, 'mort',mort(j)   
            remaining(j,k) = 1. - mort(j,k)
           
            ! print*, 'remaining', remaining(j,k), 'mort', mort(j,k), j
           
            if (remaining(j,k) .le. 0.) then
                ! print*, 'PLS dead===============================================================',j
                ! goto 10 
                dens2(j,k) = 0.
                cleaf_new(j,k) = 0.
                cwood_new(j,k) = 0.
                cheart_new(j,k) = 0.
                csap_new(j,k) = 0.
                croot_new(j,k) = 0.
                FPC_pls_2(j,k) = 0.
                ! FPC_total_accu_2(k) = 0. 
                npp_inc(j,k) = 0
                npp_inc2(j,k) = 0
                annual_npp(j,k) = 0.
                leaf_inc(j) = 0.
                root_inc(j) = 0.
                wood_inc(j) = 0.
                cl1(j,k) = 0.
                cw1(j,k) = 0.
                ch1(j,k) = 0.
                cs1(j,k) = 0.
                cr1(j,k) = 0.
            endif

            ! print*, 'testing', est_pls(j,k)

            dens2(j,k) = (dens_est(j,k) * remaining(j,k))
            ! print*, 'dens_2 (pos remaining)', dens2(j,k),j 

            !dens2(j,k) = dens2(j,k) + est_pls(j,k)


            ! print*, 'dens_2 (pos estab)', dens2(j,k)
            ! print*, '                            '
            ! print*, '                            '
            ! print*, '                            '

            cleaf_new(j,k) = cleaf_new(j,k) * remaining(j,k)
            ! print*, 'cl', cleaf_new(j,k)/1000  
            
            cwood_new(j,k) = cwood_new(j,k) * remaining(j,k)

            cheart_new(j,k) = cheart_new(j,k) * remaining(j,k)

            csap_new(j,k) = csap_new(j,k) * remaining(j,k)

            croot_new(j,k) = croot_new(j,k) * remaining(j,k)

            


            !


            !!carbon to litter
            ! carbon_litter_leaf(j) = cl2(j) - ((cl2(j))*remaining(j))
            ! print*,  'litter', carbon_litter_leaf(j), cl2(j), remaining(j)
            !add litter of previous time step
            !the litter when pls dies is equal to cl, cw, cr
            !print*,'inside',FPC_total_accu_2(k)

            
        enddo


         !----------------------------------------------------------------------------
        !updating the variables for the next year

        cl1_aux(:,k) = cleaf_new(:,k)!cl2(:,k)
        
        ! print*, 'cl1_aux', cl1_aux(:,k)/1000.

        cw1_aux(:,k) = cwood_new(:,k)

        ! print*, 'cw1 atualizado', cw1/1000.

        cs1_aux(:,k) = csap_new(:,k)

        ch1_aux(:,k) = cheart_new(:,k)
        
        cr1_aux(:,k) = croot_new(:,k)

        ! print*, 'cr1 atualizado', cr1/1000.


        dens1_aux(:,k) = dens2(:,k)


        !FPCs

        FPC_pls_1_aux(:,k) = FPC_pls_2(:,k)

        ! print*, 'FPC_atualizado', FPC_pls_1


        FPC_total_accu_1_aux(k)= FPC_total_accu_2(k)
        ! print*, 'FPC_atualizado', FPC_total_accu_1_aux(k),FPC_total_accu_2(k),k

       
        !-----------------------------------------------------------------------------
        !!!---------------------------------------------------------------------------
        !!!Fictitious allocation process in order to test the logic developed
         !------------------------------------------------------------------------------
        !NPP increment (NPPt-NPPt-1); for testing purpose a general value was defined
        !Transforms NPP increment from m2 to NPP increment for each averge individual
        
        do j=1,npls

            !!     
            ! print*, 'cl2 before alloc', cl1_aux(j,k), cw1_aux(j,k), cr1_aux(j,k),&
            !     &dwood(j,k), spec_leaf(j,k), dens1_aux(j,k), npp_inc(j,k)
            !print*,'outside alloc', cl1_aux(j,k), cw1_aux(j,k),ch1_aux(j,k),cs1_aux(j,k)

            call allocation(cl1_aux(j,k), cw1_aux(j,k),ch1_aux(j,k),cs1_aux(j,k),cr1_aux(j,k),&
                &dwood(j,k), spec_leaf(j,k), dens1_aux(j,k), npp_inc2(j,k), height(j,k),&
                &cl_inc(j,k), cw_inc(j,k),ch_inc(j,k),cs_inc(j,k), cr_inc(j,k),ctotal_inc(j,k),&
                &cl2_aux(j,k), ch2_aux(j,k), cs2_aux(j,k), cr2_aux(j,k), cw2_aux(j,k))

                !print*,'after alloc', cl2_aux(j,k), ch2_aux(j,k), cs2_aux(j,k), cr2_aux(j,k), cw2_aux(j,k)

            if(dens1_aux(j,k).le.0.) then

                cl2_aux(j,k) = 0.

                cr2_aux(j,k) = 0. 
                
                cw2_aux(j,k) = 0.

                ch2_aux(j,k) = 0.
                
                cs2_aux(j,k) = 0.

                npp_inc(j,k) = 0.

                npp_inc2(j,k) = 0.
            
                annual_npp(j,k) = 0.

                leaf_inc(j) = 0.

                root_inc(j) = 0.

                wood_inc(j) = 0.

                !sapinc/heartinc

                cl1_aux(j,k) = 0.

                cw1_aux(j,k) = 0.

                ch1_aux(j,k) = 0.

                cs1_aux(j,k) = 0.

                cr1_aux(j,k) = 0.

                ctotal_inc(j,k) = 0
                ! delta_carbon_pls(j) = 0.

                ! print*, 'delta', delta_carbon_pls(j), j
            
            else
            
            
            
            !     npp_inc(j,k) = npp_inc(j,k)/dens1_aux(j,k)

           

            ! !-------------------------------------------------------------------------------
            ! !Annual NPP available to allocation (??????? é essa NPP ou a NPP inc?)
        
            !     annual_npp(j,k) = ((npp1_initial(j,k)/dens1_aux(j,k)) + npp_inc(j,k))

            ! print*, 'annual npp', annual_npp(j)/1000.

             !-------------------------------------------------------------------------------
             ! !Increments to each compartments per individual. Here, the NPP proportions allocated
             ! to each compartment is being used for testing purpose. The actual values will be calculated
             ! in allocation routine.

                leaf_inc(j) = leaf_allocation * npp_inc2(j,k)

                root_inc(j) = root_allocation * npp_inc2(j,k) 

                wood_inc(j) = wood_allocation * npp_inc2(j,k)  

                carbon_increment(j,k) = ctotal_inc(j,k)!leaf_inc(j) + root_inc(j) + wood_inc(j)
                if(carbon_increment(j,k).le.0)then
                    ! print*,'c inc neg', carbon_increment(j,k)
                endif

                ! print*, 'final', carbon_increment(j)/1000.

                

                cl1_aux(j,k) = cl2_aux(j,k) !leaf already with allocation !cl1_aux(j,k) + leaf_inc(j)
                
                cs1_aux(j,k) = cs2_aux(j,k)
                
                ch1_aux(j,k) = ch2_aux(j,k)

                cw1_aux(j,k) = cw2_aux(j,k) !cw1_aux(j,k) + wood_inc(j)
                
                cr1_aux(j,k) = cr2_aux(j,k) !cr1_aux(j,k) + root_inc(j)

                !saving value for avg individual for outputs
                cleaf_avg_ind(j,k) = cl1_aux(j,k)
                !csap avg e cheart avg
                cwood_avg_ind(j,k) = cw1_aux(j,k)
                croot_avg_ind(j,k) = cr1_aux(j,k)

                ! print*, 'cleaf avg ind', cleaf_avg_ind(j,k)/1000
                
            endif
           
            ! print*, 'cl1 com incremento após aloca', cl1_aux(j,k)/1000., j
            ! print*, 'cl avg ind', cleaf_avg_ind(j,k)/1000., j
            ! print*, 'cw avg ind', cwood_avg_ind(j,k)/1000.
            ! print*, 'cr avg ind', croot_avg_ind(j,k)/1000.



            ! print*, ''
            ! print*, 'cw1 com incremento após aloca', cw1(j)/1000., j
            ! print*, ''
            ! print*, 'cr1 com incremento após aloca', cr1(j)/1000., j

            !print*, 'densidade p/ ano seguinte =======', dens_1(j)
            ! Loss of carbon through residence time
           

            cl1_aux(j,k) = cl1_aux(j,k) - (cl1_aux(j,k)/res_time_leaf)
            ! print*, 'cl2 after restime', cl2(j)/1000., (cl1(j)/res_time)/1000.
            ! print*, ''

            cs1_aux(j,k) = cs1_aux(j,k) - (cs1_aux(j,k)/res_time_sap)

            ! ch1_aux(j,k) = (ch1_aux(j,k) + cs1_aux(j,k)) - ch1_aux(j,k)/res_time_wood

            ch1_aux(j,k) = ch1_aux(j,k)  + (cs1_aux(j,k))

            cw1_aux(j,k) = cw1_aux(j,k) - (cw1_aux(j,k)/res_time_wood)

            cr1_aux(j,k) = cr1_aux(j,k) - (cr1_aux(j,k)/res_time_root)


            cl1_aux(j,k) = cl1_aux(j,k) * dens1_aux(j,k)
            ! print*, 'CL1', (cl1_aux(j,k)/1000)
            ! if(cl1_aux(j,k).ne.0.) then
            !     print*, 'cl * dens', cl1_aux(j,k)/1000.,j, dens1_aux(j,k)
            ! endif
            cw1_aux(j,k) = cw1_aux(j,k) * dens1_aux(j,k)
            ! print*, 'CLw', cw1_aux(j,k)/1000
            ! print*, 'cw * dens', cw1_aux(j,k)/1000, dens1_aux(j,k)

            ch1_aux(j,k) = ch1_aux(j,k) * dens1_aux(j,k)
            ! print*, 'CH1', ch1_aux(j,k)

            cs1_aux(j,k) = cs1_aux(j,k) * dens1_aux(j,k)
            ! print*, 'CS1', cs1_aux(j,k)

            cr1_aux(j,k) = cr1_aux(j,k) * dens1_aux(j,k)
            ! print*, 'CR', cr1_aux(j,k)

            npp_inc2(j,k) = npp_inc2(j,k) 

            carbon_increment(j,k) = carbon_increment(j,k)

            ! print*, 'l', cl1_aux(j,k)/1000., 's',cs1_aux(j,k)/1000., 'r', cr1_aux(j,k)/1000., 'h', ch1_aux(j,k)/1000.,&
            ! 'dens',dens1_aux(j,k), 'height',height(j,k)


            ! delta_carbon_pls(j) = delta_carbon_pls(j)

            ! print*, '==============================='
            ! print*,'cl final', cl1(j)/1000., j

        enddo


    enddo

    

    open(unit=1,file='totalFPC.csv',status='unknown')
        do k=1, time
            

            
            write(1,*) FPC_total_accu_2(k)
            
        enddo 
    close(1)

    ! open(unit=1,file='carbon_pools_time_5PLS.csv',action='write',status='replace')
    ! do k=1, time
    !    do j = 1,npls

    !         write(1,*) cl1_aux(j,k)/1000.,',',cw1_aux(j,k)/1000.,',', cr1_aux(j,k)/1000.,',','pls',j,',',k !newline
    !     enddo
    ! enddo 
    ! close(1)

    open(unit=1,file='density.csv',action='write',status='replace')
    do k=1, time
       do j = 1,npls

            write(1,*) dens1_aux(j,k),',','pls',j,',',k !newline
        enddo
    enddo 
    close(1)

    open(unit=1,file='carbons.csv',action='write',status='replace')
    do k=1, time
       do j = 1,npls

            write(1,*) cl1_aux(j,k)/1000.,',',cw1_aux(j,k)/1000.,',', cr1_aux(j,k)/1000.,',',j,',',k !newline
        enddo
    enddo 
    close(1)

end program self_thinning

! open(unit=1,file='totalFPC.csv',status='unknown')
!     ! do k=1, time
!     !     FPC(k) = FPC_total_accu_2(k)
!     !     call csv_write(1,/FPC,.true.)

!     ! enddo 

! a = (/1,2,3/)
! b = (/4,5/)

! call csv_write(1,a,.true.)
! call csv_write(1,b,.true.)  
! close(1)

! open(unit=1,file='carbon_pools_time_5PLS.csv',action='write',status='replace')
!     do k=1, time
!        do j = 1,npls

!             write(1,*) cl1_aux(j,k)/1000.,',',cw1_aux(j,k)/1000.,',', cr1_aux(j,k)/1000.,',','pls',j,',',k !newline
!         enddo
!     enddo 
! close(1)
        
!         open(unit=1,file='carbon_pools_time_5PLS_avgind.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) cleaf_avg_ind(j,k)/1000.,',',cwood_avg_ind(j,k)/1000.,',', croot_avg_ind(j,k)/1000.,',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)
        
!         open(unit=1,file='density_5PLS.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) dens1_aux(j,k),',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)

!         open(unit=2,file='FPC_total_time_5PLS.csv',status='unknown')
!         do k=1, time
        

       
!             write(2,*) FPC_total_accu_2(k),',',k, ',' , gc_area !newline
        
!         enddo    

!         close(2)
!     endif 

!     if (npls.eq.20) then
!         open(unit=1,file='carbon_pools_time_20PLS.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

       
!                 write(1,*) cl1_aux(j,k)/1000.,',',cw1_aux(j,k)/1000.,',', cr1_aux(j,k)/1000.,',','pls',j,',',k !newline
!             enddo
!         enddo
!         close(1)

!         open(unit=1,file='carbon_pools_time_20PLS_avgind.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) cleaf_avg_ind(j,k)/1000.,',',cwood_avg_ind(j,k)/1000.,',', croot_avg_ind(j,k)/1000.,',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)

!         open(unit=1,file='density_20PLS.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) dens1_aux(j,k),',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)

!         open(unit=2,file='FPC_total_time_20PLS.csv',status='unknown')
!         do k=1, time
        

       
!             write(2,*) FPC_total_accu_2(k),',',k, ',' , gc_area !newline
        
!         enddo    

!         close(2)
!     endif
    
!     if (npls.eq.50) then
!         open(unit=1,file='carbon_pools_time_50PLS.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

       
!                 write(1,*) cl1_aux(j,k)/1000.,',',cw1_aux(j,k)/1000.,',', cr1_aux(j,k)/1000.,',','pls',j,',',k !newline
!             enddo
!         enddo
!         close(1)
        
!         open(unit=1,file='carbon_pools_time_50PLS_avgind.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) cleaf_avg_ind(j,k)/1000.,',',cwood_avg_ind(j,k)/1000.,',', croot_avg_ind(j,k)/1000.,',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)

!         open(unit=1,file='density_50PLS.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) dens1_aux(j,k),',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)

!         open(unit=2,file='FPC_total_time_50PLS.csv',status='unknown')
!         do k=1, time
        

       
!             write(2,*) FPC_total_accu_2(k),',',k, ',' , gc_area !newline
        
!         enddo    

!         close(2)


!     endif

!     if (npls.eq.100) then
!         open(unit=1,file='carbon_pools_time_100PLS.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

       
!                 write(1,*) cl1_aux(j,k)/1000.,',',cw1_aux(j,k)/1000.,',', cr1_aux(j,k)/1000.,',','pls',j,',',k !newline
!             enddo
!         enddo
!         close(1)

!         open(unit=1,file='carbon_pools_time_100PLS_avgind.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) cleaf_avg_ind(j,k)/1000.,',',cwood_avg_ind(j,k)/1000.,',', croot_avg_ind(j,k)/1000.,',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)

!         open(unit=1,file='density_100PLS.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) dens1_aux(j,k),',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)

!         open(unit=2,file='FPC_total_time_100PLS.csv',status='unknown')
!         do k=1, time
        

       
!             write(2,*) FPC_total_accu_2(k),',',k, ',' , gc_area !newline
        
!         enddo    

!         close(2)
!     endif

!     if (npls.eq.500) then
!         open(unit=1,file='carbon_pools_time_500PLS.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

       
!                 write(1,*) cl1_aux(j,k)/1000.,',',cw1_aux(j,k)/1000.,',', cr1_aux(j,k)/1000.,',','pls',j,',',k !newline
!             enddo
!         enddo
!         close(1)

!         open(unit=1,file='carbon_pools_time_500PLS_avgind.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) cleaf_avg_ind(j,k)/1000.,',',cwood_avg_ind(j,k)/1000.,',', croot_avg_ind(j,k)/1000.,',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)

!         open(unit=1,file='density_500PLS.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) dens1_aux(j,k),',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)

!         open(unit=2,file='FPC_total_time_500PLS.csv',status='unknown')
!         do k=1, time
        

       
!             write(2,*) FPC_total_accu_2(k),',',k, ',' , gc_area !newline
        
!         enddo    

!         close(2)
!     endif
    

!     if (npls.eq.1000) then
!         ! print*, 'savinggggg'
!         ! open(unit=1,file='carbon_pools_time_1000PLS.csv',status='unknown')
!         ! do k=1, time
!         !     do j = 1,npls

       
!         !         write(1,*) cl1_aux(j,k)/1000.,',',cw1_aux(j,k)/1000.,',', cr1_aux(j,k)/1000.,',','pls',j,',',k !newline
!         !     enddo
!         ! enddo
!         ! close(1)

!         ! open(unit=1,file='carbon_pools_time_1000PLS_avgind.csv',status='unknown')
!         ! do k=1, time
!         !     do j = 1,npls

            
!         !         write(1,*) cleaf_avg_ind(j,k)/1000.,',',cwood_avg_ind(j,k)/1000.,',', croot_avg_ind(j,k)/1000.,',','pls',j,',',k !newline
!         !     enddo
!         ! enddo 
!         ! close(1)
!         ! open(unit=1,file='density_1000PLS.csv',status='unknown')
!         ! do k=1, time
!         !     do j = 1,npls

            
!         !         write(1,*) dens1_aux(j,k),',','pls',j,',',k !newline
!         !     enddo
!         ! enddo 
!         ! close(1)

!         open(unit=1,file='FPC_total_time_1000PLS.csv',status='unknown')
!         do k=1, time
              
!             write(1,*) FPC_total_accu_1_aux(k),',',k, ',' , gc_area !newline
        
!         enddo    

!         close(1)
!     endif


!     if (npls.eq.3000) then
!         open(unit=1,file='carbon_pools_time_3000PLS.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

       
!                 write(1,*) cl1_aux(j,k)/1000.,',',cw1_aux(j,k)/1000.,',', cr1_aux(j,k)/1000.,',','pls',j,',',k !newline
!             enddo
!         enddo
!         close(1)

!         open(unit=1,file='carbon_pools_time_3000PLS_avgind.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) cleaf_avg_ind(j,k)/1000.,',',cwood_avg_ind(j,k)/1000.,',', croot_avg_ind(j,k)/1000.,',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)

!         open(unit=1,file='density_3000PLS.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) dens1_aux(j,k),',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)

!         open(unit=2,file='FPC_total_time_3000PLS.csv',status='unknown')
!         do k=1, time
        

       
!             write(2,*) FPC_total_accu_2(k),',',k, ',' , gc_area !newline
        
!         enddo    

!         close(2)

!         open(unit=2,file='exc_area_3000PLS.csv',status='unknown')
!         do k=1, time
        

       
!             write(2,*) exc_area(k),',',k, ',' 
        
!         enddo    

!         close(2)

!         open(unit=1,file='mortality_3000PLS.csv',status='unknown')
!         do k=1, time
!             do j = 1,npls

            
!                 write(1,*) mort(j,k),',',mort_greff(j,k),',',FPC_dec(j,k),',',remaining(j,k),',','pls',j,',',k !newline
!             enddo
!         enddo 
!         close(1)
!     endif

    


! end program self_thinning 
  
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
! !   do k = 1, 3 !Loop dos anos
! !         if (k .eq. 1) then !*usando o time step t1 - INITIAL ALLOMETRY*

! !             !Carbon on tissues (wood and leaf) per average-individual (this considers the individual density [dens])
! !             cw1 = (cw1/dens)*1000. !*1000 transforma de kgC para gC - carbono
            
! !             cl1 = (cl1/dens)*1000.  !*1000 transforma de kgC para gC - carbono
           
! !             cr1 = (cr1/dens)*1000. !*1000 transforma de kgC para gC - carbono

! !             !PLS structure [diam, crown area and leaf area index]
! !             diam = ((4*(cw1))/((dwood*1000000.)*3.14*36))**(1/(2+0.22)) !nessa equação dwood deve estar em *g/m3*
! !             crown_area = k_allom1*(diam**krp)
! !             lai = (cl1*spec_leaf)/crown_area 
            
! !             !print*, crown_area, lai
! !             !Net Primary Productive logic (may be per average-individual). NPP ano atual - NPP ano anterior
! !             ! #1 NPP increment to each average-individual. Cada individuo médio de cada PLS vai ter um incremento
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