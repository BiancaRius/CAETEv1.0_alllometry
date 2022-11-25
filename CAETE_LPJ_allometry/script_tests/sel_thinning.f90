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
!A grande maioria dos cálculos aqui é feita para o indivíduo médio, portanto, a maioria das variáveis são
!divididas pela densidade (indivíduos/m²). Pools são usados em gC.
!Este código é baseado principalmente no modelo LPJ (Sitch et al 2003 e Smith et al 2001) e no código
!do LPJ-MLFire
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

    !variables for mean individual in PLS populations
    ! integer, parameter :: grassess = 0.1*npls !gramíneas não estão sendo incorporadas no momento
    real(r_8), dimension(npls,time) :: lai !Leaf Area Index (m2/m2)
    real(r_8), dimension(npls,time) :: diam !Tree diameter in m. (Smith et al., 2001 - Supplementary)
    real(r_8), dimension(npls,time) :: crown_area !Tree crown area (m2) (Sitch et al., 2003)
    real(r_8), dimension(npls,time) :: height !Height (m) (Sitch et al., 2003)
    real(r_8), dimension(npls,time) :: FPC_ind!Foliage projective cover (m²) for each average individual of a PLS (Stich et al., 2003)

    !variables for PLSs
         !variables with 'initial' in the name are used only for the 1st day
    ! real(r_8), dimension(npls) :: est_pls !establishment for a specific PLS
    real(r_8), dimension(npls,time) :: FPC_pls_1  !Total Foliage projective cover of a PLS (m²) (Stich et al., 2003)
    real(r_8), dimension(npls,time) :: FPC_pls_2  !Total Foliage projective cover of a PLS (m²) (Stich et al., 2003)
    real(r_8), dimension(npls,time) :: FPC_dec    !Quantidade de área que deve ser reduzida de um PLS em cada ano (m²)
    real(r_8), dimension(npls,time) :: nind_kill_FPC !Número de indivíduos que devem morrer para que o PLS cumpra sua redução (indivíduos)
    real(r_8), dimension(npls,time) :: nind_kill_greff !número de indivíduos a serem mortos por ineficiência de crescimento(indivíduos)
    real(r_8), dimension(npls,time) :: nind_kill_total !Número total de indivíduos que devem morrer (somando todas as fontes de mortalidade)(indivíduos)
    real(r_8), dimension(npls,time) :: mort  !% a ser descontada dos compartimentos de carbono por conta da mortalidade
    real(r_8), dimension(npls,time) :: mort_greff  ! mortality from growth efficiency takes into account mort by wd (Sakschewski et al 2015)
    real(r_8), dimension(npls,time) :: mort_wd  ! fator de mortalidade que entra na mort_greff considerando o wood density (Sakschewski et al 2015)
    real(r_8), dimension(npls,time) :: greff  ! cálculo de eficiência de crescimento (Sitch et al 2003)
    real(r_8), dimension(npls,time) :: remaining  !1 - mort (% multiplicada pelos estados dos PLSs para saber quanto do que existe deve permanecer)  real(r_8), dimension(npls,time) :: mort_wd
    real(r_8), dimension(npls,time) :: FPC_inc !incremento em área (m²) para um PLS entre um ano e outro
    real(r_8), dimension(npls,time) :: carbon_increment_initial  ! used to calculate mort greff (Sitch et al 2003)
    real(r_8), dimension(npls,time) :: carbon_increment ! (gC)used to calculate mort greff (Sitch et al 2003)
    real(r_8), dimension(npls,time) :: totalc_1
    real(r_8), dimension(npls,time) :: totalc_2
    ! real(r_8), dimension (npls) :: FPC_inc_grid
        !Variables to allocation
    real(r_8), dimension(npls,time) :: npp_inc_init  !incremento inicial anual de C para cada PLS
    real(r_8), dimension(npls,time) :: npp_inc, npp_inc2  !incremento anual de C para cada PLS
    real(r_8), dimension(npls,time) :: annual_npp !quantidade de NPP com os incrementos.
    real(r_8), dimension(npls,time) :: cl1 !gC/m2 carbon on leaves
    real(r_8), dimension(npls,time) :: cw1 !gC/m2 (Cheart + Csap)
    real(r_8), dimension(npls,time) :: cs1 !gC/m2 carbon on sapwood
    real(r_8), dimension(npls,time) :: ch1 !gC/m2 carbon on heartwood
    real(r_8), dimension(npls,time) :: cr1 !gC/m2 carbon on fine roots
    real(r_8), dimension(npls,time) :: dens1 !density of individuals (ind/m²) previous mortality and establishment
    real(r_8), dimension(npls,time) :: cl2 !carbon on leaves after allocation (gC/ind)
    real(r_8), dimension(npls,time) :: cw2 !carbon on wood after allocation (gC/ind)
    real(r_8), dimension(npls,time) :: cs2 !carbon on sapwood after allocation (gC/ind)
    real(r_8), dimension(npls,time) :: ch2 !carbon on heartwood after allocation (gC/ind)
    real(r_8), dimension(npls,time) :: cr2 !carbon on wood after allocation (gC/ind)
    real(r_8), dimension(npls,time) :: dens2 !density of individuals (ind/m²) after mortality and establishment
     !variables with initial values
    real(r_8), dimension(npls,time) :: cl1_initial
    real(r_8), dimension(npls,time) :: cw1_initial
    real(r_8), dimension(npls,time) :: cs1_initial
    real(r_8), dimension(npls,time) :: ch1_initial
    real(r_8), dimension(npls,time) :: cr1_initial
    real(r_8), dimension(npls,time) :: npp1_initial
    real(r_8), dimension(npls,time) :: FPC_pls_initial
    real(r_8), dimension(npls,time) :: dens1_initial

    !variables for grid cell scale
        !variables with 'initial' in the name are used only for the 1st day
    real(r_8), dimension(time) :: FPC_total_initial = 0.0 !sum of FPC_pls for that grid cell
    real(r_8), dimension(time) :: FPC_total_accu_initial = 0.0 !variable to accumulate the FPC_total  
    real(r_8), dimension(time) :: FPC_total_2 = 0.0 !!sum of the FPC_pls of all PLSs for that grid cell
    real(r_8), dimension(time):: FPC_total_accu_1 = 0.0 !variable to accumulate the FPC_total 
    real(r_8), dimension(time) :: FPC_total_accu_2 = 0.0 !variable to accumulate the FPC_total (used for updating)
    real(r_8) :: dead_pls, alive_pls, count_pls !variables to account for the alive PLSs
    real(r_8) :: fpc_max_tree !95% of grid-cell (in m2)
    real(r_8), dimension(time) :: gc_available !area available for individuals establishment (m²)  
    real(r_8), dimension(time) :: exc_area  !excedent in area when the sum of FPCs is higher than fpc_max_tree
    

    !Parameters and constants
    real(r_8) :: gc_area = 10000 !patch size (1ha)
    real(r_8) :: k_allom1 = 100. !allometric constant (Table 3; Sitch et al., 2003)
    real(r_8) :: k_allom2 = 40.0
    real(r_8) :: k_allom3 = 0.85
    real(r_8) :: krp = 1.6 !allometric constant (Table 3; Sitch et al., 2003)
    real(r_8) :: ltor = 0.77302587552347657 !!!ATENÇÃO ver equação 2 em sitch et al 2003 (fator de stress hídrico)
    real(r_8) :: k_est = 0.06 !establishment constant !Smith et al 2001 - Table A1
    real(r_8) :: leaf_allocation = 0.4 !% of NPP allocaed to leaves
    real(r_8) :: wood_allocation = 0.3  !% of NPP allocated to wood
    real(r_8) :: root_allocation = 0.3 !% of NPP allocated to roots
    real(r_8) :: k_mort1 = 0.01 !mortality parameter from Sitch et al 2003
    real(r_8) :: k_mort2 = 0.5
    real(r_8) :: res_time_leaf = 2 !general residence time value
    real(r_8) :: res_time_root = 2
    real(r_8) :: res_time_sap = 2
    real(r_8) :: res_time_wood = 40 !ATENÇÃO! ESSE NUMERO PRECISA SER REVISADO POIS EM SITCH ET AL 2003 APENAS O SAPWOOD É PERDIDO POR TURNOVER
    real(r_8) :: crown_area_max = 30 !m2 !number from lplmfire code (establishment.f90)
    !!!ATENÇÃO PARA CROWN AREA MAX: está limitando a quantidade folhas, porém não restringe a quantidade de tecidos lenhosos. que está ficando gigante
    real(r_8) :: pi = 3.1415

   

   

    real(r_8), dimension(npls,time) :: cleaf_avg_ind
    real(r_8), dimension(npls,time) :: csap_avg_ind
    real(r_8), dimension(npls,time) :: cheart_avg_ind
    real(r_8), dimension(npls,time) :: cwood_avg_ind
    real(r_8), dimension(npls,time) :: croot_avg_ind

    

    ! Variables with generic values for testing the logic code
    real(r_8), dimension(npls,time) :: dwood !wood density (g/cm-3) *Fearnside, 1997 - aleatory choices
    real(r_8), dimension(npls,time) :: spec_leaf !m2/gC
    real(r_8), dimension(npls) :: leaf_inc !kgC/ ind
    real(r_8), dimension(npls) :: wood_inc !kgC/ ind
    real(r_8), dimension(npls) :: root_inc !kgC/ ind

   
    

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
    real(r_8), dimension (npls,time) :: cl3_aux = 0. !leaf carbon after allocation (gC, avera_in)


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
    xmax = 40.
     
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

            FPC_pls_initial(j,k) = 0.5

            FPC_total_initial(k) = FPC_total_initial(k) + FPC_pls_initial(j,k)

        
            if (j.eq.npls) then
                FPC_total_accu_initial(k) = FPC_total_initial(k)
                ! print*, 'FPC_total_accu_initial', FPC_total_accu_initial(k),k
            endif

        enddo

    enddo

    dead_pls = 0.
    alive_pls = 0. 
    count_pls = 0.

    do k = 1, time
    !initialize variables
        
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

        totalc_1(:,k) = 0.
        totalc_2(:,k) = 0.
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
      
        if(k.eq.1)then !use initial values
            cl1(:,k) = cl1_initial(:,k)
            cw1(:,k) = cw1_initial(:,k)
            cs1(:,k) = cs1_initial(:,k)
            ch1(:,k) = ch1_initial(:,k)
            cr1(:,k) = cr1_initial(:,k)
            FPC_pls_1(:,k) = FPC_pls_initial(:,k)
            dens1(:,k) = dens1_initial(:,k)
            FPC_total_accu_1(k) = FPC_total_accu_initial(k)
            npp_inc(:,k) = npp_inc_init(:,k)
            carbon_increment(:,k) = carbon_increment_initial(:,k)
            
        else !values from previous step is being used

            cl1(:,k) = cl1_aux(:,k-1)
            cw1(:,k) = cw1_aux(:,k-1)
            cs1(:,k) = cs1_aux(:,k-1) !- provisorio
            ch1(:,k) = ch1_aux(:,k-1) !- provisorio
            cr1(:,k) = cr1_aux(:,k-1)
            dens1(:,k) = dens1_aux(:,k-1)
            FPC_pls_1(:,k) = FPC_pls_1_aux(:, k-1)
            FPC_total_accu_1(k) = FPC_total_accu_1_aux(k-1)
            carbon_increment(:,k) = carbon_increment(:,k-1)
            
            
        !cria um valor de inc de npp provisório, uma vez que não temos uma subrotina
        !de produtividade offline
        
            npp_inc(:,k)=0. !reinitializing for a new sampling
        
            xmin = 0.1
            xmax = 3.5
     
            x(:,:) = 0.
            call random_number(x)

           
            do j = 1, npls
                x(j,k) = xmin + (xmax-xmin)*x(j,k)
                
            enddo
            npp_inc(:,k) = x(:,k)*1000
            
        endif !if k=1

           

        do j = 1, npls
             
                !--------------------------------------------------------------------------
        !transforming the carbon content from gC/m2 to gc/average individual 
        !(the carbon divided by dens gives the individual carbon, as in LPJ)

            !kill PLSs with very low density or leaf =0
            if(dens1(j,k).le.1.e-1.or.cl1(j,k).le.0.) then !densidade mínima de indivíduos (no código do LPJ fire esse valor é 1e-10)
                                                           !em teoria, a dens min é quase uma mort por tamanho maximo, uma vez que quanto
                                                           !menor a dens maiores são os individuos médios
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
                
               
            else ! /densidade transforma a variável de gc/m2 em gC
               
                cl2(j,k) = (cl1(j,k)/dens1(j,k)) 
                
                ch2(j,k) = (ch1(j,k)/dens1(j,k))

                cs2(j,k) = (cs1(j,k)/dens1(j,k))

                cw2(j,k) = ch2(j,k) + cs2(j,k)

                cr2(j,k) = (cr1(j,k)/dens1(j,k)) 

                npp_inc2(j,k) = (npp_inc(j,k)/dens1(j,k))

                totalc_1(j,k) = cl2(j,k)+ cs2(j,k)+ ch2(j,k)+ cr2(j,k) !used to calculate increment
    
            !----------------------------------------------------------------------------
            !Structuring PLSs [diameter, crown area and leaf area index]
                diam(j,k) = (4*(cs2(j,k)+ch2(j,k)) / (dwood(j,k)*1000000.) / pi / k_allom2)**(1./(2. + k_allom3)) !Eqn 9                

                crown_area(j,k) = min(crown_area_max,k_allom1 * diam(j,k)**krp) !!atenção para crown_area_max (ver comentário na declaração)

                lai(j,k) = (cl2(j,k)*spec_leaf(j,k))/crown_area(j,k)
                
                height(j,k) = k_allom2 *diam(j,k)**k_allom3
                
            !------------------------------------------------------------------------------
            !---------------------------------------------------------------------------
            !Calculatin Foliage Projective Cover of average individual(FPC_ind) and Fractional
            !Projective cover for PLS (FPC_pls2)

                FPC_ind(j,k) = (1-(exp(-0.5*lai(j,k))))
            
                FPC_pls_2(j,k) = (crown_area(j,k) * dens1(j,k) * FPC_ind(j,k)) 
               
            endif
            
            
            FPC_total_2(k) = FPC_total_2(k) + (FPC_pls_2(j,k)) !accumulate the values in the variable FPC_total.
                                                        !the actual value will only be obtained when j = npls
            
            if (j.eq.npls) then   !take the value accumulated until the last pls
              
                FPC_total_accu_2(k) = FPC_total_2(k)

            endif
            
        enddo

        !--------------------------------------------------------------------------- 
        
        dead_pls = 0.
        do j=1, npls

            
            if (FPC_pls_2(j,k).le.0..or.dens1(j,k).lt.1.e-1) then !REVER ESSA MORTALIDADE POR DENSIDADE
                dead_pls = dead_pls +1 
            endif

            if (j.eq.npls) then
                dead_pls = dead_pls
                alive_pls = npls - dead_pls
            endif


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
            
        enddo

        !Verifying if FPCs together occupy more than 95% of the plot area

        fpc_max_tree = gc_area*0.95 !utilizaremos 1 ha !! 5% é destinado ao novo estabelecimento
        !---------------------------------------------------------------------------------------
                
        if (FPC_total_accu_2(k) .gt. fpc_max_tree) then
            print*, 'ULTRAPASSOU', FPC_total_accu_2(k),k, FPC_total_accu_2(k)-fpc_max_tree
                       
            est_pls(:,k) = 0.0 !if the total FPC (considering all PLS) is grater than fpc_max_tree there is no new establishment
           
            ! Excedent area
           
            exc_area(k) = FPC_total_accu_2(k) - fpc_max_tree !quanto foi o excedente em área

            !------------------------------------------------------------
            !LPJmlFire scheme ! previously it was calculated from LPJ scheme, that considers the percentage of increment to FPC of a PLS
            do j = 1, npls
                if(FPC_pls_2(j,k).gt.0) then !se a ocupação é maior q zero = PLS vivo.
                    FPC_dec(j,k) = min(FPC_pls_2(j,k), exc_area(k)*(FPC_pls_2(j,k)/FPC_total_accu_2(k)))
                else
                    FPC_dec(j,k) = 0
                endif

                if(FPC_dec(j,k).le.0.) then
                    nind_kill_FPC(j,k) = 0. !não há mortalidade nos indivíduos
                    
                else
                    nind_kill_FPC(j,k) = (dens1(j,k) * FPC_dec(j,k))/FPC_pls_2(j,k) !NIND_KILL.
                    !numero de ind. que vão morrer (ind/m2) devido ocupação maior que 95%
                endif
    
            enddo            

            do j = 1, npls
               
                if(FPC_pls_2(j,k).le.0..or.dens1(j,k).lt.1.e-1)then
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

                    mort_wd(j,k) = exp(-2.66+(0.255/dwood(j,k))) !mort by wd (from Sakschewski et al 2015)

                    mort_greff(j,k) = mort_wd(j,k)/(1+k_mort2*greff(j,k))
                    ! print*, 'mort greff dwood', mort_greff(j,k)

                    nind_kill_greff(j,k) = (dens1(j,k) * mort_greff(j,k)) !valores absolutos de ind.
                    ! print*, 'NIND_KILL', nind_kill_greff(j,k)

                    !soma nind_kill
                    nind_kill_total(j,k) = nind_kill_FPC(j,k) + nind_kill_greff(j,k) !em ind/m2

                    !mort(j,k) = ((dens1(j,k)-nind_kill_total(j,k))/dens1(j,k)) !quanto vai morrer em relação a densidade atual
                    mort(j,k) = nind_kill_total(j,k)/dens1(j,k) !% de ind que devem morrer em relação ao que tinha anteriormente

                !update variables name
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

            print*, 'n ultrapassou', FPC_total_accu_2(k), k

            !if the occupation is smaller than the stand area, the mortality is defined only by
            !the growth efficiency and the loss of carbon through turnover
            
            count_pls = 0.
           
           
            do j=1, npls
                
                
                FPC_inc(j,k) = FPC_pls_2(j,k) - FPC_pls_1(j,k)

                call establishment(j,gc_available(k),alive_pls, FPC_total_accu_2(k),gc_area, est(k),est_pls(j,k),&
            &       FPC_pls_2(j,k))

                call sapling_allometry(alive_pls,cleaf_sapl(j,k),csap_sapl(j,k),cheart_sapl(j,k),croot_sapl(j,k))
                
                call shrink(spec_leaf(j,k),dwood(j,k),cl2(j,k),ch2(j,k),cs2(j,k),cw2(j,k),cr2(j,k),est_pls(j,k),dens1(j,k),&
            &      cleaf_sapl(j,k),csap_sapl(j,k),cheart_sapl(j,k),croot_sapl(j,k),&
            &      dens_est(j,k),cleaf_new(j,k),cwood_new(j,k),cheart_new(j,k),&
            &      csap_new(j,k),croot_new(j,k))
                
                if(FPC_pls_2(j,k).le.0..or.dens_est(j,k).lt.(1.e-1)) then
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
                    mort_greff(j,k) = 1.
                    mort(j,k) = 1.
                    
                else    
                    
                    greff(j,k) = carbon_increment(j,k)/(cleaf_new(j,k)*spec_leaf(j,k))

                    mort_wd(j,k) = exp(-2.66+(0.255/dwood(j,k))) !

                    mort_greff(j,k) = mort_wd(j,k)/(1+k_mort2*greff(j,k))


                    nind_kill_greff(j,k) = (dens_est(j,k) * mort_greff(j,k))

                    nind_kill_total(j,k) = nind_kill_greff(j,k)
                
                    mort(j,k) = nind_kill_total(j,k)/dens_est(j,k)
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
            remaining(j,k) = 1. - mort(j,k) !o que deve sobrar de cada pool em %
           
           
            if (remaining(j,k) .le. 0.) then
                ! print*, 'PLS dead===============================================================',j
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

        !desconta os indivíduos mortos
            dens2(j,k) = (dens_est(j,k) * remaining(j,k))

            cleaf_new(j,k) = cleaf_new(j,k) * remaining(j,k)
            
            cwood_new(j,k) = cwood_new(j,k) * remaining(j,k)

            cheart_new(j,k) = cheart_new(j,k) * remaining(j,k)

            csap_new(j,k) = csap_new(j,k) * remaining(j,k)

            croot_new(j,k) = croot_new(j,k) * remaining(j,k)

            npp_inc2(j,k) = npp_inc2(j,k) * remaining(j,k)

            
        enddo


         !----------------------------------------------------------------------------
       
        !!!---------------------------------------------------------------------------
         !------------------------------------------------------------------------------
       
        do j=1,npls


            call allocation(cleaf_new(j,k), cwood_new(j,k),cheart_new(j,k),csap_new(j,k),croot_new(j,k),&
                &dwood(j,k), spec_leaf(j,k), dens2(j,k), npp_inc2(j,k), height(j,k),&
                &cl_inc(j,k), cw_inc(j,k),ch_inc(j,k),cs_inc(j,k), cr_inc(j,k),ctotal_inc(j,k),&
                &cl2_aux(j,k), ch2_aux(j,k), cs2_aux(j,k), cr2_aux(j,k), cw2_aux(j,k))


            if(dens2(j,k).lt.1.e-1) then

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

                carbon_increment(j,k) = 0
            
            else

                !saving value for avg individual for outputs
                cleaf_avg_ind(j,k) = cl1_aux(j,k)
                !csap avg e cheart avg
                cwood_avg_ind(j,k) = cw1_aux(j,k)
                croot_avg_ind(j,k) = cr1_aux(j,k)

                ! print*, 'cleaf avg ind', cleaf_avg_ind(j,k)/1000
                
            endif
           
            ! Loss of carbon through residence time
           
            cl2_aux(j,k) = cl2_aux(j,k) - (cl2_aux(j,k)/res_time_leaf)

            cs2_aux(j,k) = cs2_aux(j,k) - (cs2_aux(j,k)/res_time_sap)

            !the dead sapwood is added to heartwood
            ch2_aux(j,k) = (ch2_aux(j,k) - (ch2_aux(j,k)/res_time_wood)) + cs2_aux(j,k)

            cw2_aux(j,k) = cs2_aux(j,k) + ch2_aux(j,k)

            cr2_aux(j,k) = cr2_aux(j,k) - (cr2_aux(j,k)/res_time_root)

            totalc_2(j,k) = cl2_aux(j,k)+cs2_aux(j,k)+ch2_aux(j,k)+cr2_aux(j,k)

        !update to gc/m²
            
            cl2_aux(j,k) = cl2_aux(j,k) * dens2(j,k)

            cw2_aux(j,k) = cw2_aux(j,k) * dens2(j,k)

            ch2_aux(j,k) = ch2_aux(j,k) * dens2(j,k)

            cs2_aux(j,k) = cs2_aux(j,k) * dens2(j,k)

            cr2_aux(j,k) = cr2_aux(j,k) * dens2(j,k)

            npp_inc2(j,k) = npp_inc2(j,k) 

            carbon_increment(j,k) = totalc_2(j,k) - totalc_1(j,k)

 !      !updating the variables for the next year

        cl1_aux(:,k) = cl2_aux(:,k)
        
        cw1_aux(:,k) = cw2_aux(:,k)

        cs1_aux(:,k) = cs2_aux(:,k)

        ch1_aux(:,k) = ch2_aux(:,k)
        
        cr1_aux(:,k) = cr2_aux(:,k)

        dens1_aux(:,k) = dens2(:,k)

        !FPCs

        FPC_pls_1_aux(:,k) = FPC_pls_2(:,k)

        ! print*, 'FPC_atualizado', FPC_pls_1


        FPC_total_accu_1_aux(k)= FPC_total_accu_2(k)
        ! ! print*, 'FPC_atualizado', FPC_total_accu_1_aux(k),FPC_total_accu_2(k),k

       
        !-----------------------------------------------------------------------------

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

