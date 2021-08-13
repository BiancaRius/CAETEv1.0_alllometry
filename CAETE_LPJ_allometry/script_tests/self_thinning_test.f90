program self_thinning

    !Variables to use
    integer :: j,k
    integer, parameter:: npls = 20
    real, dimension(npls), allocatable :: LAI (:) !Leaf Area Index (m2/m2)
    real, dimension(npls), allocatable :: diam (:) !Tree diameter in m. (Smith et al., 2001 - Supplementary)
    real, dimension(npls), allocatable :: crown_area (:) !Tree crown area (m2) (Sitch et al., 2003)
    real, allocatable :: FPC_pls (:) !Foliage projective cover for each PLS (Sitch et al., 2003)
    real, allocatable :: cont_exc (:)
    real, allocatable :: FPC_ind_ocp (:)
    real, allocatable :: FPC_avg_ind (:) !Average idividual occupation for a PLS (Sitch et al., 2003)
    real, allocatable :: FPCgrid_perc (:) !Fractional projective cover in grid cell relative to grid cell area (Sitch et al., 2003)
    real, allocatable :: nind (:) !number of individuals per PLS (Smith, 2001, thesis)
    real, allocatable :: FPCgrid_updt (:) !Fractional projective cover in grid cell (Sitch et al., 2003) - updated to occupy 95%
    real, allocatable :: FPCgrid_perc_updt (:) !Fractional projective cover in grid cell relative to grid cell area (Sitch et al., 2003)-updated to occupy 95%
    real :: FPC_grid_total = 0.0 !Fractional projective cover in grid cell (Sitch et al., 2003)
    real :: cont_inc_perc_total = 0.0
    real :: sum_total_inc = 0.0
    real :: sum_FPCgrid_perc=0.0
    real :: sum_nind=0.0
    real :: sum_cont_inc = 0.0
    real :: gc_area = 300 !grid cell size - 300 m2 FOR TESTING PURPOSE (the real value will be 1ha or 10000 m2)
    real :: gc_area95
    real :: sum_FPCgrid_updt = 0.0 !!the new number of total PLS average-individuals after % reduction equals maximum to
                                 !! not to exceed 95% occupation.
    real :: sum_FPCgrid_perc_updt = 0.0 !the new percentage of occupation of all PLS after % reduction. 
    real :: exc_area
    real :: exc_area_perc

    !Parameters and constants
    real :: k_allom1 = 100. !allometric constant (Table 3; Sitch et al., 2003)
    real :: krp = 1.6 !allometric constant (Table 3; Sitch et al., 2003)
    real :: ltor = 0.77302587552347657

    !Variables to allocation prototype
    real, dimension(npls) :: npp1 !KgC/ano
    real, dimension(npls) :: npp_inc = 0.0 !incremento anual de C para cada PLS
    real, dimension(npls) :: annual_npp = 0.0!quantidade de NPP com os incrementos.
    real, dimension(npls) :: cl2
    real, dimension(npls) :: cw2


    !============== ALLOMETRY EQUATIONS ===============!

    ! Variables with generic values for testing the logic code

    real, dimension(npls) :: dwood !wood density (g/cm-3) *Fearnside, 1997 - aleatory choices
    real, dimension(npls) :: cw1 !KgC/m2 (Cheart + Csap)
    real, dimension(npls) :: cl1 !KgC/m2 
    real, dimension(npls) :: cr1 !KgC/m2
    real, dimension(npls) :: spec_leaf !m2/gC
    real, dimension(npls) :: leaf_inc !kgC/ ind
    real, dimension(npls) :: wood_inc !kgC/ ind
    real, dimension(npls) :: root_inc !kgC/ ind
    real, dimension(npls) :: total_inc !kgC/ ind
    real, dimension(npls) :: cont_inc!kgC/ ind
    real, dimension(npls) :: cont_inc_perc!kgC/ ind

    dwood=(/0.74,0.73,0.59,0.52,0.41,0.44,0.86,0.42,0.64,0.69,0.92,0.60,0.36,0.99,0.59,0.52,0.41,0.44,0.86,0.42/) !atenção para a unidade
    cl1=(/2.15,3.,1.18,1.6,1.5,1.8,0.3,2.,0.8,.84,0.25,1.,0.2,1.7,1.18,1.6,1.5,1.8,0.3,2./)
    spec_leaf=(/0.0153,0.0101,0.0107,0.0112,0.012,0.0141,0.0137,0.0115,0.0122,0.010,0.012,0.011,&
    &0.013,0.014,0.0112,0.012,0.0141,0.0137,0.0115,0.0122/)
    cw1=(/7.,12.,7.2,8.3,8.8,9.7,7.5,11.5,10.,8.6,7.3,10.3,6.8,9.9,5.3,9.2,15.,12.6,10.7,11.4/)
    cr1=(/0.63,0.8,0.9,1.5,1.3,0.9,0.4,1.0,0.56,0.87,0.33,0.97,0.31,0.55,0.2,0.8,0.4,0.66,0.23,1.5/)
    npp1 = (/0.5,0.8,1.5,1.2,1.9,1.3,1.7,0.8,0.6,2.0,0.7,1.1,1.9,1.85,1.96,1.77,1.33,1.54,1.62,0.55/)

    allocate (nind(1:npls))
    allocate (FPC_pls(1:npls))
    allocate (FPC_avg_ind(1:npls))
    allocate (FPC_ind_ocp(1:npls))
    allocate (cont_exc(1:npls))

    allocate (FPCgrid_perc(1:npls))
    allocate (FPCgrid_updt(1:npls))
    allocate (FPCgrid_perc_updt(1:npls))

    ! Increment of carbon on tissues per individual 
    
    do k = 1, 3
        ! Allometric Equations =================================================

        ! Grid-Cell Properties =================================================

        if (k .eq. 1) then 
            !INITIAL ALLOMETRY
            diam = ((4+(cw1))/((dwood)*3.14*40))**(1/(2+0.5))
            !print*, 'diam', diam
            
            crown_area = k_allom1*(diam**krp)
            !print*, 'crown', crown_area

            !LAI individual (Sitch et al., 2003) - Cleaf/nind
            LAI = ((cl1*1000)*spec_leaf)/crown_area !transfor SLA gC to kgC
            !print*, 'LAI', LAI
            
            annual_npp = npp1 + 0.35 !a cada ano NPP aumenta 0.35kgC/ano

            nind = diam**(-1.6) !número de individuos-médios de cada PLS

            npp_inc = annual_npp / nind !quantidade de npp pra ser alocado por individuo-médio

            !==================================================
            !INCREMENTS TO LEAF AND WOOD TISSUES PER INDIVIDUO

            leaf_inc = 0.35*npp_inc

            root_inc = 0.35*npp_inc

            wood_inc = 0.3*npp_inc

            !==================================================
            !CARBON TISSUES

            cl2 = cl1 + leaf_inc

            cw2 = cw1 + wood_inc
            !==================================================

        else
            diam = ((4+(cw2))/((dwood)*3.14*40))**(1/(2+0.5))
            !print*, 'diam', diam
            
            crown_area = k_allom1*(diam**krp)
            !print*, 'crown', crown_area

            !LAI individual (Sitch et al., 2003) - Cleaf/nind
            LAI = ((cl2*1000)*spec_leaf)/crown_area !transfor SLA gC to kgC

            annual_npp = annual_npp + 0.35

            nind = diam**(-1.6) !número de individuos-médios de cada PLS

            npp_inc = annual_npp / nind !quantidade de npp pra ser alocado por individuo-médio

            !==================================================
            !INCREMENTS TO LEAF AND WOOD TISSUES PER INDIVIDUO

            leaf_inc = 0.35*npp_inc
            root_inc = 0.35*npp_inc
            wood_inc = 0.3*npp_inc
            !==================================================
            !CARBON TISSUES

            cl2 = cl1 + leaf_inc

            cw2 = cw1 + wood_inc
            !==================================================
        endif
        ! print*, k, 'nind', nind
        ! print*, k, 'npp_inc', npp_inc

        total_inc = cleaf_inc + wood_inc + root_inc !somatoria de inc de todos os tecidos p/ um PLS
        sum_total_inc = sum(total_inc)              !somatória de inc p/ todos os PLS
        cont_inc = total_inc/sum_total_inc          !contribuição relativa de cada PLS p/ inc total
        sum_cont_inc = sum(cont_inc)
        ! print*, 'TOTAL INC=', total_inc, 'SUM DE INC', sum_total_inc, 'CONT INC', cont_inc, 'SUM CONT', sum_cont_inc


        ! cont_inc_perc = (cont_inc*100)/sum_total_inc
        ! cont_inc_perc_total = sum(cont_inc_perc)
        ! print*, 'cont_inc', cont_inc
        ! print*, 'cont_inc_perc', cont_inc_perc,'total', cont_inc_perc_total

        FPC_avg_ind = (1-exp(-0.5*LAI)) !FPC ind médio m2

        FPC_pls = crown_area*nind*FPC_avg_ind 

        FPC_ind_ocp = FPC_pls/nind                  !ocupação de cada ind em cada PLS

       ! print*, 'OCP IND', FPC_ind_ocp, 'FPC PLS', FPC_pls , 'NIND', nind, 'FPC', FPC_ind_ocp*nind
        
        FPC_grid_total = sum(FPC_pls) !FPC da célula

        ! print*, 'FPC ind medio=', FPC_avg_ind, 'FPC_grid_total', FPC_grid_total

        gc_area95 = gc_area*0.95
        print*,'gc_area95', gc_area95

        if (FPC_grid_total .gt. gc_area95) then
            exc_area = FPC_grid_total - gc_area95
            cont_exc = cont_inc*exc_area
            print*, 'teste'!'exc_area', exc_area !, 'CONT EXC',cont_exc, 'SUM CONT EXC', sum(cont_exc)        
        endif

        FPCgrid_perc = (FPC_grid_total*100)/gc_area
        ! print*, 'FPC-GRID-PERC', FPCgrid_perc(j), gc_area

    enddo

end program self_thinning