program self_thinning

    ! ================= VARIABLES TO USE DECLARATION ===================== !
    integer :: j,k
    integer, parameter:: npls = 20
    real, dimension(npls), allocatable :: lai (:) !Leaf Area Index (m2/m2)
    real, dimension(npls), allocatable :: diam (:) !Tree diameter in m. (Smith et al., 2001 - Supplementary)
    real, dimension(npls), allocatable :: crown_area (:) !Tree crown area (m2) (Sitch et al., 2003)
    real, dimension(npls) :: dens = 2  !densidade de indivíduos em cada PLS (0.5 é um número genérico)
    real, allocatable :: FPC_ind (:) !Foliage projective cover for each PLS (Stich et al., 2003)
    real, allocatable :: FPC_grid_t1 (:) !Fractional projective cover for each PLS in time 1 (Sitch et al., 2003)
    real, allocatable :: FPC_grid_t2 (:) !Fractional projective cover for each PLS in time 2 (Sitch et al., 2003)
    real, allocatable :: fpc_dec (:) !'decrease FPC' in other words: the excedend of FPC in eac year [LPJ-GUESS - Phillip's video]
    real, allocatable :: mort (:) !equivalente ao 'mort_shade' no LPJ-GUESS [Phillip's video]
    real, allocatable :: remaining (:) !taxa de redução
    real :: FPC_total_t1 = 0.0 !sum of FPC_grid in time step 1
    real :: FPC_total_t2 = 0.0 !sum of FPC_grid in time step 2
    real :: gc_area = 15 !grid cell size - 15 m2 FOR TESTING PURPOSE (the real value will be 1ha or 10000 m2)
    real :: fpc_max_tree

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
    real, dimension(npls) :: diameter

    !ARRAYS WITH GENERIC VALUES OF SOME VARIABLES
    dwood=(/0.24,0.53,0.39,0.32,0.31,0.44,0.66,0.42,0.74,0.39,0.82,0.40,0.26,0.79,0.39,0.52,0.41,0.44,0.86,0.42/) !atenção para a unidade
    cl1=(/.7,1.,0.3,1.6,1.1,1.8,0.3,0.2,0.8,.84,0.25,1.,0.2,1.7,0.4,.6,.5,.8,0.3,1.8/)
    spec_leaf=(/0.0153,0.0101,0.0107,0.0112,0.012,0.0141,0.0137,0.0115,0.0122,0.010,0.012,0.011,&
    &0.013,0.014,0.0112,0.012,0.0141,0.0137,0.0115,0.0122/)
    cw1=(/30.,22.,34.,28.3,20.2,19.7,27.5,19.5,20.,28.6,24.3,19.3,26.8,22.,18.3,22.,15.,22.6,10.7,21.4/)
    cr1=(/0.63,0.8,0.9,1.5,1.3,0.9,0.4,1.0,0.56,0.87,0.33,0.97,0.31,0.55,0.2,0.8,0.4,0.66,0.23,1.5/)
    npp1 = (/0.5,0.8,1.5,1.2,1.9,1.3,1.7,0.8,0.6,2.0,0.7,1.1,1.9,1.85,1.96,1.77,1.33,1.54,1.62,0.55/)
    diameter = (/0.16,0.45,0.17,0.25,0.34,0.4,0.23,0.49,0.37,0.5,0.53,0.12,0.75,0.22,0.63,0.31,0.41,0.63,0.52,0.15/)

    allocate (FPC_ind(1:npls))
    allocate (FPC_grid_t1(1:npls))
    allocate (FPC_grid_t2(1:npls))
    allocate (fpc_dec(1:npls))

    ! ================= END VARIABLES DECLARATION ===================== !

    ! ==================== ALLOMETRY EQUATIONS =========================!
    !        Increment of carbon on tissues per individual 

    do k = 1, 2 !Loop dos anos
        if (k .eq. 1) then !*usando o time step t1 - INITIAL ALLOMETRY*

            !Carbon on tissues (wood and leaf) per average-individual (this considers the individual density [dens])
            cw1 = (cw1/dens)*1000. !*1000 transforma de kgC para gC
            cl1 = (cl1/dens)*1000

            !PLS structure [diam, crown area and leaf area index]
            diam = ((4*(cw1))/((dwood*1000000.)*3.14*36))**(1/(2+0.22)) !nessa equação dwood deve estar em *g/m3*
            crown_area = k_allom1*(diam**krp)
            lai = (cl1*spec_leaf)/crown_area 
            
            !Net Primary Productive logic (may be per average-individual)
            ! #1 NPP increment to each average-individual
            npp_inc = 0.15/dens !this increment is fixed in 0.15kgC/yr

            ! #2 Annual NPP disponible to alloc 
            annual_npp = ((npp1/dens) + npp_inc)*1000 !*1000 transforma de kgC para gC

            !==================================================
            !INCREMENTS TO LEAF AND WOOD TISSUES PER INDIVIDUO

            leaf_inc = 0.35*npp_inc ![35% of NPP]
            root_inc = 0.35*npp_inc ![35% of NPP - Functional balance]
            wood_inc = 0.3*npp_inc  ![30% of NPP]

            !==================================================
            !CARBON TISSUES

            cl2 = cl1 + leaf_inc !cl1 e leaf_inc já está dividido peela densidade
            cw2 = cw1 + wood_inc !cw1 e wood_inc já está dividido peela densidade

            !==================================================
            !Foliage Projective Cover (FPC_ind) & Fractional Projective Cover (FPC_grid)
            
            FPC_ind = (1-exp(-0.5*lai)) !FPC ind médio [m2]
            FPC_grid_t1 = crown_area*dens*FPC_ind !FPC of pls [occupation on grid cell considers all average-individual of PLS; m2]
            FPC_total_t1 = sum(FPC_grid_t1)

        else !*usando time step 2*

            !PLS structure [diam, crown area and leaf area index]
            diam = ((4*(cw2))/((dwood*1000000.)*3.14*36))**(1/(2+0.22)) !nessa equação dwood deve estar em *g/m3*
            crown_area = k_allom1*(diam**krp)
            lai = (cl2*spec_leaf)/crown_area 

            !==================================================
            !Foliage Projective Cover (FPC_ind) & Fractional Projective Cover (FPC_grid)
            
            FPC_ind = (1-exp(-0.5*lai)) !FPC ind médio [m2]
            FPC_grid_t2 = crown_area*dens*FPC_ind !FPC of PLS [occupation on grid cell considers all average-individual of PLS; m2]
            FPC_total_t2 = sum(FPC_grid_t2)

        endif

        fpc_max_tree = gc_area*0.95
        if (FPC_total_t2 .gt. fpc_max_tree) then
            !Self-Thinning imposed when total tree cover above 'FPC MAX TREE' [max. of occupation]
            !partitioned among tree PLS in proportion to this year's FPC increment

            fpc_dec = (FPC_total_t2 - fpc_max_tree)*(FPC_grid_t2/FPC_total_t2)
            mort = 1.0-((FPC_grid_t2-fpc_dec)/FPC_grid_t2)
            remaining = 1.0-mort

            ! ==================== REDUCES INDIVIDUAL BIOMASS AND INDIVIDUAL DENSITY ==================== !
            !New individual density and leaf/wood biomass (remaining)
            dens = dens*remaining 
            cl2 = cl2*remaining 
            cw2 = cw2*remaining 
            ! print*, 'NEW DENS', dens, 'NEW CL2', cl2, 'NEW CW2', cw2
        endif
    enddo

end program self_thinning