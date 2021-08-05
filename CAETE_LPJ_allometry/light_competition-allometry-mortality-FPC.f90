program light_competition

    type :: layer_array
        real :: sum_height
        integer :: num_height !!corresponds to the number of pls
        real :: mean_height
        real :: layer_height
        real :: sum_LAI !LAI sum in a layer
        real :: mean_LAI !mean LAI in a layer
        real :: beers_law !layer's light extinction
        real :: li !layer's light incidence
        real :: lu !layer's light used (relates to light extinction - Beers Law)
        real :: la !light availability
    end type layer_array

    integer, parameter::npls=14
    real, dimension(npls), allocatable :: height (:)
    real, dimension(npls), allocatable :: Cleaf (:) !leaf increment
    real, dimension(npls), allocatable :: LA (:) !Leaf Area  (m2)
    real, dimension(npls), allocatable :: LAI (:) !Leaf Area Index (m2/m2)
    real, dimension(npls), allocatable :: diam (:) !Tree diameter in m. (Smith et al., 2001 - Supplementary)
    real, dimension(npls), allocatable :: crown_area (:) !Tree crown area (m2) (Sitch et al., 2003)
    real, dimension(npls), allocatable :: mort_WD (:) !mortality associated to WD (Sakschewski et al., 2014-SM)
    real, allocatable :: FPCind (:) !Foliage projective cover for each PLS (Sitch et al., 2003)
    real, allocatable :: FPCgrid (:) !Fractional projective cover in grid cell (Sitch et al., 2003)
    real, allocatable :: FPCgrid_perc (:) !Fractional projective cover in grid cell relative to grid cell area (Sitch et al., 2003)
    integer, allocatable :: nind (:) !number of individuals per PLS (Smith, 2001, thesis)
    real, allocatable :: nind_red (:) !reduced number of individuals per PLS according to the self-thinning rule (Smith, 2001, thesis)
    real, allocatable :: new_diam (:) !diameter updated with new number of ind (for testing purpose)
    real, allocatable :: Hcrit (:) !critical buckling height (in m) to be mechanic stable (Langam, 2017)
    real, allocatable :: FPCgrid_updt (:) !Fractional projective cover in grid cell (Sitch et al., 2003) - updated to occupy 95%
    real, allocatable :: FPCgrid_perc_updt (:) !Fractional projective cover in grid cell relative to grid cell area (Sitch et al., 2003)-updated to occupy 95%


    real :: max_height
    integer :: num_layer
    real :: layer_size
    real :: incidence_rad !Incidence radiation (relates do APAR) in J/m2/s
    real :: watt_rs = 210 !shortwave radiation in watts/m2
    real :: short_rad !shortwave radiation in joules/s
    real :: k_allom1 = 100. !allometric constant (Table 3; Sitch et al., 2003)
    real :: k_allom2 = 36.0 !allometric constant Value from Philip's code (Allocation.py)
    real :: k_allom3 = 0.22 !allometric constant Value from Philip's code (Allocation.py)
    real :: krp = 1.6 !allometric constant (Table 3; Sitch et al., 2003)
    real :: klatosa = 6000.0 !leaf_area:sapwood_area. Value from Philip's code (Allocation.py)
    real :: ltor = 0.77302587552347657 !leaf:root  Value from Philip's code (Allocation.py)
    !real :: spec_leaf = 15.365607091853349 !specific leaf area (variant trait) generic value to calculate leaf area index (LAI) value from Philip's code (Allocation.py)
    real :: sum_FPCgrid=0.0
    real :: sum_FPCgrid_perc=0.0
    real :: sum_nind=0.0
    real :: mort_occupation=0.0
    real :: gc_area = 300 !grid cell size - 300 m2 FOR TESTING PURPOSE (the real value will be 1ha or 10000 m2)
    real :: gc_area_95 = 0. !95% of grid cell size 
    real :: sum_FPCgrid_updt = 0.0 !!the new number of total PLS average-individuals after % reduction equals maximum to
                                 !! not to exceed 95% occupation.
    real :: sum_FPCgrid_perc_updt = 0.0 !the new percentage of occupation of all PLS after % reduction. 
    real :: mort_LPJ = 0.0 !testing the base-area Mortality
    real :: mort_max=0.01
    real :: par_min=4.05 !valor chutado
    real :: k_par=4.05
    real :: log_teste=0.0

    integer::i,j

    integer :: last_with_pls

    type(layer_array), allocatable :: layer(:)

! Variables with generic values for testing the logic code

    real, dimension(npls) :: dwood !wood density (g/cm-3) *Fearnside, 1997 - aleatory choices
    real, dimension(npls) :: carbon_stem !KgC/m2 (Cheart + Csap)
    real, dimension(npls) :: carbon_leaf !KgC/m2 
    real, dimension(npls) :: spec_leaf

    dwood=(/0.74,0.73,0.59,0.52,0.41,0.44,0.86,0.42,0.64,0.69,0.92,0.60,0.36,0.99/) !atenção para a unidade
    carbon_leaf=(/2.15,3.,1.18,1.6,1.5,1.8,0.3,2.,0.8,.84,0.25,1.,0.2,1.7/)
    spec_leaf=(/15.35,10.1,10.7,11.2,12.,14.1,13.7,11.5,12.2,10.,12.,11.,13.,14./)

! Allometric Equations
    
    diam = ((4+(carbon_stem))/((dwood)*3.14*40))**(1/(2+0.5))
    print*, 'diam', diam

    height = k_allom2*(diam**k_allom3)
    print*, 'height', height
    
    crown_area = k_allom1*(diam**krp)
    print*, 'crown', crown_area

    LAI = (carbon_leaf*spec_leaf)/crown_area
    print*, 'LAI', LAI

 

! Grid-Cell Properties

    allocate (nind(1:npls))
    allocate (nind_red(1:npls))
    allocate (new_diam(1:npls))
    allocate (FPCind(1:npls))
    allocate (FPCgrid(1:npls))
    allocate (FPCgrid_perc(1:npls))
    allocate (FPCgrid_updt(1:npls))
    allocate (FPCgrid_perc_updt(1:npls))
    allocate (mort_WD(1:npls))



!!!!===========		QUESTIONS TO BE SOLVED      ===========!!!!
    !!!to be verified (what is the size of the gc? and also will we deal with bare soil?
    !!!what about the grasses?
!!!!=======================================================!!!!

    do j=1,npls
        nind(j) = diam(j)**(-1.6)
        print*, 'Nind', nind(j)

        FPCind(j) = (1-exp(-0.5*LAI(j)))
        !print*, 'FPC', FPCind(j)

        FPCgrid(j) = crown_area(j)*nind(j)*FPCind(j)
        !print*, 'FPC-GRID', FPCgrid(j)

        FPCgrid_perc(j) = (FPCgrid(j)*100)/gc_area
        !print*, 'FPC-GRID-PERC', FPCgrid_perc(j), gc_area
    enddo

! Mortality Dynamic
	!!!!Mortality associated to WD
	do j=1,npls
		mort_WD(j)=(0.255/dwood(j))
		!mort_WD(j)=
		!print*,'WD mortality',mort_WD(j)
	enddo




    !Mortality relates to height critic (Langam, 2017 - aDGVM)

    allocate (Hcrit(1:npls))

    do j=1,npls
        !dwood*1000 -> transform wood density to g cm-3 to kg m-3
        !the symbol 'rho' is diam/height -> equal to resistivity (physical concepts)

        Hcrit(j) = (0.79*(((11.852*(diam(j)/height(j))+37)/9.81)*&
        &((dwood(j)*1000)**1/3))*(diam(j)**2/3))
        !print*, 'Hcrit', Hcrit(j)
    enddo

    !Mortality relates to Light Competition 

    !The result of 'mligh' formulation is the % of reduction of PLS ocupation. 
    !This value indicate the reduction percent of tree FPC due to the light (or space) competition.

    do j=1,npls
        sum_FPCgrid = sum_FPCgrid+FPCgrid(j)
        sum_FPCgrid_perc = sum_FPCgrid_perc+FPCgrid_perc(j)
        sum_nind = sum_nind+nind(j)
    enddo
    
    !print*, 'NIND_SUM', sum_nind
    !print*, 'SUM-FPC-GRID-PERC', sum_FPCgrid_perc
    !print*, 'SUM-FPC-GRID', sum_FPCgrid

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
    
    !Percentage of tree population reduction in all area: (IAP-DGVM; Zeng et al., 2014) - !ATTENTION!
    mort_occupation = (1.-(0.95/sum_FPCgrid))*sum_nind
    !print*, '*******mort_light', mort_occupation

    !!!!testing the reduction in individuals number!!!
    do j=1,npls

		FPCgrid_updt(j)=FPCgrid(j)-FPCgrid(j)*(mort_occupation/100)
		!print*, 'FPC with mort_occupation', FPCgrid_updt(j)

        FPCgrid_perc_updt(j) = (FPCgrid_updt(j)*100)/gc_area
        !print*, 'FPC-GRID-PERC with mort_occupation', FPCgrid_perc_updt(j), gc_area
    enddo	

    

    do j=1,npls
        sum_FPCgrid_updt = sum_FPCgrid_updt+FPCgrid_updt(j)
        sum_FPCgrid_perc_updt = sum_FPCgrid_perc_updt+FPCgrid_perc_updt(j)
    enddo
    !print*, 'sumfpc updated', sum_FPCgrid_updt
    !print*, 'sumfpc_perc updated', sum_FPCgrid_perc_updt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!testing LPJ equation!!!!!
	mort_LPJ=((mort_max)*(par_min/k_par)*(exp(-5*(1-(sum_FPCgrid_perc/100)))))*100
	!print*, '*************mort_LPJ', mort_LPJ, sum_FPCgrid_perc

! Layer's dynamics

    max_height = maxval(height)
    ! print*, 'max_height',max_height
    
    num_layer = nint(max_height/5)
    ! print*, 'num_layer',num_layer

    allocate(layer(1:num_layer))
    
    last_with_pls=num_layer

    layer_size = max_height/num_layer
    ! print*, 'layer_size', layer_size

    layer(i)%layer_height=0

    do i=1,num_layer
        layer(i)%layer_height=layer_size*i
        ! print*, 'layer_height',layer(i)%layer_height, i
    enddo

    layer(i)%num_height=0
    layer(i)%sum_height=0
    layer(i)%mean_height=0
    layer(i)%sum_LAI=0
   

    do i=1, num_layer
        do j=1,npls
            if ((layer(i)%layer_height.ge.height(j)).and.&
                &(layer(i-1)%layer_height.lt.height(j))) then

                layer(i)%sum_height=&
                &layer(i)%sum_height+height(j)

                layer(i)%num_height=&
                &layer(i)%num_height+1

                layer(i)%sum_LAI=&
                &layer(i)%sum_LAI+LAI(j)
                                      
            endif
        enddo

        layer(i)%mean_height=layer(i)%sum_height/&
            &layer(i)%num_height

        if(layer(i)%sum_height.eq.0.) then
            layer(i)%mean_height=0.
        endif

        layer(i)%mean_LAI=layer(i)%sum_LAI/&
            &layer(i)%num_height
             !print*,'mean_LAI',layer(i)%mean_LAI

        if(layer(i)%sum_LAI.eq.0.) then
            layer(i)%mean_LAI=0.
            !  print*, 'mean_LAI2', layer(i)%mean_LAI
        endif
        
        ! print*,'lyr',i,'mean_height',&
        !     &layer(i)%mean_height,'lai',&
        !     &layer(i)%mean_LAI
    enddo

    layer(i)%li = 0

    layer(i)%la = 0

    layer(i)%lu = 0

!Light Extinction

    short_rad = watt_rs*1000.
    ! print*,'short_rad',short_rad

    incidence_rad = 0.5*short_rad
    ! print*,'APAR', incidence_rad

!=================== TEST ====================
    do i=num_layer,1,-1
        layer(i)%beers_law = incidence_rad*&
            &(1-exp(-0.5*layer(i)%mean_LAI))
        !  print*,'law',layer(i)%beers_law
    enddo
!=============================================


    do i=num_layer,1,-1
        if(i.eq.num_layer)then
            layer(i)%li = incidence_rad
        else
            if(layer(i)%mean_height.gt.0.)then

                layer(i)%li = layer(last_with_pls)%la
                last_with_pls=i

            else
                continue
            endif
        endif

        layer(i)%lu = layer(i)%li * (1-exp(-0.5*layer(i)%mean_LAI))
        
        layer(i)%la = layer(i)%li - layer(i)%lu

        ! print*,i, 'inc', layer(i)%li, 'used', layer(i)%lu,& 
        !     & 'avai', layer(i)%la

    enddo
    

end program light_competition