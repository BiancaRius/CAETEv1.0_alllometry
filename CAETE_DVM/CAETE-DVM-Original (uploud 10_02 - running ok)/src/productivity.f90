! Copyright 2017- LabTerra

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

! contacts :: David Montenegro Lapola <lapoladm ( at ) gmail.com>
! Author: JP Darela
! This program is based on the work of those that gave us the INPE-CPTEC-PVM2 model

module productivity
  implicit none
  private

  public :: prod,light_compet


contains

  subroutine prod(dt,light_limit,catm,temp,ts,p0,w,ipar,rh,emax,cl1_prod,&
       & ca1_prod,cf1_prod,beta_leaf,beta_awood,beta_froot,wmax,ph,ar,&
       & nppa,laia,f5,vpd,rm,rg,rc,wue,c_defcit,vm_out,sla, e)

    use types
    use global_par
    use allometry_par
    use photo_par
    use photo
    use water

!Input
!-----
    integer(i_4), parameter :: npft = npls
    real(r_8),dimension(ntraits),intent(in) :: dt ! PLS data
    real(r_4), intent(in) :: temp, ts                 !Mean monthly temperature (oC)
    real(r_4), intent(in) :: p0                   !Mean surface pressure (hPa)
    real(r_8), intent(in) :: w                    !Soil moisture kg m-2
    real(r_4), intent(in) :: ipar                 !Incident photosynthetic active radiation (w/m2)
    real(r_4), intent(in) :: rh,emax !Relative humidity/MAXIMUM EVAPOTRANSPIRATION
    real(r_8), intent(in) :: catm, cl1_prod, cf1_prod, ca1_prod        !Carbon in plant tissues (kg/m2)
    real(r_8), intent(in) :: beta_leaf            !npp allocation to carbon pools (kg/m2/day)
    real(r_8), intent(in) :: beta_awood
    real(r_8), intent(in) :: beta_froot, wmax
    logical(l_1), intent(in) :: light_limit                !True for no ligth limitation

!     Output
!     ------
    real(r_4), intent(out) :: ph                   !Canopy gross photosynthesis (kgC/m2/yr)
    real(r_4), intent(out) :: rc                   !Stomatal resistence (not scaled to canopy!) (s/m)
    real(r_8), intent(out) :: laia                 !Autotrophic respiration (kgC/m2/yr)
    real(r_4), intent(out) :: ar                   !Leaf area index (m2 leaf/m2 area)
    real(r_4), intent(out) :: nppa                 !Net primary productivity (kgC/m2/yr)
    real(r_4), intent(out) :: vpd
    real(r_8), intent(out) :: f5                   !Water stress response modifier (unitless)
    real(r_4), intent(out) :: rm                   !autothrophic respiration (kgC/m2/day)
    real(r_4), intent(out) :: rg
    real(r_4), intent(out) :: wue
    real(r_4), intent(out) :: c_defcit     ! Carbon deficit gm-2 if it is positive, aresp was greater than npp + sto2(1)
    real(r_8), intent(out) :: sla, e        !specific leaf area (m2/kg)
    real(r_8), intent(out) :: vm_out
!     Internal
!     --------

    real(r_8) :: tleaf,awood            !leaf/wood turnover time (yr)
    real(r_8) :: g1
    real(r_8) :: c4

    real(r_8) :: n2cl
    real(r_8) :: n2cl_resp
    real(r_8) :: n2cw_resp
    real(r_8) :: n2cf_resp
    real(r_8) :: p2cl
    integer(i_4) :: c4_int
    real(r_8) :: jl_out

    real(r_8) :: f1       !Leaf level gross photosynthesis (molCO2/m2/s)
    real(r_8) :: f1a      !auxiliar_f1
    ! real(r_8) :: diam1    !test print to diameter function (in m)
    ! real(r_8) :: area_crown !test print to crown area function (in m2)
    ! real(r_8) :: height_tree !test print to height function (in m)
    ! real(r_8) :: LAI_calc !test print to LAI function (in m2 m-2)
    ! real(r_8) :: max_height_tree 
    ! integer(i_4) :: p

!getting pls parameters


    g1  = dt(1)
    tleaf = dt(3)
    awood = dt(7)
    c4  = dt(9)
    n2cl = dt(10)
    n2cl_resp = n2cl
    n2cw_resp = dt(11)
    n2cf_resp = dt(12)
    p2cl = dt(13)


    n2cl = n2cl * (cl1_prod * 1D3) ! N in leaf g m-2
    p2cl = p2cl * (cl1_prod * 1D3) ! P in leaf g m-2

    c4_int = idnint(c4)
    
    ! diam1 = diameter(ca1_prod)
    !print*, 'diameter =', diam1

    ! area_crown = crownarea(diam1)
    ! !print*, 'crown area =', area_crown

    ! height_tree = tree_height(diam1)
    ! print*, 'height =', height_tree

!     ==============
!     Photosynthesis
!     ==============
! rate (molCO2/m2/s)

    call photosynthesis_rate(catm,temp,p0,ipar,light_limit,c4_int,n2cl,&
         & p2cl,tleaf,f1a,vm_out,jl_out)

    ! print*, 'jl: ',jl_out, catm

!    _____     ____   _     _____    _    _____
!    | __ )_ _/ ___| | |   | ____|  / \  |  ___|
!    |  _ \| | |  _  | |   |  _|   / _ \ | |_
!    | |_)|| | |_| | | |___| |___ / ___ \|  _|
!    |____/___\____| |_____|_____/_/   \_\_|
!     Leaf area index (m2/m2)
    sla = spec_leaf_area(tleaf)
    ! laia = leaf_area_index(cl1_prod, sla)

    laia = 0.2D0 * dexp((2.5D0 * f1a)/p25)
! VPD
!========
    vpd = vapor_p_defcit(temp,rh)

!Stomatal resistence
!===================
    rc = canopy_resistence(vpd, f1a, g1, catm) * real(laia, kind=r_4)!s m-1

! Novo calculo da WUE

    wue = water_ue(f1a, rc, p0, vpd)

!     calcula a transpiração em mm/s

    e = transpiration(rc, p0, vpd, 2)

!     Water stress response modifier (dimensionless)
!     ----------------------------------------------
    ! print*,cf1_prod, 'CF in F5'
    ! print*, w, 'w'
    ! print*, rc, 'rc'
    ! print*, emax, 'emax'

    ! wsoil + h +
    f5 =  water_stress_modifier(w, cf1_prod, rc, emax, wmax)


!     Photosysthesis minimum and maximum temperature
!     ----------------------------------------------

    if ((temp.ge.-10.0).and.(temp.le.50.0)) then
       f1 = f1a * f5 ! :water stress factor ! Ancient floating-point underflow spring (from CPTEC-PVM2)
    else
       f1 = 0.0      !Temperature above/below photosynthesis windown
    endif



!     Canopy gross photosynthesis (kgC/m2/yr)
!     =======================================x
    ph =  gross_ph(f1,cl1_prod,sla)       ! kg m-2 year-1

!     Autothrophic respiration
!     ========================
!     Maintenance respiration (kgC/m2/yr) (based in Ryan 1991)
    rm = m_resp(temp,ts,cl1_prod,cf1_prod,ca1_prod &
         &,n2cl_resp,n2cw_resp,n2cf_resp,awood)

! c     Growth respiration (KgC/m2/yr)(based in Ryan 1991; Sitch et al.
! c     2003; Levis et al. 2004)
    rg = g_resp(beta_leaf,beta_awood, beta_froot,awood)

    if (rg.lt.0) then
       rg = 0.0
    endif

!     c Autotrophic (plant) respiration -ar- (kgC/m2/yr)
!     Respiration minimum and maximum temperature
!     -------------------------------------------
    if ((temp.ge.-10.0).and.(temp.le.50.0)) then
       ar = rm + rg
    else
       ar = 0.0               !Temperature above/below respiration windown
    endif
!     Net primary productivity(kgC/m2/yr)
!     ====================================
    nppa = ph - ar
! this operation affects the model mass balance
! If ar is bigger than ph, what is the source or respired C?

    if(ar .gt. ph) then
       c_defcit = ((ar - ph) * 2.73791) ! tranform kg m-2 year-1 in  g m-2 day-1
       nppa = 0.0
    else
       c_defcit = 0.0
    endif

  end subroutine prod

    subroutine light_compet(cwood, cleaf)
        use types
        use photo
        use global_par, only: npls

        type :: layer_array
            real(r_8) :: sum_height
            integer(r_8) :: num_height !!corresponds to the number of pls
            real(r_8) :: mean_height
            real(r_8) :: layer_height
            real :: sum_LAI !LAI sum in a layer
            real :: mean_LAI !mean LAI in a layer
            real :: beers_law !layer's light extinction
            real :: linc !layer's light incidence
            real :: lused !layer's light used (relates to light extinction - Beers Law)
            real :: lavai !light availability
        end type layer_array

        integer(kind=i_4),parameter :: npft = npls
        real(r_8),dimension(npft), intent(in) :: cwood !C in wood tissues (Kg/m2)
        real(r_8),dimension(npft), intent(in) :: cleaf !C in leaf (kg/m2)
        real(r_8),dimension(npft) :: sla 
        real(r_8),dimension(npft) :: diam_aux, crown_aux, height_aux, lai_aux
        real(r_8) :: max_height !maximum height in m. in each grid-cell
        integer(kind=i_4) :: num_layer !number of layers according to max height in each grid-cell
        real(r_8) :: layer_size !size of each layer in m. in each grid-cell
        integer(kind=i_4) :: p, i, j
        integer(kind=i_4) :: last_with_pls
        type(layer_array), allocatable :: layer(:)


        do p=1,npft
            diam_aux(p) = diameter(cwood(p))

            crown_aux(p) = crownarea(diam_aux(p))

            height_aux(p) = tree_height(diam_aux(p))

            lai_aux(p) = leaf_area_index(cleaf(p), sla(p))

            !print*,diam_aux(p),crown_aux(p),height_aux(p), lai_aux(p)
        enddo

        max_height = maxval(height_aux(:))
        print*, 'max_height', max_height

        num_layer = nint(max_height/5)
        print*, 'num_layer', num_layer

        layer_size = max_height/num_layer
        print*, 'layer_size', layer_size

        last_with_pls=num_layer

        allocate(layer(1:num_layer))

        layer(i)%layer_height=0.0D0

        do i=1,num_layer
            layer(i)%layer_height=layer_size*i
            !print*, 'layer_height',layer(i)%layer_height, i
        enddo

        do i=1,num_layer
            layer(i)%num_height=0.0D0
            layer(i)%sum_height=0.0D0
            layer(i)%mean_height=0.0D0
            layer(i)%sum_LAI=0.0D0
        enddo
        
        do i=1, num_layer
            do j=1,npft
                
                if ((layer(i)%layer_height .ge. height_aux(j)).and.&
                    &(layer(i-1)%layer_height .lt. height_aux(j))) then
    
                    layer(i)%sum_height=&
                    &layer(i)%sum_height + height_aux(j)
                    print*, 'sum_height =', layer(i)%sum_height
    
                    layer(i)%num_height=&
                    &layer(i)%num_height+1
                    print*, 'num_height=', layer(i)%num_height
    
                    layer(i)%sum_LAI =&
                    &layer(i)%sum_LAI + lai_aux(j)
                    print*, 'sum_LAI =', layer(i)%sum_LAI
    
                endif
            enddo
            
        enddo

    end subroutine light_compet

end module productivity
