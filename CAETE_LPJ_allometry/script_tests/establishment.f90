! Copyright 2017- LabTerra

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR 2PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have recei(i_4)am.  If not, see <http://www.gnu.org/licenses/>.

! AUTHOR: Bianca Fazio Rius, Bárbara Cardeli, Carolina Blanco



!-------------------------------------------------------------------------------------------------
    !Module defining functions related to the establishment of average individuals in
    !a PLS if the area of occupation of all PLSs is greater than 95% of the considered area
    !It informs both the new density of individuals and the "shrink" process (see section 
    !4.6.1.Establishment in Smith et al. 2001)
    !Most of the equations are based on LPJ population mode (as described in Smith et al 2001) 
    !and can be found in Sitch et al., 2003 and Smith et al., 2001.
    !In this mode, each PLS (i_4)is represented by a population of average individuals.
    !All individuals of a PLS present the same features and the population is defined by
    !the density of individuals per m2.
!--------------------------------------------------------------------------------------------------   


module establish
    
    use types

    implicit none
    public ::                    &      
        establishment           ,&
        shrink

    contains

    subroutine establishment(j,gc_available,npls_alive, FPC_total_accu_2, gc_area, est, est_pls, FPC_pls)
    implicit none

    !input variables
    integer(i_4), intent(in) :: j
    real(r_8), intent(in) :: npls_alive
    real(r_8), intent(in) :: FPC_total_accu_2
    real(r_8), intent(in) :: gc_area
    real(r_8), intent(in) :: gc_available
    real(r_8), intent(in) :: FPC_pls
    
    

    !output variables
    real(r_8), intent(out):: est
    real(r_8), intent(out):: est_pls
    
    !internal variables
    real(r_8) :: est_max  !2 individuals m -2 yr -1 - reference: Levis et al 2004 (Eq 53)
    real(r_8) :: FPC_total_perc
    ! print*, 'alive pls in est', npls
    real(r_8), parameter :: dens_min = 1.e-1      !minimum individual density for persistence of PFT (indiv/m2)
    real(r_8), parameter :: aprec_min_estab = 100. !miminum annual prec for establishment (mm) !reference: code LPL-LMfire
    real(r_8), parameter :: eps = 1.e-6

    ! print*, npls_alive
    ! if (dens.le.dens_min) then
    !     print*, 'DENS MIN', dens, j
    ! endif

    est_max = 2 !based on Sitch et al 2003
    

    FPC_total_perc = FPC_total_accu_2/gc_area !The overall establishment rate for trees (esttree) is proportional to the fractional ground area not covered by trees.
    ! print*, 'fpc perc', FPC_total_perc

    ! if(FPC_total_perc.le.0.9) then
    !     est = 0.06*(1-FPC_total_perc)
    !     ! print*, '<0.9',est
    ! else
    !     est = 0.06*(1-exp(-5*(1-FPC_total_perc)))*(1-FPC_total_perc)
    !     ! print*, '>0.9',est
    ! endif

    ! est_pls = est*(est_max/est_max*npls_alive)*FPC_pls*(1-FPC_total_perc)
    ! print*, 'est_pls smith', est_pls  

    est = est_max * (1. - exp(-5. * (1-FPC_total_perc)) / npls_alive)
    est_pls = max(est * (1. - FPC_total_perc),0.)
    
    ! est_pls = ((est_max*(1-exp(-5*(1-FPC_total_perc))))*(1-FPC_total_perc))/npls_alive
    ! print*, 'est pls sitch', est_pls

    end subroutine

    subroutine shrink(sla,wooddens,cl_old,ch_old,cs_old,cw_old,cr_old,est_pls,dens_old,cleaf_sapl_npls,csap_sapl_npls,&
    &           cheart_sapl_npls,croot_sapl_npls, dens_new, cleaf_new,& 
    &           cwood_new,cheart_new, csap_new, croot_new)
    implicit none

        ! !input variables
        ! integer, intent(in) :: npls
        real(r_8), intent(in) :: cl_old
        real(r_8), intent(in) :: cw_old
        real(r_8), intent(in) :: ch_old
        real(r_8), intent(in) :: cs_old
        real(r_8), intent(in) :: cr_old
        real(r_8), intent(in) :: est_pls
        real(r_8), intent(in) :: dens_old
        real(r_8), intent(in) :: cleaf_sapl_npls
        real(r_8), intent(in) :: csap_sapl_npls
        real(r_8), intent(in) :: cheart_sapl_npls
        real(r_8), intent(in) :: croot_sapl_npls
        real(r_8), intent(in) :: wooddens
        real(r_8), intent(in) :: sla

        ! !output variables
        real(r_8), intent(out):: dens_new
        real(r_8), intent(out):: cleaf_new
        real(r_8), intent(out):: cwood_new
        real(r_8), intent(out):: cheart_new
        real(r_8), intent(out):: csap_new
        real(r_8), intent(out):: croot_new
        
        ! !internal variables
        ! real(r_8) :: dens_est_pls 
        
        real(r_8) :: cwood_sapl_npls
        real(r_8) :: csap_new_tmp
        real(r_8) :: allom1 = 100. !allometric constant (Table 3; Sitch et al., 2003)
        real(r_8) :: allom2 = 40.0
        real(r_8) :: allom3 = 0.85
        real(r_8) :: pi = 3.1415
        real(r_8) :: stemdiam
        real(r_8) :: height
        real(r_8) :: crown_area
        real(r_8) :: crown_area_max = 30 !m2 !number from lplmfire code (establishment.f90)
        real(r_8) :: reinickerp = 1.6
        real(r_8) :: latosa = 8000.0
        
        
        
        if(cl_old.le.0.)then
            dens_new = 0.
            cleaf_new = 0.
            cheart_new = 0.
            csap_new = 0.
            csap_new_tmp = 0.
            cwood_sapl_npls = 0.
            cwood_new = 0.
            croot_new = 0.
        else
            dens_new = dens_old + est_pls


        ! print*, dens_new, dens_old 

        ! print*, 'dens_new', dens_new, 'dens_old', dens_old

            ! csap_new_tmp = ((cs_old*dens_old)+(csap_sapl_npls*est_pls))/dens_new !temporary since part of it goes to heartwood
            csap_new = ((cs_old*dens_old)+(csap_sapl_npls*est_pls))/dens_new !temporary since part of it goes to heartwood


            cleaf_new = ((cl_old*dens_old)+(cleaf_sapl_npls*est_pls))/dens_new
        ! print*,'cleaf_new',cleaf_new,'cl_old', cl_old, 'dens_olds', dens_old,'cleaf_sapl',cleaf_sapl_npls

            cheart_new = ((ch_old*dens_old)+(cheart_sapl_npls*est_pls))/dens_new

            cwood_new = cheart_new + csap_new
            ! print*,'cw_new',cwood_new/1000.

            croot_new = ((cr_old*dens_old)+(croot_sapl_npls*est_pls))/dens_new
        ! print*, 'cr new', croot_new/1000., 'crold', cr_old/1000.

            ! stemdiam = (4. * (csap_new_tmp + cheart_new) / wooddens*1000000. / pi / allom2)**(1./(2. + allom3)) !Eqn 9

            ! height = allom2 * stemdiam**allom3                           !Eqn C

            ! crown_area = min(crown_area_max,allom1 * stemdiam**reinickerp) !Eqn D


            !Recalculate sapwood mass, transferring excess sapwood to heartwood compartment, if necessary to satisfy Eqn A

            ! csap_new = cleaf_new * height * wooddens * sla / latosa 

            ! cheart_new = max(cheart_new + (csap_new_tmp - csap_new),0.)
            ! print*, 'ch new', cheart_new, csap_new_tmp, csap_new
        endif
        
    end subroutine shrink

    subroutine sapling_allometry(npls_alive,cleaf_sapl_npls, csap_sapl_npls, cheart_sapl_npls,croot_sapl_npls)
        !subroutine to calculate carbon in saplings compartments

        implicit none
        !input variables
        real(r_8), intent(in) :: npls_alive
        
        ! !output variables
        real(r_8), intent(out) :: cleaf_sapl_npls     
        real(r_8), intent(out) :: csap_sapl_npls
        real(r_8), intent(out) :: cheart_sapl_npls
        real(r_8), intent(out) :: croot_sapl_npls

        real(r_8) :: cleaf_sapl !gC
        real(r_8) :: csap_sapl  !gC
        real(r_8) :: cheart_sapl !gC
        real(r_8) :: croot_sapl   !gC


        !internal variables
       

        real(r_8) :: sla_sapl, diam_sapl, lai_sapl, height_sapl,dwood_sapl

        real(r_8) :: klatosa_sapl = 8000.

        real(r_8) :: pi = 3.14159

        real(r_8) :: k_allom1_sapl = 100.

        real(r_8) :: k_allom2_sapl = 40.

        real(r_8) :: k_allom3_sapl = 0.5

        real(r_8) :: x_sapl = 3. !from lpjlmfire (pftparametersmod.f90)

        real(r_8) :: reinickerp = 1.6
        

       
        sla_sapl = 0.021 !from lpjlmfire (pftparametersmod.f90, line 229) m2/gC
        lai_sapl = 4 !lpjmlfire (pft parameter)
        dwood_sapl = 2e5 !gc/m3

        cleaf_sapl = (lai_sapl * k_allom1_sapl * x_sapl**reinickerp * (4. *sla_sapl / pi / klatosa_sapl)**(reinickerp * 0.5) / & 
                      sla_sapl)**(1. - 1. / reinickerp)  
        !print*, 'cleafsapl, lpjmlfire', cleaf_sapl

        diam_sapl = x_sapl * (4. * cleaf_sapl * sla_sapl / pi / klatosa_sapl)**0.5

        !print*, 'diamsapl, lpjmlfire', diam_sapl*100

        height_sapl = k_allom2_sapl*diam_sapl**k_allom3_sapl
        !print*, 'height, lpjmlfire', height_sapl

        csap_sapl = dwood_sapl * height_sapl * cleaf_sapl * sla_sapl / klatosa_sapl
        !print*, 'csapl, lpjmlfire', csap_sapl

        cheart_sapl = (x_sapl - 1.) * csap_sapl
        ! print*, 'cheart, lpjmlfire', cheart_sapl
        

        !update to gC and to distribution to the pls's

        !cleaf_sapl_npls = (cleaf_sapl*1000)/npls_alive
        !csap_sapl_npls = (csap_sapl*1000)/npls_alive
        !cheart_sapl_npls = (cheart_sapl*1000)/npls_alive
        !croot_sapl_npls = (croot_sapl*1000)/npls_alive

        cleaf_sapl_npls  = cleaf_sapl
        csap_sapl_npls   = csap_sapl
        cheart_sapl_npls = cheart_sapl
        croot_sapl_npls  = croot_sapl



        ! print*, cleaf_sapl_npls, csap_sapl_npls, cheart_sapl_npls, croot_sapl_npls


    end subroutine sapling_allometry

end module establish