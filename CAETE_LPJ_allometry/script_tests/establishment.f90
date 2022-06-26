! Copyright 2017- LabTerra

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR 2PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

! AUTHOR: Bianca Fazio Rius, Bárbara Cardeli, Carolina Blanco



!-------------------------------------------------------------------------------------------------
    !Module defining functions related to the establishment of average individuals in
    !a PLS if the area of occupation of all PLSs is greater than 95% of the considered area
    !It informs the both the new density of individuals and the "shrink" process (see section 
    !4.6.1.Establishment in Smith et al. 2001)
    !Most of the equations are based on LPJ population mode (as described in Smith et al 2001) 
    !and can be found in Sitch et al., 2003 and Smith et al., 2001.
    !In this mode, each PLS is represented by a population of average individuals.
    !All individuals of a PLS present the same features and the population is defined by
    !the density of individuals per m2.
!--------------------------------------------------------------------------------------------------   


module establish
    
    implicit none
    public ::                    &      
        establishment           ,&
        shrink

    contains

    subroutine establishment(gc_available,npls, FPC_total_accu_2, gc_area, est, est_pls)
    
    !input variables
    real, intent(in) :: npls
    real, intent(in) :: FPC_total_accu_2
    real, intent(in) :: gc_area
    real, intent(in) :: gc_available
    

    !output variables
    real, intent(out):: est
    real, intent(out):: est_pls
    
    !internal variables
    real :: est_max  !2 individuals m -2 yr -1 - reference: Levis et al 2004 (Eq 53)
    real :: FPC_total_perc
    ! print*, 'alive pls in est', npls
    
    est_max = 2*(gc_available)
    ! est_max = 1*(gc_area)
    FPC_total_perc = FPC_total_accu_2/gc_area
        
        ! print*, 'fpc perc', FPC_total_perc

    if(FPC_total_perc.lt.0.9) then

        est = est_max*(1 - FPC_total_perc)
        ! print*, 'lt 0.9'
    else

        est = (1 - exp(-5. * (1 - FPC_total_perc)))*(1 - FPC_total_perc)
        ! print*, 'gt 0.9'

    endif

        est_pls = est/npls

        ! print*, 'estab', est, est_pls, npls
    
    end subroutine

    subroutine shrink(cl_old,cw_old,cr_old,est_pls,dens_old,cleaf_sapl_npls,csap_sapl_npls,&
    &           cheart_sapl_npls,croot_sapl_npls, dens_new, cleaf_new,& 
    &           cwood_new,croot_new)
    
        ! !input variables
        ! integer, intent(in) :: npls
        real, intent(in) :: cl_old
        real, intent(in) :: cw_old
        real, intent(in) :: cr_old
        real, intent(in) :: est_pls
        real, intent(in) :: dens_old
        real, intent(in) :: cleaf_sapl_npls
        real, intent(in) :: csap_sapl_npls
        real, intent(in) :: cheart_sapl_npls
        real, intent(in) :: croot_sapl_npls

        ! !output variables
        real, intent(out):: dens_new
        real, intent(out):: cleaf_new
        real, intent(out):: cwood_new
        real, intent(out):: croot_new
        
        ! !internal variables
        ! real :: dens_est_pls 
        
        real :: cwood_sapl_npls
        
        
        dens_new = dens_old + est_pls


        ! print*, dens_new, dens_old 

        ! print*, 'dens_new', dens_new, 'dens_old', dens_old

        cleaf_new = ((cl_old*dens_old)+(cleaf_sapl_npls*est_pls))/dens_new
        ! print*,'cleaf_new',cleaf_new,'cl_old', cl_old, 'dens_olds', dens_old,'cleaf_sapl',cleaf_sapl_npls

        cwood_sapl_npls = csap_sapl_npls + cheart_sapl_npls
        ! print*, 'cwood_sapl', cwood_sapl_npls, 'est_pls',est_pls, cwood_sapl_npls*est_pls

        cwood_new = ((cw_old*dens_old)+(cwood_sapl_npls*est_pls))/dens_new
        ! print*,'cw_new',cwood_new/1000.,'cw_old', cw_old/1000.


        croot_new = ((cr_old*dens_old)+(croot_sapl_npls*est_pls))/dens_new
        ! print*, 'cr new', croot_new/1000., 'crold', cr_old/1000.

    end subroutine shrink

    subroutine sapling_allometry(npls,cleaf_sapl_npls, csap_sapl_npls, cheart_sapl_npls,croot_sapl_npls)
    
        !input variables
        real, intent(in) :: npls
        
        ! !output variables
        real, intent(out) :: cleaf_sapl_npls     
        real, intent(out) :: csap_sapl_npls
        real, intent(out) :: cheart_sapl_npls
        real, intent(out) :: croot_sapl_npls

        real :: cleaf_sapl 
        real :: csap_sapl
        real :: cheart_sapl
        real :: croot_sapl


        !internal variables
       
        real :: aux, aux1, aux2, aux3, aux4, sla ! auxiliary variables to calculate cleaf_sapl

        real :: aux5, aux6, aux7, aux8, aux9, aux10, dwood ! auxiliary variables to calculate cleaf_sapl


        real :: klatosa_sapl = 8000.

        real :: pi = 3.14159

        real :: k_rp = 1.6

        real :: sapl_hw = 0.2

        real :: lai_sapl = 1.5

        real :: k_allom1_sapl = 100.

        real :: k_allom2_sapl = 40.

        real :: k_allom3_sapl = 0.5

        

        sla = 20 !m2/kgC
        lai_sapl = 1.5 !m2/m2
        dwood = 200 !kgC/m3
        
        
        ! sla_sapl2 = (2*(0.0001))*((exp(6.15)/((12.*2.)**0.46)))
        ! print*, 'sla_sapl', sla_sapl, sla_sapl2

        ! sla = sla_sapl/1000. !transforms from m2/gC to m2/KgC
       
        ! dwood = dwood_sapl*1000 !transforms from gC/cm3 to KgC/m3
        

        ! aux4 = (4*sla)/(pi/klatosa_sapl)
  

        ! ! aux3 = ((4*sla_sapl)/(klatosa_sapl*pi))**(k_rp*0.5)
        ! aux3 = (aux4)**(k_rp*0.5)
        
        ! aux2 = (1.0 + sapl_hw)**k_rp

        ! aux1 = lai_sapl*k_allom1_sapl

        ! aux = aux1 * aux2 * aux3
              
       
        ! ! cleaf_sapl = (aux/(sla**(2.0/(2.0-k_rp))))
        ! cleaf_sapl = (aux/sla)**(2.0/(2.0-k_rp))

        ! ! print*,'cleaf_sapl',cleaf_sapl

        ! ! aux5 = dwood * k_allom2_sapl

        ! ! aux6 = (1. + sapl_hw)

        ! ! aux7 = (4*cleaf_sapl*sla)/(pi)

        ! ! aux8 = sqrt(aux7/klatosa_sapl)

        ! ! aux9 = (cleaf_sapl*sla)/klatosa_sapl

        ! ! aux10 = (aux6*aux8)**k_allom3_sapl

        ! ! csap_sapl = aux5 * aux10 * aux9

        ! ! ! csap_sapl = (dwood*(1+sapl_hw)*cleaf_sapl*sla)/klatosa_sapl

        ! ! print*, 'csap_sapl', csap_sapl*1000

        ! ! cheart_sapl = (1+sapl_hw) * csap_sapl

        ! ! print*, 'cheart_sapl', cheart_sapl*1000. + csap_sapl*1000
        
        ! aux1 = 4*sla_sapl
! 
        ! aux2 = pi*klatosa_sapl

        ! aux3 = aux1/aux2

        ! aux4 = aux3**(0.5*k_rp)
! 
        ! ! ! aux5 = 1.2**k_rp
! 
        ! ! ! aux6 = 1.5*k_allom1_sapl
! 
        ! ! aux7 = 2./(2-k_rp)
! 
        ! ! aux8 = sla_sapl**aux7
! 
        ! aux9 = aux6*aux5*aux4
        
        ! cleaf_sapl = aux9/aux8
        
        aux3=lai_sapl*k_allom1_sapl
        aux4=(1.+sapl_hw)**k_rp
        aux7=((4*sla)/pi)/klatosa_sapl
        aux6=k_rp*0.5
        aux5 = aux7**aux6
        aux1 = (aux3*aux4*aux5)/sla
        aux2 = 2./(2.-k_rp)
        aux = aux1**aux2
        cleaf_sapl = aux !*1000
        
        
        ! print*,'cleaf', cleaf_sapl
        
        aux10 = (cleaf_sapl*sla)/klatosa_sapl
        aux9 = (1+sapl_hw)*sqrt(((4*cleaf_sapl*sla)/pi)/klatosa_sapl)
        aux8 = dwood*k_allom2_sapl
        csap_sapl = (aux8 * (aux9**k_allom3_sapl)*aux10) !*1000.
        ! print*, csap_sapl

        cheart_sapl = (sapl_hw*csap_sapl)!*1000.
        ! print*, cheart_sapl

        croot_sapl = (1./1.*cleaf_sapl)!*1000.

        

        !update to gC and to distribution to the pls's

        cleaf_sapl_npls = (cleaf_sapl*1000)/npls
        csap_sapl_npls = (csap_sapl*1000)/npls
        cheart_sapl_npls = (cheart_sapl*1000)/npls
        croot_sapl_npls = (croot_sapl*1000)/npls


        ! print*, cleaf_sapl_npls, csap_sapl_npls, cheart_sapl_npls, croot_sapl_npls


    end subroutine sapling_allometry

end module establish