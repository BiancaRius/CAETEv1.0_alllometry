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

    subroutine establishment (npls, FPC_total_accu_2, gc_area, est, est_pls)
    
    !input variables
    integer, intent(in) :: npls
    real, intent(in) :: FPC_total_accu_2
    real, intent(in) :: gc_area
  
    !output variables
    real, intent(out):: est
    real, intent(out):: est_pls
    
    !internal variables
    real :: est_max = 0.24 !individuals m -2 yr -1 - reference: Levis et al 2004 (Eq 53)
    real :: FPC_total_perc
       
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

        ! print*, 'estab', est, est_pls
    
    end subroutine

    subroutine shrink()
    
        ! !input variables
        ! integer, intent(in) :: npls
        ! real, intent(in) :: FPC_total_accu_2
        ! real, intent(in) :: gc_area
      
        ! !output variables
        ! real, intent(out):: est
        ! real, intent(out):: est_pls
        
        ! !internal variables
        ! real :: est_max = 0.24 !individuals m -2 yr -1 - reference: Levis et al 2004 (Eq 53)
        ! real :: FPC_total_perc
    
    end subroutine shrink

    subroutine sapling_allometry(sla_sapl,dwood_sapl, cleaf_sapl, csap_sapl)
    
        !input variables
        ! integer, intent(in) :: npls
        real, intent(in) :: sla_sapl 
        real, intent(in) :: dwood_sapl
      
        ! !output variables
        real, intent(out) :: cleaf_sapl 
        real, intent(out) :: csap_sapl
        !internal variables
       
        real :: aux, aux1, aux2, aux3
        
        real :: klatosa_sapl = 8000.

        real :: pi = 3.14159

        real :: k_rp = 1.6

        real :: sapl_hw = 0.2

        real :: lai_sapl = 1.5

        real :: k_allom1_sapl = 100.

        real :: k_allom2_sapl = 40.
        
        ! sla_sapl = (2*(0.0001))*((exp(6.15)/((12.*2.)**0.46)))
        ! print*, 'sla_sapl', sla_sapl

        aux3 = ((4*sla_sapl)/(klatosa_sapl*pi))**(k_rp*0.5)

        aux2 = (1.0 + sapl_hw)**k_rp

        aux1 = lai_sapl*k_allom1_sapl

        aux = aux1 * aux2 * aux3
              
       
        cleaf_sapl = (aux/(sla_sapl**(2.0/(2.0-k_rp))))

        print*,'cleaf_sapl',cleaf_sapl

        csap_sapl = 
        
        



    end subroutine sapling_allometry

end module establish