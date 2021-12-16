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


! module allocation

!variable names

! REMAINING FILE | CAETE MODEL
!

! end module allocation


!---------------------------------


! module allometry

!Variable names (inputs)
! REMAINING FILE | CAETE MODEL
!   dwood        |     dt (valor da posição 18 no array de traits) --! wood density
!   cw1          |     cawood1
! end module allometry

module FPC
    
    implicit none
    public
!-------------------------------------------------------------------------------------------------
    !Module defining functions related to the Foliar projective cover of PLS.
    !It informs the portion of the gridcell/plot that is occupied by each PLS.
    !Most of the equations are based on LPJ population mode (as described in Smith et al 2001) 
    !and can be found in Sitch et al., 2003 and Smith et al., 2001.
    !In this mode, each PLS is represented by a population of average individuals.
    !All individuals of a PLS present the same features and the population is defined by
    !the density of individual per m2.
!--------------------------------------------------------------------------------------------------    
 
    integer, parameter:: npls = 20      !this number is being used for testing purpose
  

contains

     
    subroutine show_var(npls,dwood,FPC_ind)          
        
        integer, intent(in):: npls
       
        real, dimension(npls),intent(in):: dwood
       
        real, dimension(npls),intent(out):: FPC_ind

        
            FPC_ind = 1.*dwood
       
        print*, "npls =  ", npls 
        print*,'FPC_ind = ', FPC_ind
        print*, 'ok'            
    end subroutine show_var

end module FPC

program test_FPC
       
    use FPC

    real, dimension(npls) :: dwood      !wood density (g/cm-3) *Fearnside, 1997 - aleatory choices
    real, dimension(npls):: FPC_ind
    
    ! ! !Arrays with generic values for start and or test the code logic
    dwood=(/0.24,0.53,0.39,0.32,0.31,0.44,0.66,0.42,0.74,0.39,0.82,0.40,0.26,0.79,0.39,0.52,0.41,0.44,0.86,0.42/) !atenção para a unidade
    
    ! ! cl1=(/.7,1.,0.3,1.6,1.1,1.8,0.3,0.2,0.8,.84,0.25,1.,0.2,1.7,0.4,.6,.5,.8,0.3,1.8/)
    ! ! spec_leaf=(/0.0153,0.0101,0.0107,0.0112,0.012,0.0141,0.0137,0.0115,0.0122,0.010,0.012,0.011,&
    ! ! &0.013,0.014,0.0112,0.012,0.0141,0.0137,0.0115,0.0122/)
    ! ! cw1=(/30.,22.,34.,28.3,20.2,19.7,27.5,19.5,20.,28.6,24.3,19.3,26.8,22.,18.3,22.,15.,22.6,10.7,21.4/)
    ! ! cr1=(/0.63,0.8,0.9,1.5,1.3,0.9,0.4,1.0,0.56,0.87,0.33,0.97,0.31,0.55,0.2,0.8,0.4,0.66,0.23,1.5/)
    ! ! npp1 = (/0.5,0.8,1.5,1.2,1.9,1.3,1.7,0.8,0.6,2.0,0.7,1.1,1.9,1.85,1.96,1.77,1.33,1.54,1.62,0.55/)
    ! ! diameter = (/0.16,0.45,0.17,0.25,0.34,0.4,0.23,0.49,0.37,0.5,0.53,0.12,0.75,0.22,0.63,0.31,0.41,0.63,0.52,0.15/)
    ! ! dens = (/1.,2.,5., 0.3, 0.1, 7., 0.8, 3.,0.5,0.7,0.6,9.,4.,0.45,0.27,4.6,8.2,0.29,3.,0.8/)

    call show_var(npls,dwood,FPC_ind)

end program test_FPC

