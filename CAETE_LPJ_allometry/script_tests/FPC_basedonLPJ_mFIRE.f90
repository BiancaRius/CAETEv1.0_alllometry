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
!Module defining functions related to occupation of the area simulated (FPC - fractional projective cover)
!Basically it checks if the area occupied by all PLS is higher than 95%. If so, there is a mortality
!factor that decrease the number of individuals per PLS
!Most of the equations are based on LPJ population mode (as described in Smith et al 2001) 
!and can be found in Sitch et al., 2003 and Smith et al., 2001.
!However, the code follows as reference the code found in https://github.com/ARVE-Research/LPJ-LMfire
!In this mode, each PLS is represented by a population of average individuals.
!All individuals of a PLS present the same features and the population is defined by
!the density of individuals per m2.
!--------------------------------------------------------------------------------------------------

!initially it will be done without time

program FPC
    !parameters
    integer, parameter:: npls = 20
    real :: k_allom1 = 100. !allometric constant (Table 3; Sitch et al., 2003)
    real :: krp = 1.6 !allometric constant (Table 3; Sitch et al., 2003)

    !average individual properties of an specific PLS
    real, dimension(npls) :: lai_ind !Leaf Area Index (m2/m2)
    real, dimension(npls) :: crown_area_ind !Tree crown area (m2) (Sitch et al., 2003)
    real, dimension(npls) :: fpc_ind !Foliage projective cover for each average individual of a PLS (Stich et al., 2003)
    real, dimension(npls) :: cl_2_ind !amount of leaf carbon in an average individual !gc/average individual 
    real, dimension(npls) :: cw_2_ind !amount of wood carbon in an average individual !gc/average individual 
    real, dimension(npls) :: diam_ind !amount of carbon in an average individual !gc/average individual 
 

    !PLS properties
    real, dimension(npls) :: fpc_grid !Foliage projective cover of a PLS (Stich et al., 2003)
    real, dimension(npls) :: pls_excess
    real, dimension(npls) :: nind_kill
    real, dimension(npls) :: nind_2
    
    !grid propertie
    real :: fpc_tree_total !total FPC considering trees
    real :: fpc_tree_max_perc = 0.95 !maximum of occupation = 95% of simulation area
    real :: fpc_tree_max !maximum of occupation in m2
    real :: cell_area = 15 !m2 only for testing purpose
    real :: excess !excess of area occupied by the sum of PLS's FPC in comparison to the fpc_tree_max
    

    !variable traits
    real, dimension(npls) :: spec_leaf !m2/gC
    real, dimension(npls) :: dwood !wood density (g/cm-3) *Fearnside, 1997 - aleatory choices

    !variables with initial values to start the logic
    real, dimension(npls) :: cl_1 !carbon leaf to start the logic !KgC/m2 
    real, dimension(npls) :: nind_1 !initial number of individuals per PLS !ind/m2
    real, dimension(npls) :: cw_1 !carbon wood to start the logic !KgC/m2 

    cl_1 = (/.7,1.,0.3,1.6,1.10,1.8,0.3,0.2,0.8,0.84,0.25,1.,0.2,1.7,0.4,.6,.5,.8,0.3,1.8/)

    nind_1 = (/1.,2.,5.,3.3, 1.3, 7., 2.8, 3.,4.5,1.7,3.6,9.,4.,2.45,5.27,4.6,8.2,9.29,3.,4.8/)
    
    cw_1 = (/30.,40.,34.,28.3,20.2,26.7,27.5,19.5,20.,28.6,24.3,19.3,26.8,22.,25.,22.,15.,22.6,10.7,21.4/)

    dwood = (/0.24,0.53,0.39,0.32,0.31,0.44,0.66,0.42,0.74,0.39,0.82,0.40,0.26,0.79,0.39,0.52,0.41,0.44,0.86,0.42/)
    
    spec_leaf = (/0.0153,0.0101,0.0107,0.0112,0.012,0.0141,0.0137,0.0115,0.0122,0.010,0.012,0.011,&
    &0.013,0.014,0.0112,0.012,0.0141,0.0137,0.0115,0.0122/)

    !!!------------------------------------------------------
        !Transforms from kgC to gC (as in LPJ)

    do j = 1, npls 

       cl_1(j) = cl_1(j)*1000.
        
       cw_1(j) = cw_1(j)*1000.

    enddo

    !!!-------------------------------------------------------

    !transforming the carbon content from gC/m2 to gc/average individual 
    !(the carbon divided by nind gives the individual carbon, as in LPJ)

    do j = 1, npls 

        cl_2_ind(j) = cl_1(j)/nind_1(j)

        cw_2_ind(j) = cw_1(j)/nind_1(j)
        
    enddo

!!!!---------------------------------------------------------
     !Structuring PLSs [crown area and diameter]
    
    do j = 1, npls

        diam_ind(j) = ((4*(cw_2_ind(j)))/((dwood(j)*1000000.)*3.14*36))**(1/(2+0.22)) !nessa equação dwood deve estar em *g/m3*
      
        crown_area_ind(j) = k_allom1*(diam_ind(j)**krp)

         
    enddo 

    
    !calculates FPC
    do j = 1, npls
        if (crown_area_ind(j) .gt. 0.) then 

            lai_ind(j) = (cl_2_ind(j)*spec_leaf(j))/crown_area_ind(j) 

            fpc_ind(j)  = 1. - exp(-0.5 * lai_ind(j))

            fpc_grid(j) = fpc_ind(j) * nind_1(j) * crown_area_ind(j)
            
        else

            lai_ind(j)  = 0.

            fpc_ind(j)  = 0.

            fpc_grid(j) = 0.

        endif

        
    enddo
    
    fpc_tree_total = sum(fpc_grid(:))
   
!-------------------------------
!light competition, woody plants

pls_excess(:) = 0.

fpc_tree_max = fpc_tree_max_perc * cell_area

if(fpc_tree_total.gt.fpc_tree_max) then !reduce tree cover
 
    excess = fpc_tree_total - fpc_tree_max

    do j = 1, npls

        if (fpc_grid(j) .gt. 0.) then

            !this formulation ensures equal competition (precludes total dominance by one PFT)

            pls_excess(j) = min(fpc_grid(j), excess *  fpc_grid(j) / fpc_tree_total)
        
        else
        
            pls_excess(j) = 0.
        
        endif

    end do

    do j = 1,npls
        
        if (pls_excess(j) .gt. 0.) then
      
            !Reduce individual density (and thereby gridcell-level biomass) so that total tree FPC reduced to 'fpc_tree_max'
      
            nind_kill(j) = nind_1(j) * (pls_excess(j) / fpc_grid(j))
            print*, nind_1(j), (pls_excess(j) / fpc_grid(j)), pls_excess(j), fpc_grid(j), nind_kill(j)
           
      
            nind_2(j) = nind_1(j) - nind_kill(j)
            ! print*,cl_2_ind(j)- (nind_kill(j)*cl_2_ind(j))
            

        endif

    enddo
else 
    
    do j = 1,npls
        
    enddo

endif


    

    
end program FPC