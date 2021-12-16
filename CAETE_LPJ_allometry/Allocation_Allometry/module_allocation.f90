module allocation
    use ISO_FORTRAN_ENV, only: REAL32, REAL64, REAL128
    use constants
    implicit none

    ! real(REAL64) :: H = 0.0
    ! real(REAL64) :: L = 0.0
    ! real(REAL64) :: R = 0.0
    ! real(REAL64) :: S = 0.0

    !the 1st thing to do is to calculate the potential npp, verify the nutrients availability
    !!!!!!

    
    
    
    
    ! Variables with generic values for testing the logic code


    ! real(REAL64) :: H = 18.91909828977032 !in kgC/m2 --> /densindiv (numindv/m2) !HEARTWOOD - SOMENTE PARA TESTES (Valores: Cod. Phillip)
    ! real(REAL64) :: L = 1.2279169651518438 !LEAF BIOMASS - SOMENTE PARA TESTES (Valores: Cod. Philipe)
    ! real(REAL64) :: S = 9.790591253578555 !SAPWOOD - SOMENTE PARA TESTES (Valores: Cod. Philipe)
    ! real(REAL64) :: R = 0.88026193814051623 !ROOT BIOMASS - SOMENTE PARA TESTES (Valores: Cod. Philipe)

    real(REAL64) :: H = 108.91909828977032 !in kgC/m2 --> /densindiv (numindv/m2) !HEARTWOOD - SOMENTE PARA TESTES (Valores: Cod. Phillip)
    real(REAL64) :: L = 1.2279169651518438 !LEAF BIOMASS - SOMENTE PARA TESTES (Valores: Cod. Philipe)
    real(REAL64) :: S = 29.790591253578555 !SAPWOOD - SOMENTE PARA TESTES (Valores: Cod. Philipe)
    real(REAL64) :: R = 0.88026193814051623 !ROOT BIOMASS - SOMENTE PARA TESTES (Valores: Cod. Philipe)
    
    
    real(REAL64) :: stem = 0.0   !stem (heartwood + sapwood) pool update 
    real(REAL64) :: L_updt = 0.0 !Leaf pool update (L year before + delta_heartwood)
    real(REAL64) :: S_updt = 0.0 !Sapwood pool update (S year before + delta_heartwood)
    real(REAL64) :: R_updt = 0.0 !Root pool update (R year before + delta_heartwood)
    real(REAL64) :: H_updt = 0.0 !heartwood pool update (R year before + delta_heartwood)
    real(REAL64) :: turnover_sap = 0.0 !Sapwood turnover to heartwood (yr-1; Sitch et al 2003)
    real(REAL64) :: turnover_leaf = 0.0 !Leaf turnover (yr-1; Sitch et al 2003) 
    real(REAL64) :: turnover_root = 0.0 !Root turnover (yr-1; Sitch et al 2003)
    real(REAL64) :: funcs_calc_tau1 
    
    contains


    
    !==============================!
	!= Subrotines
	!==============================!
    ! Just to test comparisons with "tol" value
    ! It can be deleted later
    ! subroutine show_consts()
    !     implicit none
    !     real(REAL128) :: x = 0.0000001
        
    !     print*, 'show_consts()'
    !     print '(F11.9)', tol

    !     if(x .eq. tol) then
    !         print*, 'eq'
    !     elseif(x .gt. tol) then
    !         print*, 'gt'
    !     else
    !         print*, 'lt'
    !     end if

    ! end subroutine show_consts

    ! Use the bisection method to solve the leaf mass increment



    subroutine leaf_increment(delta_leaf)
        real(REAL64) :: delta_leaf

        delta_leaf = bisection_method(0.0, 10.0)
        
        return
    end subroutine leaf_increment

	! Once we have the leaf mass increment we can cant get 
    ! root mass increment based on the LTOR constant
    subroutine root_increment(delta_leaf, delta_root)
        real(REAL64) :: delta_leaf
        real(REAL64) :: delta_root
        
        delta_root = (delta_leaf + L) / ltor - R
        
        return
    end subroutine root_increment

    ! Finally using the cmass_increment mass conservation we can calculate sapwood increment
    subroutine sapwood_increment(delta_leaf, delta_root, delta_sapwood)
        real(REAL64) :: delta_leaf
        real(REAL64) :: delta_root
        real(REAL64) :: delta_sapwood
        
        delta_sapwood = bminc - delta_leaf - delta_root
        
        return
    end subroutine sapwood_increment

    !Updating carbon pools in each compartments with the deltas and turnover (this part go to allometry dynamic)

    subroutine updating_pool_leaf(delta_leaf,L,L_updt,turnover_leaf,turnover_rate_leaf)
        real(REAL64) :: delta_leaf
        real(REAL64) :: L 
        real(REAL64) :: L_updt
        real(REAL64) :: turnover_leaf 
        real(REAL64) :: turnover_rate_leaf
        
        turnover_leaf = L*turnover_rate_leaf

        L_updt = delta_leaf  + L - turnover_leaf

        return
    end subroutine updating_pool_leaf

  
    subroutine updating_pool_root(delta_root,R,R_updt,turnover_root,turnover_rate_root)
        real(REAL64) :: delta_root
        real(REAL64) :: R
        real(REAL64) :: R_updt
        real(REAL64) :: turnover_root 
        real(REAL64) :: turnover_rate_root

        turnover_root = R*turnover_rate_root

        R_updt = delta_root + R - turnover_root

        return
    end subroutine updating_pool_root

    subroutine updating_pool_sapwood(delta_sapwood,S,S_updt,turnover_sap,turnover_rate_sap)
        real(REAL64) :: delta_sapwood
        real(REAL64) :: S
        real(REAL64) :: S_updt
        real(REAL64) :: turnover_sap 
        real(REAL64) :: turnover_rate_sap 

        turnover_sap = S*turnover_rate_sap

        S_updt = delta_sapwood + S - turnover_sap

        return
    end subroutine updating_pool_sapwood

    subroutine updating_pool_heartwood(H,turnover_sap,H_updt)

        real(REAL64) :: H
        real(REAL64) :: turnover_sap
        real(REAL64) :: H_updt
       
        H_updt = H + turnover_sap

        return
    end subroutine updating_pool_heartwood

    subroutine updating_pool_stem(S_updt, H_updt, stem)
        real(REAL64) :: S_updt
        real(REAL64) :: H_updt
        real(REAL64) :: stem

        stem = S_updt + H_updt

        return
    end subroutine updating_pool_stem

    subroutine initial_c_pool (cl1,cw1,cr1,cleaf_init,cheart_init,csap_init,croot_init)
        real,dimension(npls),intent(in):: cl1,cw1,cr1
    
        real,dimension(npls),intent(out):: cleaf_init,cheart_init,csap_init,croot_init
      
        cleaf_init = cl1
        cheart_init = cw1*0.95
        csap_init = cw1*0.05
        croot_init = cr1

    end subroutine initial_c_pool


	!==============================!
	!= Functions
	!==============================!
    function bisection_method(a, b) result(midpoint)
        implicit none
        real(REAL32) :: a, b
        real(REAL64) :: aux_a, aux_b
        real(REAL64) :: midpoint
        
        aux_a = a
        aux_b = b

        if((f(aux_a) * f(aux_b)) .gt. 0) then
            midpoint = -2.0
            return
        endif
        
        do while((aux_b - aux_a) / 2.0 .gt. tol)
            midpoint = (aux_a + aux_b) / 2
            
            if(f(midpoint) .eq. 0.0) then
                exit            
            elseif(f(aux_a) * f(midpoint) .lt. 0) then
                aux_b = midpoint
            else
                aux_a = midpoint
            endif
        end do
    end function bisection_method

    function f(x) result(searched_x)
        implicit none
        real(REAL64) :: x
        real(REAL64) :: searched_x
        
        searched_x = & 
            calc_tau1() * &
            (sapwood() - x - x / ltor + H) - &
            ( &
                (sapwood() - x - x / ltor) / &
                (L + x) * calc_tau3() &
            ) ** calc_tau2()
    end function f

    function calc_tau1() result(tau1)
        implicit none
        real(REAL64) :: tau1
        
        tau1 = k_allom2 ** (2.0 / k_allom3) * 4.0 / 3.14159 / dw

    end function calc_tau1

    function calc_tau2() result(tau2)
        implicit none
        real(REAL64) :: tau2 
        
        tau2 = 1.0 + 2.0 / k_allom3

    end function calc_tau2

    function calc_tau3() result(tau3)
        implicit none
        real(REAL64) :: tau3
        
        tau3 = klatosa / dw / spec_leaf
        
    end function calc_tau3

    ! function sapwood (csap_init,cleaf_init,croot_init,bminc) result (SS)
    function sapwood () result (SS)
         implicit none
         real(REAL64):: SS
        !  real, dimension(npls) :: SS
        !  real, dimension(npls) :: csap_init
        !  real, dimension(npls) :: cleaf_init
        !  real, dimension(npls) :: croot_init
        !  real :: bminc
        
         SS = S + bminc - L/ ltor + R
        
        !  SS = csap_init + bminc - cleaf_init / ltor + croot_init
     end function sapwood

    

     

end module allocation
