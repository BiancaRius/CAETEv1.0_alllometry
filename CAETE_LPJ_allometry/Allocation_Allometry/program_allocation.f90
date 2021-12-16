program caete
    use ISO_FORTRAN_ENV, only: REAL32, REAL64, REAL128
    use allocation
    use constants
    implicit none

    real(REAL64) :: delta_leaf
    real(REAL64) :: delta_root
    real(REAL64) :: delta_sapwood
    
    
    integer :: j, k
    real, dimension(npls) :: dwood !wood density (g/cm-3) *Fearnside, 1997 - aleatory choices
    real, dimension(npls) :: cw1,cl1,cr1 !KgC/m2 (Cheart + Csap)
    real, dimension(npls) :: cleaf_init, cheart_init, csap_init, croot_init, SS
    real, dimension(npls):: diam !Tree diameter in m. (Smith et al., 2001 - Supplementary)
    real, dimension(npls) :: nind !number of individuals per PLS (Smith, 2001, thesis)

    dwood =(/0.74,0.73,0.59,0.52,0.41,0.44,0.86,0.42,0.64,0.69,0.92,0.60,0.36,0.99,0.59,0.52,0.41,0.44,0.86,0.42/) !atenção para a unidade
    cw1 =(/7.,12.,7.2,8.3,8.8,9.7,7.5,11.5,10.,8.6,7.3,10.3,6.8,9.9,5.3,9.2,15.,12.6,10.7,11.4/)
    cl1=(/2.15,3.,1.18,1.6,1.5,1.8,0.3,2.,0.8,.84,0.25,1.,0.2,1.7,1.18,1.6,1.5,1.8,0.3,2./)
    cr1=(/0.63,0.8,0.9,1.5,1.3,0.9,0.4,1.0,0.56,0.87,0.33,0.97,0.31,0.55,0.2,0.8,0.4,0.66,0.23,1.5/)

   
    
    call initial_c_pool(cl1,cw1,cr1,cleaf_init,cheart_init,csap_init,croot_init)
    ! print*, cleaf_init,cheart_init,csap_init,croot_init

    ! SS = sapwood(csap_init,cleaf_init,croot_init,bminc)
    ! ! print*,'function SS', SS

    call leaf_increment(delta_leaf)
    print*, 'CARBON LEAF INCREMENT =', delta_leaf

    call root_increment(delta_leaf, delta_root)
    print*, 'CARBON ROOT INCREMENT =', delta_root
    
    call sapwood_increment(delta_leaf, delta_root, delta_sapwood)
    print*, 'CARBON SAPWOOD INCREMENT =', delta_sapwood

    call updating_pool_leaf(delta_leaf,L,L_updt,turnover_leaf,turnover_rate_leaf)
    print*, 'LEAF POOL UPDATED =', L_updt

    call updating_pool_root(delta_root,R,R_updt,turnover_root,turnover_rate_root)
    print*, 'ROOT POOL UPDATED =', R_updt, R, delta_root

    call updating_pool_sapwood(delta_sapwood,S,S_updt,turnover_sap,turnover_rate_sap)
    print*, 'SAPWOOD POOL UPDATED =', S_updt, S, delta_sapwood

    call updating_pool_heartwood(H,turnover_sap,H_updt)
    print*, 'HEARTWOOD POOL UPDATED =', H_updt, H, turnover_sap

    call updating_pool_stem(S_updt, H_updt, stem)
    print*, 'STEM POOL UPDATED =', stem !H, S_updt

   
    
   



end program caete