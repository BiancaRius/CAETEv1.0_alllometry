module constants
    use ISO_FORTRAN_ENV, only: REAL32, REAL64, REAL128
    implicit none

    !!!!PRECISAMOS REVER AS UNIDADES.

    !==============================!
    != Constants
    !==============================!
    integer, parameter:: npls = 20
    real(REAL64), parameter :: klatosa = 6000.0
    real(REAL64), parameter :: dw = 200 !(*1000) converte de g/cm3 para kg/m3
    real(REAL64), parameter :: ltor = 0.77302587552347657
    real(REAL64), parameter :: k_allom1 = 100.0 !allometric constant (Table 3; Sitch et al., 2003)
    real(REAL64), parameter :: k_allom2 = 36.0
    real(REAL64), parameter :: k_allom3 = 0.22
    real(REAL64), parameter :: spec_leaf = 15.365607091853349 
    real(REAL64), parameter :: bminc = 2.33 !3.5 (valor medio de NPP)/12 (p/ calcular NPP mensal) * 8(densidade de individuos )
    real(REAL64), parameter :: dens = 1.
    real(REAL64), parameter :: tol = 0.0000001
    real(REAL64), parameter :: pi = 3.1415926536
    real(REAL64), parameter :: krp = 1.6 !allometric constant (Table 3; Sitch et al., 2003)
    real(REAL64), parameter :: turnover_rate_sap = 0.004166667!!value/month!!0.05!turnover rate for sapwood(Table 1; Sitch et al., 2003)
    real(REAL64), parameter :: turnover_rate_leaf = 0.041666667!value/month!!!0.5!turnover rate for leaves(Table 1; Sitch et al., 2003)
    real(REAL64), parameter :: turnover_rate_root = 0.041666667!value/month!!!0.5!turnover rate for roots(Table 1; Sitch et al., 2003)
   

    !==============================!

end module constants