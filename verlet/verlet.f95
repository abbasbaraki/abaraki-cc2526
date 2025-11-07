PROGRAM verlet
  IMPLICIT NONE ! inline comment
  INTEGER, PARAMETER :: wp =SELECTED_REAL_KIND (p=13, r=300)

  ! PARAMETERS & STATE
  REAL (KIND=wp) :: tau
  REAL (KIND=wp) :: m
  REAL (KIND=wp) :: x, y, z, xn, yn, zn, v_x, v_y, v_z, v_xn, v_yn, v_zn, f_x, f_y, f_z, f_xn, f_yn, f_zn
  INTEGER :: k, i

  tau = 0.2_wp
  m = 1_wp
  x = 0_wp
  y = 0_wp
  z = 0_wp
  v_x = 0_wp
  v_y = 0_wp
  v_z = 0_wp
  f_x = 0_wp
  f_y = 0.1_wp
  f_z = 0_wp
  f_xn = 0_wp
  f_yn = 0.1_wp
  f_zn = 0_wp
  i = 100

  DO k = 1, i
    xn = x + tau * v_x + (tau * tau) * f_x / (2 * m)
    yn = y + tau * v_y + (tau * tau) * f_y / (2 * m)
    zn = z + tau * v_z + (tau * tau) * f_z / (2 * m)
    v_xn = v_x + tau / (2 * m) * (f_x + f_xn)
    v_yn = v_y + tau / (2 * m) * (f_y + f_yn)
    v_zn = v_z + tau / (2 * m) * (f_z + f_zn)
    v_x = v_xn
    v_y = v_yn
    v_z = v_zn
    x = xn
    y = yn
    z = zn
    f_x = f_xn
    f_y = f_yn
    f_z = f_zn

    PRINT *, "The positions of the particle at the", k, "-th iteration at time =", tau * k, ":"
    PRINT *, xn
    PRINT *, yn
    PRINT *, zn
    PRINT *, v_xn
    PRINT *, v_yn
    PRINT *, v_zn

  END DO

END PROGRAM VERLET