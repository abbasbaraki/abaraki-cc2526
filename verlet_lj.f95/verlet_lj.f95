PROGRAM verlet_LJ
    
    IMPLICIT NONE                                                               ! inlilne comment
    INTEGER, PARAMETER :: wp =SELECTED_REAL_KIND (p=13, r=300)

    ! PARAMETERS & SATATE
    INTEGER :: i, k, number                                                     ! number of atoms, iteration
    REAL (KIND=wp) :: tau                                                       ! each time interval between two iterations
    REAL (KIND=wp) :: sigma, epsilon                                            ! obviously sigma and epsilon
    REAL (KIND=wp) :: m, r                                                      ! mass and distance between atoms
    REAL (KIND=wp) :: v_prime                                                   ! derivative of Lennard-Jones potential with respect to distance "r"
    REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: dx                             ! displacement vector [dx, dy, dz] between the two atoms
    REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: xmat, x_nmat, vmat, v_nmat, fmat, f_nmat     
    ! matrix of position, velocity and force in different arrays

    ALLOCATE ( dx (3))                                                          ! "dx" is a length-3 vector (x,y,z components)
    ALLOCATE ( xmat (2,3))                                                      ! "x" positions of 2 atoms in 3D (row = atom, col = x,y,z)
    ALLOCATE ( x_nmat (2,3))                                                    ! "x_n" new positions after one time step
    ALLOCATE ( vmat (2,3))                                                      ! "v" velocities of 2 atoms (x,y,z)
    ALLOCATE ( v_nmat(2,3))                                                     ! "v_n" new velocities after one time step
    ALLOCATE ( fmat (2,3))                                                      ! "f" current forces on each atom (x,y,z)
    ALLOCATE (f_nmat (2,3))                                                     ! "f_n" new forces (after update of positions)

    xmat(1,:) = (/ 0.0_wp, 0.0_wp, 8.0_wp /)                                    ! atom 1 at (0,0,8) or {xmat(1,1) = 0.0_wp; xmat(1,2) = 0.0_wp; xmat(1,3) = 8.0_wp}
    xmat(2,:) = (/ 0.0_wp, 0.0_wp, 0.0_wp /)                                    ! given by question atom 2 at (0,0,0)
    
    vmat(1,:) = (/ 0.0_wp, 0.0_wp, 0.0_wp /)                                    ! initial velocities set to 0
    vmat(2,:) = (/ 0.0_wp, 0.0_wp, 0.0_wp /)                                    ! given by question
   
    epsilon = 0.000112991_wp                                                    ! given by question
    sigma = 5.2186_wp                                                           ! given by question
    m = 20.1797_wp * 1822.888486_wp                                             ! given by question
    !number = 2
    tau = 1.0_wp
    i = 6000                                                                    ! just to chek
     
    OPEN (UNIT=12, FILE="Neon.xyz", STATUS="replace", ACTION="write")
    
    DO k = 1, i

        dx = xmat(2,:) - xmat(1,:)                                                        ! displacement vector r_2 - r_1 = [dx,dy,dz]

        r = sqrt(sum( dx ** 2 ))                                                          ! distance "r" between the two atoms

        v_prime = 4_wp*epsilon * ( ( -12_wp * (sigma**12) / (r**13) ) + ( 6_wp * (sigma**6) / (r**7) ) )    ! v_prime = dV/dr for Lennard-Jones

        fmat(1,:) = -v_prime * ( (dx(:)) / r )                                            ! new force on first atom in 3 directions
        fmat(2,:) = -fmat(1,:)                                                            ! new force on second atom due to the NEWTON's 3rd law

        x_nmat(1,:) = xmat(1,:) + tau * vmat(1,:) + (tau**2) * fmat(1,:) / (2_wp * m)     ! new positions of first atom in 3 directions
        x_nmat(2,:) = xmat(2,:) + tau * vmat(2,:) + (tau**2) * fmat(2,:) / (2_wp * m)     ! new positions of second atom in 3 directions

        dx = x_nmat(2,:) - x_nmat(1,:)                                                    ! "dx" in updated position
        r = sqrt(sum( dx ** 2 ))                                                          ! "r" in updated position
        v_prime = 4_wp * epsilon * ( ( -12_wp * (sigma**12) / (r**13) ) + ( 6_wp * (sigma**6) / (r**7) ) )

        f_nmat(1,:) = -v_prime * ( (dx(:)) / r )                                          ! new force on first atom in 3 directions
        f_nmat(2,:) = -f_nmat(1,:)                                                        ! new force on second atom due to the NEWTON's 3rd law

        v_nmat(1,:) = vmat(1,:) + (tau / (2.0_wp*m)) * (fmat(1,:) + f_nmat(1,:))            ! new velocitiy of first atom in "x,y,z"
        v_nmat(2,:) = vmat(2,:) + (tau / (2.0_wp*m)) * (fmat(2,:) + f_nmat(2,:))            ! new velocitiy of second atom

        xmat(:,:) = x_nmat(:,:)                                                           ! new positions become current
        vmat(:,:) = v_nmat(:,:)                                                           ! new velocities become current
        fmat(:,:) = f_nmat(:,:)                                                           ! new forces become current


        WRITE (UNIT=12, FMT=*) xmat(1,:) 
        WRITE (UNIT=12, FMT=*) xmat(2,:)
        WRITE (UNIT=12, FMT=*) f_nmat(1,:)
        WRITE (UNIT=12, FMT=*) f_nmat(2,:)
        WRITE (UNIT=12, FMT=*) v_nmat(1,:)
        WRITE (UNIT=12, FMT=*) v_nmat(2,:)

        PRINT *, "The positions of the particle at the", k, "-th iteration at time =", tau * k, ":"
        PRINT *, x_nmat(1,:)
        PRINT *, x_nmat(2,:)
        PRINT *, "---------------X------------------"
        PRINT *, xmat(1,:)
        PRINT *, xmat(2,:)
        PRINT *, "---------------V_N-----------------"
        PRINT *, v_nmat(1,:)
        PRINT *, v_nmat(2,:)
        PRINT *, "----------------F-----------------"
        PRINT *, fmat(1,:)
        PRINT *, fmat(2,:)
        PRINT *, "----------------F_N-----------------"
        PRINT *, f_nmat(1,:)
        PRINT *, f_nmat(2,:)
        

    END DO
    CLOSE (12)
END PROGRAM VERLET_LJ