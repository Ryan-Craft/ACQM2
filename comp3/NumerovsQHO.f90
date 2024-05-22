
subroutine NumerovForwards(psi_L, V, rgrid, nr, E)

    implicit none
    
    ! initialise local and inbound variables
    integer :: i, j
    real *8, intent(in) :: nr

    ! array init

    real* 8, dimension(nr), intent(in) :: rgrid
    real* 8, dimension(nr), intent(in) :: V
    real* 8, dimension(nr), intent(in) :: psi_L
    real*8, dimension(nr) :: g


    g = 2*(V-E)
    
    do i=1,nr
        Print *, g(i)
    end do
    
    do i=1,nr
        

    end do
    
    

    !generate local variables








end NumerovForwards








program NumerovsQHO
    implicit none

    integer*8 :: i, j, nr
    real*8 :: dr, rmax, n, E, E_min, E_max

    !ARRAY INITIALISATION
    real*8, dimension(:), allocatable :: rgrid, V
    real*8, dimension(:), allocatable :: psi_L, psi_R, psi    



    
    !!Read in params

    open(unit=1, file="NumerovsQHOParams.txt", action="read")
    read(1,*) dr, rmax, n
    Print *, dr, rmax, n 
    close(1)




    !create nr and allocate rgrid
    ! allocate wavefunction vectors

    nr = (2*rmax)/dr + 1

    allocate(rgrid(nr))
    allocate(V(nr))
    allocate(psi_L(nr), psi_R(nr), psi(nr))
    
    !set wavefunctions to initial zero
    psi_L=0.0d0 
    psi_R=0.0d0
    psi=0.0d0


    ! generate rgrid and use rgrid to generate V(r)

    Print *, "Generating rgrid from -rmax to rmax"
    do i =1, nr
        rgrid(i) = -rmax + (i-1.0d0)*dr
    end do

    Print *, "Generating V array from -rmax to rmax"
    do i=1, nr
        V(i) = 0.5*rgrid(i)**2
    end do

    !WRITE potentia out an output file:

    open(unit=1, file="potentialout.txt", action="write")
    do i=1,nr
        write(1, '(*(f12.8))') rgrid(i), V(i)
    end do

    
    !Initial guess E, based on n
    ! I find this a bit weird because we know the energies of the QHO exactly; En = (n+1/2)*w*h_bar
    ! so considering that this is evenly spaced, my first energy is going to be around 0.5. 
    ! my guesses for En; E_min = n and E_max = n+0.75 so that we can get the iteration scheme to work

    E_min = n
    E_max = n+0.75
    E = (E_min + E_max)/2 
    





















end program
