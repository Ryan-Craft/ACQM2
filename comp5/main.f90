program main

    implicit none
    
    ! Variables
    real*8 :: alpha, l, dr, rmax, dk, kmax, &
              nr, nk
    integer :: i,j,ier, N
    integer*8 :: p

    ! 1D allocatable arrays
    real*8, dimension(:), allocatable :: rgrid, kgrid, rweights, A1, A2

    ! 2D allocatable arrays
    real*8, allocatable :: &
              Vdirect(:,:), &
              Exchange(:,:), &
              wf(:,:)                
 

    ! READ IN ALL THE VARIABLES 
    open(1, file="LaguerreParams.txt", action="read")
    read(1,*) alpha, N, l, dr, rmax, dk, kmax
    close(1) 
    Print *, "alpha, N, l, dr, rmax, dk, kmax"
    Print *, alpha, N, l, dr, rmax, dk, kmax

    ! determine nr and nk
    nr = rmax/dr 
    nk = kmax/dk

    ! MAKE MEMORY ALLOCATIONS
    allocate(rgrid(nr))
    allocate(kgrid(nk))
    allocate(rweights(nr))
    allocate(Vdirect(nk,nk))
    allocate(Exchange(nk,nk))
    allocate(wf(nr,N))  




    ! DEFINE R AND K GRIDS, MAKE SURE TO CONVERT TO eV, there is no 0 energy element
    do i=1,nr
       rgrid(i) = i*dr  
    end do
    
    do i=1,nk
       kgrid(i) = i*dk/27.21136 ! dividing by 27 here to immediately convert to hartrees
    end do

    ! DEFINE RWEIGHTS
    rweights = 0.0d0
    do i=1,nr
        rweights(i) = 4.0d0 - 2.0d0*mod(i+1,2)
    end do   

    ! PRINT TO SCREEN FOR VALIDATION OF GRIDS AND WEIGHTS
    Print *, "In order of appearance: rgrid, kgrid and rweights first 5 elements" 
    Print *, rgrid(1:5)
    Print *, kgrid(1:5)
    Print *, rweights(1:5)
 

    ! CREATE RADIAL WAVEFUNCTION OF HYDROGEN
    
    

    ! DEALLOCATE MEMORY 



end program























