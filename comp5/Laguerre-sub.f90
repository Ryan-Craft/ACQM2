

program LaguerreSub
        implicit none

        ! need int N, real alpha, real l

        real :: alpha, l
        real :: dr, rmax
        integer :: N, nr
        integer :: i,j,k
        real, dimension(:), allocatable :: rgrid
        real, dimension (:,:), allocatable :: basis


        

        open(unit=1, file="LaguerreParams.txt", action="read")

        read(1,*) alpha, N, l, dr, rmax
       
        Print *, alpha, N, l, dr, rmax


        !calculate rgrid params
        nr = rmax/dr + 1
        Print *, nr

        !allocate rgrid and basis
        allocate(rgrid(nr))
        allocate(basis(nr,N))

        
        do i = 1, nr
                rgrid(i) = (i-1)*dr
        end do     

 
        basis(:,1) = (2.0d0*alpha*rgrid(:))**(l+1) *exp(-alpha*rgrid(:))
        basis(:,2) = 2.0d0*(l+1-alpha*rgrid(:)) * (2.0d0*alpha*rgrid(:))**(l+1) *exp(-alpha*rgrid(:))      

        !generate N laguerre basis using recurrence relation
        do i = 3, N
                basis(:,i) = (2*(i-1+l-alpha*rgrid(:))*basis(:,i-1) - (i+2*l-1)*basis(:,i-2) )  / (i-1)

        end do
 
        Print *, "BASIS 1: \n" 
        do j = 1, nr
                Print *, rgrid(j), basis(j,1), basis(j,2), basis(j,3), basis(j,4)
        end do


        !generate N laguerre basis using recurrence relation






end program
