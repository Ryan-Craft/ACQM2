
subroutine LaguerreSub(alpha, l, nr, N, rgrid, basis)
        implicit none

        ! initialise alpha, l, dr, rmax, N and others
        integer :: i,j
        integer, intent(in) :: nr
        integer, INTENT(IN) :: N
        real*8, INTENT(IN) :: alpha, l
        real*8, dimension(nr), INTENT(IN) :: rgrid
        real*8, dimension(nr,N) :: basis
        
        basis(:,1) = (2.0d0*alpha*rgrid(:))**(l+1) *exp(-alpha*rgrid(:))
        basis(:,2) = 2.0d0*(l+1-alpha*rgrid(:)) * (2.0d0*alpha*rgrid(:))**(l+1) *exp(-alpha*rgrid(:))

        !generate N laguerre basis using recurrence relation
        do i = 3, N
                basis(:,i) = (2*(i-1+l-alpha*rgrid(:))*basis(:,i-1) - (i+2*l-1)*basis(:,i-2) )  / (i-1)
        end do
        
        return
end subroutine LaguerreSub


program main
         implicit none


         real*8 :: normalise
         real*8 :: alpha, l
         real*8 :: dr, rmax
         integer*8 :: p
         integer :: N, nr, ier
         integer :: i,j
         real*8, dimension(:), allocatable :: rgrid
         real*8, dimension(:,:), allocatable :: basis
         
         !add in overlap and hamiltonian
         real*8, dimension(:,:), allocatable :: H
         real*8, dimension(:,:), allocatable :: B
         real*8, dimension(:,:), allocatable :: K
         ! energies and expansion coefficients
         real*8, dimension(:,:), allocatable :: w
         real*8, dimension(:,:), allocatable :: z
         real*8, dimension(:,:), allocatable :: V
         
         ! create array for wavefunctions
         real,dimension(:,:), allocatable :: wf

         ! fortran constant
         real, parameter :: pi = 4.0d0*atan(1.0d0)
         real, parameter :: radtodeg = 180.0d0/pi

         ! additions for assignment 5, energy grid
         real*8, dimension(:), allocatable :: kgrid, rweights
         integer :: nk, ii, jj
         real*8 :: kmax, dk, deltafi, H_init, H_final
         real*8, dimension(:), allocatable :: A1
         real*8, dimension(:), allocatable :: A2
         real*8, dimension(:), allocatable :: f
         real*8, dimension(:), allocatable :: g
         real*8, dimension(:,:), allocatable :: Vdirect
         real*8, dimension(:), allocatable :: foverlap
         real*8, dimension(:), allocatable :: ioverlap
         real*8, dimension(:,:), allocatable :: V1
         real*8, dimension(:,:), allocatable :: V2
         real*8, dimension(:,:), allocatable :: V12
         !open file location: hard coded for now but could become flexible
         !read stored values into relevent variables
         
         open(unit=1, file="LaguerreParams.txt", action="read")
         read(1,*) alpha, N, l, dr, rmax, dk, kmax, H_init, H_final
         Print *, alpha, N, l, dr, rmax, dk, kmax, H_init, H_final

         !calculate rgrid params
         nr = rmax/dr
         Print *, "nr ::", nr
        
         ! calculate kgrid params
         nk = kmax/dk +1
         Print *, "nk ::", nk



         ! based on options from file, allocate appropriate memory to rgrid and the basis array
         allocate(rgrid(nr))
         allocate(basis(nr,N))
         allocate(H(N,N))
         allocate(B(N,N))
         allocate(V(N,N))
         allocate(K(N,N))
         allocate(w(N,1)) 
         allocate(z(N,N))
         allocate(wf(nr,N))
         allocate(kgrid(nk))         
         allocate(rweights(nr))
         allocate(Vdirect(nk,nk))
         allocate(ioverlap(nk), foverlap(nk), V1(nk,nk), V2(nk,nk), V12(nk,nk))

!!! YOU HAVE CHANGED WHAT RANGE THE RGRID GOES OVER DO NOT FORGOR
         !allocate values to the rgrid and kgrid
         do i = 1, nr
                rgrid(i) = (i)*dr
         end do 
         
         do i =1,nk
                kgrid(i) = dk*(i-1)
         end do

         ! generate weights for r dimension integration process on an even grid
         rweights = 0.0d0
         do i =1,nr
             rweights(i) = 4.0d0 - 2.0d0*mod(i+1,2)
         end do
         rweights = rweights * dr / 3.0d0


         Print *, "rgrid and kgrid first 5 values"
         Print *, rgrid(1:5)
         Print *, kgrid(1:5)
         Print *, "rweights"
         Print *, rweights(1:5)


         !use recurrence relation to compute the basis functions
         CALL LaguerreSub(alpha, l, nr, N, rgrid, basis)
        
         !implement normalisation condition using simplified factorial
         do i = 1, N
                 p=1.0
                 do j = 0, 2*l
                         p = p*(i+2*l-j)
                 !       Print *, p, j
                 end do
                 normalise = sqrt(alpha /((i+l)* p))
                 !Print *, "Norm:: ", normalise
                 basis(:,i) = normalise*basis(:,i)  
         end do

 
 
         !write basis to file for plotting 
         
         open(1, file='basisout.txt', action='write')
         do i =1,nr
                         write(1, '(*(f12.8))'), rgrid(i), basis(i,:)
         end do
         close(1)
        
         !calculate overlap matrix
         B = 0.0d0
         do i =1, N-1
                 B(i,i) = 1.0d0
                 ! basis function matrix elements are real valued, thus dot product is the same as inner prouct here
                 B(i,i+1) = -0.5d0*sqrt(1.0d0-( (l*(l+1.0d0)) / ((i+l)*(i+l+1.0d0))))
                 B(i+1,i) = B(i,i+1)
         end do
         B(N,N) = 1.0d0 
         Print *, "B ARRAY ::" 
         do i=1,N
                 !Print *, B(i,:)
         end do
 
         !calculate V matrix:
         V = 0.0d0
         do i=1,N
                 V(i,i) = -(alpha/(i+l))
         end do

         !create K to check
         K = (-alpha**2/2.0) * B
         do i=1,N
                 K(i,i) = K(i,i) + alpha**2
         end do
         
         !compute H-matrix Elements
         H = (-alpha**2/2.0) * B
         do i =1,N
                 H(i,i) = H(i,i) + alpha**2 - (alpha/(i+l)) 
         end do
         
         Print *, "H MATRIX::"
         do i=1,N
                 !Print *, H(i,:)
         end do         
 
         Print *, "V MATRIX::"
         do i=1,N
                 !Print *, V(i,:)
         end do  
         
         Print *, "K MATRIX::"
         do i=1,N
                 !Print *, K(i,:)
         end do 
         
         

         CALL rsg(N,N,H,B,w,1,z,ier)


         !recover wavefunctions:
         wf = 0.0d0
         do i =1,N
                 do j = 1,N
                         wf(:,i) = z(j,i)*basis(:,j) + wf(:,i)
                         !Print *, wf(:,i)
                 end do
         end do
         Print *, "W::"
         !Print *, w(1,:)
         Print *, "Z::"
         !Print *, z(1,:)

         open(1, file='wfout.txt', action='write')
         do i =1,nr
                         write(1, '(*(f12.8))'), rgrid(i), wf(i,:)
         end do
         close(1)
         
         open(1, file='wout.txt',action='write', access='append')
         do i =1,N
                 write(1, '(*(f12.8))'), real(N), w(i,1)
         end do
 



         !!! Assignment 5 additional

         
         ! memory allocations

         allocate(A1(nr)) 
         allocate(A2(nr))
         allocate(f(nr)) 
         allocate(g(nr))

!        G is a function of the initial and final wavefunctions of H
         g = wf(:,H_init)*wf(:,H_final) 

! inner sums
         A1(1) = rweights(1)*g(1) 
         do i=2,nr
             A1(i) = A1(i-1) + rweights(i)*g(i)
         end do 

         A2(nr) = (1/rgrid(nr))*rweights(nr)*g(nr) 
         do i=nr-1,1,-1
             A2(i) = A2(i+1) + (1/rgrid(i))*rweights(i)*g(i)
         end do

         Print *, "A1 and A2"
         Print *, A1(1:5)
         Print *, A2(1:5)

         Vdirect=0.0d0
         

        ! generating direct V matrix
         do i=1,nk
             do j=1,nk
                 f = sin(kgrid(i)*rgrid(:)) * sin(kgrid(j)*rgrid(:))
                 if(H_init==H_final) then
                     deltafi = 1.0d0
                 else
                     deltafi = 0.0d0
                 end if                
  
                 do ii=1,nr
                     Vdirect(j,i) = Vdirect(j,i) + rweights(ii)*f(ii)*((1/rgrid(ii) * A1(ii) + A2(ii)) - deltafi/rgrid(ii))
                 end do
             end do
         end do
         Vdirect = Vdirect * (2.0d0/pi) 
 
         Print *, "f and g"
         Print *, f(1:5)
         Print *, g(1:5)
 
         open(1, file="onshellVdirect.txt", action="write")
         do i=2,nk
             write(1, *) kgrid(i), ((-0.25)*(kgrid(i)**2/(kgrid(i)**2+1)) - (0.25)*log(1+kgrid(i)**2))*(2/pi)!, (16.0*kgrid(i)*(4.0*kgrid(i)**2.0 + 3.0)*sqrt(8.0*kgrid(i)**2.0 - 6.0))!/(81.0*(4.0*kgrid(i)**2.0 +1.0)**2.0) 
         end do
         close(1)

         open(1, file="Vdirectout.txt", action="write")
         do i=1,nk
             do j=1,nk
                 write(1, *) kgrid(i), kgrid(j), Vdirect(i,j)
             end do
                 write(1, *) ""
         end do  

         close(1)



         ! Overlap Matricies
        
         do i=1,nk
             foverlap(i) = sum(rweights(:)*wf(:,H_final)*sin(rgrid(:)*kgrid(i)))
             ioverlap(i) = sum(rweights(:)*wf(:,H_init)*sin(rgrid(:)*kgrid(i)))
         end do
 
         

! i = initial, j = final
         do i=1,nk
            !g = wf(:,H_final)*sin(kgrid(i)*rgrid(:)
            !A1(1) = rweights(1)*g(1)

            !do jj=1,nr
            !   A1(jj) = A1(jj-1) + rweights(jj)*g(jj) 
            !end do

            !A2(nr) = rweights(nr)*g(nr)*(1/rgrid(nr))
            !do jj=nr-1,1,-1
            !    A2(jj) = A2(jj+1) + rweights(jj)*g(jj)*(1/rgrid(jj))
            !end do

             

            !do j=1,nk
            !    f = sin(kgrid(j)*rgrid(:))*wf(:,H_init) 
            !    do ii=1,nr
            !        V12(i,j) = V12(i,j) + rweights(ii)*f(ii)*(1/rgrid(ii) * A1(ii) + A2(ii)) 
            !    end do
 
            !end do
         end do
         






         ! we are going to reuse f and g her. also A1 and A2
         g = 0.0d0 
         f= 0.0d0
         A1 = 0.0d0 
         A2 = 0.0d0
         
        


         do i=1,nk
             do j=1,nk
                 !f = sin(kgrid(i)*rgrid(:))*wf(:,)
                 do ii=1,nr
                     V12(i,j) = V12(i,j) !+          
                 end do
             end do
         end do








         deallocate(rgrid)
         deallocate(basis)
         if (N>1) then
            deallocate(H)
            deallocate(B)
            deallocate(V)
            deallocate(w)
            deallocate(z)
            deallocate(wf)
         endif

       

end program























