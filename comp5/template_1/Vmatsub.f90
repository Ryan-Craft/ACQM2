
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


subroutine Vmatsub(kgrid, Vtotal, nk, rmax, dr, alpha, l, S, theta, H_init, H_final, N)
         implicit none

         !! INTENT IN
         integer, intent(in) :: nk, H_init, H_final, N
         real*8, dimension(nk), intent(in) :: kgrid
         real*8, intent(in) :: dr, alpha, rmax, S, theta, l
 
         !! INTENT INOUT
         real*8, dimension(nk,nk), intent(out) :: Vtotal

         !! LOCAL

         real*8 :: normalise
         integer*8 :: p
         integer :: nr, ier
         integer :: i,j
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

         ! additions for assignment 5, energy grid
         integer :: ii, jj
         real*8 :: deltafi, energy, kp
         real*8, dimension(:), allocatable :: A1
         real*8, dimension(:), allocatable :: A2
         real*8, dimension(:), allocatable :: f
         real*8, dimension(:), allocatable :: g
         real*8, dimension(:), allocatable :: rgrid
         real*8, dimension(:), allocatable :: rweights
         real*8, dimension(:,:), allocatable :: Vdirect
         real*8, dimension(:), allocatable :: foverlap
         real*8, dimension(:), allocatable :: ioverlap
         real*8, dimension(:), allocatable :: V1
         real*8, dimension(:), allocatable :: V2
         real*8, dimension(:,:), allocatable :: V1mat
         real*8, dimension(:,:), allocatable :: V2mat
         real*8, dimension(:,:), allocatable :: V12
         real*8, dimension(:,:), allocatable :: Exchange
         real*8, dimension(:,:), allocatable :: Eoverlap
         real*8, dimension(:,:), allocatable :: analyticalEx
         real*8, dimension(:,:), allocatable :: analyticalV
         real*8, dimension(:,:), allocatable :: totalV
         real*8, dimension(:), allocatable :: onshellv, onshellEx

         Print *, "SUBROUTINE RUNNING: THETA, S"
         Print *, theta, S


         !open file location: hard coded for now but could become flexible
         !read stored values into relevent variables

         nr = rmax/dr
         Print *, "nr ::", nr
         if (mod(nr,2)==0) nr=nr+1
        

         ! based on options from file, allocate appropriate memory to rgrid and the basis array
         allocate(basis(nr,N), rgrid(nr), rweights(nr))
         allocate(H(N,N))
         allocate(B(N,N))
         allocate(V(N,N))
         allocate(K(N,N))
         allocate(w(N,1)) 
         allocate(z(N,N))
         allocate(wf(nr,N))
         allocate(Vdirect(nk,nk))
         allocate(ioverlap(nk), foverlap(nk), V1(nk), V2(nk), V12(nk,nk), V2mat(nk,nk), &
                  V1mat(nk,nk), Exchange(nk,nk), Eoverlap(nk,nk))
         allocate(analyticalEx(nk,3)) 
         allocate(analyticalV(nk,3))

         do i = 1, nr
                rgrid(i) = (i)*dr
         end do
         rweights = 1.d0
         do i =2,nr-1
             rweights(i) = 2.0d0 + 2.0d0*mod(i+1,2)
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
         
         CALL rsg(N,N,H,B,w,1,z,ier)


         !recover wavefunctions:
         wf = 0.0d0
         do i =1,N
                 do j = 1,N
                         wf(:,i) = z(j,i)*basis(:,j) + wf(:,i)
                 end do
                 if (wf(1, i) < 0) wf(:, i) = -wf(:, i)
         end do

         open(1, file='wfout.txt', action='write')
         do i =1,nr
                         write(1, '(*(f12.8))'), rgrid(i), wf(i,:)
         end do
         close(1)
         

         !!! Assignment 5 additional
         ! memory allocations

         allocate(A1(nr)) 
         allocate(A2(nr))
         allocate(f(nr)) 
         allocate(g(nr))
         allocate(onshellV(nk), onshellEx(nk)) 
         !g is a function of the initial and final wavefunctions of H
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

         open(1, file="Vdirectout.txt", action="write")
         do i=1,nk
             do j=1,nk
                 write(1, *) kgrid(i), kgrid(j), Vdirect(i,j)
             end do
                 write(1, *) ""
         end do
         close(1)



!!! SETUP EXCHANGE ELEMENTS
!!!! hard coded projectile energy
         energy = 54.4232/27.21136 - 0.5
         energy = energy * (1 - theta + theta*(-1)**S)
         ! V1 and V2
         do i =1,nk
             V1(i) = -sum(rweights * sin(kgrid(i)*rgrid) * wf(:,H_init)/rgrid) 
             V2(i) = -sum(rweights * sin(kgrid(i)*rgrid) * wf(:,H_final)/rgrid)
             ioverlap(i) = sum(rweights*sin(kgrid(i)*rgrid)*wf(:,H_final))
             foverlap(i) = sum(rweights*sin(kgrid(i)*rgrid)*wf(:,H_init))
         end do  
       
         ! V1mat V2mat calculation from overlap matrix
         do i = 1,nk
             do j = 1,nk
                 V1mat(j,i) = V1(j) * ioverlap(i)
                 V2mat(j,i) = V2(i) * foverlap(j) 
                 Eoverlap(j,i) = (energy - kgrid(j)**2/2.0d0 - kgrid(i)**2/2.0d0) * ioverlap(i)* foverlap(j)
             end do
         end do


         ! we are going to reuse f and g her. also A1 and A2
         g = 0.0d0 
         f = 0.0d0
         A1 = 0.0d0 
         A2 = 0.0d0


!!! CALCULATING EXCHANGE
         do i=1,nk
             g = sin(kgrid(i)*rgrid(:)) * wf(:,H_final)

             A1(1) = rweights(1)*g(1)
             do ii=2,nr
                 A1(ii) = A1(ii-1) + rweights(ii)*g(ii)
             end do

             A2(nr) = (1/rgrid(nr))*rweights(nr)*g(nr)
             do ii=nr-1,1,-1
                 A2(ii) = A2(ii+1) + (1/rgrid(ii))*rweights(ii)*g(ii)
             end do

             do j=1,nk
                 f = sin(kgrid(j)*rgrid(:))*wf(:,H_init)
                 V12(j,i) = sum(rweights*f*(A1/rgrid + A2))
             end do
         end do
         Exchange = 0.0d0
         Exchange = (Eoverlap - V1mat - V2mat - V12)*(2.0d0/pi)
 
         open(1, file="Exchange.txt", action="write")
         do i=1,nk
             do j=1,nk
                 write(1, *) kgrid(i), kgrid(j), Exchange(i,j)
             end do
                 write(1, *) ""
         end do
         close(1)

 
         allocate(totalV(nk,nk))
         Vtotal = Vdirect - (-1.0d0**S)*Exchange

         open(1, file="Vtotal.txt", action="write")
         do i=1,nk
             do j=1,nk
                 write(1, *) kgrid(i), kgrid(j), Vtotal(i,j)
             end do
                 write(1, *) ""
         end do
         close(1)



         deallocate(basis)
         if (N>1) then
            deallocate(H)
            deallocate(B)
            deallocate(V)
            deallocate(w)
            deallocate(z)
            deallocate(wf)
         endif

       

end subroutine Vmatsub























