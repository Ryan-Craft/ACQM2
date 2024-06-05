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


subroutine hwfsubroutine(alpha, l, N, nr, dr, rmax, rgrid, wf)
         implicit none

         ! LOCAL
         real*8 :: normalise
         integer :: i,j, ier
	 integer*8 :: p
         real*8, dimension(nr,N) :: basis
         real*8, dimension(N,N) :: H
         real*8, dimension(N,N) :: B
         real*8, dimension(N,N) :: K

         ! energies and expansion coefficients
         real*8, dimension(N,1) :: w
         real*8, dimension(N,N) :: z
         real*8, dimension(N,N) :: V



         ! IN
         real*8, intent(in) :: alpha, l, dr, rmax
         integer, intent(in) :: N,nr 
         real, dimension(nr) :: rgrid
         

         ! INOUT
         real*8, dimension(nr,N), intent(inout) :: wf
       
         !use recurrence relation to compute the basis functions
         CALL LaguerreSub(alpha, l, nr, N, rgrid, basis)
        
         !implement normalisation condition using simplified factorial
         do i = 1, N
                 p=1.0
                 do j = 0, 2*l
                         p = p*(i+2*l-j)
                         Print *, p, j
                 end do
                 normalise = sqrt(alpha /((i+l)* p))
                 Print *, "Norm:: ", normalise
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
                 Print *, B(i,:)
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
         do i =1,N
                 do j = 1,N
                         !wf(:,i) = z(j,i)*basis(:,j) + wf(:,i)
                         !Print *, wf(:,i)
                 end do
         end do
         !Print *, "W::"
         !Print *, w(1,:)
         !Print *, "Z::"
         !Print *, z(1,:)

         open(1, file='wfout.txt', action='write')
         do i =1,nr
!                         write(1, '(*(f12.8))'), rgrid(i), wf(i,:)
         end do
         close(1)
         
         open(1, file='wout.txt',action='write', access='append')
         do i =1,N
 !                write(1, '(*(f12.8))'), real(N), w(i,1)
         end do
 
       

end subroutine hwfsubroutine























