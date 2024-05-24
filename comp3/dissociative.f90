
subroutine LaguerreSub(alpha, l, nr, N, rgrid, basis)
        implicit none

        ! initialise alpha, l, dr, rmax, N and others
        integer :: i
        integer*8, intent(in) :: nr
        integer*8, INTENT(IN) :: N
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


program vibe
         implicit none

         real*8 :: normalise
         real*8 :: alpha, l
         real*8 :: dr, rmax, mu
         integer*8 :: p
         integer*8 :: N, nr, ier
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
         real*8,dimension(:,:), allocatable :: wf

         !create an array to hold the interpolated ssg potential
         real*8, dimension(:,:), allocatable :: V_interp

         !an array to hold potential (called ssg but really could be any of them, big readability issue)
         real*8, dimension(:,:), allocatable :: ssg, psu
         real*8, dimension(:), allocatable :: weights
         integer :: lines
 
         !!! We are going to need more declarations for our dissociative states now
         real*8 :: Emax, dE, E 
         real*8, dimension(:), allocatable :: g, psi_L, psi_R, psi, Ekgrid, potent_Numerov
         real*8, dimension(:,:), allocatable :: Ekwf, frank
         integer, dimension(4) :: idx = [1,4,7,10]
         integer*8 :: nk, x_m, nodes
         logical :: pass_condition

         !open file location: hard coded for now but could become flexible
         !read stored values into relevent variables
         open(unit=1, file="dissociativeParams.txt", action="read")
         read(1,*) alpha, N, l, dr, rmax, mu, Emax, dE
         Print *, alpha, N, l, dr, rmax, mu, Emax, dE
         close(1)
         !calculate rgrid params
         nr = rmax/dr + 1.0
         Print *, nr


         !!! Numerov allocations
         nk = Emax/dE + 1
         Print *, nk

         allocate(psi_L(nr), psi_R(nr), psi(nr))
         allocate(Ekgrid(nk))
         allocate(Ekwf(nr,nk))
         allocate(g(nr))
         !allocate(potent_Numerov)

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
         allocate(V_interp(nr-1,2))
         allocate(weights(nr))
         allocate(potent_Numerov(nr)) 
         

         !allocate values to the rgrid: allows arbitrary size of PEC.1ssg, though its a bit clumsy
         do i = 1, nr
                rgrid(i) = (i-1.0)*dr
         end do 
         
         !read in PEC.1ssg potential into array: generalised in case we get bigger potentials later
                ! check how many lines are in PEC.1ssg
         lines = 0
         open(1, file='PEC_good/PEC.1ssg', action='read')
         do
                read(1,*, END=10)
                lines = lines +1
         end do
         10 close(1)
         Print *, "lines"
         Print *, lines

         allocate(ssg(lines, 2))
                ! read in PEC
         open(1, file='PEC_good/PEC.1ssg', action='read')
         do i=1,lines
                read(1,*), ssg(i,1), ssg(i,2)
         end do
         close(1)


         ! interpolate the electron potential:
         ! intrpl (number of data points, x array for input, yarray for input, number of points
         ! to interpolate, x values of desired points)
         V_interp = 0.0d0
         call intrpl(lines, ssg(:,1), ssg(:,2), nr-1, rgrid(2:), V_interp(:,2))

         V_interp(:,1) = rgrid(2:)
         do i=1,nr-1
                !Print *, V_interp(i,:)
         end do

         !!! create integration weights
         weights = 0.0d0
         weights(1) = 1.0d0
         do i=2, nr-1
                 weights(i) = 2.0d0 + 2.0d0*mod(i+1,2)
         end do
         weights(nr) = 1.0d0
         weights(:) = weights(:)*dr/3.0d0


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
                 !Print *, "Norm:: ", normalise
                 basis(:,i) = normalise*basis(:,i)  
         end do
         !write basis to file for plotting 
         
         open(1, file='basisout.txt', action='write')
         do i =1,nr
                         write(1, '(*(f12.8))') rgrid(i), basis(i,:)
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
         !do i=1,N
                 !Print *, B(i,:)
         !end do
 
         !calculate V matrix:
         V = 0.0d0
         do i=1,N
                do j=1,N
                        V(i,j) = sum(basis(2:,i) * V_interp(:,2) * basis(2:,j) * weights(2:))
                end do
         end do

         !create K to create H
         K = (-alpha**2/2.0) * B
         do i=1,N
                 K(i,i) = K(i,i) + alpha**2
         end do
         
         !compute H-matrix Elements
         H = (1/mu)*K + V
         
         ! optional prints in case we need debugging the matricies
         Print *, "H MATRIX::"
         !do i=1,N
                 !Print *, H(i,:)
         !end do         
 
         !Print *, "V MATRIX::"
         !do i=1,N
                 !Print *, V(i,:)
         !end do  
         
         !Print *, "K MATRIX::"
         !do i=1,N
                 !Print *, K(i,:)
         !end do 
         
         

         CALL rsg(N,N,H,B,w,1,z,ier)
         !recover wavefunctions:
         wf = 0.0d0
         do i =1,N
                 do j = 1,N
                         wf(:,i) = z(j,i)*basis(:,j) + wf(:,i)
                 end do
         end do
         Print *, "W::"
         Print *, w(:,1)

         ! write wavefunctions to a text file (caution, file sizes a pretty big)
         open(1, file='boundwfout.txt', action='write')
         do i =1,nr
                         write(1, '(*(f12.8))') rgrid(i), wf(i,:)**2
         end do
         close(1) 



         !!! the code for vibe.f90 stopped here and deallocated everthing, but we have a few more things to do
         !!! all declarations are at the top, so we dont make a complete mess, but this section onwards is dedicated to 
         !!! applying NUMEROVs iteration scheme to get the dissociative wavefunctions in the PEC.2psu potential. 


         ! need to calculate nk; the number of steps along Ek (to make our Ek grid)

         Print *, "NUMEROV METHOD :::::"
         
         Ekgrid = 0.0d0 
         Ekwf = 0.0d0 
         psi_L = 0.0d0 
         psi_R =0.0d0 
         psi=0.0d0
         g = 0.0d0  

         do i=1,nk
             Ekgrid(i) = (i-1.0)*dE
         end do
         
         Print *, Ekgrid(:)
         !!! CONVERT TO HARTREES
         Ekgrid = Ekgrid
 
         

         !for i in ekgrid
             !call numerovForwards, and backwards
             !match them up (choose a good xm using some kind of scheme)
             ! normalise them 
             ! send to file called something appropriate
         lines = 0
         open(1, file='PEC_good/PEC.2psu', action='read')
         do
                read(1,*, END=20)
                lines = lines +1
         end do
         20 close(1)

         Print *, "lines"
         Print *, lines

         allocate(psu(lines, 2))
         open(1, file='PEC_good/PEC.2psu', action='read')
         do i=1,lines
                read(1,*), psu(i,1), psu(i,2)
         end do
         close(1)

               
         call intrpl(lines, psu(:,1), psu(:,2), nr, rgrid(:), potent_Numerov)
 
         do i =1, nk
             
             E = Ekgrid(i)/27.21136 - 0.5             
             g = 2*mu*(potent_Numerov - E)
             
             Print *, "G size"
              
             
             call NumerovForwardsMu(psi_L, nr, 0, 0.0001, dr, mu, g)
             
             ! finding x_m
             !do j=1,
             j=nr
             pass_condition = .false.
             nodes=0
             do while(pass_condition .eqv. .false.)
                 
                  if(psi_L(j) < 0 .AND. 0 < psi_L(j-1)) then
                      nodes = nodes + 1
                  elseif(psi_L(j) > 0 .AND. 0 > psi_L(j-1)) then
                      nodes = nodes + 1
                  end if
                  Print *, "nodes"
                  Print *, nodes

                  if (nodes > 1) then
                      pass_condition =.true. 
                      Print *, "Index"
                       
                      psi_L = psi_L / maxval(abs(psi_L(j:nr)))
 
                  endif

             j = j-1
             end do            

 
             ! we can normalise here as well
             normalise = sqrt(2*mu/(4.0d0*datan(1.d0)*sqrt(2*mu*Ekgrid(i)/27.21136)))
             Print *, "NORMALISE:::", normalise
             Ekwf(:,i) = psi_L(:) 

         end do 

         open(1, file='unboundout.txt', action='write')
         do i =1,nr
                         write(1, *) rgrid(i), Ekwf(i,:)
         end do
         close(1)         

                  
         ! okay now lets frank condon:

         allocate(frank(nk,4))
         Print *, size(Ekgrid(:)), size(frank(:,1))
         do i=1,4
             do j=1,nk
                 if (j==1) then
                     frank(j,i)=0.0d0                 
                 else
                     frank(j,i) = sum(weights(:)*Ekwf(:,j)*wf(:,idx(i)))**2
                 end if
             end do
         end do  

         open(1, file='condonout.txt', action='write')
         do i =1,nk
                         write(1, *) Ekgrid(i), frank(i,:)
         end do
         close(1)



















         ! deallocate memory
         deallocate(V_interp)
         deallocate(weights)
         deallocate(rgrid)
         deallocate(basis)
         if (N>1) then
            deallocate(Ekwf)
            deallocate(Ekgrid)
            deallocate(H)
            deallocate(B)
            deallocate(V)
            deallocate(K)
            deallocate(w)
            deallocate(z)
            deallocate(wf)
         endif

       

end program























