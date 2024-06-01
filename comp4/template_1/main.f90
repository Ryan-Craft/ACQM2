!  This file is an outline of a code to perform a calculation of a
!    charged particle scattered by a spherically-symmetric short-range local
!    potential.
!  There are a number of comments throughout the file beginning with !>>> which
!    indicate sections of code you will need to complete.

! Last modified: May 25 2024
!                Liam Scarlett

program main

  use constants
  implicit none

  real*8, allocatable :: &
    Vmat(:,:),       & !V-matrix elements Vmat(kf,ki)
    V(:),            & !radial potential function V(r)
    rgrid(:),        & !radial grid
    rweights(:),     & !radial integration weights
    kgrid(:),        & !momentum-space grid
    kweights(:),     & !momentum-space integration weights (and Green's func)
    contwaves(:,:),  & !projectile radial continuum waves contwaves(k,r)
    DCS(:),          & !array to hold differential cross section - DCS(theta)
    theta(:),        & !array to hold values of theta - in degrees
    ICS(:),          & !integrated cross section per l
    psi(:)             !defined so that numerov has a temp vector to work with


  real*8 :: &
    rmax,   & !max value of radial grid
    dr,     & !radial grid step size
    energy, & !projectile energy
    k,      & !projectile momentum
    kg_A, kg_B, kg_P !some parameters for setting up kgrid

  complex*16, allocatable :: Ton(:) !on-shell T-matrix element per l

  integer :: &
    nrmax,      & !number of rgrid points
    nkmax,      & !number of kgrid points
    zproj,      & !projectile charge 
    l,          & !partial-wave angular momentum
    lmin, lmax, & !min and max values of l
    iounit,     & !a unit number for input/ouput
    ntheta,     & !an index to iterate over theta
    nthetamax,  & !max number of theta
    kg_Na, kg_Nb, kg_Np, &!some more parameters for setting up kgrid
    i, j !RC: I need some iteration params so they are going here as well


  !set kgrid parameters - leave this as is
    kg_Na = 30; kg_Nb = 30; kg_Np = 10
    kg_a = 0.85; kg_b = 2.5; kg_p = 4.0
    nkmax=kg_Na+kg_Nb+kg_Np+1

  !>>> open data.in file and read input parameters
  !    note: energy should be read in electron volts
  !      and grid parameters in atomic units

!RC: I have read in the energy directly as eV, which is what it is in the data.in file I presume
    open(1, file="data.in", action="read")
    read(1,*) energy, rmax, dr, zproj, lmin, lmax
    Print *, "INPUT PARAMS FROM data.in::"
    Print *, energy, rmax, dr, zproj, lmin, lmax
    close(1)

  !>>> do any input validation you think is necessary here
!RC: might come back to this


  !>>> convert the energy to atomic units and calculate the
  !      projectile momentum
    energy = energy/eV
    k = sqrt(energy*2)
    Print *, "Energy in Hartrees"
    Print*, energy


  !>>> determine number of rgrid points nrmax
  !    note: nrmax should be even for simpson's integration
  !          to take into account that the r=0 point has been omitted 

!RC: for the values of dr and rmax in the data.in, this is an even number, so we follow the rules on the 
!    lecture slides for creating the rgrid
    nrmax = rmax/dr
    Print *, nrmax 



  !allocate memory
    allocate(rgrid(nrmax),rweights(nrmax))
    allocate(kgrid(nkmax),kweights(nkmax))
    allocate(contwaves(nkmax,nrmax))
    allocate(V(nrmax))
    allocate(Ton(lmin:lmax),ICS(lmin:lmax))
    allocate(Vmat(nkmax,nkmax))
    allocate(psi(nrmax))

  !setup grids
    call setup_rgrid(nrmax, dr, rgrid, rweights)
    call setup_kgrid(k, nkmax, kg_Na, kg_a, kg_Nb, kg_b, kg_Np, kg_p, kgrid, kweights)

    Print *, "RGRID and RWEIGHTS" 
    Print *, rgrid(1:5)
    Print *, rweights(1:5)
    Print *, "KGRID FIRST 5 ELEMENTS and the end one ::"
    Print *, kgrid(1:6), kgrid(nkmax)
    Print *, "KGRID SIZE ::", size(kgrid)  

  
  !>>> define short-range potential V(r)
!RC: I am also writing this potential to file because, A: im probably going to want to plot with it
!                                                      B: I want to inspect it so I know its not broken
!RC: I checked it, it looks how youd expect.
  open(1, file="Scattering-Potential.txt", action="write")
  do i=1,nrmax
    V(i) = zproj*(1.0d0+(1.0d0/rgrid(i)))*exp(-2.0d0*rgrid(i))
    write(1, *) rgrid(i), V(i)
  end do
  close(1)
  


  !begin loop over angular momenta
  do l=lmin, lmax
    !populate contwaves matrix with a continuum wave for each off-shell k
!RC: implemented but all cont waves normalised to unity
      contwaves = 0.0d0
      call setup_contwaves(nkmax,kgrid,l,nrmax,rgrid,contwaves)
 
    !evaluate the V-matrix elements  
!RC :: Implemented, checking now. Checked it, looks pretty good
      call calculate_Vmatrix(nkmax,kgrid,contwaves,nrmax,rgrid,rweights,V,Vmat)
      
      open(1, file="Vmat-halfonshell.txt", action="write")
      do i=1,nkmax
        write(1, *) kgrid(i), Vmat(i,1)
      end do
      close(1)




    !solve the Lippman-Schwinger equation for the on-shell T-matrix
!RC ::     
      call tmatrix_solver(nkmax,kgrid,kweights,Vmat,Ton(l))




  enddo

  !populate theta - theta = 0, 1,..., 180 degrees
  nthetamax = 181
  allocate(theta(nthetamax),DCS(nthetamax))
  do ntheta = 1, nthetamax
    theta(ntheta) = dble(ntheta-1)
  enddo

  !call subroutines to calculate DCS and ICS
    call compute_dcs(nthetamax, theta, lmin, lmax, Ton, k, DCS)
    call compute_ics(lmin, lmax, Ton, k, ICS)

  !>>> output the DCS and ICS to files
  !    to easily study the convergence you can write the ics to file as a function of l
  !    along with a running total so that the running total in the final line is your
  !    total ICS summed over l


end program main

subroutine compute_ics(lmin, lmax, Ton, k, ICS)
  use constants
  implicit none
  integer, intent(in) :: lmin, lmax
  complex*16, intent(in) :: Ton(lmin:lmax)
  real*8, intent(in) :: k
  real*8, intent(out) :: ICS(lmin:lmax)
  integer :: l

  !>>> populate the ICS array with the partial-wave ICS per l
end subroutine compute_ics

subroutine compute_dcs(nthetamax, theta, lmin, lmax, Ton, k, DCS)
  use constants
  implicit none
  integer, intent(in) :: nthetamax, lmin, lmax
  real*8, intent(in) :: theta(nthetamax), k
  complex*16, intent(in) :: Ton(0:lmax)
  real*8, intent(out) :: DCS(nthetamax)
  integer :: l, ntheta !loop indices
  real*8:: PL !Legendre polynomials - from file plql.f
  real*8 :: costheta !use this to store cos(theta in radians)
  complex*16 :: f(nthetamax) !scattering amplitude

  !>>> calculate the scattering amplitude f(theta) for each theta
  !    by iterating over l and using the partial-wave
  !    expansion of f

  !>>> obtain the DCS from the scattering amplitude

end subroutine compute_dcs

subroutine setup_rgrid(nrmax, dr, rgrid, rweights)
  implicit none
  integer, intent(in) :: nrmax
  real*8, intent(in) :: dr
  real*8, intent(out) :: rgrid(nrmax), rweights(nrmax)
  integer :: ir !index to iterate over r

  !>>> iterate over r and populate the rgrid and rweights arrays
  !      - rweights should contain Simpson's integration weights:
  !        (4, 2, 4, 2, ..., 2, 4, 2) * dr / 3.0
  !      - you can make use of the intrinsic MOD function for the 
  !        alternating 4, 2 terms
  !      - note we have neglected the terms with a coefficient of 1 (rather than 4 or 2) since
  !        the first term (r=0) is skipped and the last term corresponds to the end of the
  !        radial grid where we assume all functions should be zero (and if not the grid is not large enough)

  do ir=1, nrmax
    rgrid(ir) = ir*dr    
    !RC: started at dr instead of 0
  end do
  
  rweights = 0.0d0
  do ir =1,nrmax
    rweights(ir) = 4.0d0 - 2.0d0*mod(ir+1,2)
  end do
  rweights = rweights * dr / 3.0d0


end subroutine setup_rgrid

subroutine setup_contwaves(nkmax, kgrid, l, nrmax, rgrid, contwaves)
  implicit none
  integer, intent(in) :: nkmax, l, nrmax
  real*8, intent(in) :: kgrid(nkmax), rgrid(nrmax)
  real*8, intent(out) :: contwaves(nkmax,nrmax)
  real*8 :: ncontwaves(nkmax,nrmax)
  real*8 :: g(nrmax) ! RC: I added this one
  integer :: nk, nr, nodes !indices to loop over k and r
  real*8 :: E
  logical :: pass_condition 
  !>>> iterate over k, populating the contwaves matrix                                 
!RC : We need to use forwards numerov to get the continuum waves

  do nk=1,nkmax
      E = kgrid(nk)**2/2
      g = 2*( l*(l+1)/(2*rgrid**2)  - E)
      call NumerovForwards(nrmax, rgrid, nkmax, kgrid(nk), contwaves(nk,:), g, l) 
  
    
    nr=nrmax
    pass_condition = .false.
    nodes=0
    do while(pass_condition .eqv. .false.)
      if(contwaves(nk,nr) < 0 .AND. 0 <= contwaves(nk,nr-1)) then
        nodes = nodes + 1
      elseif(contwaves(nk,nr) >= 0 .AND. 0 > contwaves(nk,nr-1)) then
        nodes = nodes + 1
      end if 

      if (nodes > 1) then
        pass_condition =.true.
        !Print *, "Index"
        contwaves(nk,:) = contwaves(nk,:) / maxval(abs(contwaves(nk,nr:nrmax)))
      endif

      nr = nr-1
    end do 
  end do

  
   
! These need to made unit valued at their asymptotic region and then normalised


  Print *, contwaves(1,1:5)
  Print *, contwaves(2,1:5)
  Print *, contwaves(3,1:5)

  Print *, "Writing Numerovl0vsSin.txt, contains rgrid, first continuum wave and sin(kr)" 
  open(1, file="Numerovl0vsSin.txt", action="write")
  do nr=1,nrmax
    write(1, *) rgrid(nr), contwaves(20,nr), sin(kgrid(20)*rgrid(nr))
  end do
  close(1)




end subroutine setup_contwaves

   
!>>> your forwards Numerov subroutine can go here
subroutine NumerovForwards(nrmax, rgrid, nkmax, kval, psi, g, l)
    implicit none 
    integer, intent(in) :: nrmax, l, nkmax
    real*8, intent(in) :: g(nrmax), rgrid(nrmax), kval 
    real*8, intent(inout):: psi(nrmax)
    real*8 :: psi_ip1, psi_ip2, denom, dr
    integer*8 :: dfactorial, i, j ! RC : High precision integers becaues of the factorial

    dr = rgrid(1)

!RC : Well we need to make a code to calculate a double factorial now
    dfactorial=1
    do i = 2*(l+1), 0, -2
        if(i==0 .or. i==1) then
              dfactorial=dfactorial
        else
            dfactorial = dfactorial * i
        end if
    end do

    psi(1) = (rgrid(1)*kval)**(l+1) / dfactorial !RC : because of page 72 of the lectures
    psi(2) = (rgrid(2)*kval)**(l+1) / dfactorial
    !Print *, "NUMEROV LEFT BOUNDARY ::"
    !Print *, psi(1), psi(2)

    do i=3, nrmax 
        denom = 1-(dr**2/12)*g(i)
        psi_ip1 = (1 + (5*dr**2/12)*g(i-1) )*psi(i-1)
        psi_ip2 = (1 - (dr**2/12)*g(i-2) )*psi(i-2)
        psi(i) = (1/denom) * ( 2*psi_ip1 - psi_ip2)
    end do
end subroutine NumerovForwards



subroutine calculate_Vmatrix(nkmax,kgrid,contwaves,nrmax,rgrid,rweights,V,Vmat)
  use constants
  implicit none
  integer, intent(in) :: nkmax, nrmax
  real*8, intent(in) :: kgrid(nkmax), contwaves(nkmax,nrmax), rgrid(nrmax), rweights(nrmax), V(nrmax)
  real*8, intent(out) :: Vmat(nkmax,nkmax)
  integer :: nkf,nki !indices for looping over on- and off-shell k

  !>>> evaluate the V-matrix elements and store in the Vmat matrix
  !    note: the V-matrix is symmetric, make use of this fact to reduce the 
  !          amount of time spent in this subroutine
  nkf=0.0d0
  do nkf =1, nkmax
    do nki =nkf,nkmax
      
      Vmat(nkf,nki) = (2/pi)*sum(contwaves(nkf,:)*V(:)*contwaves(nki,:)*rweights(:))
      Vmat(nki,nkf) = Vmat(nkf,nki) 
    end do
  end do
  Print *, "V Matrix top left corner ::"
  Print *, Vmat(1,1:4)
  Print *, Vmat(2,1:4)
  Print *, Vmat(3,1:4)
  Print *, Vmat(4,1:4)

end subroutine calculate_Vmatrix
    
subroutine tmatrix_solver(nkmax,kgrid,kweights,Vmat,Ton)
  use constants
  implicit none
  integer, intent(in) :: nkmax
  real*8, intent(in) :: kgrid(nkmax), kweights(nkmax), Vmat(nkmax,nkmax)
  complex*16, intent(out) :: Ton !on-shell T-matrix element
  complex*16 :: denom
  real*8 :: &
    Koff(nkmax-1), & !half-off-shell K-matrix elements
    Kon,           & !on-shell K-matrix element
    Von,           & !on-shell V-matrix element
    A(nkmax-1,nkmax-1) !Coefficient matrix for the linear system Ax=b
  integer :: f,n,j, ipiv(nkmax-1), info

  !>>> store the on-shell V-matrix element in Von
  Von = Vmat(1,1)

  !>>> populate the matrix A according to Eq (142) in the slides
 
  do f=1, nkmax-1
    do n=1, nkmax-1
      if(f==n) then
        A(f,n) = 1 - kweights(n+1)*Vmat(f+1,n+1)
      else
        A(f,n) = - kweights(n+1)*Vmat(f+1,n+1)
      end if
    end do
  end do

 
  !>>> populate the vector Koff with the half-on-shell V-matrix elements (RHS of Eq (141)) 
  do n=1, nkmax-1
    Koff(n) = Vmat(n+1,1) 
  end do
 
  !Here is the call to DGESV
  call dgesv(nkmax-1, 1, A, nkmax-1, ipiv, Koff, nkmax-1, info )
  if(info /= 0) then
    print*, 'ERROR in dgesv: info = ', info
  endif

  !>>> Now use the half-on-shell K matrix which has been stored in Koff to get the on-shell K-matrix element Kon
  
  Kon = Vmat(1,1) + sum(kweights*Koff*Vmat(2:,1))
 
  !>>> And then use Kon to get the on-shell T-matrix element Ton
  
! RC :: POSSIBLE BUG !!!!!! idk what k_f was, page 110 says the magnitudes of ki and kf are the same so... is it just ki???
  denom = complex(1, (pi/kgrid(1))*Kon)

  Ton = Kon/denom
 
  Print *, "Complex declaration :::"
  Print *, "Kon", Kon
  Print *, Ton

end subroutine tmatrix_solver

!A subroutine provided for you to set up the kgrid and kweights
!note: the kgrid is setup with the on-shell point in the first element
!      and the corresponding kweights include the integration weights
!      AND the Green's function
subroutine setup_kgrid(k,nkmax,Na,a,Nb,b,Np,p,kgrid,kweights)
  implicit none
  real*8, intent(in) :: k
  integer, intent(in) :: Na, Nb, Np
  integer, intent(in) :: nkmax
  real*8, intent(out) :: kgrid(nkmax), kweights(nkmax)
  integer :: nk
  real*8, intent(in) ::  a, b, p
  real*8 :: grid1(nkmax-1), weight1(nkmax-1)

  call kgrid_igor(0.0,k,a,Na,b,Nb,p,Np,nkmax-1,grid1,weight1)
  
  kgrid(1) = k
  kgrid(2:nkmax) = grid1
  kweights(1) = 0.0d0
  kweights(2:nkmax) = weight1
  do nk=2, nkmax
    kweights(nk) = 2.0d0* kweights(nk) / (k**2 - kgrid(nk)**2)
  enddo
end subroutine setup_kgrid
