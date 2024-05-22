
subroutine NumerovForwards(psi_L, V, rgrid, nr, E, n, s, dr)

    implicit none
    
    ! initialise local and inbound variables
    integer :: i, j
    real *8, intent(in) :: s, E, dr, n
    integer, intent(in) :: nr
    ! array init

    real* 8, dimension(nr), intent(in) :: rgrid
    real* 8, dimension(nr), intent(in) :: V
    real* 8, dimension(nr):: psi_L
    real*8, dimension(nr) :: g


    g = 2*(V-E)
    
    do i=1,nr
    !    Print *, g(i)
    end do

   
    ! for the recurrence relation we need to initialise two values of psi
    psi_L(1) = 0.0d0
    psi_L(2) = (-1)**n * s

    ! this is essentially a recurrence relation like the Laguerre situation 
    do i=3,nr
        psi_L(i) = ( 2*(1+(5*dr**2/12)*g(i-1))*psi_L(i-1) - (1 - (dr**2/12)*g(i-2))*psi_L(i-2) )  / (1 - (dr**2/12)*g(i))     
    end do
   

end subroutine NumerovForwards

subroutine NumerovBackwards(psi_R, V, rgrid, nr, E, n, s, dr)

    implicit none

    ! initialise local and inbound variables
    integer :: i, j
    real *8, intent(in) :: s, E, dr, n
    integer, intent(in) :: nr

    ! array init

    real* 8, dimension(nr), intent(in) :: rgrid
    real* 8, dimension(nr), intent(in) :: V
    real* 8, dimension(nr) :: psi_R
    real*8, dimension(nr) :: g

    g = 2*(V-E)

    do i=1,nr
     !   Print *, g(i)
    end do


    ! for the recurrence relation we need to initialise two values of psi
    psi_R(nr) = 0.0d0
    psi_R(nr-1) = s

    !numerov but hes backwards
    do i= nr-2,1, -1 
        psi_R(i) =  (psi_R(i+2)*(1-(dr**2/12)*g(i+2)) - 2*(1 + (5*dr**2 / 12)*g(i+1))*psi_R(i+1) )  / (1 - (dr**2/12)*g(i))
    end do 
    
    psi_R = psi_R*(-1)



end subroutine NumerovBackwards




program NumerovsQHO
    implicit none

    integer :: i, j, nr, nodes, x_m
    real*8 :: dr, rmax, E, E_min, E_max, n, cooley_correct, e_lim

    !ARRAY INITIALISATION
    real*8, dimension(:), allocatable :: rgrid, V, g
    real*8, dimension(:), allocatable :: psi_L, psi_R, psi    

    ! some stuff for the while loops
    logical :: pass_condition
    
    
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
    close(1)
    
    
    ! need to choose an x_m that is suitable. Looks like every second function is zero at zero. So for every second function ill displace 
    ! xm by some value to the right.

    if (mod(n,2.0) > 0 ) then
        !odd functions have a node at zero
        ! we actually want to reference an array element here for x_m so we can get the corresponding element of psi_~
        x_m = nr/2 + 1/dr

    else
        x_m = nr/2

    end if
    Print *, "x_m set to", x_m
 

    !Initial guess E, based on n
    ! I find this a bit weird because we know the energies of the QHO exactly; En = (n+1/2)*w*h_bar
    ! so considering that this is evenly spaced, my first energy is going to be 0.5, the next 1.5 and so on
    ! my guesses for En; E_min = n and E_max = n+0.75 so that we can get the iteration scheme to work

    
    E_min = n
    E_max = n+0.75



    !! what is going on here?
    !! we actually use the method of bisections and the Numerov function to make initial guess
    !! for E and 

    pass_condition = .false.
    do while(pass_condition .eqv. .false.)   
 
        E = (E_min + E_max)/2        

        call NumerovForwards(psi_L, V, rgrid, nr, E, n, 0.0001, dr)
        call NumerovBackwards(psi_R, V, rgrid, nr, E, n, 0.0001, dr)

        psi_L = psi_L / psi_L(x_m)
        psi_R = psi_R / psi_R(x_m)
   
        do i=1,nr
            if(rgrid(i) .le. x_m) then
                psi(i) = psi_L(i)
            else
                psi(i) = psi_R(i)
            end if 
       end do

   ! count the number of nodes
       do i=1,nr-1
           if(psi(i) < 0 .AND. 0 <= psi(i+1)) then
               nodes = nodes + 1
           elseif(psi(i) >= 0 .AND. 0 > psi(i+1)) then       
               nodes = nodes + 1
           end if
       end do
       Print *, nodes
      
       if(nodes==n) then
           pass_condition = .true.
       else if (nodes < n) then 
           E_min = E
       else if (nodes > n) then 
           E_max = E
       end if
       
       Print *, E

   end do
   
   e_lim = 0.0001    
   pass_condition = .false.
   g = 2*(V-E)
   do while(pass_condition .eqv. .false.)
       
       !compute numerov cooley
      
       ! there is a bug right here in the following line somewhere 
       cooley_correct = (psi(x_m)/sum(psi**2)) * (-0.5 * ( (1.0- (dr**2/12.0)*g(x_m +1))*psi(x_m+1) - 2.0*( (1.0-(dr**2/12.0)*g(x_m))*psi(x_m)) + (1.0- (dr**2/12.0)*g(x_m-1))*psi(x_m-1))/dr**2 + (V(x_m)-E)*psi(x_m)  )
       
       if(abs(cooley_correct) <= e_lim ) then
           pass_condition = .true.

       else if (abs(cooley_correct) > e_lim) then
           E = E + cooley_correct
           call NumerovForwards(psi_L, V, rgrid, nr, E, n, 0.0001, dr)
           call NumerovBackwards(psi_R, V, rgrid, nr, E, n, 0.0001, dr)
 
       end if 
       Print *, "Cooley Correction" 
       Print *, cooley_correct

       do i=1,nr
            if(rgrid(i) .le. x_m) then
                psi(i) = psi_L(i)
            else
                psi(i) = psi_R(i)
            end if
       end do

   end do

   Print *, "Corrected Energy"
   Print *, E


    open(unit=1, file="psi.txt", action="write")
    do i=1,nr
        write(1, '(*(f12.4))') rgrid(i), psi_L(i), psi_R(i), psi(i)
    end do
    close(1)





















end program
