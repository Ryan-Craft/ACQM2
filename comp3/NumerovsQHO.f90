
subroutine NumerovForwards(psi_L, V, nr, E, n, s, dr)

    implicit none
    
    ! initialise local and inbound variables
    integer*8 :: i, j
    real *8, intent(in) :: s, E, dr, n
    integer*8, intent(in) :: nr
    ! array init

    real* 8, dimension(nr), intent(in) :: V
    real* 8, dimension(nr):: psi_L
    real*8, dimension(nr) :: g
    real*8 :: psi_ip1, psi_ip2, denom


    g = 2*(V-E)
    
    ! for the recurrence relation we need to initialise two values of psi
    psi_L(1) = 0.0d0
    psi_L(2) = (-1)**n * s
    Print *,  psi_L(2)

    ! this is essentially a recurrence relation like the Laguerre situation 
    do i=3,nr
        !psi_L(i) = ( 2*(1+(5*dr**2/12)*g(i-1))*psi_L(i-1) - (1 - (dr**2/12)*g(i-2))*psi_L(i-2) )  / (1 - (dr**2/12)*g(i))     
    
        denom = 1-(dr**2/12)*g(i)
        psi_ip1 = (1 + (5*dr**2/12)*g(i-1) )*psi_L(i-1)
        psi_ip2 = (1 - (dr**2/12)*g(i-2) )*psi_L(i-2)
        psi_L(i) = (1/denom) * ( 2*psi_ip1 - psi_ip2)

    end do
    
    Print *,  psi_L(2)   

end subroutine NumerovForwards

subroutine NumerovBackwards(psi_R, V, nr, E, n, s, dr)

    implicit none

    ! initialise local and inbound variables
    integer*8 :: i, j
    real *8, intent(in) :: s, E, dr, n
    integer*8, intent(in) :: nr

    real*8 :: psi_ip1, psi_ip2, denom
    
    ! array init
 
    real* 8, dimension(nr), intent(in) :: V
    real* 8, dimension(nr) :: psi_R
    real*8, dimension(nr) :: g

    g = 2*(V-E)


    ! for the recurrence relation we need to initialise two values of psi
    psi_R(nr) = 0.0d0
    psi_R(nr-1) = s

    !numerov but hes backwards
    do i= nr-2,1,-1
        denom = (1-((dr**2.0)/12.0)*g(i))
        psi_ip1 = (1.0 + (5.0*(dr**2.0)/12.0)*g(i+1) )*psi_R(i+1)
        psi_ip2 = (1.0 - ((dr**2.0)/12.0)*g(i+2) )*psi_R(i+2)
    !    Print *, denom, psi_ip1, psi_ip2
        psi_R(i) = (1.0/denom) * (2.0*psi_ip1-psi_ip2)

    end do 
    


end subroutine NumerovBackwards



program NumerovsQHO
    implicit none

    integer*8 :: i, j, nr, nodes, x_m
    real*8 :: dr, rmax, E, E_min, E_max, n, cooley_correct, e_lim
    real*8 :: fract, Yx, Yx1, Yxm1, Ecorr

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
    E_max = n + 0.75



    !! what is going on here?
    !! we actually use the method of bisections and the Numerov function to make initial guess
    !! for E and 

    pass_condition = .false.
    do while(pass_condition .eqv. .false.)   
        nodes = 0 
        E = (E_min + E_max)/2        

        call NumerovForwards(psi_L, V, nr, E, n, 0.00001, dr)
        call NumerovBackwards(psi_R, V, nr, E, n, 0.00001, dr)

        psi_L = psi_L / psi_L(x_m)
        psi_R = psi_R / psi_R(x_m)
   
        do i=1,nr
            if(i .le. x_m) then
                psi(i) = psi_L(i)
            else 
                psi(i) = psi_R(i)
            end if 
       end do

   ! count the number of nodes
       do i=1,nr-1
           if(psi(i) < 0 .AND. 0 < psi(i+1)) then
               nodes = nodes + 1
           elseif(psi(i) > 0 .AND. 0 > psi(i+1)) then       
               nodes = nodes + 1
           end if
       end do
       Print *, "Number of Nodes"
       Print *, nodes
      
       if(nodes==n) then
           pass_condition = .true.
       else if (nodes < n) then 
           E_min = E
       else if (nodes > n) then 
           E_max = E
       end if
       Print *, "Energy" 
       Print *, E

   end do

   e_lim = 0.0000001
   pass_condition = .false.
   Ecorr = E
   do while (pass_condition .eqv. .false.)
      
       g = 2*(V-Ecorr) 
       call NumerovForwards(psi_L, V, nr, Ecorr, n, 0.00001, dr)
       call NumerovBackwards(psi_R, V, nr, Ecorr, n, 0.00001, dr)

       psi_L = psi_L / psi_L(x_m)
       psi_R = psi_R / psi_R(x_m)

       do i=1,nr
           if(i .le. x_m) then
               psi(i) = psi_L(i)
           else
               psi(i) = psi_R(i)
           end if
       end do

     
       Yx = (1-(dr**2/12)*g(x_m))*psi(x_m)
       Yx1 = (1-(dr**2/12)*g(x_m+1))*psi(x_m+1)
       Yxm1 = (1-(dr**2/12)*g(x_m-1))*psi(x_m-1)
       fract = (psi(x_m)/sum(psi**2))

       cooley_correct = fract* ( -(0.5/dr**2)*(Yx1 - 2.0*Yx + Yxm1) + (V(x_m)-Ecorr)*psi(x_m)   )

       if(abs(cooley_correct) > e_lim) then
           Ecorr = Ecorr + cooley_correct

       else if (abs(cooley_correct) <= e_lim) then
           pass_condition = .true.
       end if 


       Print *, "cooley_correct:"
       Print *, cooley_correct
       Print *, "Corrected E"
       Print *, Ecorr

   end do
   


    open(unit=1, file="psi.txt", action="write")
    do i=1,nr
        write(1, *) rgrid(i), psi_L(i), psi_R(i), psi(i)
    end do
    close(1)



    g = 2*(V-E)
    open(unit=1, file="g.txt", action="write")
    do i=1,nr
        write(1, '(*(f12.5))') rgrid(i), g(i)
    end do
    close(1)


















end program
