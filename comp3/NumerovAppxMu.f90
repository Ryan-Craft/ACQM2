subroutine NumerovForwardsMu(psi_L, V, nr, E, n, s, dr, mu)

    implicit none

    ! initialise local and inbound variables
    integer*8 :: i, j
    real *8, intent(in) :: s, E, dr, n, mu
    integer*8, intent(in) :: nr
    ! array init

    real* 8, dimension(nr), intent(in) :: V
    real* 8, dimension(nr):: psi_L
    real*8, dimension(nr) :: g
    real*8 :: psi_ip1, psi_ip2, denom


    g = 2*mu*(V-E)

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

end subroutine NumerovForwardsMu

subroutine NumerovBackwardsMu(psi_R, V, nr, E, n, s, dr, mu)

    implicit none

    ! initialise local and inbound variables
    integer*8 :: i, j
    real *8, intent(in) :: s, E, dr, n, mu
    integer*8, intent(in) :: nr

    real*8 :: psi_ip1, psi_ip2, denom

    ! array init

    real* 8, dimension(nr), intent(in) :: V
    real* 8, dimension(nr) :: psi_R
    real*8, dimension(nr) :: g

    g = 2*mu*(V-E)


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



end subroutine NumerovBackwardsMu
