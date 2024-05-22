subroutine sortenergies(energies, nstates)
    implicit none
    logical :: sorted
    integer, intent(in) :: nstates
    real*8, dimension(nstates), intent(inout) :: energies
    real*8 :: E1,E2,temp
    integer :: i

    sorted = .false.
    do while(.not. sorted)
        sorted = .true.
        do i=1,nstates-1
            E1 = energies(i)
            E2 = energies(i+1)
            
            if(E2<E1) then
                sorted = .false.
                temp = energies(i)
                energies(i) = energies(i+1)
                energies(i+1) = temp
            endif
        end do  
    end do

    return
end subroutine