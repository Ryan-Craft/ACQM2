


program signtest
    
    integer :: i, j, nodes
    real :: a(8) 


    nodes = 0

    a = [-1.0,-0.5, -0.5, 1.0, 0.5, 1.0, 0.0, -1.0]

    do i=1,7
       if(  a(i) <= 0 .AND. 0 < a(i+1)) then
           nodes = nodes + 1
       else if ( a(i) > 0 .AND. 0 >= a(i+1)  ) then
           nodes = nodes + 1
       end if
    end do

    Print *, nodes



end program
