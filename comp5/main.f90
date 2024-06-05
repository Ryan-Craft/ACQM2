program main

    implicit none
    
    real*8, dimension(2001,10) :: wf
    real*8, dimension(2001) :: rgrid
    integer :: i,j, N
    wf=0.0d0
    rgrid = 0.0d0
  
    do i=2,2001
        rgrid(i) = (i-1)*0.01
    end do 
    
    N=50
         
    
    

    call hwfsubroutine(1.0d0,1.0d0,N,2001,0.01,20.0d0,rgrid,wf)
                      !(alpha, l, N, nr, dr, rmax, rgrid, wf)
    



end program























