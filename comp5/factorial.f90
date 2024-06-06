
program factorial

        integer :: i,t
        integer*8 :: p
        real :: a, l,k
        real*8 :: f1, f2, s
       
        k= 2097150

        a=1.0
        l=1.0
        
        f1 = sqrt((a*gamma(k))/((k+l)*gamma(k+2*l+1))) 
        

        Print *, "Gamma factorial: ", f1

        !implement my new method
        p=1.0
        do i = 0, 2*l
                p = p*(k+2*l-i)
                Print *, p, i
        end do

        f2 = sqrt(a /((k+l)* p))

        Print *, "Simple factorial: ", f2

end program

