program readtest

    real :: a, b, c, d, e,f 

    open(1, file="data.in", action="read")
    read(1,*) a, b, c, d, e, f
    Print *, a, b, c, d, e, f
    close(1)




end program readtest
