FC := gfortran

CFLAGS := -O3 -ffree-line-length-512 -fdefault-real-8 -fdefault-double-8

#UNCOMMENT THE BELOW LINE FOR DEBUGGING:
CFLAGS := $(CFLAGS) -g -Og -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow

#Write your program in main.f90

main : main.o rsg.o
	$(FC) $(CFLAGS) rsg.o main.o -o main

rsg.o: rsg.f
	$(FC) $(CFLAGS) -c rsg.f


main.o: main.f90 
	$(FC) $(CFLAGS) -c main.f90
                          	



clean : 
	rm *.o main
