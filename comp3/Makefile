FC := gfortran

CFLAGS := -O3 -ffree-line-length-512 -fdefault-real-8 -fdefault-double-8

#UNCOMMENT THE BELOW LINE FOR DEBUGGING:
CFLAGS := $(CFLAGS) -g -Og -fcheck=all -fbacktrace

#Write your program in main.f90

main : NumerovsQHO.f90
	$(FC) $(CFLAGS) NumerovsQHO.f90  -o NumerovsQHO

clean :
	rm NumerovsQHO
