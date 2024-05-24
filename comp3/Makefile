FC := gfortran

CFLAGS := -O3 -ffree-line-length-512 -fdefault-real-8 -fdefault-double-8

#UNCOMMENT THE BELOW LINE FOR DEBUGGING:
CFLAGS := $(CFLAGS) -g -Og -fcheck=all -fbacktrace

#Write your program in main.f90

Numerov : NumerovsQHO.f90
	$(FC) $(CFLAGS) NumerovsQHO.f90  -o NumerovsQHO

dissociative : dissociative.f90
	$(FC) $(CFLAGS) rsg.f intp.f NumerovAppxMu.f90 dissociative.f90  -o dissociative


clean :
	rm NumerovsQHO
	rm dissociative 
