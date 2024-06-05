FC := gfortran

CFLAGS := -O3 -ffree-line-length-512 -fdefault-real-8 -fdefault-double-8

#UNCOMMENT THE BELOW LINE FOR DEBUGGING:
CFLAGS := $(CFLAGS) -g -Og -fcheck=all -fbacktrace

#Write your program in main.f90

main : 
	$(FC) $(CFLAGS) rsg.f H-wf-subroutines.f90 main.f90 -o main

clean : 
	rm main
