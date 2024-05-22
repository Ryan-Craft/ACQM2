FC := gfortran

CFLAGS := -O3 -ffree-line-length-512 -fdefault-real-8 -fdefault-double-8

#UNCOMMENT THE BELOW LINE FOR DEBUGGING:
CFLAGS := $(CFLAGS) -g -Og -fcheck=all -fbacktrace

#Write your program in main.f90

main : main.f90
	$(FC) $(CFLAGS) sort.f90 wigner.f rsg.f main.f90 -o main

vibe : vibe.f90
	$(FC) $(CFLAGS) rsg.f intp.f vibe.f90 -o vibe

clean : 
	rm main
