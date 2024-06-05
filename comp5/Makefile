FC := gfortran

CFLAGS := -O3 -ffree-line-length-512 -fdefault-real-8 -fdefault-double-8

#UNCOMMENT THE BELOW LINE FOR DEBUGGING:
CFLAGS := $(CFLAGS) -g -Og -fcheck=all -fbacktrace

#Write your program in main.f90

main : main.o H-wf-subroutines.o rsg.o 
	$(FC) $(CFLAGS) rsg.o H-wf-subroutines.o main.o -o main

main.o: main.f90
	$(FC) $(CFLAGS) -c main.f90

H-wf-subroutines.o: H-wf-subroutines.f90
	$(FC) $(CFLAGS) -c H-wf-subroutines.f90

rsg.o: rsg.f
	$(FC) $(CFLAGS) -c rsg.f 


clean : 
	rm *.o main
