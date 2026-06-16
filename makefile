all: .o .mod mesa_mod_demo mesa_mod_summary

.o:
	-mkdir .o
	
.mod:
	-mkdir .mod

.o/mesa_iso_mass.o: mesa_iso_mass.f90
	gfortran -J.mod -I.mod -I.o -c mesa_iso_mass.f90 -o .o/mesa_iso_mass.o

.o/mesa_mod_reader.o: mesa_mod_reader.f90
	gfortran -J.mod -I.mod -I.o -c mesa_mod_reader.f90 -o .o/mesa_mod_reader.o

.o/mesa_mod_demo.o: mesa_mod_demo.f90
	gfortran -J.mod -I.mod -I.o -c mesa_mod_demo.f90 -o .o/mesa_mod_demo.o

.o/mesa_mod_summary.o: mesa_mod_summary.f90
	gfortran -J.mod -I.mod -I.o -c mesa_mod_summary.f90 -o .o/mesa_mod_summary.o

mesa_mod_demo: .o/mesa_mod_reader.o .o/mesa_mod_demo.o
	gfortran -J.mod -I.mod -I.o .o/mesa_mod_reader.o .o/mesa_mod_demo.o -o mesa_mod_demo

mesa_mod_summary: .o/mesa_mod_reader.o .o/mesa_iso_mass.o .o/mesa_mod_summary.o 
	gfortran -J.mod -I.mod -I.o .o/mesa_mod_reader.o .o/mesa_iso_mass.o .o/mesa_mod_summary.o -o mesa_mod_summary

clean:
	-rm .o/*.*
	-rm .mod/*.*
