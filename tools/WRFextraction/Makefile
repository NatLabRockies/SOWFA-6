.SUFFIXES: .F .o .f90

############ COMPILER SPECIFIC SETTINGS
#
############ uncomment this section for Intel Fortran
#FC      = ifort

############ uncomment this section for gfortran
#FC      = gfortran

############ uncomment this section for Intel OneAPI Fortran
FC      = ifx

# Also uncomment these next two lines and set ICCLIBPATH to correct path on your 
# system if using gfortran but with a netCDF library that was compiled with Intel C 
# This will get rid of undefined refs to intel_fast_memset, etc. The path shown
# is for the Redrock system at NREL, as of 20101230.
#ICCLIBPATH = /projects/global/x86_64/compilers/intel/intel-11.1-cproc-064/lib/intel64
#LIB = -L$(ICCLIBPATH) -lirc
# Leave the lines above commented out for NREL's Eagle
############ LIBRARY SETTINGS
#
# Set to the path that contains the netecdf lib and include directories on your system
#
#NETCDFPATH      =    /projects/nrel/apps/netcdf/4.0.1
#NETCDFPATH      =    /nopt/nrel/apps/netcdf/dist/netcdf
#NETCDFPATH      =    /usr
#NETCDFPATH       =    /nopt/nrel/apps/netcdf-c/4.7.4-intel/
# for fortran-compiled netcdf.
# Make sure the appropriate module is loaded in your environment
#    module use /nopt/nrel/apps/base/2020-05-12/spack/share/spack/modules/linux-centos7-x86_64/gcc-8.4.0
#    module load netcdf-fortran/4.5.2
NETCDFPATH       = /home/mchurchf/software-2024-03/spack/opt/spack/linux-rhel8-sapphirerapids/oneapi-2023.2.0/netcdf-fortran-4.5.2-ebrxi6jgymlvskxabuc5x4wa46ixxch6
#
############
#
LD      = $(FC)
CPP     = /lib/cpp -C -P
#FCFLAGS = -g -C
#FCFLAGS = -O3 -g
FCFLAGS = -O3 -cpp
LDFLAGS = $(FCFLAGS)

wrftoof : wrftoof.o module_dm.o module_openfoam_bc.o
	$(LD) -o wrftoof $(LDFLAGS) wrftoof.o module_dm.o module_openfoam_bc.o \
		-L$(NETCDFPATH)/lib -L$(NETCDFPATH)/lib64 -lnetcdff $(LIB)
	@ if [ -n "$(ICCLIBPATH)" ] ; then \
            echo ------------------------------------------------------------------------- ; \
            echo A D V I S O R Y : ;\
            echo "  " You may need to add: ;\
            echo "     " $(ICCLIBPATH) ;\
            echo "  " to LD_LIBRARY_PATH in your shell environment. ; \
            echo ------------------------------------------------------------------------- ; \
        fi

clean :
	/bin/rm -f *.o *.f90 *.mod

superclean : clean
	/bin/rm -f wrftoof


.F.o:
	$(CPP) $*.F > $*.f90
	$(FC) -o $@ -c $(FCFLAGS) -I$(NETCDFPATH)/include $*.f90

#
wrftoof.o : module_openfoam_bc.o module_dm.o

module_openfoam_bc.o : module_dm.o

