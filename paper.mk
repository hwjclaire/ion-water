 CXX=g++
 LD=$(CXX)

 INCLUDE =  -I/usr/include -I../lib -I../ewald

 OPT= -O3 -Wunused -DASSERT

 CXXFLAGS= $(OPT)  $(INCLUDE) $(PLTFLAGS) $(DFLAGS)
 GCCDIR= /share/apps/gcc/gcc-4.5.2-build1

 LIBPATH = -L$(GCCDIR)/lib64

 LIBS =  -lmkl_core -lmkl_sequential -lmkl_gf_lp64 \
         -lm -lc -lgfortran -lgcc

 PLAT=Linux_x8664
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = /share/apps/netlib/BLACS-OMPI-gcc4.5/LIB
 BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-LINUX-0.a
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-LINUX-0.a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-LINUX-0.a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = /share/apps/netlib/scalapack_mkl_blacs_gcc4.5/lib
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

# Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)
