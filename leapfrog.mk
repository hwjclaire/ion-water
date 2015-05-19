 CXX=g++
 LD=$(CXX)

 INCLUDE =  -I/usr/include -I../lib -I../ewald

 OPT= -O3 -Wunused #-D_assert

 CXXFLAGS= $(OPT)  $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L/usr/lib64 

 LIBS =  -lmkl_core -lmkl_sequential -lmkl_gf_lp64 \
         -lm -lc -lgfortran -lgcc

 PLAT=Linux_x8664
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = /opt/netlib/blacs-ompi143-intel12
 BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-LINUX-0.a
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-LINUX-0.a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-LINUX-0.a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = /opt/netlib/scalapack-ompi143-intel12/lib
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

# Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

