 CXX=icpc
 LD=$(CXX)

 INCLUDE =  -I/usr/include -I../lib -I../ewald -I$(MKLROOT)/include -I~/apps/fftw3/include

 OPT= -O3  -DASSERT -DADD_ -openmp -qopt-report=5 -qopt-report-phase=vec

 CXXFLAGS= $(OPT)  $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm
         

