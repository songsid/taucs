# These are for a Pentium4 version of ATLAS (not bundled with taucs)
LDFLAGS   = -Wl,-rpath /specific/a/home/cc/cs/stoledo/Public

#LIBBLAS   =   -L /specific/a/home/cc/cs/stoledo/Public -lgoto_p4_512-r0.9 -L external/lib/linux/ -lcblas 

#LIBBLAS   = -L /home/stoledo/Public -lgoto_p3_256-r0.9 \
#            
#LIBLAPACK =  external/lib/linux/blas_aux.o /specific/a/home/cc/cs/stoledo/Public/lapack_linux.a

LIBLAPACK = -L external/lib/linux -llapackStoledo
LIBBLAS = -L external/lib/linux  -lf77blas -lcblas -latlas -lf2c -lm
CFLAGS = -g
#########################################################







