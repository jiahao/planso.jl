#
# Definitions of commands and arguements used in generating the library
#
# CRAY T3E flags
#CC=cc
#FC=f90
#MPIFC=f90
#CFLAGS=-O3
#FFLAGS=-O3 -dp
#NOOPT=-O1 -dp
#LIBS=-L$(HOME)/lib  -lpspars -lskit
#LDFLAGS=-Wl"-D equiv(DCOPY)=SCOPY;equiv(DDOT)=SDOT;equiv(DAXPY)=SAXPY;equiv(IDAMAX)=ISAMAX;equiv(DSCAL)=SSCAL;equiv(DSQRT)=SQRT;equiv(DNRM2)=SNRM2;equiv(DSWAP)=SSWAP;equiv(DGEMV)=SGEMV;equiv(DGEMM)=SGEMM;equiv(DLAMCH)=SLAMCH;equiv(DSYEV)=SSYEV"
#
# Sparc/Solaris flags
#
CC=cc
FC=f77
MPIFC=mpif77
CFLAGS=-xO3
FFLAGS=-O4
NOOPT=-O2
#CFLAGS=-g
#FFLAGS=-g
#NOOPT=-g -O0
LIBS=-L/afs/nersc.gov/projects/PEIG/lib/sol2  -lpspars -lskit -lblas
LDFLAGS=-z muldefs
#
# SGI r1000 SMP floags
#
#CC=cc
#FC=f77
#MPIFC=mpif77
#CFLAGS=-64 -O
#FFLAGS=-64 -O
#NOOPT=-64 -O1
#LIBS=/project/darpa2/kewu/Shop/PSP2_15/libpspars.a /project/darpa/lib/irix64/libskit.a -lcomplib.sgimath -lmalloc
#LDFLAGS=
#
# Sun SMP HPC environment
#
#CC=cc
#FC=f77
#MPIFC=tmf77
#CFLAGS=-O
#FFLAGS=-O
#NOOPT=-O1
#LIBS=-L${HOME}/lib  -lpspars -lskit -lblas
#LDFLAGS=-z muldefs -lmpi
#
AR=ar
RANLIB=ranlib
#
CD=cd
RM=rm
ECHO=echo
#
