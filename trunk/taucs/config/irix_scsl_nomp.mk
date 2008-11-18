#########################################################
# Irix                                                  #
#########################################################
OBJEXT=.o
LIBEXT=.a
EXEEXT= 
F2CEXT=.f
PATHSEP=/
DEFFLG=-D

# Strange: I couldn't get SCSL to work with CC... 
# So I had to use GNU's compilers (Haim)

FC        = g77
FFLAGS    = -O3 -fno-second-underscore -Wall
FOUTFLG   = -o ./

CC        = gcc
COUTFLG   = -o ./
CFLAGS    = -O3 -Wall -Werror -std=c89

LD        = $(CC) 
LDFLAGS   = $(CFLAGS)
LOUTFLG   = $(COUTFLG)

AR        = ar cr
AOUTFLG   =

RANLIB    = echo
RM        = /bin/rm -rf

# We use SCSL without it's multiprocessing ability
LIBBLAS   = -lscs

LIBLAPACK = 
LIBMETIS  = -L external/lib/irix -lmetis32

LIBF77 = -L /usr/local/gcc3/lib/gcc-lib/mips-sgi-irix6.5/3.0.4/ -lg2c -lgcc
LIBC = -L /usr/lib32/mips3 -lc -lm 


#########################################################
