ARB_INC=/data/local/include
ARB_LIB=/data/local/lib64
PS_INC=/home/madjp/primesieve-5.4.2/include
PS_LIB=/home/madjp/primesieve-5.4.2/.libs
CC=gcc
CFLAGS=-O2 -c -fPIC -I${ARB_INC} -I ../include -I ${PS_INC}
DEPS=../include/glfunc.h ../include/glfunc_internals.h
OBJ=glfunc.o g.o acb_fft.o error.o coeff.o buthe.o compute.o upsample.o zeros.o rank.o io.o special_values.o clear.o
all: lib

lib: $(OBJ)
	${CC} -shared -o ../libs/libglfunc.so ${OBJ} -L ${PS_LIB} -lprimesieve -L ${ARB_LIB} -larb

%.o: %.c ${DEPS}
	${CC} ${CFLAGS} $<


clean:
	rm -f *.o
	rm -f ../libs/libglfunc.so
