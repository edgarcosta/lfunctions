ARB_INC=/data/local/include
ARB_LIB=/data/local/lib64
PS_INC=/home/madjp/primesieve-5.4.2/include
PS_LIB=/home/madjp/primesieve-5.4.2/.libs
L_INC=../include
L_LIB=../libs
all: dir_test dir_test1

dir_test: dir_test.c
	gcc -O2 -o dir_test dir_test.c -I ${ARB_INC} -I ${L_INC} -L ${L_LIB} -lglfunc -L ${ARB_LIB} -larb

dir_test1: dir_test1.c
	gcc -O2 -o dir_test1 dir_test1.c -I ${ARB_INC} -I ${PS_INC} -I ${L_INC} -L ${L_LIB} -lglfunc -L ${ARB_LIB} -larb -L ${PS_LIB} -lprimesieve

clean:
	rm -f dir_test dir_test1
