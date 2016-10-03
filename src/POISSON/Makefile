all:
	${CC} ${CL} ${CFLAGS} -o ${ODIR}/poisson_openmp poisson_openmp.c

clean:
	rm -rf $(shell ls | grep -v Makefile | grep -v poisson_openmp.c)
