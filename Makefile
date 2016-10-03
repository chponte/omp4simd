include config.make
export

all: mxm_openmp poisson_openmp md_openmp

mxm_openmp :
	cd src/MATMUL && ${MAKE}

poisson_openmp :
	cd src/POISSON && ${MAKE}

md_openmp :
	cd src/MD && ${MAKE}

clean:
	cd src/MATMUL && ${MAKE} clean
	cd src/POISSON && ${MAKE} clean
	cd src/MD && ${MAKE} clean
