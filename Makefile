include config.make
export

all: mxm_openmp poisson_openmp md_openmp

outdir:
	mkdir -p ${ODIR}

mxm_openmp: outdir
	cd src/MATMUL && ${MAKE}

poisson_openmp: outdir
	cd src/POISSON && ${MAKE}

md_openmp: outdir
	cd src/MD && ${MAKE}

clean:
	cd src/MATMUL && ${MAKE} clean
	cd src/POISSON && ${MAKE} clean
	cd src/MD && ${MAKE} clean
	rm -rf ${ODIR}