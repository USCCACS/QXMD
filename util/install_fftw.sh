fftw_version=3.3.8
wget https://fftw.org/pub/fftw/fftw-${fftw_version}.tar.gz
tar -xzvf fftw-${fftw_version}.tar.gz
cd ./fftw-${fftw_version}/
mkdir build
./configure --prefix=/content/QXMD/fftw-${fftw_version}/build CC=gcc CFLAGS=-O3 MPICC=mpicc F77=mpif77 F90=mpif90 FFLAGS=-O3 --enable-shared=yes --enable-sse2 --enable-avx --enable-avx2 --disable-avx512 --disable-avx-128-fma --disable-kcvi --disable-altivec --disable-vsx --disable-neon --disable-generic-simd128 --disable-generic-simd256 --disable-fma --enable-static=yes --enable-mpi --enable-openmp --enable-threads
make -j 8
make install