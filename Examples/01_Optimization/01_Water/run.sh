cp ../../../qxmd .
export PYTHONPATH=../../../util/
export LD_LIBRARY_PATH=../../../fftw-3.3.8/build/lib:${LD_LIBRARY_PATH}
./qxmd | tee log
# fsdfgsdgsdh
