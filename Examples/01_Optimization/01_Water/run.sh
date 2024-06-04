cp ../../../qxmd .
cp ${SOURCE_DIR}/util/qxREAD.py ${SOURCE_DIR}/Examples/01_Optimization/01_Water/analysis/DOS/.
cp ${SOURCE_DIR}/util/qxREAD.py ${SOURCE_DIR}/Examples/01_Optimization/01_Water/analysis/eng/.
cp ${SOURCE_DIR}/util/qxREAD.py ${SOURCE_DIR}/Examples/01_Optimization/01_Water/analysis/pdb/.
export PYTHONPATH=../../../util/
export LD_LIBRARY_PATH=../../../fftw-3.3.8/build/lib:${LD_LIBRARY_PATH}
mkdir data
./qxmd | tee log
