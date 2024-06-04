cp ../../../qxmd .
cp ${SOURCE_DIR}/util/qxREAD.py ${SOURCE_DIR}/Examples/03_NAQMD/01_Water_NVE/analysis/bond-angle/.
cp ${SOURCE_DIR}/util/qxREAD.py ${SOURCE_DIR}/Examples/03_NAQMD/01_Water_NVE/analysis/eig/.
cp ${SOURCE_DIR}/util/qxREAD.py ${SOURCE_DIR}/Examples/03_NAQMD/01_Water_NVE/analysis/rate/.
mkdir data
export PYTHONPATH=../../../util/
export LD_LIBRARY_PATH=../../../fftw-3.3.8/build/lib:${LD_LIBRARY_PATH}
#Get data for restart from optimized water molecule
cp ../../01_Optimization/01_Water/data/QM_* ./data/
cp ../../01_Optimization/01_Water/data/MD_* ./data/
./qxmd | tee log
