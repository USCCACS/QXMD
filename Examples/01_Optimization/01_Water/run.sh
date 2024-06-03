cp ../../../qxmd .
cp /content/QXMD/util/qxREAD.py /content/QXMD/Examples/01_Optimization/01_Water/analysis/DOS/.
cp /content/QXMD/util/qxREAD.py /content/QXMD/Examples/01_Optimization/01_Water/analysis/eng/.
cp /content/QXMD/util/qxREAD.py /content/QXMD/Examples/01_Optimization/01_Water/analysis/pdb/.
export PYTHONPATH=../../../util/
export LD_LIBRARY_PATH=../../../fftw-3.3.8/build/lib:${LD_LIBRARY_PATH}
mkdir data
./qxmd | tee log
