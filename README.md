# RateStudies

The package contains a framework to produce ntuples with L1 information out of miniAOD events and codes for rate studies of L1 trigger seeds.

```
cmsrel CMSSW_8_0_20
cd CMSSW_8_0_20/src
cmsenv
git clone https://github.com/camendola/RateStudies
cd RateStudies
git checkout master
scram b
```

## Usage
### L1 Ntuples
The stored variables are set in L1ntuples/plugins/L1ntuples.cc.

To run interactively, define the input and output in L1ntuples/python/L1ntuples.py, then launch
```
cmsRun L1ntuples/python/L1ntuples.py
```
To submit on the T3:
```
cd L1ntuples/test
source MakeFileListDAS.sh -t "<fileListName>" -o <fileListName>.py -p <datasetName>
voms-proxy-init -voms cms   
./submitOnTier3_L1Ntuplizer.py -o <outputDirectory> -t <tag> -s <fileListName>.py
```
### Rate studies
All the information needs to be modified in the code EvalRate.cpp
```
c++ -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`
./EvalRate
```