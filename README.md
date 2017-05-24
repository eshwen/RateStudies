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
source /opt/exp_soft/cms/t3/t3setup
voms-proxy-init -voms cms   
source MakeFileListDAS.sh -t "<fileListName>" -o <fileListName>.py -p <datasetName>
./submitOnTier3_L1Ntuplizer.py -o <outputDirectory> -t <tag> -s <fileListName>.py
```
### Rate studies
All the information needs to be modified in the code EvalRate.cpp
```
make eval
./EvalRate
```

### Acceptance studies
All the information needs to be modified in the code OfflineSelection.cpp
```
#compile
make off
#launch through AcceptanceTest.py script
python AcceptanceTest.py <ptTauL1BoostedDiSeed> <ptTauPairL1BoostedDiSeed> % options
```

usage:
```
positional arguments:
 ptSingleL1            pt_tau L1 cut for boosted ditau seed
 ptPairL1              pt_tautau L1 cut for boosted ditau seed

optional arguments:
  -h, --help                 show this help message and exit
  --ditaupt PTSINGLEL1DITAU  pt_tau L1 cut for ditau seed
  --offptpair PTPAIROFFBOOST pt_tautau offline cut for boosted ditau seed
  --VBFsubL1 VBFSUBL1        subleading jet L1 cut for VBF seed
  --VBFleadL1 VBFLEADL1	     leading jet L1 cut for VBF seed
  --VBFMjjL1 VBFMJJL1        Mjj L1 cut for VBF seed
signal sample (mandatory argument; exclusive group)
    --VBF                      VBF signal sample
    --ggH                      ggH signal sample
    --Ztt                      Ztt signal sample
    --HHbbtt                   HHbbtt signal sample
 --ditauptOR PTSINGLEL1DITAUOR	 pt_tau L1 cut for ditau seed complementary to boosted	
offline analyis selections (optional argument; exclusive group)
    --VBFtag                   offline selection: VBF-tag
    --boosted                  offline selection: 1jet boosted category

    --emu                 bool: L1 emulated										

```


positional arguments:
  ptSingleL1            pt_tau L1 cut for boosted ditau seed
    ptPairL1              pt_tautau L1 cut for boosted ditau seed

optional arguments:
  -h, --help            show this help message and exit
   
   			
