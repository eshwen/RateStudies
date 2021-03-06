#fill 6245 triggersID: https://cmswbm.cern.ch/cmsdb/servlet/L1Summary?RUN=303948&KEY=l1_trg_collisions2017/v99

#### VBF INCLUSIVE
### unpacked
#python scripts/submitEvalRateVBF.py --njobs 200 --mjj 620 --lead 90 --sub 30 --filelist fileLists/ZeroBias2017_fill6245.list -o /data_CMS/cms/amendola/Rate2017_VBF/HighPUFill --run 303948 --fill 6245 --trID 276 --tag VBF --pu utils/PU_per_LS_fill6245_2017.txt --VBFincl
#python scripts/submitEvalRateVBF.py --njobs 20 --mjj 620 --lead 100 --sub 35 --filelist fileLists/ZeroBias2017_fill6245.list -o /data_CMS/cms/amendola/Rate2017_VBF/HighPUFill --run 303948 --fill 6245 --trID 278 --tag VBF2 --pu utils/PU_per_LS_fill6245_2017.txt --VBFincl
python scripts/submitEvalRateVBF.py --njobs 100 --mjj 620 --lead 110 --sub 35 --filelist fileLists/ZeroBias2017_fill6245.list -o /data_CMS/cms/amendola/Rate2017_VBF/HighPUFill --run 303948 --fill 6245 --trID 279 --tag VBF3 --pu utils/PU_per_LS_fill6245_2017.txt --VBFincl --lmin 71 --lmax 132
#python scripts/submitEvalRateVBF.py --njobs 20 --mjj 620 --lead 115 --sub 35 --filelist fileLists/ZeroBias2017_fill6245.list -o /data_CMS/cms/amendola/Rate2017_VBF/HighPUFill --run 303948 --fill 6245 --trID 281 --tag VBF4 --pu utils/PU_per_LS_fill6245_2017.txt --VBFincl
#python scripts/submitEvalRateVBF.py --njobs 20 --mjj 620 --lead 115 --sub 40 --filelist fileLists/ZeroBias2017_fill6245.list -o /data_CMS/cms/amendola/Rate2017_VBF/HighPUFill --run 303948 --fill 6245 --trID 282 --tag VBF5 --pu utils/PU_per_LS_fill6245_2017.txt --VBFincl


#### VBF+TAU
#python scripts/submitEvalRateVBF.py --njobs 100 --mjj 450 --sub 35 --tau 45 --filelist fileLists/ZeroBias2017_fill6245.list -o /data_CMS/cms/amendola/Rate2017_VBF/HighPUFill --run 303948 --fill 6245 --trID 286 --tag VBF --pu utils/PU_per_LS_fill6245_2017.txt --VBFtau
#python scripts/submitEvalRateVBF.py --njobs 100 --mjj 450 --sub 40 --tau 45 --filelist fileLists/ZeroBias2017_fill6245.list -o /data_CMS/cms/amendola/Rate2017_VBF/HighPUFill --run 303948 --fill 6245  --tag VBFtau1 --pu utils/PU_per_LS_fill6245_2017.txt --VBFtau 
#python scripts/submitEvalRateVBF.py --njobs 100 --mjj 450 --sub 45 --tau 45 --filelist fileLists/ZeroBias2017_fill6245.list -o /data_CMS/cms/amendola/Rate2017_VBF/HighPUFill --run 303948 --fill 6245  --tag VBFtau2 --pu utils/PU_per_LS_fill6245_2017.txt --VBFtau
#python scripts/submitEvalRateVBF.py --njobs 100 --mjj 450 --sub 45 --tau 50 --filelist fileLists/ZeroBias2017_fill6245.list -o /data_CMS/cms/amendola/Rate2017_VBF/HighPUFill --run 303948 --fill 6245  --tag VBFtau3 --pu utils/PU_per_LS_fill6245_2017.txt --VBFtau
