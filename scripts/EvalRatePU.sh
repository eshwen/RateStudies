#./EvalRatePU 620 90 30
#./EvalRatePU 620 100 30
#../test/EvalRatePU 620 110 35 2017 /data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesHighPU_2017_fill6245/ ../fileLists/ZeroBias2017_fill6245.list ../utils/PU_per_LS_fill6245_2017.txt

./EvalRatePU.py 620 110 35 -y 2017 -o /data_CMS/cms/amendola/RateStudiesL1Ntuples/L1NtuplesHighPU_2017_fill6245/ -i ../fileLists/ZeroBias2017_fill6245.list --pu ../utils/PU_per_LS_fill6245_2017.txt
#./EvalRatePU 620 110 40

#./EvalRatePU 620 115 40   
