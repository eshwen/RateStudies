#!/usr/bin/env python2
import os
import sys
from subprocess import check_output

commonStr = "/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/bundocka/ZeroBias/zbF_MET/180209_101449/0000/"
nFiles = int( check_output("ls -1 {0} | wc -l".format(commonStr), shell=True) ) - 1

fileName = "../fileLists/ZeroBias_fill6255_new.list"

for i in xrange(1, nFiles):
    fileIn = open(fileName, "w+")

    fileIn.write("{0}/L1Ntuple_{1}.root\n".format(commonStr, str(i)) )

    fileIn.close()
    
    print "Written entry", i
    
