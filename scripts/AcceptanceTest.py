#!/usr/bin/env python
import os,sys
import argparse


if __name__ == "__main__":

    usage = 'usage: %prog [arguments]'
    parser = argparse.ArgumentParser(usage)
    parser.add_argument ('ptSingleL1', type=int, help ='pt_tau L1 cut for boosted ditau seed')
    parser.add_argument ('ptPairL1', type=int, help ='pt_tautau L1 cut for boosted ditau seed')
    parser.add_argument ('--ditaupt', dest='ptSingleL1ditau', type=int, help ='pt_tau L1 cut for ditau seed', default = 30)
    parser.add_argument ('--offptpair', dest='ptPairOffBoost', type=int, help ='pt_tautau offline cut for boosted ditau seed', default = 0)
    signal = parser.add_mutually_exclusive_group(required=True)
    signal.add_argument ('--VBF' , dest='VBFSignal' , help='VBF signal sample' , default=False, action="store_true")
    signal.add_argument ('--ggH' , dest='ggHSignal' , help='ggH signal sample' , default=False, action="store_true")
    signal.add_argument ('--Ztt' , dest='ZttSignal' , help='Ztt signal sample' , default=False, action="store_true" )
#    signal.add_argument ('--HHbbtt' , dest='HHbbttSignal' , help='HHbbtt signal sample' , default=False ,action="store_true")
    signal.add_argument ('--VBFHHbbtt' , dest='VBFHHbbttSignal' , help='VBFHHbbtt signal sample' , default=False ,action="store_true")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument ('--VBFtag', dest='VBFtag' , help='offline selection: VBF-tag' , default=False,action="store_true" )
    selection.add_argument ('--boosted', dest='boosted' , help='offline selection: 1jet boosted category' , default=False,action="store_true")
    parser.add_argument ('--ditauptOR', dest='ptSingleL1ditauOR', type=int, help ='pt_tau L1 cut for ditau seed complementary to boosted', default = 30)
    parser.add_argument ('--emu', dest='emu' , help='L1 emulated' , default=False,action="store_true")
    parser.add_argument ('--VBFsubL1', dest = 'VBFsubL1',type=int, help ='subleading jet L1 cut for VBF seed',default = 30)
    parser.add_argument ('--VBFleadL1',dest ='VBFleadL1', type=int, help ='leading jet L1 cut for VBF seed',default =90)
    parser.add_argument ('--VBFMjjL1', dest = 'VBFMjjL1',type=int, help ='Mjj L1 cut for VBF seed',default = 620)
    args = parser.parse_args()
    
####emu
line='./../test/OfflineSelectionL1_bbtt '+str(args.ptSingleL1)+' '+str(args.ptPairL1)+' '+str(args.ptSingleL1ditau)

if (args.VBFSignal==True):
    line+=' 1'
else:
    line+=' 0'


if (args.ggHSignal==True):
    line+=' 1'
else:
    line+=' 0'

    
if (args.ZttSignal==True):
    line+=' 1'
else:
    line+=' 0'


if (args.VBFHHbbttSignal==True):
    line+=' 1'
else:
    line+=' 0'


if (args.VBFtag==True):
    line+=' 500 3.5'
else:
    line+=' 0 0'

if (args.boosted==True or args.VBFtag==True):
    line+=' 1'
else:
    line+=' 0'

line+=' '+str(args.ptPairOffBoost)+' '+str(args.ptSingleL1ditauOR)

#### new for emu

if(args.emu==True):
    line+=' 1'
else: 
    line+=' 0'

line+=' '+str(args.VBFsubL1)+' '+str(args.VBFleadL1)+' '+str(args.VBFMjjL1)
    
print '{0}'.format(line)

os.system(line)
 
