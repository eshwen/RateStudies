#!/usr/bin/env python
import os,sys
import argparse


if __name__ == "__main__":

    usage = 'usage: %prog [arguments]'
    parser = argparse.ArgumentParser(usage)
    parser.add_argument ('ptSingleL1', type=int, help ='pt_tau L1 cut for boosted ditau seed')
    parser.add_argument ('ptPairL1', type=int, help ='pt_tautau L1 cut for boosted ditau seed')
    parser.add_argument ('--ditaupt', dest='ptSingleL1ditau', type=int, help ='pt_tau L1 cut for ditau seed', default = 30)
    signal = parser.add_mutually_exclusive_group(required=True)
    signal.add_argument ('--VBF' , dest='VBFSignal' , help='VBF signal sample' , default=False, action="store_true")
    signal.add_argument ('--ggH' , dest='ggHSignal' , help='ggH signal sample' , default=False, action="store_true")
    signal.add_argument ('--Ztt' , dest='ZttSignal' , help='Ztt signal sample' , default=False, action="store_true" )
    signal.add_argument ('--HHbbtt' , dest='HHbbttSignal' , help='HHbbtt signal sample' , default=False ,action="store_true")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument ('--VBFtag', dest='VBFtag' , help='offline selection: VBF-tag' , default=False,action="store_true" )
    selection.add_argument ('--boosted', dest='boosted' , help='offline selection: 1jet boosted category' , default=False,action="store_true")

    args = parser.parse_args()
    

line='./OfflineSelection '+str(args.ptSingleL1)+' '+str(args.ptPairL1)+' '+str(args.ptSingleL1ditau)

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


if (args.HHbbttSignal==True):
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

    
print '{0}'.format(line)

os.system(line)
 
