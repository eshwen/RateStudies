#!/usr/bin/env python
import os,sys
import argparse


if __name__ == "__main__":

    usage = 'usage: %prog [arguments]'
    parser = argparse.ArgumentParser(usage)
    parser.add_argument ('mjj', type=int, help ='mjj threshold')
    parser.add_argument ('leadjet', type=int, help ='leadjet threshold')
    parser.add_argument ('subjet',  type=int, help ='subjet threshold')
    args = parser.parse_args()
    optparser = optparse.OptionParser(usage)
    optparser.add_option ('-y',   '--year',   help='year' , default='2017', action="store_true")
    optparser.add_option ('-o',   '--output', help='output directory' , default =  None, action="store_true")
    optparser.add_option ('-i',   '--inp',    help='input filelist' , default =  None, action="store_true")
    optparser.add_option ('--pu',             help='pu per lumi file' , default =  None, action="store_true")

    
line='../test/EvalRatePU '+str(args.mjj)+' '+str(args.leadjet)+' '+str(args.subjet)+' '+str(args.year)+' '+str(args.output)+' '+str(args.inp)+' '+str(args.pu)     
 
