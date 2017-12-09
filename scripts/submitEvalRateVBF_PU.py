import os
import sys
import argparse
import datetime



def parseInputFileList (fileName) :
    filelist = []
    with open (fileName) as fIn:
        for line in fIn:
            line = (line.split("#")[0]).strip()
            if line:
                filelist.append(line)
    return filelist

usage = 'usage: %prog [arguments]'
parser = argparse.ArgumentParser(usage)
parser.add_argument('--njobs',  dest='njobs', type=int, help='number of jobs for parallelization',  default=10)
parser.add_argument('--fill',  dest='fill', type=int, help='fill number', default=None)
parser.add_argument('--run',  dest='run', type=int, help='run number', default=None)
parser.add_argument('--mjj',  dest='mjj', type=int, help='mjj threshold', default=None)
parser.add_argument('--lead',  dest='lead', type=int, help='lead threshold', default=None)
parser.add_argument('--sub',  dest='sub', type=int, help='sub threshold', default=None)
parser.add_argument('--tau',  dest='tau', type=int, help='tau threshold', default=None)
parser.add_argument('--trID',  dest='trID', type=int, help='triggerID', default=0)
parser.add_argument('--emu',  dest='emu', type=int, help='emulated events', default=0)
parser.add_argument('--filelist','-i',  dest='filelist', help='filelist',  default=False)
parser.add_argument('--output','-o',  dest='output', help='output directory', default=False)
parser.add_argument('--pu',  dest='pu', help='PU per LS file', default=False)
parser.add_argument('--maxevents', '-m', dest='m', type=int, help='maximum events', default=-999)
parser.add_argument('--tag',  dest='tag', help='tag',  default=None)

process = parser.add_mutually_exclusive_group(required=True)
process.add_argument('--VBFincl',  dest='VBFincl', help='inclusive VBF seed',  default=False, action="store_true")
process.add_argument('--VBFtau',  dest='VBFtau', help='VBF+tau seed',  default=False, action="store_true")

args = parser.parse_args()


# datetime.datetime.now()
outDir = datetime.datetime.now().strftime('%Y.%m.%d_%H.%M.%S')
if args.tag:
    outDir = args.tag
outDir = "jobs_"+outDir


here = os.getcwd()

if(args.VBFincl):
    program = './EvalRatePU'
elif(args.VBFtau):
    program = './EvalRatePU_tau'
    
proto = 'job_' ## job .sh fie name
logproto = 'log_' ## job .sh fie name

os.system('mkdir ' + outDir)

if(args.VBFincl):
    outFolder = args.output +'/VBFincl_' + str(args.mjj) + '_' +str(args.lead) + '_' + str(args.sub)    
elif(args.VBFtau):
    outFolder = args.output +'/VBFtau_' + str(args.mjj) + '_' +str(args.sub) + '_' + str(args.tau)    

if args.emu == 1:
    outFolder = outFolder + '_emu'

inputfiles = parseInputFileList (args.filelist)
if args.njobs > len (inputfiles) : args.njobs = len (inputfiles)
nfiles = (len (inputfiles) + len (inputfiles) % args.njobs) / args.njobs
inputlists = [inputfiles[x:x+nfiles] for x in xrange (0, len (inputfiles), nfiles)]

os.system('mkdir ' + outFolder)
nj = 0

for listname in inputlists:

    listFileName = "filelist_%i.txt" % nj
    thisinputlistFile = open(outDir + "/" + listFileName, 'w')
    for line in listname:
        thisinputlistFile.write(line+"\n")
    thisinputlistFile.close()
    scriptName = proto + str(nj) + '.sh'
    logName    = logproto + str(nj) + '.txt'
    scriptFile = open (outDir + '/' + scriptName, 'w')
    scriptFile.write ('#!/bin/bash\n')
    scriptFile.write ('export X509_USER_PROXY=~/.t3/proxy.cert\n')
    scriptFile.write ('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
    scriptFile.write ('cd %s\n' % here)
    scriptFile.write ('export SCRAM_ARCH=slc6_amd64_gcc491\n')
    scriptFile.write ('eval `scram r -sh`\n')
    if (args.VBFincl):
        command = program + ' ' +str(args.mjj) + ' ' +str(args.lead) + ' ' + str(args.sub)
    elif(args.VBFtau):
        command = program + ' ' +str(args.mjj) + ' ' +str(args.sub) + ' ' + str(args.tau)
    command = command + ' ' + outFolder + ' ' + outDir+'/'+listFileName  + ' ' + str(args.pu) + ' ' + str(args.trID) + ' ' +str(args.emu) + ' '+ str(args.run) + ' '+ str(args.fill) +' '+  str(args.m)
    command = command + ' ' + str(nj) + ' ' + str(args.njobs) + ' 2>&1 | tee ' + outDir + '/' + logName
    scriptFile.write(command)
    scriptFile.close()
    
    os.system ('chmod u+rwx ' + outDir + '/' + scriptName)
    launchcommand = ('/opt/exp_soft/cms/t3/t3submit -short \'' + outDir + '/' + scriptName +"\'")
    print launchcommand
    os.system (launchcommand)
    nj = nj + 1
   # command = '/opt/exp_soft/cms/t3/t3submit -short ' + outDir + '/' + proto + str (nj) + '.sh'
