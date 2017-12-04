CC=c++
TEST = test
all:
	$(CC) -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`

	$(CC) -lm -o OfflineSelection OfflineSelection.cpp `root-config --glibs --cflags`


off:

	$(CC) -lm -o test/OfflineSelectionL1 test/OfflineSelectionL1.cpp `root-config --glibs --cflags`
offbbtt:

	$(CC) -lm -o test/OfflineSelectionL1_bbtt test/OfflineSelectionL1_bbtt.cpp `root-config --glibs --cflags`

eval:
	$(CC) -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`

ploff:
	$(CC) -lm -o MakePlotsOff MakePlotsOff.cpp `root-config --glibs --cflags`

pl:
	$(CC) -lm -o MakePlots MakePlots.cpp `root-config --glibs --cflags`
evaltak:
	$(CC) -lm -o EvalRateL1 EvalRateL1Ntuples.cpp `root-config --glibs --cflags`
ratepu:
	$(CC) -lm -o EvalRatePU test/EvalRatePU.cpp `root-config --glibs --cflags`
VBF:
	$(CC) -lm -o EvalRateVBF test/EvalRateVBF.cpp `root-config --glibs --cflags`
