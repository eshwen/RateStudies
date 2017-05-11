CC=c++
all:
	$(CC) -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`

	$(CC) -lm -o OfflineSelection OfflineSelection.cpp `root-config --glibs --cflags`


off:

	$(CC) -lm -o OfflineSelection OfflineSelection.cpp `root-config --glibs --cflags`

eval:
	$(CC) -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`

pl:
	$(CC) -lm -o MakePlots MakePlots.cpp `root-config --glibs --cflags`

evaltak:
	$(CC) -lm -o EvalRateL1 EvalRateL1Ntuples.cpp `root-config --glibs --cflags`

test:
	$(CC) -lm -o RateTest RateTest.cpp `root-config --glibs --cflags`
