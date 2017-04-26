CC=c++
all:
	$(CC) -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`

	$(CC) -lm -o OfflineSelection OfflineSelection.cpp `root-config --glibs --cflags`


off:

	$(CC) -lm -o OfflineSelection OfflineSelection.cpp `root-config --glibs --cflags`

eval:
	$(CC) -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`

old:
	$(CC) -lm -o EvalRate EvalRateOld.cpp `root-config --glibs --cflags`
