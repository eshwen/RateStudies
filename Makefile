CC=c++
all:
	$(CC) -lm -o EvalRate EvalRate.cpp `root-config --glibs --cflags`

	$(CC) -lm -o OfflineSelection OfflineSelection.cpp `root-config --glibs --cflags`
