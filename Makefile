CC = gcc
CXX = g++
CCFLAGS = -D_XOPEN_SOURCE -Wall -c -Wno-unused-variable -Wno-unused-function -Wno-unused-but-set-variable -std=c99 -fPIE
LDFLAGS = -lm

ifeq ($(IS),1)
    CCFLAGS += -D__USE_IS
endif

ifeq ($(ES),1)
    CCFLAGS += -D__USE_ES
endif

ifeq ($(RR),1)
    CCFLAGS += -D__USE_RR
endif

all: CCFLAGS += -O3 -fopenmp
all: LDFLAGS += -fopenmp
all: prepare PathTracer

debug: CCFLAGS += -g -O0 -D__DEBUG
debug: prepare PathTracer


prepare:
	@mkdir -p build
	
clean:
	rm -r build

PathTracer: build/PathTracer.o build/utils_path.o build/magic.o build/svdDynamic.o build/meshes.o
	$(CC) $^ $(LDFLAGS) -o $@

build/PathTracer.o: PathTracer.c buildScene.c utils_path.h
	$(CC) $(CCFLAGS) PathTracer.c -o $@

build/utils_path.o: utils_path.c utils_path.h
	$(CC) $(CCFLAGS) utils_path.c -o $@
	
build/magic.o: magic.c
	$(CC) $(CCFLAGS) $^ -o $@

build/svdDynamic.o: svdDynamic.c
	$(CC) $(CCFLAGS) $^ -o $@
	
build/meshes.o: meshes.c
	$(CC) $(CCFLAGS) $^ -o $@

.PHONY: all prepare clean debug
