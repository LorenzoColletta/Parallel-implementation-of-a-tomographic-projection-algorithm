ifeq ($(FILE_BINARY),yes)
	BINARY = -DBINARY
else
	BINARY =
endif

CFLAGS=-Wall -Wpedantic -std=c99 -fopenmp -I./source/

all: inputgeneration projector

inputgeneration: inputgeneration.c ./bin/voxel.o ./source/voxel.h
	gcc ${CFLAGS} inputgeneration.c ./bin/voxel.o -lm -o inputgeneration

./bin/voxel.o: ./source/voxel.c ./source/voxel.h bin
	gcc ${CFLAGS} -c ./source/voxel.c -o ./bin/voxel.o

projector: projector.c ./bin/projection.o ./source/projection.h
	gcc ${CFLAGS} $(BINARY) projector.c ./bin/projection.o -lm -o projector

./bin/projection.o: ./source/projection.h ./source/projection.c bin
	gcc ${CFLAGS} -I./source/ -c ./source/projection.c -o ./bin/projection.o

bin:
	-mkdir ./bin

.PHONY:	clean

clean:
	-rm projector inputgeneration ./bin/*.o
