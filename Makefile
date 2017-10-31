TARGET = diagmc
LIBS = -lncurses -lgsl -lm
CC = gcc-mp-6
CC = clang-mp-5.0
#CC = /opt/local/libexec/llvm-5.0/libexec/ccc-analyzer
CFLAGS = -DNDEBUG -O2 -Wall -fopenmp -std=c11 -I/opt/local/include/
CFLAGS =  -Wall -pedantic -O0 -g -I/opt/local/include/ -fopenmp
LDFLAGS = -L/opt/local/lib/

# Linux-specific flags
#CC=gcc
#CFLAGS+=-std=gnu99
#LDFLAGS+=-lcblas

.PHONY: clean all default

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c libprogressbar/*.c inih/*.c njsummat/*.c murmurhash3/*.c cubature/hcubature.c gnuplot_i/*.c))
HEADERS = $(wildcard *.h  libprogressbar/*.h inih/*.h njsummat/*.h murmurhash3/*.h gnuplot_i/*.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o libprogressbar/*.o inih/*.o njsummat/*.o murmurhash3/*.o cubature/*.o gnuplot_i/*.o
	-rm -f $(TARGET)
