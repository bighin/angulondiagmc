TARGET = diagmc
LIBS = -lncurses -lgsl -lm
#CC = gcc
CC = clang-mp-5.0
#CC = /opt/local/libexec/llvm-5.0/libexec/ccc-analyzer
CFLAGS = -DNDEBUG -O2 -Wall -std=gnu11 -I/opt/local/include/
#CFLAGS =  -Wall -std=gnu11 -O0 -g -I/opt/local/include/
LDFLAGS = -L/opt/local/lib/

# Linux-specific flags
#CC=gcc
#CFLAGS+=-std=gnu99
#LDFLAGS+=-lcblas

.PHONY: clean all default

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c libprogressbar/*.c inih/*.c njsummat/*.c murmurhash3/*.c cubature/hcubature.c gnuplot_i/*.c rbitree/*.c))
HEADERS = $(wildcard *.h  libprogressbar/*.h inih/*.h njsummat/*.h murmurhash3/*.h gnuplot_i/*.h rbitree/*.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -Wall $(LIBS) -o $@

libangulon: libangulon.a

libangulon.a: $(OBJECTS)
	ar rcs $@ $(OBJECTS)

clean:
	-rm -f *.o libprogressbar/*.o inih/*.o njsummat/*.o murmurhash3/*.o cubature/*.o gnuplot_i/*.o rbitree/*.o
	-rm -f $(TARGET)
