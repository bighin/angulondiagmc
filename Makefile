TARGET = diagmc
LIBS = -lncurses -lgsl -lm
CC = gcc-mp-6
#CC = clang-mp-4.0
#CC = /opt/local/libexec/llvm-5.0/libexec/ccc-analyzer
#CFLAGS = -DNDEBUG -O2 -Wall -fopenmp -std=c11 -I/opt/local/include/
CFLAGS = -Wall -pedantic -O2 -DNDEBUG -I/opt/local/include/
LDFLAGS = -L/opt/local/lib/

.PHONY: clean all default

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c libprogressbar/*.c inih/*.c njsummat/*.c murmurhash3/*.c cubature/hcubature.c))
HEADERS = $(wildcard *.h  libprogressbar/*.h inih/*.h njsummat/*.h murmurhash3/*.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o libprogressbar/*.o inih/*.o njsummat/*.o murmurhash3/*.o cubature/*.o
	-rm -f $(TARGET)
