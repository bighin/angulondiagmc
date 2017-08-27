TARGET = diagmc
LIBS = 
CC = gcc-mp-6
#CC = clang-mp-4.0
#CFLAGS = -DNDEBUG -O2 -Wall -fopenmp -std=c11 -I/opt/local/include/
CFLAGS = -Wall -pedantic -O2 -g -pg -DNDEBUG -I/opt/local/include/
LDFLAGS = -lncurses -lgsl -L/opt/local/lib/

.PHONY: clean all default

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c libprogressbar/*.c inih/*.c njsummat/*.c))
HEADERS = $(wildcard *.h  libprogressbar/*.h inih/*.h njsummat/*.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o libprogressbar/*.o inih/*.o njsummat/*.o
	-rm -f $(TARGET)
