TARGET = diagmc
LIBS = 
CC = gcc-mp-6
#CC = clang-mp-3.8
#CFLAGS = -DNDEBUG -O2 -Wall -fopenmp -std=c11 -I/opt/local/include/
CFLAGS = -Wall -g -pg -I/opt/local/include/
LDFLAGS = -lncurses -lgsl -lf2c -L/opt/local/lib/

.PHONY: clean all default

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c libprogressbar/*.c inih/*.c njgraf/*.c))
HEADERS = $(wildcard *.h  libprogressbar/*.h inih/*.h njgraf/*.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o libprogressbar/*.o inih/*.o njgraf/*.o
	-rm -f $(TARGET)
