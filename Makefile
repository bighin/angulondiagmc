TARGET = diagmc
LIBS = 
CC = gcc-mp-6
#CC = clang-mp-3.8
CFLAGS = -O2 -Wall -fopenmp -std=c11 -I/opt/local/include/
#CFLAGS = -Wall -g -pg
LDFLAGS = -lgsl -L/opt/local/lib/

.PHONY: clean all default

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c inih/*.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)
