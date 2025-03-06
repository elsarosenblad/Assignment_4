CC = gcc
LD = gcc
CFLAGS = -O3 -g -Wall -Werror -march=native -funroll-loops -ffast-math -fopt-info-vec-optimized 
LDFLAGS = -lm -lpthread 
RM = /bin/rm -f
OBJS = galsim.o
EXECUTABLE = galsim

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(OBJS) -o $(EXECUTABLE) $(LDFLAGS)

galsim.o: galsim.c 
	$(CC) $(CFLAGS) -c galsim.c 

clean:
	$(RM) $(EXECUTABLE) $(OBJS) result.gal

