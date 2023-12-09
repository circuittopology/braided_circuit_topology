CC = gcc
CFLAGS = -O0 -g -Wall -std=c99
LIBS = -lm -lgsl -lgslcblas

all: CT_analysis.exe initial_condition.exe reshape.exe CT_analysis_batch.exe

%.exe: %.c
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.exe *.o core.*
