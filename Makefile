CC=icc
CFLAGS=-std=c++11	

all:
		$(CC) $(CFLAGS) -o  neighbor main.cpp

clean:
		rm -rf neighbor
