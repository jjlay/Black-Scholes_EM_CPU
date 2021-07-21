CC = g++-10
CFLAGS = -std=c++20 -g

sde : sde.o
	$(CC) $(CFLAGS) -o sde.out sde.o

sde.o : sde.cpp
	$(CC) $(CFLAGS) -c sde.cpp

clean :
	rm -f *.o
	rm -f sde.out
