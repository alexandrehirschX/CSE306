CC = g++
CFLAGS  = -g -Wall


render: renderer.o  
	$(CC) $(CFLAGS) -o renderer renderer.o


renderer.o: renderer.cpp
	$(CC) $(CFLAGS) -O3 -c renderer.cpp 
