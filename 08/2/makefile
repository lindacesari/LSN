CFLAGS = -Wall -O3 --std=c++11

main08.exe : main08.o random.h random.o func.h
	g++ random.o main08.o -o main08.exe
main08.o : main08.cpp func.h random.h
	g++ -c main08.cpp -o main08.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main08.exe
