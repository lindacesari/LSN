CFLAGS = -Wall -O3 --std=c++11

main05.1.exe : main05.1.o random.h random.o func.h
	g++ random.o main05.1.o -o main05.1.exe
main05.1.o : main05.1.cpp func.h random.h
	g++ -c main05.1.cpp -o main05.1.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main05.1.exe
