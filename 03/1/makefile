CFLAGS = -Wall -O3 --std=c++11

main03.1.exe : main03.1.o random.h random.o func.h
	g++ random.o main03.1.o -o main03.1.exe
main03.1.o : main03.1.cpp func.h random.h
	g++ -c main03.1.cpp -o main03.1.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main03.1.exe
