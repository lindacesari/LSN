CFLAGS = -Wall -O3 --std=c++11

main01.1.exe : main01.1.o random.o func.h
	g++ random.o main01.1.o -o main01.1.exe
main01.1.o : main01.1.cpp func.h
	g++ -c main01.1.cpp -o main01.1.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main01.1.exe
