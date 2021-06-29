CFLAGS = -Wall -O3 --std=c++11

main01.2.exe : main01.2.o random.o func.h
	g++ random.o main01.2.o -o main01.2.exe
main01.2.o : main01.2.cpp func.h
	g++ -c main01.2.cpp -o main01.2.o $(CFLAGS)
random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main01.2.exe seed.out
