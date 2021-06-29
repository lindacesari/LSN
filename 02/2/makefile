CFLAGS = -Wall -O3 --std=c++11

main02.2.exe : main02.2.o random.h random.o func.h  
	g++ random.o main02.2.o -o main02.2.exe
main02.2.o : main02.2.cpp func.h random.h  
	g++ -c main02.2.cpp -o main02.2.o $(CFLAGS)
random.o : random.cpp random.h  
	g++ -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main02.2.exe
