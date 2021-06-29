CFLAGS = -Wall -O3 --std=c++11

main10.1.exe : main10.1.o popolazione.o popolazione.h func.h
	g++ popolazione.o main10.1.o -o main10.1.exe -larmadillo
main10.1.o : main10.1.cpp func.h popolazione.h
	g++ -c main10.1.cpp -o main10.1.o $(CFLAGS)
popolazione.o : popolazione.cpp popolazione.h
	g++ -c popolazione.cpp -o popolazione.o $(CFLAGS)
clean :
	rm *.o main10.1.exe
