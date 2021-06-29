CFLAGS = -Wall -O3 --std=c++11

main10.2.exe : main10.2.o popolazione.o func.h popolazione.h
	mpicxx popolazione.o main10.2.o  -o main10.2.exe  $(CFLAGS) -larmadillo
main10.2.o : main10.2.cpp
	mpicxx -c main10.2.cpp -o main10.2.o $(CFLAGS)
popolazione.o : popolazione.cpp popolazione.h
	mpicxx -c popolazione.cpp -o popolazione.o $(CFLAGS)
clean :
	rm *.o main10.2.exe
