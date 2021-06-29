CFLAGS = -Wall -O3 --std=c++11

main09.1.exe : main09.1.o popolazione.o popolazione.h func.h
	g++ popolazione.o main09.1.o -o main09.1.exe -larmadillo
main09.1.o : main09.1.cpp func.h popolazione.h
	g++ -c main09.1.cpp -o main09.1.o $(CFLAGS)
popolazione.o : popolazione.cpp popolazione.h
	g++ -c popolazione.cpp -o popolazione.o $(CFLAGS)
clean :
	rm *.o main09.1.exe
