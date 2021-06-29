CFLAGS = -Wall -O3 --std=c++11

main02.1.exe : main02.1.o random.h random.o func.h funzionebase.h
	g++ random.o main02.1.o -o main02.1.exe
main02.1.o : main02.1.cpp func.h random.h funzionebase.h
	g++ -c main02.1.cpp -o main02.1.o $(CFLAGS)
random.o : random.cpp random.h funzionebase.h
	g++ -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main02.1.exe
