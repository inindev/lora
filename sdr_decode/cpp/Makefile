
# Copyright (c) 2024, John Clark <inindev@gmail.com>

CC = clang++
#CC = g++

#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#  -std=  c++17, c++20, c++23, c++2c
CFLAGS = -Wall -g -std=c++17

TARGET = loraphy
 
$(TARGET): main.o LoraPhy.o #kfft.o
	$(CC) $(CFLAGS) -o $(TARGET) main.o LoraPhy.o #kfft.o

main.o: main.cpp LoraPhy.h sdr_stream.h types.h
	$(CC) $(CFLAGS) -c main.cpp

LoraPhy.o: LoraPhy.cpp LoraPhy.h sdr_stream.h types.h fft.h
	$(CC) $(CFLAGS) -c LoraPhy.cpp

kfft.o: kfft.cpp kfft.h types.h
	$(CC) $(CFLAGS) -c kfft.cpp

.PHONY: clean
clean:
	rm -f *.o $(TARGET)
