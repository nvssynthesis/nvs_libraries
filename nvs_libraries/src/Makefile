CC=gcc
CFLAGS=-I../include -Wall -std=c++17

HEADERS = nvs_LFO.h nvs_filters.h

all: nvs_LFO.o nvs_filters.o 

# %.o: %.cpp $(HEADERS)
# 	$(CC) -c -o $@ $< $(CFLAGS)

nvs_LFO.o: nvs_LFO.cpp #$(HEADERS)
	$(CC) -c -o ../obj/$@ $^ $(CFLAGS)

nvs_filters.o: nvs_filters.cpp
	$(CC) -c -o ../obj/$@ $^ $(CFLAGS)