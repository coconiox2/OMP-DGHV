run: build 
	./main

build: 
	g++ -g -std=c++11 Code/main.cpp Code/Flat_DGHV.cpp Code/Params.cpp Code/utilities.cpp -o main -O3 -lntl -lgmp -fopenmp -lpthread
clean:
	rm -rf main
