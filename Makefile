CC = g++
CCFLAGS = -std=c++20 -g

LD = g++
LDFLAGS = -g

#INC = -IMESTypes

main: main.o MESTypes.o
	${LD} ${LDFLAGS} main.o MESTypes.o -o main ${LIB}

main.o: main.cpp MESTypes.hpp
	${CC} -c ${CCFLAGS} main.cpp ${INC}

MESTypes.o: MESTypes.cpp MESTypes.hpp
	${CC} -c ${CCFLAGS} MESTypes.cpp ${INC}

.PHONY: run
run: main
	./main
