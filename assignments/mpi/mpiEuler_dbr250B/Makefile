

default:
	mpicc -Wall -o eulerSolve eulerSolve.c config.c euler.c -lm
	mpicc -Wall -o test test.c config.c euler.c -lm	
	mpicc -Wall -o test2 test2.c config.c euler.c -lm	
clean:
	rm -r test
