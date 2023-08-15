#make file

all:
	g++ -g main.c forward_problem.c kalman.c math_utils.c data_processing.c -o main

