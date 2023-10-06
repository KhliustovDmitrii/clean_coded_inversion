#make file

all:
	gcc -g main.c forward_problem.c kalman.c math_utils.c data_processing.c -o main.exe

#all:
#	gcc -g synt_data_generation.c forward_problem.c math_utils.c -o main.exe
