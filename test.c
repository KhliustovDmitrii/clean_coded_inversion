#include <stdio.h>
#include <string.h>

int main(int argc, char **argv)
{
        printf("%s\n", strcat(argv[0]+2*sizeof(char), ".conf"));

	return 0;
}
