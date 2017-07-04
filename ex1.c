#include <stdio.h>

/* This is a comment. */
int main(int argc, char *argv[])
{
	int distance = 100, distance2 = 200, empty;
	
	// this is also a comment
	printf("You are %d miles away. \n", distance);
	printf("Now you are %d miles away. \n", distance2);
	empty = distance+distance2;
	printf("And now you are %d miles away. \n", empty);
	empty = distance*distance2;
	printf("Even still now %d miles away. \n", empty);
	empty = distance2/distance;
	printf("Now %d miles away. \n", empty);
	
	
	return 0;
}