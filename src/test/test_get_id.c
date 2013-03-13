#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

int main(void)
{
	char test_str_true[]   = "toto_003",
	     test_str_true2[]  = "toto0.5_003",
	     test_str_false[]  = "toto0.5.dat",
	     test_str_false2[] = "toto.dat";
	
	int val = get_id(test_str_true);
	if( val == -1 )
	{
		fprintf(stderr, "\033[31mFAILED :\033[00m %s hasn't give any result!\n", test_str_true);
		return EXIT_FAILURE;
	}
	else if( val != 3 )
	{
		fprintf(stderr, "\033[31mFAILED :\033[00m Found %d instead of 3!\n", val);
		return EXIT_FAILURE;
	}
	printf("\033[33mSUCCESS :\033[00m Found the correct value : %d (3)!\n", val);

	val = get_id(test_str_true2);
	if( val == -1 )
	{
		fprintf(stderr, "\033[31mFAILED :\033[00m %s hasn't give any result!\n", test_str_true2);
		return EXIT_FAILURE;
	}
	else if( val != 3 )
	{
		fprintf(stderr, "\033[31mFAILED :\033[00m Found %d instead of 3!\n", val);
		return EXIT_FAILURE;
	}
	printf("\033[33mSUCCESS :\033[00m Found the correct value : %d (3)!\n", val);

	val = get_id(test_str_false);
	if( val != -1 )
	{
		fprintf(stderr, "\033[31mFAILED :\033[00m Found something with %s :: %d!\n", test_str_false, val);
		return EXIT_FAILURE;
	}
	printf("\033[33mSUCCESS :\033[00m Found the correct value : nothing!\n");
	
	val = get_id(test_str_false2);
	if( val != -1 )
	{
		fprintf(stderr, "\033[31mFAILED :\033[00m Found something with %s :: %d!\n", test_str_false2, val);
		return EXIT_FAILURE;
	}
	printf("\033[33mSUCCESS :\033[00m Found the correct value : nothing!\n");
	
	return EXIT_SUCCESS;
}
