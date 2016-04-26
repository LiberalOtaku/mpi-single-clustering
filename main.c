// File: main.c
// Authors: Edward Ly, Eamon Roosa, George Crowson
// Last Updated: 26 April 2016
// Single-linkage clustering program on strings using MPI implementation.
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>
#include <string.h>

void q_sort(int *start, int *end);
int get_tuples(char* strings, int* counts, int max_string_count);
bool is_match(char *a, char *b, int len, int max);
int is_match_char(char a, char b);

// Main program, by Edward
int main(int argc, char **argv)
{
	// MPI variables
	int numprocs, myrank, namelen;	
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status;
	
	// Program variables
	bool debug = true;
	int i, j; // counters
	double start = 0.0, end = 0.0, wall_server, wall_client, wall_comp;
	
	// Main variables
	int string_len, max_num_strings = 1e6, num_strings = 0, threshold = 2;
	int diff, num_matches, local_distance, next_target = 0, current_string;
	char *a, *local_a;
	int *count;
	bool *merged;
	size_t count_size = sizeof(int) * max_num_strings;
	size_t merged_size = sizeof(bool) * max_num_strings;
	
	// Initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Get_processor_name(processor_name, &namelen);

	// Setup input
	if (argc < 2) {
		fprintf(stderr, "Usage: %s string-length\n", argv[0]);
		exit(1);
	}
	string_len = atoi(argv[1]);
	
	// Variables dependent on string length
	int array_len = string_len + 1;
	int a_char_size = num_strings * string_len;
	int num_local_strings = (num_strings + numprocs - 1) / numprocs; // round up
	int local_a_char_size = num_local_strings * string_len;
	size_t size = sizeof(char) * max_num_strings * array_len;
	size_t local_size = sizeof(char) * num_local_strings * array_len;
	stream = (char *)malloc(sizeof(char) * array_len);
	char target[array_len], other[array_len];
	
	// Malloc a, count, and merged arrays
	a = (char *)malloc(size);
	if (!a) {
		if (myrank == 0)
			perror("unable to allocate array a: "); 
		exit(-1);
	}
	count = (int *)malloc(count_size);
	if (!count) {
		if (myrank == 0)
			perror("unable to allocate array count: "); 
		exit(-1);
	}
	merged = (bool *)malloc(merged_size);
	if (!merged) {
		if (myrank == 0)
			perror("unable to allocate array merged: "); 
		exit(-1);
	}	
	
	

	// Start clock
	start = MPI_Wtime();
	
	// Read data file from memory
	num_strings = get_tuples(a, count);

	// Stop clock and record read time
	end = MPI_Wtime();
	wall_server = end - start;
	

	// Sort strings from highest count to lowest count here (WIP)
	int *start = &count[0];
	int *end = &count[num_strings - 1];
	q_sort(start, end);
	
	
	
	// Start clock
	start = MPI_Wtime();
	
	// Malloc local copy of a
	local_a = (char *)malloc(local_size);
	if (!local_a) {
		if (myrank == 0)
			perror("unable to allocate array local_a: "); 
		exit(-1);
	}
	
	// Copy relevant strings to local process
	MPI_Scatter(a, local_a_char_size, MPI_CHAR, local_a, local_a_char_size, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	// Stop clock and record copy time
	end = MPI_Wtime();
	wall_client = end - start;
	
	// Print wall times
	if (debug && myrank == 0) {
		fprintf(stderr, "Tot. read/copy wall time: server = %f, clients = %f\n", wall_server, wall_client);
		fprintf(stderr, "Avg. read/copy wall time: server = %f, clients = %f\n", wall_server / (double)num_strings, wall_client / (double)num_strings);
	}
	
	
	// Now for the actual work!
	while (true)
	{
		// Find next target string
		if (myrank == 0) {
			while (merged[next_target] == true && next_target < num_strings)
				next_target++;
			if ( next_target < num_strings )
				strncpy(target, &a[string_len * next_target], string_len);
			else break; // we are done
		}		
		
		// Start clock
		MPI_Barrier();
		start = MPI_Wtime();
			
		// Do the actual work (find the strings to merge and merge them)
		num_matches = 0;
		for (i = 0; i < local_a_char_size; i += string_len)
		{
			current_string = (myrank * num_local_strings) + (i / string_len);
			if ((next_target != current_string) && (current_string < num_strings) && (merged[current_string] == false))
			{
				strncpy(other, &local_a[i], string_len);
				if (is_match(target, other, string_len, threshold))
				{
					// Merge current string and add count
					num_matches++;
					merged[(myrank * num_local_strings) + (i / string_len)] = true;
					count[next_target] += count[current_string];
				}
			}
		}
		merged[next_target] = true;
		
		// Stop clock and record/print comp time
		end = MPI_Wtime();
		wall_comp = end - start;
		if (debug)
			fprintf(stdout, "num_strings = %d, target = %d, rank = %d, num_matches = %d, wall_time = %.6f\n", num_strings, next_target, myrank, num_matches, wall_comp);
		MPI_Barrier();
	}
	
	// Do something to output list here (optional?)
	
	
	// Cleanup
	free(a); free(local_a);
	MPI_Finalize();
	exit(0);
}


// Read strings and counts from file, by George
int get_tuples(char* strings, int* counts, int max_string_count)
{
	char str[length + 1];
	int current_string_count = 0;

	{
		event("opening file");
		FILE* f = fopen("/cluster/home/charliep/courses/cs360/single-linkage-clustering/Iceland2014.trim.contigs.good.unique.good.filter.unique.count.fasta", "r");
		if(f == NULL) error("file was unable to be opened");

		event("reading file");
		//@note #this could be more robust but the input is well-formatted so meh
		while (fscanf(f, "%s", str) != EOF) {
			strcpy(&strings[length * current_string_count], str);
			//printf("str: %s\n", str);

			//get the count
			//@note #could probably use fscanf instead
			fscanf(f, "%s", str);
			counts[current_string_count] = atoi(str);

			current_string_count++;

			//@hack #debugging
			if (current_string_count >= 100) break;

			if (current_string_count > max_string_count)
				error("number of tuples in file is greater than maximum allowed\n");
		}

		//debug("loaded %d tuples", current_string_count);

		event("closing file");
		fclose(f);
		free(f);
	}

	return current_string_count;
}


// The actual string comparison (working in chunks with early exit), by Edward
bool is_match(char *a, char *b, int len, int max)
{
	int diff = 0;
	for (int i = 0; i < len; i += 6)
	{
		diff += is_match_char(*(a + i    ), *(b + i    ));
		diff += is_match_char(*(a + i + 1), *(b + i + 1));
		diff += is_match_char(*(a + i + 2), *(b + i + 2));
		diff += is_match_char(*(a + i + 3), *(b + i + 3));
		diff += is_match_char(*(a + i + 4), *(b + i + 4));
		diff += is_match_char(*(a + i + 5), *(b + i + 5));
		if (diff > max)
			return false;
	}
	return true;
}


// The individual character comparison, by Eamon
int is_match_char(char a, char b)
{
	char is_matched;
	
	//Check whether the twho characters match, but don't
	//use a branch. ASM is needed to do a left rotate
	//so that we can interpret the sign bit separately.

	//If a == b, then a ^ b == 0. We then subtract one and
	//rotate the byte left by 1 bit so that the sign bit
	//can be masked off. The result will be 0 if the characters
	//match and 1 if they don't. We can then maintain a
	//running count of the distance between the two strings
	//without using a branch. This is probably unnecessary, but
	//is way too cool to do without.

	//NOTE: this only works if a and b are less than 128.
	//Otherwise the sign bits will
	asm ("mov %1, %%al\n\t"
		"xor %2, %%al\n\t"
		"sub $1, %%al\n\t"
		"rol $1, %%al\n\t"
		"and $0x01, %%al\n\t"
		"xor $0x01, %%al\n\t"
		"mov %%al, %0"
		: "=r" (is_matched)
		: "r" (a), "r" (b)
		: "%al");

	//is_match will be 0 when the chars are equal and will
	//be one when they are different
	
	return (int)is_matched;
}