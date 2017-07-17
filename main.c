// File: main.c
// Group Members: Edward Ly, Eamon Roosa, George Crowson
// Last Updated: 29 April 2016
// Single-linkage clustering program on strings using MPI implementation.
// Compile command: mpicc -O0 -g -Wall -lm -o main main.c

/* Basic algorithm, copied from George:
    master:
        load tuples into arrays
        sort tuples
        iterate tuples:
            send target to client
            wait for client to finish
            receive indices of merged tuples from client
            update tuples to reflect merges

    client:
        define stride of the decomposition
        load tuples into arrays
        sort tuples
        iterate tuples:
            receive target from master
            find first relevant tuple
            iterate through tuples doing merging
            send master the indices of merged tuples
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

// The individual character comparison, by Eamon
int is_match_char(char a, char b) {
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

// The actual string comparison (working in chunks with early exit), by Edward
// Assumed that string length is divisible by 6
bool is_match(char *a, char *b, int len, int max) {
	int diff = 0;
	for (int i = 0; i < len; i += 6) {
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

// Functions to sort counts and strings simultaneously from highest to lowest count, by George
void swap(int *file_index, char *strings, int *count, int left, int right, const int string_len) {
    int tmp_index     = file_index[left];
    file_index[left]  = file_index[right];
    file_index[right] = tmp_index;

    char tmp_string[string_len + 1];
    strncpy(tmp_string,                   &strings[string_len * left],  string_len);
    strncpy(&strings[string_len * left],  &strings[string_len * right], string_len);
    strncpy(&strings[string_len * right], tmp_string,                   string_len);

    int tmp_count;
    tmp_count    = count[left];
    count[left]  = count[right];
    count[right] = tmp_count;
}

int partition(int *file_index, char *strings, int *count, int low, int high, const int string_len) {
    int left  = low;
    int right = high;
    // int pivot = left; // unused
    int pivot_item = count[left];

    while (left < right) {
        while (count[left] <= pivot_item) left++;
        while (count[right] > pivot_item) right--;
        if (left < right) swap(file_index, strings, count, left, right, string_len);
    }

    swap(file_index, strings, count, low, right, string_len);

    return right;
}

void quicksort(int *file_index, char *strings, int *count, int low, int high, const int string_len) {
    if (high > low) {
        int pivot = partition(file_index, strings, count, low, high, string_len);
        quicksort(file_index, strings, count, low, pivot - 1, string_len);
        quicksort(file_index, strings, count, pivot + 1, high, string_len);
    }
}

void sort_tuples(int *file_index, char *strings, int *count, int low, int high, const int string_len) {
    quicksort(file_index, strings, count, low, high, string_len);

    int middle = ceil((double)(high/2.0f));
    for (int i = low; i < middle; i++)
        swap(file_index, strings, count, i, high - i, string_len);
}

// ---------------------------------------------------------------------------

// Main program, by Edward
int main(int argc, char **argv) {
	// Constant variables
	const bool debug_stderr  = true;
	const bool debug_event   = true;
	const bool debug_stream  = false;
	const bool enable_output = false;
	const int string_len = 594, MAX_NUM_STRINGS = 100, THRESHOLD = 2, CLIENT_STOP = -1;
	const char* filename = "/cluster/home/charliep/courses/cs360/single-linkage-clustering/Iceland2014.trim.contigs.good.unique.good.filter.unique.count.fasta";

	// MPI variables
	int numprocs, myrank, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status;

	// Initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Get_processor_name(processor_name, &namelen);

	// Check if we have more than one process
    if (numprocs < 2) {
		fprintf(stderr, "Unable to run program with less than two processes.\n");
		MPI_Finalize();
		exit(2);
	}

	// Program variables
	int num_clients = numprocs - 1;
	double start = 0.0, end = 0.0, wall_t, max_wall_t;
	int i, j, k; // counters

	// Main variables
	int num_strings = 0; // num_matches, next_target = 0, current_string;
	int match_count = 0, match_count_sum;
	char *stream, *a, *other;
	// char *local_a;
	int *string_index, *match_index, *count;
	bool *merged;

	// Variables dependent on array length
	int array_len = string_len + 1;
	int num_local_strings = (MAX_NUM_STRINGS + numprocs - 1) / numprocs; // round up
	// int local_a_char_size = num_local_strings * string_len;
	size_t stream_size = sizeof(char) * array_len;
	size_t other_size = stream_size;
	size_t a_size = sizeof(char) * MAX_NUM_STRINGS * array_len;
	// size_t local_a_size = sizeof(char) * num_local_strings * array_len;
	size_t count_size = sizeof(int) * MAX_NUM_STRINGS;
	size_t merged_size = sizeof(bool) * MAX_NUM_STRINGS;
	size_t str_index_size = sizeof(int) * MAX_NUM_STRINGS;
	size_t match_index_size = sizeof(int) * num_local_strings;

	// Malloc all global arrays
	stream = (char *)malloc(stream_size);
	if (!stream) {
		if (myrank == 0) perror("Unable to allocate array stream: ");
		MPI_Finalize();
		exit(-1);
	}
	other = (char *)malloc(other_size);
	if (!other) {
		if (myrank == 0) perror("Unable to allocate array other: ");
		MPI_Finalize();
		exit(-1);
	}
	a = (char *)malloc(a_size);
	if (!a) {
		if (myrank == 0) perror("Unable to allocate array a: ");
		MPI_Finalize();
		exit(-1);
	}
	count = (int *)malloc(count_size);
	if (!count) {
		if (myrank == 0) perror("Unable to allocate array count: ");
		MPI_Finalize();
		exit(-1);
	}
	merged = (bool *)malloc(merged_size);
	if (!merged) {
		if (myrank == 0) perror("Unable to allocate array merged: ");
		MPI_Finalize();
		exit(-1);
	}
	string_index = (int *)malloc(str_index_size);
	if (!string_index) {
		if (myrank == 0) perror("Unable to allocate array string_index: ");
		MPI_Finalize();
		exit(-1);
	}
	match_index = (int *)malloc(match_index_size);
	if (!match_index) {
		if (myrank == 0) perror("Unable to allocate array match_index: ");
		MPI_Finalize();
		exit(-1);
	}

	// Set malloc'd arrays to initial values
	for (i = 0; i < array_len; i++)
		stream[i] = other[i] = ' ';
	for (i = 0; i < MAX_NUM_STRINGS * array_len; i++)
		a[i] = ' ';
	for (i = 0; i < MAX_NUM_STRINGS; i++) {
		count[i] = 0;
		merged[i] = false;
		string_index[i] = i;
	}

	if (debug_event && myrank == 0)
		fprintf(stderr, "Setup of variables successful.\n");


	// Open and read strings and counts from file, based on George's code
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	int stream_count;
	FILE* fstream = fopen(filename, "r");
	if (!fstream) {
		if (myrank == 0) perror("Unable to open file: ");
		MPI_Finalize();
		exit(1);
	}

	while (fscanf(fstream, "%s %d", stream, &stream_count) != EOF) {
		if (debug_stream && myrank == 0)
			fprintf(stderr, "%.30s, %d, %d\n", stream, stream_count, num_strings);
		strncpy(&a[string_len * num_strings], stream, string_len);
		count[num_strings] = stream_count;
		num_strings++;

		if (num_strings >= MAX_NUM_STRINGS) break;
	}

	fclose(fstream);
	end = MPI_Wtime();
	wall_t = end - start;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&wall_t, &max_wall_t, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (debug_stderr && myrank == 0) fprintf(stderr, "File successfully read with %d strings in %.6f seconds.\n", num_strings, max_wall_t);

	// Sort strings from highest count to lowest count
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	sort_tuples(string_index, a, count, 0, num_strings - 1, string_len);
	end = MPI_Wtime();
	wall_t = end - start;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&wall_t, &max_wall_t, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (debug_event && myrank == 0)
		fprintf(stderr, "Sorting of strings completed successfully in %.6f seconds.\n", max_wall_t);

/* 	Edward's obsolete implementation (collapse to ignore):
	// Malloc local arrays
	local_a = (char *)malloc(local_a_size);
	if (!local_a) {
		if (myrank == 0) perror("Unable to allocate array local_a: ");
		MPI_Finalize();
		exit(-1);
	}


	// Copy relevant strings to local processes
	start = MPI_Wtime();
	MPI_Scatter(a, num_local_strings, MPI_CHAR, local_a, num_local_strings, MPI_CHAR, 0, MPI_COMM_WORLD);
	end = MPI_Wtime();
	wall_client = end - start;

	// Print wall times
	if (debug_stderr && myrank == 0) {
		fprintf(stderr, "Tot. read/copy wall time: server = %f, clients = %f\n", wall_server, wall_client);
		fprintf(stderr, "Avg. read/copy wall time: server = %f, clients = %f\n", wall_server / (double)num_strings, wall_client / (double)num_strings);
	} */


	// Now for the actual work (based on George's implementation)!
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	if (myrank == 0) {
		for (i = 0; i < num_strings; i++) {
			if (merged[i]) continue;
			merged[i] = true;

			// send target to clients
			for (j = 1; j < numprocs; j++) {
				MPI_Send(&i, sizeof(i), MPI_INT, j, 1, MPI_COMM_WORLD);
			}

			// get merged indices from clients
			match_count_sum = count[i];
			for (j = 1; j < numprocs; j++) {
				MPI_Recv(&match_count, sizeof(match_count), MPI_INT, j, 1, MPI_COMM_WORLD, &status);
				MPI_Recv(&match_index, sizeof(int) * match_count, MPI_INT, j, 1, MPI_COMM_WORLD, &status);

				if (match_count > 0) {
					// resolve matches
					for (k = 0; k < match_count; k++) {
						int index = match_index[k];
						match_count_sum += count[index];
						merged[index] = true;
					}
				}
			}
			count[i] = match_count_sum;
		}

		for (i = 1; i < numprocs; i++)
			MPI_Send(&CLIENT_STOP, sizeof(CLIENT_STOP), MPI_INT, i, 1, MPI_COMM_WORLD);
	}
	else {
		while (true) {
			// get target
			int t;
			MPI_Recv(&t, sizeof(t), MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

			if (t == CLIENT_STOP) break;

			merged[t] = true;

			// find first tuple that needs to be checked
			int group_index = t - (t % num_clients); //(int) floor((double) (t / num_clients)) * num_clients;
			int start_index = group_index + myrank - 1;

			// iterate through tuples
			// @note #tuples with the same count that were searched previously will have already been matched
			match_count = 0;
			for (i = start_index; i < num_strings; i += num_clients) {
				if (merged[i]) continue;

				// check if index is a match
				if (is_match(&a[string_len * t], &a[string_len * i], string_len, THRESHOLD)) {
					match_index[match_count] = i;
					match_count++;
					merged[i] = true;
				}
			}

			// returning matched indices
			MPI_Send(&match_count, sizeof(match_count), MPI_INT, 0, 1, MPI_COMM_WORLD);
			MPI_Send(&match_index, sizeof(int) * match_count, MPI_INT, 0, 1, MPI_COMM_WORLD);
		}
	}


/* 	Edward's failed attempt at a better implementation (collapse to ignore):
	for (i = 0; i < num_strings; i++)
	{
		if (merged[i]) continue;
		merged[i] = true;
		strncpy(target, &a[string_len * next_target], string_len);

		// Find the strings to merge and merge them
		num_matches = 0;
		MPI_Barrier(MPI_COMM_WORLD);
		start = MPI_Wtime();
		for (j = 0; j < num_local_strings; j++)
		{
			current_string = (myrank * num_local_strings) + i;
			if ((next_target != current_string)
				&& (current_string < num_strings)
				&& (merged[current_string] == false))
			{
				strncpy(other, &local_a[i * string_len], string_len);
				if (is_match(target, other, string_len, THRESHOLD))
				{
					// Merge current string and add count
					num_matches++;
					merged[current_string] = true;
					count[next_target] += count[current_string];
				}
			}
		}
		end = MPI_Wtime();
		wall_comp = end - start;
		if (debug_stderr)
			fprintf(stderr, "num_strings = %d, target = %d, rank = %d, num_matches = %d, wall_time = %.6f\n", num_strings, next_target, myrank, num_matches, wall_comp);
	} */

	end = MPI_Wtime();
	wall_t = end - start;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&wall_t, &max_wall_t, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (debug_event && myrank == 0)
		fprintf(stderr, "Single-linkage clustering completed in %.6f seconds.\n", max_wall_t);

	// Print (to terminal or file) sorted list of strings and final counts
	MPI_Barrier(MPI_COMM_WORLD);
	if (enable_output && myrank == 0) {
		for (i = 0; i < num_strings; i++) {
			strncpy(other, &a[i * string_len], string_len);
			fprintf(stdout, "%s %d\n", other, count[i]);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (debug_event)
		fprintf(stderr, "Program execution successful on rank %d of %d.\n", myrank, numprocs);

	// Cleanup, we are done
	MPI_Finalize();
	// free(a); free(other); free(stream);
	// free(count); free(merged);
	// free(match_index); free(string_index);
	exit(0);
}
