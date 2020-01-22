#include "headers.h"

/// INFO ///////////////////////////////////////////////////////////////////////
//
// CLM_utils.c: General Algorithms
//
// Author:  Carlos Luna-Mota
// Created: 2015 Feb 12
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Monte-Carlo unbiased random integer generation in the interval [ini, end]
int randint(int ini, int end) {
	int x;
	int range     = ABS(end-ini) + 1;
	int big_range = RAND_MAX - (RAND_MAX % range);
    do { x = rand(); } while (x >= big_range);
    return ini + (x % range);
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
double list_max(int N, double *list) {
	int i, best = 0;
	for (i=1; i<N; i++) { if (list[best] < list[i]) { best = i; } }
	return list[best];
}
double list_min(int N, double *list) {
	int i, best = 0;
	for (i=1; i<N; i++) { if (list[best] > list[i]) { best = i; } }
	return list[best];
}
double list_sum(int N, double *list) {
	int i;
	double result = 0.0;
	for (i=0; i<N; i++) { result += list[i]; }
	return result;
}
double list_prod(int N, double *list) {
	int i;
	double result = 1.0;
	for (i=0; i<N; i++) { result *= list[i]; }
	return result;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Finds the smallest k'th element of an array of size N.
// Note: It makes a copy of list to avoid modifiying it!
double find_smallest(int N, double *list, int k) {
    int l,h;
	int low  = 0;
	int high = N-1;
    double pivot;
	double aux = list[0];
	double *aux_list = NULL;
	
	if (k >= N-1)      { for (l=1; l<N; l++) { aux = MAX(aux, list[l]); } } // Special case: FIND MAX
	else if (k <=  0 ) { for (l=1; l<N; l++) { aux = MIN(aux, list[l]); } }	// Special case: FIND MIN
	
	// GENERAL CASE:
	else {				

		// Copy the whole list to avoid modifiying 'list'
		aux_list = create_double_vector(N);
		for (l=0; l<N; l++) { aux_list[l] = list[l]; }

		// Now apply the algorithm over the 
		while (low<high) {
			pivot = aux_list[k];		
			l=low;
			h=high;
			do {
				while (aux_list[l]<pivot) { l++; }
				while (pivot<aux_list[h]) { h--; }
				if (l<=h) {
					aux         = aux_list[l];	// Swap: 
					aux_list[l] = aux_list[h];	//  aux_list[l] 
					aux_list[h] = aux;		    //  aux_list[h]
					l++; h--;
				}
			} while (l<=h);
		
			// Now we have:
			// [__smaller than elem__] [__bigger than elem__]
			// ^low                 h^ ^l               high^

			// So it is easy to decide where the k'th element is:
			if (h<k) { low  = l; }
			if (k<l) { high = h; }
		}

		aux = aux_list[k];				// get result
		free_double_vector(&aux_list);	// free memory
	}

	return aux;						// return result
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Compares the first N cells of two integer arrays.
int ARE_EQUAL(int N, int *arr1, int *arr2) {
	int i;
	for (i=0; i<N; i++) { if (arr1[i] != arr2[i]) { return NO; } }
	return YES;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Returns the elapsed time (in seconds) since the last "init_time = clock();"
double elapsed(clock_t init_time) {
	return ((double) (clock() - init_time)) / CLOCKS_PER_SEC;
}
////////////////////////////////////////////////////////////////////////////////
