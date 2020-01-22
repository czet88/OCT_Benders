#include "headers.h"

/// INFO ///////////////////////////////////////////////////////////////////////
//
// main.c: Main function
//
// Author:  Carlos Luna-Mota
// Created: 2015 Jan 21
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
	
	CONFIGURATION param;		// Global Paramenters of the algorithm
	INSTANCE *instances;		// Array of instance data
	RESULTS	 *results;			// Array of results
	int n_of_instances;			// Number of instances to solve
	int i;						// Index variable
	FILE *output = NULL;		// Output file
	char filename[50] = "";		// Output filename
	time_t now = time(NULL);			// Clock
	struct tm tm = *localtime(&now);	// Timestamp holder
	int status = 0;

	// User must provide a filename with the experimental settings
	// and another with the instances list
	if (argc < 3) {	show_help(argv[0]); }
	else {
		// read experimental settings:
		status = read_CONFIGURATION(argv[1], &param);
		assert(status == 0);

		// read the instance list
		instances = create_INSTANCE_vector(&n_of_instances, argv[2]);
		results	  = create_RESULTS_vector(n_of_instances);
		
		// create an output file to print results:
		sprintf(filename, "RESULTS - %04d%02d%02d%02d%02d%02d.txt", tm.tm_year+1900, tm.tm_mon+1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
		output_settings(filename, param, n_of_instances);
		
		for (i=0; i<n_of_instances; i++) {
			countfails=0;
			// Read the instance data:
			status = read_INSTANCE(&(instances[i]));
			assert(status == 0);

			// HEURISTIC METHODS ///////////////////////////////////////////////
			status = Heuristics(param, instances[i], &(results[i]));
			assert(status == 0);
			////////////////////////////////////////////////////////////////////
			
			// EXACT METHODS ///////////////////////////////////////////////////
			if (param.ALGORITHM != NONE) {
				if (param.ALGORITHM == 0 || param.ALGORITHM == 1) {		// Using Benders
					status = Benders(param, instances[i], &(results[i]));
					assert(status == 0);	
				} else if (param.ALGORITHM == 2) {		// Using the 2-index formulation
					status = TwoIndex(param, instances[i], &(results[i]));
					assert(status == 0);
				} else if (param.ALGORITHM == 3) {		// Using the 3-index formulation
					status = ThreeIndex(param, instances[i], &(results[i]));
					assert(status == 0);	
				} else if (param.ALGORITHM == 4) {		// Using the 4-index formulation
					status = FourIndex(param, instances[i], &(results[i]));
					assert(status == 0);	
				} 
			}
			////////////////////////////////////////////////////////////////////

			output_RESULTS(filename, param, instances[i], results[i]);

			// Free instance data:
			free_INSTANCE(&(instances[i]));
		}

		// Output average times & gaps...
		//output_averaged_RESULTS(filename, param, n_of_instances, results);

		// clean:
		free_RESULTS_vector(&results, n_of_instances);
		free_INSTANCE_vector(&instances);
	}

	return status; // 0 meaning success
}
////////////////////////////////////////////////////////////////////////////////
