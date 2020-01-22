#include "headers.h"

/// INFO ///////////////////////////////////////////////////////////////////////
//
// output.c: Funtions for writting data
//
// Author:  Carlos Luna-Mota
// Created: 2015 Jan 21
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Shows a command line message to help the user.
void show_help(char *program_name){
	fprintf(stderr, "\n\n\tHELP:\n\nYou must provide 2 valid filenames:\n\n\t%s configuration.txt instancelist.txt\n", program_name);
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void output_settings(char *filename, CONFIGURATION param, int n_of_instances) {

	FILE *output = NULL;
	
	// Open file:
	output = open_file(filename, "w+");
	assert(output != NULL);
	
	// Write settings to file:
	fprintf(output, "# %s\n#\n", filename);
	fprintf(output, "# SETTINGS:\n");
	if (param.TO_INTEGRALITY == NO) { fprintf(output, "# Solve the Linear Relaxation with "); }
	else                            { fprintf(output, "# Solve the Integer Problem with ");   }
	if (param.ALGORITHM == 0)      { fprintf(output, "Iterative Benders and...\n"); }
	else if (param.ALGORITHM == 1) { fprintf(output, "Branch&Cut Benders and...\n");       }
	else if (param.ALGORITHM == 2) { fprintf(output, "the 2-index formulation and...\n"); }
	else if (param.ALGORITHM == 3) { fprintf(output, "the 3-index formulation ...\n");       }
	else                           { fprintf(output, "the 4-index formulation and...\n");      }
	fprintf(output, "#\t-> A maximum of %.0f seconds of CPU time.\n", param.MAX_CPU_TIME);	

	if (param.ALGORITHM < 2) {
		if (param.FEAS_CUT_TYPE == 0) { fprintf(output, "#\t-> CutSet inequalities as Feasibility cuts"); }
		else                          { fprintf(output, "#\t-> Subtour elimination constraints as Feasibility cuts"); }
		fprintf(output, ".\n");

		if (param.OPT_CUT_TYPE == 0)      { fprintf(output, "#\t-> Regular Benders Optimality Cuts"); }
		else if (param.OPT_CUT_TYPE == 1) { fprintf(output, "#\t-> Magnanti-Wong Optimality Cuts"); }
		else if (param.OPT_CUT_TYPE == 2) { fprintf(output, "#\t-> Papadakos Optimality Cuts"); }
		if      (param.ADD_NON_VIOLATED_CUTS == NO) { fprintf(output, " (adding only violated cuts)"); }
		fprintf(output, ".\n");

		if (param.HEUR_LIFTING == YES) { fprintf(output, "#\t-> Lifting Integer Cuts (Heuristically).\n"); }
		if (param.ALGORITHM == 1) {
			if (param.ADD_FRAC_OPT_CUTS == YES)	{ fprintf(output, "#\t-> Adding Fractional Optimality Cuts (if depth is multiple of %d).\n", param.FRAC_CUTS_FREQ); }
			if (param.ADD_FRAC_FEAS_CUTS == YES)	{ fprintf(output, "#\t-> Adding Fractional Feasibility Cuts (if depth is multiple of %d).\n", param.FRAC_CUTS_FREQ); }
			if (param.ADD_LOCAL_CUTS == YES)	{ fprintf(output, "#\t-> Adding Local Cuts.\n"); }
		}
		if (param.USE_ROUNDING_CALLBACK == YES) { fprintf(output, "#\t-> Providing a Rounding Callback to CPLEX.\n"); }
		if (param.UPPER_CUTOFF == YES)          { fprintf(output, "#\t-> Providing a Combinatorial Upper Cutoff.\n"); }

		fprintf(output, "#\t-> Providing the following additional constraints at root node:\n");
		if (param.GLOBAL_LOWER_BOUND == YES)    { fprintf(output, "#\t\t-> Global Combinatorial Lower Bound.\n"); }
		if (param.COMMODITY_LOWER_BOUND == YES) { fprintf(output, "#\t\t-> Basic Commodity Lower Bounds.\n"); }
		if (param.TREE_3_LOWER_BOUND == YES)    { fprintf(output, "#\t\t-> 3-Trees Lower Bounds.\n"); }
		if (param.TREE_4_LOWER_BOUND == YES)    { fprintf(output, "#\t\t-> 4-Trees Lower Bounds.\n"); }
		if (param.MST_OPT_CUTS == YES)			{ fprintf(output, "#\t\t-> MST Optimality Cuts.\n"); }
		if (param.AMT_OPT_CUTS == YES)			{ fprintf(output, "#\t\t-> AMT Optimality Cuts.\n"); }
		if (param.GHT_OPT_CUTS == YES)			{ fprintf(output, "#\t\t-> GHT Optimality Cuts.\n"); }
		if (param.MCMCT_OPT_CUTS == YES)		{ fprintf(output, "#\t\t-> MCMCT Optimality Cuts.\n"); }
		if (param.STARS_OPT_CUTS == YES)		{ fprintf(output, "#\t\t-> Star-Trees Optimality Cuts.\n"); }
		if (param.ROOT_CLEANUP == 0.0)			{ fprintf(output, "#\t\t-> ...and then delenting all non-binding root-node constraints.\n"); }

	}
	
	if (param.FILE_OUTPUT >= 1) {
		fprintf(output, "\n\n# Inst.\t\t\tSize\t K/E\t| Initial /   Root  / Final GAP\t| \t\tTotal\t/ \t\tRoot\t/ Master/Total / Subprob/Total / Feas/Subprob / Opt/Subprob TIME | B&B Nodes\t");
		if (param.FILE_OUTPUT >= 2 && param.ALGORITHM < 2) { 
			fprintf(output, "|  IntOpt / IntFeas / FracOpt / FracFeas CUTS\t| IntOpt / IntFeas / FracOpt / FracFeas CUTS/SUBPROB\t"); 
		}
		fprintf(output, "| Best_LwrBound / Best_UprBound | ");
		/*fprintf(output, "Reinelt LwrBound       / ");
		fprintf(output, "Elena   LwrBound       / ");
		fprintf(output, "Best Star   (val&time) / ");
		fprintf(output, "MST         (val&time) / ");
		fprintf(output, "GHT         (val+time) / ");
		fprintf(output, "AMT         (val+time) / ");
		fprintf(output, "AMT+AMLS    (val+time) / ");
		fprintf(output, "D&C         (val&time) / ");
		fprintf(output, "D&C+AMLS    (val+time) / ");
		fprintf(output, "Best Star + AMLS       / ");
		fprintf(output, "Best Star + DLS        / ");*/
        fprintf(output, "\n");
	}		
	// Close file:
	fclose(output);

}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void output_RESULTS(char *filename, CONFIGURATION param, INSTANCE data, RESULTS results) {
	int e, ee, new_line;
	int N = data.N;
	int K = data.K;
	int E = data.E;
	FILE *output = NULL;


	int i,j,k,d,size, diameter, sum_of_dist;
	int *degree = create_int_vector(N);
	int **neighbor = create_int_matrix(N,N);
	int *stack = create_int_vector(N);
	int **dist = create_int_matrix(N,N);

	if (param.SCREEN_OUTPUT >= YES) {

		printf("\nRESULTS:\n");

		printf("\tInstance Name: %s\n\n", data.filename);

		printf("\tNumber of vertices:     %d\n", data.N);
		printf("\tNumber of edges:        %d\n", data.E);
		printf("\tNumber of commodities:  %d (%.2f%%)\n\n", data.K, 100.0*((double)(data.K))/((double)(data.E)) );

		printf("\tBase Model:      %d constraints\n",   data.N+1);
		printf("\tImproved Model:  %d constraints\n",   results.initial_cons);
		printf("\tFinal Model:     %d constraints\n\n", results.root_cons);
		
		printf("\tRoot  time:   %.2f sec\n", results.root_time);
		printf("\tTotal time:   %.2f sec\n\n", results.total_time);

		printf("\t\tMaster Problem: %6.2f%%\n\t\tSub Problems:   %6.2f%%\n", 100.0*results.master_time/results.total_time, 100.0*results.sub_time/results.total_time);
		printf("\t\tFeasibility Subproblems: %6.2f%%\n\t\tOptimality Subproblems:  %6.2f%%\n\n", 100.0*results.feas_time/(results.feas_time+results.opt_time), 100.0*results.opt_time/(results.feas_time+results.opt_time));


		printf("\tInitial GAP:  %6.2f%%\n",   100.0*results.initial_gap);
		printf("\tRoot  GAP:    %6.2f%%\n",   100.0*results.root_gap);
		printf("\tFinal GAP:    %6.2f%%\n\n", 100.0*results.final_gap);
		
		printf("\t\tFinal UpperBound:  %.2f\n",   results.UpperBound);
		printf("\t\tFinal LowerBound:  %.2f\n\n", results.LowerBound);
		
		printf("\tBest Solution Found:\n\n\t");

		ee = 0;
		new_line = YES;
		for (e=0; e<E; e++) {
			if (results.best_solution[e] > EPSILON) {
				printf("\t(%d,%d)", data.index_i[e], data.index_j[e]); 
				ee++;
				new_line = NO;
			}	
			if (ee%7 == 0 && new_line == NO && ee+1 <= N) { 
				new_line = YES;
				printf("\n\t");
			}
		} printf("\n");	
	}

	if (param.FILE_OUTPUT >= 1) {

		// Open file
		output = open_file(filename, "a+");
		assert(output != NULL);

		// Append a row to the file
		fprintf(output, "%14s;%4d;%3.0f%% ; %6.2f%% ; %6.2f%% ; %6.2f%%; %10.2f ; %10.2f ; %8.2f%%    ; %8.2f%%     ; %8.2f%%    ; %8.2f%%     ;%9d",
			&((data.filename)[12]), data.N, 100.0*((double)(data.K))/((double)(data.E)),
			100*results.initial_gap,100*results.root_gap,100*results.final_gap,
			results.total_time, results.root_time, 
			100*results.master_time/(results.total_time+EPSILON), 100*results.sub_time/(results.total_time+EPSILON),
			100*results.feas_time/(results.sub_time+EPSILON), 100*results.opt_time/(results.sub_time+EPSILON),
			results.n_explored_nodes);
		if (param.FILE_OUTPUT >= 2 && param.ALGORITHM != 0) {
			fprintf(output, "; %6d ; %7d ; %7d ; %8d; %6d ; %7d ; %7d ;%8d;",
				results.n_int_opt_subproblems, results.n_int_feas_subproblems, 
				results.n_frac_opt_subproblems, results.n_frac_feas_subproblems, 
				results.n_int_opt_cuts, results.n_int_feas_cuts, 
				results.n_frac_opt_cuts, results.n_frac_feas_cuts);
		}
		fprintf(output, "%13.2f; %13.2f ; %13.2f; %d ", results.LowerBound, results.UpperBound, FinalLP, countfails);

		/*for (e=0; e<11; e++) { fprintf(output, "%13.2f %8.2f / ", 	results.heuristic_values[e], results.heuristic_times[e]); }*/
		fprintf(output, "\n");

		// Close file
		fclose(output);

	}

	free_int_vector(&degree);
	free_int_matrix(&neighbor, N);
	free_int_vector(&stack);
	free_int_matrix(&dist,N);

}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void output_averaged_RESULTS(char *filename, CONFIGURATION param, int n_of_instances, RESULTS *results) {

	int i;
	double ave_time = 0.0;
	double ave_gap  = 0.0;
	FILE *output = NULL;

	// Compute average times, gaps, etc...
	for (i=0; i<n_of_instances; i++) {
		ave_time += results[i].total_time;
		ave_gap  += results[i].final_gap;
	}
	ave_time /= ((double)(n_of_instances));
	ave_gap  /= ((double)(n_of_instances));

	if (param.FILE_OUTPUT >= YES) {
		// OPEN a file with a timestamp in the name (to avoid overwritting and confussion)
		output = open_file(filename, "a+");
		assert(output != NULL);

		fprintf(output, "\n# Average time: %.2f", ave_time);
		fprintf(output, "\n# Average gap:  %.2f%%\n", 100*ave_gap);

		// Clean:
		fclose(output);
	}

	if (param.SCREEN_OUTPUT >= YES) {
		printf("\n# Average: %.2f\t%.2f%%\n", ave_time, 100*ave_gap);
	}
}
////////////////////////////////////////////////////////////////////////////////
