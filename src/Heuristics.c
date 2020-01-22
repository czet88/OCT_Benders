#include "headers.h"

/// INFO ///////////////////////////////////////////////////////////////////////
//
// Heuristics.c: Optimum Communication Spanning Tree Problem Heuristic Phase
//
// Author:  Carlos Luna-Mota
// Created: 2015 Jun 21
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Applies some OCSTP heuristics to set up an initial lower bound, an initial
// upper bound and an initial best known solution.
int Heuristics(CONFIGURATION param, INSTANCE data, RESULTS *results) {
	int		i,j,k;
	double	lwr, upp, imp;
	clock_t	start;
	double  total_heuristic_time = 0.0;
	
	int		N = data.N;
	int		E = data.E;
	int		K = data.K;
	COMMODITY *Com = data.Com;
	double	**C = data.C;
	double	**W = data.W;
	double	**C1 = data.C1;
	double	**C2 = data.C2;
	int     **index_e = data.index_e;
	int		status = 0;

	double	**dist  = create_double_matrix(N,N);
	double	**Y	    = create_double_matrix(N,N);

	if (param.SCREEN_OUTPUT >= 1) { 
		printf("\n#############################################################\n");
		printf("\nAPPLYING HEURISTIC PHASE: %s\n", data.filename); 
	}


	//////////////////
	// LOWER BOUNDS //
	//////////////////																			
	
	// METHOD 0: Assume that the full graph is available //////////////////////
	results->LowerBound = 0.0;		
	start = clock();
	for (k=0; k<K; k++) { results->LowerBound += Com[k].w * C1[Com[k].o][Com[k].d]; }
	results->heuristic_values[0] = results->LowerBound;
	results->heuristic_times[0]  = elapsed(start);
	///////////////////////////////////////////////////////////////////////////
	
	// METHOD 1: REINELT'S LwrBound ///////////////////////////////////////////
	start = clock();
	if (data.IS_EUCLIDEAN) { lwr = Reinelt_lower_bound(N, W, C1, C2); }
    else                   { lwr = 0.0;                               }
	results->LowerBound = MAX(lwr, results->LowerBound);
	results->heuristic_values[1] = lwr;
	results->heuristic_times[1]  = elapsed(start);	
	///////////////////////////////////////////////////////////////////////////

	// METHOD 2: ELENA'S LwrBound /////////////////////////////////////////////
	start = clock();
	lwr = Elena_lower_bound(N, W, C);
	results->LowerBound = MAX(lwr, results->LowerBound);
	results->heuristic_values[2] = lwr;
	results->heuristic_times[2]  = elapsed(start);
	///////////////////////////////////////////////////////////////////////////
	
	//////////////////
	// UPPER BOUNDS //
	//////////////////

	// METHOD 0: random tree //////////////////////////////////////////////////
	start = clock();	
	random_tree(N, Y);
	tree_dist(N, C, Y, dist);
	upp = 0.0;
	for (k=0; k<K; k++) { upp += Com[k].w * dist[Com[k].o][Com[k].d]; }
	results->heuristic_times[29]  = elapsed(start);
	results->heuristic_values[29] = upp;
	results->UpperBound = 1000*upp;
	results->best_solution = create_double_vector(E+K);
	for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { results->best_solution[index_e[i][j]] = Y[i][j]; } }
	for (k=0; k<K; k++) { results->best_solution[E+k] = dist[Com[k].o][Com[k].d]; }
	total_heuristic_time += results->heuristic_times[29];
	///////////////////////////////////////////////////////////////////////////

	
	// METHOD 1: MST //////////////////////////////////////////////////////////
	if (param.DO_MST_HEUR == YES) {
		start = clock();	
		MST(N, C, Y, NULL);
		tree_dist(N, C, Y, dist);
		if (param.ADD_LOCAL_SEARCH_HEUR == YES) { imp = AM_LS(N, C, W, Y, dist, NULL); }
		upp = 0.0;
		for (k=0; k<K; k++) { upp += Com[k].w * dist[Com[k].o][Com[k].d]; }
		results->heuristic_times[4]  = elapsed(start);
		results->heuristic_values[4] = upp;
		if (upp < results->UpperBound) {
			results->UpperBound = upp;
			for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { 
				results->best_solution[index_e[i][j]] = Y[i][j]; 
			} }
			for (k=0; k<K; k++) { results->best_solution[E+k] = dist[Com[k].o][Com[k].d]; }
		}
		total_heuristic_time += results->heuristic_times[4];		
	}
	///////////////////////////////////////////////////////////////////////////


	// METHOD 2: GHT //////////////////////////////////////////////////////////
	if (param.DO_GHT_HEUR == YES) {
		start = clock();	
		GHT(N, W, Y);
		tree_dist(N, C, Y, dist);
		if (param.ADD_LOCAL_SEARCH_HEUR == YES) { imp = AM_LS(N, C, W, Y, dist, NULL); }
		upp = 0.0;
		for (k=0; k<K; k++) { upp += Com[k].w * dist[Com[k].o][Com[k].d]; }
		results->heuristic_times[5]  = elapsed(start);
		results->heuristic_values[5] = upp;
		if (upp < results->UpperBound) {
			results->UpperBound = upp;
			for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { 
				results->best_solution[index_e[i][j]] = Y[i][j]; 
			} }
			for (k=0; k<K; k++) { results->best_solution[E+k] = dist[Com[k].o][Com[k].d]; }
		}
		total_heuristic_time += results->heuristic_times[5];
	}
	///////////////////////////////////////////////////////////////////////////
	

	// METHOD 3: STARS //////////////////////////////////////////////////////////
	if (param.DO_STAR_HEUR == YES) {
		start = clock();	
		STARS(N, C, W, Y);
		tree_dist(N, C, Y, dist);
		if (param.ADD_LOCAL_SEARCH_HEUR == YES) { imp = AM_LS(N, C, W, Y, dist, NULL); }
		upp = 0.0;
		for (k=0; k<K; k++) { upp += Com[k].w * dist[Com[k].o][Com[k].d]; }
		results->heuristic_times[3]  = elapsed(start);
		results->heuristic_values[3] = upp;

		// AMLS:
		start = clock();	
		imp = AM_LS(N, C, W, Y, dist, NULL);
		results->heuristic_times[8]  = elapsed(start);
		upp = 0.0;
		for (k=0; k<K; k++) { upp += Com[k].w * dist[Com[k].o][Com[k].d]; }
		results->heuristic_values[8] = upp;
	
		// Dandelion
		STARS(N, C, W, Y);
		start = clock();	
		imp = DLS2(N, C, W, Y);
		results->heuristic_times[9]  = elapsed(start);
		tree_dist(N, C, Y, dist);
		upp = 0.0;
		for (k=0; k<K; k++) { upp += Com[k].w * dist[Com[k].o][Com[k].d]; }
		results->heuristic_values[9] = upp;
		
		if (upp < results->UpperBound) {
			results->UpperBound = upp;
			for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { 
				results->best_solution[index_e[i][j]] = Y[i][j]; 
			} }
			for (k=0; k<K; k++) { results->best_solution[E+k] = dist[Com[k].o][Com[k].d]; }
		}
		total_heuristic_time += results->heuristic_times[3];
		total_heuristic_time += results->heuristic_times[8];
		total_heuristic_time += results->heuristic_times[9];
	}
	///////////////////////////////////////////////////////////////////////////


	// METHOD 4: AMT //////////////////////////////////////////////////////////
	if (param.DO_AMT_HEUR == YES) {
		start = clock();	
		AMT(N, C, W, Y, NULL);
		tree_dist(N, C, Y, dist);
		if (param.ADD_LOCAL_SEARCH_HEUR == YES) { imp = AM_LS(N, C, W, Y, dist, NULL); }
		upp = 0.0;
		for (k=0; k<K; k++) { upp += Com[k].w * dist[Com[k].o][Com[k].d]; }
		results->heuristic_times[6]  = elapsed(start);
		results->heuristic_values[6] = upp;

		// AMLS:
		start = clock();
		AMT(N, C, W, Y, NULL);
		tree_dist(N, C, Y, dist);
		imp = AM_LS(N, C, W, Y, dist, NULL);
		upp = 0.0;
		for (k=0; k<K; k++) { upp += Com[k].w * dist[Com[k].o][Com[k].d]; }
		results->heuristic_times[12]  = elapsed(start);		
		results->heuristic_values[12] = upp;

		if (upp < results->UpperBound) {
			results->UpperBound = upp;
			for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { 
				results->best_solution[index_e[i][j]] = Y[i][j]; 
			} }
			for (k=0; k<K; k++) { results->best_solution[E+k] = dist[Com[k].o][Com[k].d]; }
		}
		total_heuristic_time += results->heuristic_times[6];
		total_heuristic_time += results->heuristic_times[12];
	}
	///////////////////////////////////////////////////////////////////////////


	// METHOD 5: BUDaC ////////////////////////////////////////////////////////
	if (param.DO_BUDaC_HEUR == YES) {
		start = clock();	
		upp = OCSTP_BUDaC(N, C, W, Y, NO, NO);
		results->heuristic_times[7]  = elapsed(start);
		results->heuristic_values[7] = upp;

		start = clock();	
		upp = OCSTP_BUDaC(N, C, W, Y, YES, YES);
		results->heuristic_times[13]  = elapsed(start);
		results->heuristic_values[13] = upp;

		if (upp < results->UpperBound) {
			results->UpperBound = upp;
			for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { 
				results->best_solution[index_e[i][j]] = Y[i][j]; 
			} }
			tree_dist(N, C, Y, dist);
			for (k=0; k<K; k++) { results->best_solution[E+k] = dist[Com[k].o][Com[k].d]; }
		} 
		total_heuristic_time += results->heuristic_times[7];
		total_heuristic_time += results->heuristic_times[13];
	}
	///////////////////////////////////////////////////////////////////////////

	// METHOD 6: Best + LS ////////////////////////////////////////////////////
	if (param.ADD_LOCAL_SEARCH_HEUR == YES) {
		// Dandelion
		for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { Y[i][j] = 0.0; } }
		for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { Y[i][j] = Y[j][i] = results->best_solution[index_e[i][j]]; } }
		tree_dist(N, C, Y, dist);
		start = clock();	
		imp = DLS2(N, C, W, Y);
		results->heuristic_times[11]  = elapsed(start);
		tree_dist(N, C, Y, dist);
		upp = 0.0;
		for (k=0; k<K; k++) { upp += Com[k].w * dist[Com[k].o][Com[k].d]; }
		results->heuristic_values[11] = upp;

		// AMLS:
		for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { Y[i][j] = 0.0; } }
		for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { Y[i][j] = Y[j][i] = results->best_solution[index_e[i][j]]; } }
		tree_dist(N, C, Y, dist);
		start = clock();	
		imp = AM_LS(N, C, W, Y, dist, NULL);
		results->heuristic_times[10]  = elapsed(start);
		upp = 0.0;
		for (k=0; k<K; k++) { upp += Com[k].w * dist[Com[k].o][Com[k].d]; }
		results->heuristic_values[10] = upp;	
		
		if (upp < results->UpperBound) {
			results->UpperBound = upp;
			for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { 
				results->best_solution[index_e[i][j]] = Y[i][j]; 
			} }
			for (k=0; k<K; k++) { results->best_solution[E+k] = dist[Com[k].o][Com[k].d]; }
		}

		total_heuristic_time += results->heuristic_times[10];
		total_heuristic_time += results->heuristic_times[11];
		//results->total_time=total_heuristic_time;
		
	}
	///////////////////////////////////////////////////////////////////////////



	// COMPUTE INITIAL GAP: ////////////////////////////////////////////////////
	results->final_gap   = (results->UpperBound-results->LowerBound) / results->UpperBound;	// A priori GAP
	results->initial_gap = results->final_gap;					// Remember that value!
	if (param.SCREEN_OUTPUT >= 1) {
		printf("HeurPhase\t\t%.2f <- %.2f%% -> %.2f  \tTime: %.2f\"\n", results->LowerBound, 100*results->final_gap, results->UpperBound, total_heuristic_time);
	}
	////////////////////////////////////////////////////////////////////////////

	// Clean
	free_double_matrix(&dist, N);
	free_double_matrix(&Y, N);

	return status;
}
////////////////////////////////////////////////////////////////////////////////

