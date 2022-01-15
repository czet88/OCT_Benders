#include "headers.h"


/// INFO ///////////////////////////////////////////////////////////////////////
//
// master.c: Master Problem Functions
//
// Author:  Carlos Luna-Mota 
// Created: 2015 Jan 21
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int create_CPLEX_master_enviroment(CPXENVptr *env, CONFIGURATION param, GLOBAL_INFO global) {
	int status = 0;

	// Create enviroment:
	*env = CPXopenCPLEX(&status);
	assert(*env != NULL);

	// Set parameters:
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_THREADS,    1));						// Threads usados
    status = MAX(status, CPXsetdblparam(*env, CPX_PARAM_TILIM,		param.MAX_CPU_TIME));		// Time Limit
	status = MAX(status, CPXsetdblparam(*env, CPX_PARAM_EPGAP,		EPSILON*EPSILON));			// Gap de Epsilon Optimalidad
    status = MAX(status, CPXsetdblparam(*env, CPX_PARAM_CUTSFACTOR, 1.0));						// <= 1.0 No cuts will be generated, >1: Limits the number of cuts that can be added
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_MIPSEARCH,  CPX_MIPSEARCH_TRADITIONAL));// Turn on traditional search for use with control callbacks 
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_MIPCBREDLP, CPX_OFF));					// Let MIP callbacks work on the original model 
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_HEURFREQ,	-1));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_PRELINEAR, CPX_OFF));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_MIPINTERVAL, 1));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_PRESLVND, -1));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_PREIND, 0));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_REPEATPRESOLVE, 0));
	//status = MAX(status, CPXsetdblparam(*env, CPX_PARAM_CUTUP, (global.results->UpperBound)+0.01));
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_MIPEMPHASIS, 2));		//	0	CPX_MIPEMPHASIS_BALANCED		Balance optimality and feasibility; default
																				//	1	CPX_MIPEMPHASIS_FEASIBILITY		Emphasize feasibility over optimality
																				//	2	CPX_MIPEMPHASIS_OPTIMALITY		Emphasize optimality over feasibility
																				//	3	CPX_MIPEMPHASIS_BESTBOUND		Emphasize moving best bound
																				//	4	CPX_MIPEMPHASIS_HIDDENFEAS		Emphasize finding hidden feasible solutions	
	// The Lazy Constraints Callback will switch off these anyway:
    status = MAX(status, CPXsetintparam(*env, CPX_PARAM_REDUCE,		CPX_OFF));	// 0: No primal or dual reductions,   1: Only primal reductions ???
	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_PRELINEAR,	CPX_OFF));	// Assure linear mappings between the presolved and original models 
	

	status = MAX(status, CPXsetintparam(*env, CPX_PARAM_VARSEL,  3)); // 0: default, 1: max_infeas, 2: pseudo_cost, 3: strong_branching


	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int create_CPLEX_master_lp(CPXENVptr env, CPXLPptr *lp_ptr, GLOBAL_INFO global) {

	clock_t	start;
	int i, j, k, l, c, e, t;
	int added_cuts, status = 0;

	// Just for lazyness:
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int		**index_e = global.data->index_e;
	int		*index_i  = global.data->index_i;
	int		*index_j  = global.data->index_j;
	int		**index_k = global.data->index_k;
	double	**W = global.data->W;
	double	**C = global.data->C;
    double	**C1 = global.data->C1;
    double	**C2 = global.data->C2;
	COMMODITY *Com = global.data->Com;
	double longest_path_ub = LongestHamiltonianPathUB(N, C);

	// New data:
	double **Y = create_double_matrix(N,N);

	// 4-tree data structures:
	double best_tree_value, tree_value;
	double COST[6];
	int TREE_4[16][6] = {{1,1,1,0, 0,0}, {0,1,1,1, 0,0}, {1,0,1,1, 0,0}, {1,1,0,1, 0,0}, 
					     {1,0,1,0, 1,0}, {0,1,0,1, 0,1}, {1,0,1,0, 0,1}, {0,1,0,1, 1,0},
					     {1,0,0,0, 1,1}, {0,1,0,0, 1,1}, {0,0,1,0, 1,1}, {0,0,0,1, 1,1},
					     {1,1,0,0, 0,1}, {0,1,1,0, 1,0}, {0,0,1,1, 0,1}, {1,0,0,1, 1,0}};

    /// CPLEX MASTER POBLEM ////////////////////////////////////////////////////
    int       numcols; // número de variables ..................................
    int       numrows; // número de restricciones sin contar las cotas .........
    int       numnz;   // número de elementos no nulos de la matriz ............
    int       objsen;  // sentido de la optimizacion (min:1, max:-1 ) ..........
    double    *obj;    // coeficientes de la función objetivo ..................
    double    *rhs;    // términos independientes de las restricciones .........
    char      *sense;  // sentido de las restricciones (<=: 'L', =:'E', >=:'G'). 
    int       *matbeg; // índice del primer coeficiente no nulo de cada columna.
    int       *matcnt; // número de elementos no nulos de cada columna .........
    int       *matind; // fila a la que corresponde cada elemento no nulo  .....
    double    *matval; // valores de los coef no nulos en las restricciones ....
    double    *lb;     // cotas inferiores de las variables ....................
    double    *ub;     // cotas superiores de las variables ....................
    ////////////////////////////////////////////////////////////////////////////
	
	objsen  = CPX_MIN;		// Dirección de la optimización (minimizar!)
	numcols = E + K;		// Number of variables of the Master Prob (Y's & Z's)
	numrows = N + 1;		// Number of initial Constrains of the Master Problem (single_node cutsets & sum y = N-1)
	numnz   = N*(N-1) + E;	// Number of Non-Zero coefs (single_node cutsets & sum y = N-1)
		
	obj = create_double_vector(numcols);
    rhs = create_double_vector(numrows);
    sense = create_char_vector(numrows);
    matbeg = create_int_vector(numcols);
    matcnt = create_int_vector(numcols); 
    matind = create_int_vector(numnz); 
    matval = create_double_vector(numnz);
    lb = create_double_vector(numcols);
    ub = create_double_vector(numcols);
	
	numnz = 0;	// index of the current nonzero element
	for (i=0; i<N-1; i++) {             // Y_ij variables...
        for (j=i+1; j<N; j++) {	     // ...with i<j
            lb[index_e[i][j]]     = 0.0;	// lower bound
            ub[index_e[i][j]]     = 1.0;	// upper bound
			obj[index_e[i][j]]    = 0.0;	// Objective Function Coef.
            matcnt[index_e[i][j]] = 1 + 2;		// Number of constraints it appears
			matbeg[index_e[i][j]] = numnz;	// first NZ coeff related to this variable
            matind[numnz]      = 0;			// Cardinality constraint (first row)
            matval[numnz]      = 1.0;		// Coef
            numnz++;
            matind[numnz]      = i+1;      // Basic cut-set inequalities (i-row)
            matval[numnz]      = 1.0;      // Coef
            numnz++;
		    matind[numnz]      = j+1;      // Basic cut-set inequalities (j-row)
            matval[numnz]      = 1.0;      // Coef
            numnz++;
	    }
	}
	
	for(k=0; k<K; k++) {			// Z_k variables
		lb[E+k]     = 0.0;				// Lower Bound
		ub[E+k]     = longest_path_ub;	// Upper Bound (heuristic global upper bound on longest path)
		obj[E+k]    = Com[k].w;			// Obj. Fun. coeff
		matcnt[E+k] = 0;				// Number of constraints in which appear
		matbeg[E+k] = numnz;			// first NZ coeff related to this variable
	}
    
	// Cardinality constraint (first row)
	sense[0] = 'E';             // Sense of constraints (Equality)
    rhs[0]   = ((double)(N-1));  // Right Hand Side      (N-1)
		
	for (i=0; i<N; i++) {		// Basic cut-set inequalities (next N rows)
	    sense[i+1] = 'G';         // Sense of constraints (Greater or Equal)
        rhs[i+1]   = 1.0;         // Right Hand Side      (1)
	}

	assert(numnz == N*(N-1) + E);

	// Define the linear program
	*lp_ptr = CPXcreateprob(env, &status, "MASTER_PROBLEM");
	assert(*lp_ptr != NULL);

	// Load the linear program to Cplex enviroment
    status = CPXcopylp(env, *lp_ptr, numcols, numrows, objsen, obj, rhs, sense, matbeg, matcnt, matind, matval, lb, ub, NULL);
	assert(status == 0);

	////////////////////////////////////////////////////////////////////////////

	// Now we are going to add some additional cuts
	
	//	Add a Global Lower Bound:
	if (global.param->GLOBAL_LOWER_BOUND == YES) {
		numnz	 = 0;
		rhs[0]	 = global.results->LowerBound;
		sense[0] = 'G';		
		matbeg[0] = 0;
		for (k=0; k<K; k++) {			// Sum W_k Z_k >= LowerBound  OR  Z >= LowerBound
			matval[numnz] = Com[k].w;
			matind[numnz] = E+k;
			numnz++;
		}
		status = CPXaddrows(env, *lp_ptr, 0, 1, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
		assert(status == 0);
	}

	//	Add a Commodity Lower Bound:
	if (global.param->COMMODITY_LOWER_BOUND == YES && global.data->IS_EUCLIDEAN == YES) {
		for (k=0; k<K; k++) {			// Z_k + (C2_k - C1_k) Y_k >= C2_k
			numnz	 = 2;
			rhs[0]	 = C2[Com[k].o][Com[k].d];
			sense[0] = 'G';	
			matbeg[0] = 0;
			matval[0] = C2[Com[k].o][Com[k].d] - C1[Com[k].o][Com[k].d];
			matind[0] = index_e[Com[k].o][Com[k].d];
			matval[1] = 1.0;
			matind[1] = E+k;
			status = CPXaddrows(env, *lp_ptr, 0, 1, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
			assert(status == 0);
		}		
	}

	//	Add 3-trees Lower Bounds:
	if (global.param->TREE_3_LOWER_BOUND == YES && global.data->IS_EUCLIDEAN == YES) {
		for (c=0; c<K; c++) {
			i = Com[c].o;
			k = Com[c].d;
			for (j=0; j<N; j++) {
				if (j!= i && j!= k && W[i][j]+W[j][i] > 0 && W[j][k]+W[k][j] > 0) {
					// Add a 3-tree lower bound for each {i,j,k} such that {i,j}, {j,k} & {k,i} have communication requests
					numnz	  = 3;
					rhs[0]	  = C1[i][j] + C1[j][k] + C1[k][i] + MIN(C2[i][j]-C1[i][j], MIN(C2[j][k]-C1[j][k], C2[k][i]-C1[k][i]));
					sense[0]  = 'G';
					matbeg[0] = 0;
					matval[0] = 1.0;
					matval[1] = 1.0;
					matval[2] = 1.0;
					matind[0] = E+index_k[i][j];
					matind[1] = E+index_k[j][k];
					matind[2] = E+index_k[k][i];
					status = CPXaddrows(env, *lp_ptr, 0, 1, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
					assert(status == 0);
				}
			}
		}
	}

	//	Add 4-trees Lower Bounds:
	if (global.param->TREE_4_LOWER_BOUND == YES && global.data->IS_EUCLIDEAN == YES) {
		for (c=0; c<K; c++) {
			i = Com[c].o;
			l = Com[c].d;
			for (j=0; j<N-1; j++) {
				if (j!=i && j!=l && index_k[i][j]>=0 && index_k[j][l]>=0) {
					for (k=j+1; k<N; k++) {
						if (k!=i && k!=l && index_k[i][k]>=0 && index_k[j][k]>=0 && index_k[k][l]>=0) {
						
							// COMPUTE_ T^4_ijkl
							COST[0] = C2[i][j]-C1[i][j]; 	COST[1] = C2[j][k]-C1[j][k]; 	COST[2] = C2[k][l]-C1[k][l];
							COST[3] = C2[l][i]-C1[l][i]; 	COST[4] = C2[i][k]-C1[i][k]; 	COST[5] = C2[j][l]-C1[j][l];
							best_tree_value = 0.0;
							for (e=0; e<6; e++) { best_tree_value += ((1-TREE_4[0][e]) * COST[e]); } 
							for (t=1; t<16; t++) {
								tree_value = 0.0;
								for (e=0; e<6; e++) { tree_value += ((1-TREE_4[t][e]) * COST[e]); }
								best_tree_value = MIN(best_tree_value, tree_value);
							}

							// Add the corresponding 4-tree lower bound
							numnz	  = 6;
							rhs[0]	  = (C1[i][j] + C1[j][k] + C1[k][l] + C1[l][i] + C1[j][l] + C1[k][i] + best_tree_value);
							sense[0]  = 'G';
							matbeg[0] = 0;

							matval[0] = 1.0;
							matval[1] = 1.0;
							matval[2] = 1.0;
							matval[3] = 1.0;
							matval[4] = 1.0;
							matval[5] = 1.0;
							matind[0] = E+index_k[i][j];
							matind[1] = E+index_k[j][k];
							matind[2] = E+index_k[k][i];
							matind[3] = E+index_k[i][l];
							matind[4] = E+index_k[j][l];
							matind[5] = E+index_k[k][l];
							status = CPXaddrows(env, *lp_ptr, 0, 1, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
							assert(status == 0);
						}
					}
				}
			}
		}
	}

	// OPT CUTS FROM HEURISTIC SOLUTIONS
	if (global.param->MST_OPT_CUTS == YES) {		// Minimum Spanning Tree (w.r.t. matrix C)
		MST(N, global.data->C, Y, NULL);
		for (e=0; e<E; e++) { global.var->master[e]   = Y[index_i[e]][index_j[e]]; }
		for (k=0; k<K; k++) { global.var->master[E+k] = 0.0;                       }
		start  = clock();
		status = solve_optimality_subproblems(global, YES);
		assert(status == 0);
		global.results->opt_time += elapsed(start);
		global.results->sub_time += elapsed(start);
		added_cuts = add_OPT_root_cuts (global, env, *lp_ptr, YES);
		assert(added_cuts > 0);
	}

	if (global.param->GHT_OPT_CUTS == YES) {		// Gomory-Hu Tree (w.r.t. matrix W)
		GHT(N, global.data->W, Y);
		for (e=0; e<E; e++) { global.var->master[e]   = Y[index_i[e]][index_j[e]]; }
		for (k=0; k<K; k++) { global.var->master[E+k] = 0.0;                       }
		start  = clock();
		status = solve_optimality_subproblems(global, YES);
		assert(status == 0);
		global.results->opt_time += elapsed(start);
		global.results->sub_time += elapsed(start);
		added_cuts = add_OPT_root_cuts (global, env, *lp_ptr, YES);
		assert(added_cuts > 0);
	}
	if (global.param->AMT_OPT_CUTS == YES) {		// Ahuja-Murty Heuristic Tree (w.r.t. matrices C & W)
		AMT(N, global.data->C, global.data->W, Y, NULL);
		for (e=0; e<E; e++) { global.var->master[e]   = Y[index_i[e]][index_j[e]]; }
		for (k=0; k<K; k++) { global.var->master[E+k] = 0.0;                       }
		start  = clock();
		status = solve_optimality_subproblems(global, YES);
		assert(status == 0);
		global.results->opt_time += elapsed(start);
		global.results->sub_time += elapsed(start);
		added_cuts = add_OPT_root_cuts (global, env, *lp_ptr, YES);
		assert(added_cuts > 0);
	}
	if (global.param->MCMCT_OPT_CUTS == YES) {	// MinCut-MinCost Heuristic Tree (w.r.t. matrices MIN_CUT & MIN_CUT_MIN_COST)
		MCMCT(N, global.data->MIN_CUT, global.data->MIN_CUT_MIN_COST, Y, NULL);
		for (e=0; e<E; e++) { global.var->master[e]   = Y[index_i[e]][index_j[e]]; }
		for (k=0; k<K; k++) { global.var->master[E+k] = 0.0;                       }
		start  = clock();
		status = solve_optimality_subproblems(global, YES);
		assert(status == 0);
		global.results->opt_time += elapsed(start);
		global.results->sub_time += elapsed(start);
		added_cuts = add_OPT_root_cuts (global, env, *lp_ptr, YES);
		assert(added_cuts > 0);
	}
	if (global.param->STARS_OPT_CUTS == YES) {	// All Star-Shapped Trees
		for (c=0; c<N; c++) {
			for (e=0; e<E; e++) {				// Get the i-centered Star:
				if (index_i[e]==c || index_j[e]==c) { global.var->master[e] = 1.0; }
				else                                { global.var->master[e] = 0.0; }
			}
			for (k=0; k<K; k++) { global.var->master[E+k] = 0.0; }
			start  = clock();
			status = solve_optimality_subproblems(global, YES);
			assert(status == 0);
			global.results->opt_time += elapsed(start);
			global.results->sub_time += elapsed(start);
			added_cuts = add_OPT_root_cuts (global, env, *lp_ptr, YES);
			assert(added_cuts > 0);
		}
	}
	////////////////////////////////////////////////////////////////////////////

	// Remember initial number of constraints:
	global.results->initial_cons = CPXgetnumrows(env, *lp_ptr);
 
	// Clean:
	free_double_matrix(&Y,N);
	free_double_vector(&obj);
    free_double_vector(&rhs);
    free_char_vector(&sense);
    free_int_vector(&matbeg);
    free_int_vector(&matcnt); 
    free_int_vector(&matind); 
    free_double_vector(&matval);
    free_double_vector(&lb);
    free_double_vector(&ub);
	
	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int solve_root_node(CPXENVptr env, CPXLPptr lp, GLOBAL_INFO global) {
	clock_t start;
	int		status = 0;
	int		cuts_found = 0;
	int     imp_it=0;
	int		no_imp=0;
	double	old_LowerBound;
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	
	if (global.param->SCREEN_OUTPUT >= 2) { printf("\t\tNumRows: %d\n", CPXgetnumrows(env, lp)); }

	// Solve (Linear-Relaxed) Master Problem
	start = clock();
	status = CPXlpopt(env, lp);
	assert(status == 0);
	global.results->master_time += elapsed(start);
	
	// Obtain a valid Lower Bound for the integer problem
	old_LowerBound = 0.0;
	if (elapsed(global.var->init_time) < global.param->MAX_CPU_TIME) { status = CPXgetobjval(env, lp, &(global.results->LowerBound)); }
	assert(status == 0);

	// Update GAP
	global.results->final_gap = (global.results->UpperBound-global.results->LowerBound) / global.results->UpperBound;
	if (global.param->SCREEN_OUTPUT >= 1) { printf("Root\t\t\t%.2f <- %.2f%% -> %.2f  \tTime: %.2f\"\n", global.results->LowerBound, 100*global.results->final_gap, global.results->UpperBound, elapsed(global.var->init_time)); }
	
	// while we are improving 'enough':
	while (elapsed(global.var->init_time) < global.param->MAX_CPU_TIME && (ABS(global.results->LowerBound - old_LowerBound) > 0.0000001 || imp_it<1 || no_imp<3 )) {
	
		// Check if the OPT_GAP has beeen closed:
		if (global.results->UpperBound - global.results->LowerBound < global.param->MIN_ABS_GAP) { break; }
		
		// If we are using iterative Benders we shouldn't add any fractional cut!!!:
		if (global.param->ALGORITHM == 0) { break; }
		
		// Solve the subproblem and add some cuts:
		status = CPXgetx(env, lp, global.var->master, 0, E+K-1);
		
		// First look for feasibility cuts:
		start  = clock();
		status = solve_feasibility_subproblem(global, NO);
		assert(status == 0);
		global.results->feas_time += elapsed(start); 
		global.results->sub_time  += elapsed(start);
		
		// Then try to add them to the master problem:
		cuts_found = add_FEAS_cuts(global, env, lp, NO);
		if (global.param->SCREEN_OUTPUT >= 2) { printf("\t\t%d violated feasibility cuts found\n", cuts_found); }

		// If the problem is feasible:
		if (cuts_found == 0) {
		
			// Look for optimality cuts:
			global.var->flag_generate_cuts = YES;
			start  = clock();
			status = solve_optimality_subproblems(global, NO);
			assert(status == 0);
			global.results->opt_time += elapsed(start);
			global.results->sub_time += elapsed(start);
		
			// And try to add them to the master problem:
			cuts_found += add_OPT_root_cuts(global, env, lp, NO);
			if (global.param->SCREEN_OUTPUT >= 2) { printf("\t\t%d violated optimality cuts found\n", cuts_found); }
		}

		if (cuts_found == 0) { break; }
		
		// Solve the new (Linear-Relaxed) Master Problem
		start = clock();
		status = CPXlpopt(env, lp);
		assert(status == 0);
		global.results->master_time += elapsed(start);

		// Obtain a valid Lower Bound for the integer problem
		old_LowerBound = global.results->LowerBound;
		if (elapsed(global.var->init_time) < global.param->MAX_CPU_TIME) { status = CPXgetobjval(env, lp, &(global.results->LowerBound)); }
		if(global.results->LowerBound-old_LowerBound>0.1){imp_it++;} //We know now that we've improved one more.
		if(imp_it>=1 && global.results->LowerBound-old_LowerBound<0.1) no_imp++;
		assert(status == 0);
		
		// Update GAP
		global.results->final_gap = (global.results->UpperBound-global.results->LowerBound) / global.results->UpperBound;
		if (global.param->SCREEN_OUTPUT >= 1) { printf("Root\t\t\t%.2f <- %.2f%% -> %.2f  \tTime: %.2f\"\n", global.results->LowerBound, 100*global.results->final_gap, global.results->UpperBound, elapsed(global.var->init_time)); }
	} 	
	FinalLP=global.results->LowerBound;
	// Remember GAP
	global.results->root_gap = global.results->final_gap;
	
	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int clean_root_node(CPXENVptr env, CPXLPptr lp, GLOBAL_INFO global) {

	int     row, status = 0;
	int		N = global.data->N;
	double  slack_threshold, old_LowerBound = 0.0;

	int		numrows;
	double  *slack;
	char	*sense;
	
	numrows = CPXgetnumrows(env, lp);
	slack = create_double_vector(numrows);
	sense = create_char_vector(numrows);

	// Get Constraint Senses:
	status = CPXgetsense(env, lp, sense, 0, numrows-1);
	assert(status==0);

	// Get Constraint Slacks:
	status = CPXgetslack(env, lp, slack, 0, numrows-1);
	assert(status==0);

	// Amend slack vector so a positive slack means a non-binding constraint
	for (row=0; row<numrows; row++) { if (sense[row] == 'G') { slack[row] *= -1.0; } }

	// Get max acceptable slack:
	slack_threshold = MAX(EPSILON, find_smallest(numrows, slack, ((int)(((double)(numrows))*global.param->ROOT_CLEANUP))));
	
	// We are going to check all rows that don't form part of the base model (= the first N+1 rows)
	for (row=numrows-1; row>=N+1; row--) {		// We must erase backwards to avoid loosing the indices!
		if (slack[row] > slack_threshold) {		// We only erase non-binding constraints !
			 status = CPXdelrows(env, lp, row, row);	// One at a time.
			 assert(status == 0);
		}
	}
	
	// Inform:
	global.results->root_cons = CPXgetnumrows(env, lp);
	if (global.param->SCREEN_OUTPUT >= 1) { 	printf("\t%d rows have been deleted from master problem\n", numrows-global.results->root_cons); }	

	// Free data structures:
	free_double_vector(&slack);
	free_char_vector(&sense);
	
	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int create_CPLEX_master_milp(CPXENVptr env, CPXLPptr lp, GLOBAL_INFO global) {

	int  e, k;
	int  status = 0;
	int	 E = global.data->E;
	int	 K = global.data->K;
	char *ctype	= create_char_vector(E+K);
	int *beg = create_int_vector(1);
	int *varindices = create_int_vector(E+K);
	int *effortlevel = create_int_vector(1);

	beg[0] = 0;
	effortlevel[0] = 4;

	// Y_e
	for (e=0; e<E; e++) { 
		ctype[e] = 'B'; // type of variable (Binary)
	}	

	for (e = 0; e < E + K; e++) {
		varindices[e] = e;
	}
            
	// Z_k
	for (k=0; k<K; k++) { ctype[E+k] = 'C'; }	// type of variable (Continous)

	// Now tell CPLEX what kind of variables do we have:
	status = CPXcopyctype(env, lp, ctype);               
	assert(status == 0);

	// Unless we are using CPLEX as a black box for the Master Problem...
	if (global.param->ALGORITHM == 1) {		// Branch & Cut Benders
		
		printf("Current upper bound %.lf\n", global.results->UpperBound);
		// Provide the initial solution found so far
		//status = CPXaddmipstarts(env, lp, 1, E, beg, varindices, global.results->best_solution, effortlevel, NULL);
		// Provide a Lazy Constraint Callback Function
		status = CPXsetlazyconstraintcallbackfunc(env, lazy_cut_callback, &global);
		assert(status == 0);
		
		// Provide a User Cut Callback Function
		status = CPXsetusercutcallbackfunc(env, user_cut_callback, &global);
		assert(status == 0);

		if (global.param->USE_BRANCHING_CALLBACK == YES) {
			// Provide a User Cut Callback Function
			status = CPXsetbranchcallbackfunc(env, set_branch_callback, &global);
			assert(status == 0);
		}
	}

	// Provide heuristic callback for rounding fractional solutions:
	if (global.param->USE_ROUNDING_CALLBACK == YES) {
		status = CPXsetheuristiccallbackfunc(env, rounding_callback, &global);
		assert(status == 0);
	}	
		
	// clean:
	free_char_vector(&ctype);
	
	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int solve_to_integrality(CPXENVptr env, CPXLPptr lp, GLOBAL_INFO global) {
	int		status = 0;
	int		E = global.data->E;
	int		K = global.data->K;
	int		cuts_found = 0;
	clock_t	start;

	// Update GAP
	global.results->final_gap = (global.results->UpperBound-global.results->LowerBound) / global.results->UpperBound;
	if (global.param->SCREEN_OUTPUT >= 1) { printf("Master\t\t\t%.2f <- %.2f%% -> %.2f  \tTime: %.2f\"\n", global.results->LowerBound, 100*global.results->final_gap, global.results->UpperBound, elapsed(global.var->init_time)); }

	// while we haven't closed the GAP:
	while (elapsed(global.var->init_time) < global.param->MAX_CPU_TIME 
		   && global.results->UpperBound - global.results->LowerBound > global.param->MIN_ABS_GAP) {

		// Solve Master Problem
	    CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);
		CPXsetintparam(env,CPX_PARAM_MIPDISPLAY,3); //different levels of output display
		// New stuff added
		/*CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_AbsMIPGap, global.param->MIN_ABS_GAP);
		CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_Integrality, 0.0005);
		CPXsetintparam(env, CPXPARAM_MIP_Strategy_LBHeur, CPX_ON);
		CPXsetintparam(env, CPXPARAM_MIP_Strategy_HeuristicFreq, 5);
		CPXsetintparam(env, CPXPARAM_MIP_Strategy_FPHeur, 1);*/

		
		start = clock();
		status = CPXmipopt(env, lp);
		assert(status == 0);
		global.results->master_time += elapsed(start);
		
		if (CPXgetstat(env, lp) == 103) {				// Infeasible Master
			if (global.param->UPPER_CUTOFF == NO) {		// Maybe because of the Upper CutOff?
				printf("\nERROR: Infeasible Master Problem\n");
				return 1;
			} else {
				printf("\nWARNING: Infeasible Master Problem. Assuming UpperCutOff as a cause.\n");
				global.results->LowerBound = global.results->UpperBound;
			}
		} else {

			// Obtain a valid Lower Bound for the integer problem
			status = CPXgetbestobjval(env, lp, &(global.results->LowerBound));
			assert(status == 0);

			global.results->n_explored_nodes = CPXgetnodecnt(env, lp);
			// Update GAP
			global.results->final_gap = (global.results->UpperBound-global.results->LowerBound) / global.results->UpperBound;

		}

		// Report current state:
		if (global.param->SCREEN_OUTPUT >= 1) { 	printf("Master\t\t\t%.2f <- %.2f%% -> %.2f  \tTime: %.2f\"\n", global.results->LowerBound, 100*global.results->final_gap, global.results->UpperBound, elapsed(global.var->init_time)); }

		// Remember final constraints and explored nodes
		global.results->final_cons = CPXgetnumrows(env, lp);

		// Now we must distinguish between iterative and B&C benders variants:
		if      (elapsed(global.var->init_time) > global.param->MAX_CPU_TIME)                          { break; }
		else if (global.results->UpperBound - global.results->LowerBound <= global.param->MIN_ABS_GAP) { break; }
		else {

			if (global.param->ALGORITHM == 1) { global.var->flag_fist_master = NO; }

			// Obtain the current master problem solution:
			status = CPXgetx(env, lp, global.var->master, 0, E+K-1);
			assert(status == 0);
			 
			// First look for feasibility cuts:
			start = clock();
			status = solve_feasibility_subproblem(global, YES);
			assert(status == 0);
			global.results->feas_time += elapsed(start); 
			global.results->sub_time += elapsed(start);

			// Then try to add them to the master problem:
			cuts_found = add_FEAS_cuts(global, env, lp, YES);
	
			if (global.param->SCREEN_OUTPUT >= 2) { printf("\t\t%d violated feasibility cuts found\n", cuts_found); }

			// If the problem is feasible:
			if (cuts_found == 0) {

				// Look for optimality cuts:
				start = clock();
				status = solve_optimality_subproblems(global, YES);
				assert(status == 0);
				global.results->opt_time += elapsed(start);
				global.results->sub_time += elapsed(start);
		
				// And try to add them to the master problem:
				if (global.var->flag_generate_cuts == YES) { cuts_found += add_OPT_root_cuts(global, env, lp, YES); }

				if (global.param->SCREEN_OUTPUT >= 2) { printf("\t\t%d violated optimality cuts found\n", cuts_found); }
			}

			if (global.param->SCREEN_OUTPUT >= 2 && cuts_found == 0) { printf("WARNING: 0 cuts have been added this iteration\n"); break; }		
		}
	}

	return status;
}
////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////
int CPXPUBLIC lazy_cut_callback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p) {
	
	clock_t start;
	GLOBAL_INFO global = *((GLOBAL_INFO_ptr) cbhandle);
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int		**index_e = global.data->index_e;
	double	*ub;
	double	*lb;
	int		previous_cuts, cuts_found = 0;
	int		status = 0;
	int		NODES_LEFT, SEQNUM, depth;

	/// GET INFORMATION: ///////////////////////////////////////////////////
	// Sometimes the cut_callback is called in a situation where: CPX_CALLBACK_INFO_NODES_LEFT == 0.
	// That makes the program crash since you cannot call CPXgetcallbacknodeinfo if NODES_LEFT == 0.
	// For those rare cases this function just returns gracefully without doing anything: 
	status = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODES_LEFT, &NODES_LEFT);
	assert(status == 0);
	if (NODES_LEFT == 0) { *useraction_p = CPX_CALLBACK_DEFAULT; return 0; }

	// We must be sure that we are in a integer solution:
	if (wherefrom != CPX_CALLBACK_MIP_CUT_FEAS) { printf("ERROR: UNBOUNDED NODE FOUND\n"); return 1; }		

	// Get callback depth
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
	assert(status == 0);
		
	// Get the Node Id
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &SEQNUM);
	assert(status == 0);

	// Remember the Node Id
	if (global.var->last_id == SEQNUM) { global.var->last_id_visits++; }
	else { global.var->last_id = SEQNUM; global.var->last_id_visits = 0; }
	////////////////////////////////////////////////////////////////////////


	if (elapsed(global.var->init_time) < global.param->MAX_CPU_TIME) {
			
		/// GET MORE INFORMATION: //////////////////////////////////////////////
			
		// Get current Master Problem Solution:
		status = CPXgetcallbacknodex (env, cbdata, wherefrom, global.var->master, 0, E+K-1);
		assert(status == 0);

		// Get global LowerBound
		status = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_REMAINING, &(global.results->LowerBound));
		assert(status == 0);
		///////////////////////////////////////////////////////////////////////


		// SOLVE FEASIBILITY SUBPROBLEM ///////////////////////////////////////
		// First look for feasibility cuts:
		start = clock();
		status = solve_feasibility_subproblem(global, YES);
		assert(status == 0);
		global.results->feas_time += elapsed(start); 
		global.results->sub_time  += elapsed(start);
		global.results->master_time -= elapsed(start);
		
		// Then try to add them to the master problem:
		cuts_found = add_FEAS_callback_cuts(global, env, cbdata, wherefrom, YES);

		// Report:
		if (global.param->SCREEN_OUTPUT >= 2) { 	printf("\t\t%d violated feasibility cuts found\n", cuts_found); }
		///////////////////////////////////////////////////////////////////////
		
		// SOLVE OPTIMALITY SUBPROBLEM ////////////////////////////////////////
		// If the problem is feasible:
		if (cuts_found == 0 && elapsed(global.var->init_time) < global.param->MAX_CPU_TIME) {

			// Look for regular optimality cuts:
			global.var->flag_generate_cuts = YES;
			start = clock();
			status = solve_optimality_subproblems(global, YES);
			assert(status == 0);
			global.results->opt_time += elapsed(start);
			global.results->sub_time += elapsed(start);
			global.results->master_time -= elapsed(start);
		
			// And try to add them to the master problem:
			previous_cuts = cuts_found;
			if (global.var->flag_generate_cuts == YES){ cuts_found += add_OPT_callback_cuts(global, env, cbdata, wherefrom, YES); }
			else { update_core_point(global, NULL); }

			// Report:
			if (global.param->SCREEN_OUTPUT >= 2) { 	printf("\t\t%d violated optimality cuts found\n", cuts_found-previous_cuts); }
		}
		///////////////////////////////////////////////////////////////////////


		// LOCAL CUTS /////////////////////////////////////////////////////////
		if (global.param->ADD_LOCAL_CUTS == YES ){ // && global.param->USE_BRANCHING_CALLBACK == NO && elapsed(global.var->init_time) < global.param->MAX_CPU_TIME) {
			
			// Allocate Memory:
			lb = create_double_vector(E+K);
			ub = create_double_vector(E+K);
		
			// Get Information:
			status = CPXgetcallbacknodelb(env, cbdata, wherefrom, lb, 0, E+K-1);
			assert(status == 0);
			status = CPXgetcallbacknodeub(env, cbdata, wherefrom, ub, 0, E+K-1);
			assert(status == 0);

			// Add Local Cuts:
			previous_cuts = cuts_found;
			start = clock();
			cuts_found += add_local_cuts(global, env, cbdata, wherefrom, lb, ub);
			global.results->sub_time += elapsed(start);
			global.results->master_time -= elapsed(start);

			// Report:
			if (global.param->SCREEN_OUTPUT >= 2) { 	printf("\t\t%d violated local cuts found\n", cuts_found-previous_cuts); }

			// Clean:
			free_double_vector(&lb);
			free_double_vector(&ub);
		}
		///////////////////////////////////////////////////////////////////////

		if (cuts_found == 0) { printf("WARNING: 0 cuts have been added in this cut_callback call\n"); }
		else if (global.param->SCREEN_OUTPUT >= 2) { printf("\t\t%d cuts have been added during this lazy_cut_callback call\n", cuts_found); }

		// Update GAP
		global.results->final_gap = (global.results->UpperBound-global.results->LowerBound) / global.results->UpperBound;
		if (global.param->SCREEN_OUTPUT >= 1) { 	printf("Lazy #%d depth:%d \t%.2f <- %.2f%% -> %.2f  \tTime: %.2f\"\n", SEQNUM, depth, global.results->LowerBound, 100*global.results->final_gap, global.results->UpperBound, elapsed(global.var->init_time)); }
	}	

	// Tell CPLEX if we had found new cuts:
	if (cuts_found == 0) { *useraction_p = CPX_CALLBACK_DEFAULT; }
	else                 { *useraction_p = CPX_CALLBACK_SET;     }

	return status;
}
////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////
int CPXPUBLIC user_cut_callback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p) {
	
	clock_t start;
	GLOBAL_INFO global = *((GLOBAL_INFO_ptr) cbhandle);
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int		**index_e = global.data->index_e;
	
	double	*lb;
	double	*ub; 

	int		cuts_found = 0;
	int		previous_cuts;
	int		status = 0;
	
	int		NODES_LEFT, SEQNUM, depth;
	
	/// GET INFORMATION: ///////////////////////////////////////////////////
	// Sometimes the cut_callback is called in a situation where: CPX_CALLBACK_INFO_NODES_LEFT == 0.
	// That makes the program crash since you cannot call CPXgetcallbacknodeinfo if NODES_LEFT == 0.
	// For those rare cases this function just returns gracefully without doing anything: 
	status = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODES_LEFT, &NODES_LEFT);
	assert(status == 0);
	if (NODES_LEFT == 0) { *useraction_p = CPX_CALLBACK_DEFAULT; return 0; }
	
	// Make sure that the node is not unbounded:
	if (wherefrom == CPX_CALLBACK_MIP_CUT_UNBD) { printf("ERROR: UNBOUNDED NODE FOUND\n"); return 1; }		

	// Get callback depth
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
	assert(status == 0);
		
	// Get the Node Id
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &SEQNUM);
	assert(status == 0);
	////////////////////////////////////////////////////////////////////////

	if (elapsed(global.var->init_time) < global.param->MAX_CPU_TIME 
		&& (wherefrom == CPX_CALLBACK_MIP_CUT_LOOP || wherefrom == CPX_CALLBACK_MIP_CUT_LAST) 
		&& (global.param->ADD_FRAC_FEAS_CUTS == YES || global.param->ADD_FRAC_OPT_CUTS == YES || global.param->ADD_LOCAL_CUTS == YES)
		&& global.var->last_id_visits < 3 && (depth%global.param->FRAC_CUTS_FREQ==0 /*SEQNUM%100==0|| (depth>10 && depth<12) */) ) { 

		// Remember the Node Id ////////
		if(global.var->last_id == SEQNUM) { 	global.var->last_id_visits++; }
		else { global.var->last_id = SEQNUM; global.var->last_id_visits = 0; }
		
			
		/// GET MORE INFORMATION: //////////////////////////////////////////////
			
		// Get current Master Problem Solution:
		status = CPXgetcallbacknodex (env, cbdata, wherefrom, global.var->master, 0, E+K-1);
		assert(status == 0);

		// Get global LowerBound
		status = CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_REMAINING, &(global.results->LowerBound));
		assert(status == 0);
		///////////////////////////////////////////////////////////////////////

		
		// SOLVE FEASIBILITY SUBPROBLEM ///////////////////////////////////////
		if (global.param->ADD_FRAC_FEAS_CUTS == YES) {
			// First look for feasibility cuts:
			start = clock();
			status = solve_feasibility_subproblem(global, NO);
			assert(status == 0);
			global.results->feas_time += elapsed(start); 
			global.results->sub_time  += elapsed(start);
			global.results->master_time -= elapsed(start);
		
			// Then try to add them to the master problem:
			cuts_found = add_FEAS_callback_cuts(global, env, cbdata, wherefrom, NO);

			// Report:
			if (global.param->SCREEN_OUTPUT >= 2) { 	printf("\t\t%d violated feasibility cuts found\n", cuts_found); }
		}
		///////////////////////////////////////////////////////////////////////



		// SOLVE OPTIMALITY SUBPROBLEM ////////////////////////////////////////
		// If the problem is feasible:
		if (global.param->ADD_FRAC_OPT_CUTS == YES && cuts_found == 0 
			&& elapsed(global.var->init_time) < global.param->MAX_CPU_TIME) {

			// Look for regular optimality cuts:
			global.var->flag_generate_cuts = YES;
			start = clock();
			status = solve_optimality_subproblems(global, NO);
			assert(status == 0);
			global.results->opt_time += elapsed(start);
			global.results->sub_time += elapsed(start);
			global.results->master_time -= elapsed(start);
		
			// And try to add them to the master problem:
			previous_cuts = cuts_found;
			if (global.var->flag_generate_cuts == YES){ cuts_found += add_OPT_callback_cuts(global, env, cbdata, wherefrom, NO); }
			else { update_core_point(global, NULL); }
			
			// Report:
			if (global.param->SCREEN_OUTPUT >= 2) { 	printf("\t\t%d violated optimality cuts found\n", cuts_found-previous_cuts); }
		}
		///////////////////////////////////////////////////////////////////////

		// LOCAL CUTS /////////////////////////////////////////////////////////
		if (global.param->ADD_LOCAL_CUTS == YES){ //&& global.param->USE_BRANCHING_CALLBACK == NO) {
			
			// Allocate Memory:
			lb = create_double_vector(E+K);
			ub = create_double_vector(E+K);
		
			// Get Information:
			status = CPXgetcallbacknodelb(env, cbdata, wherefrom, lb, 0, E+K-1);
			assert(status == 0);
			status = CPXgetcallbacknodeub(env, cbdata, wherefrom, ub, 0, E+K-1);
			assert(status == 0);

			// Add Local Cuts:
			previous_cuts = cuts_found;
			start = clock();
			cuts_found += add_local_cuts(global, env, cbdata, wherefrom, lb, ub);
			global.results->sub_time += elapsed(start);
			global.results->master_time -= elapsed(start);

			// Report:
			if (global.param->SCREEN_OUTPUT >= 2) { 	printf("\t\t%d violated local cuts found\n", cuts_found-previous_cuts); }

			// Clean:
			free_double_vector(&lb);
			free_double_vector(&ub);
		}
		///////////////////////////////////////////////////////////////////////

		if (global.param->SCREEN_OUTPUT >= 2) { 	printf("\t\t%d cuts have been added during this user_cut_callback call\n", cuts_found); }

		// Update GAP
		global.results->final_gap = (global.results->UpperBound-global.results->LowerBound) / global.results->UpperBound;
		if (global.param->SCREEN_OUTPUT >= 1) {	printf("User #%d depth:%d \t%.2f <- %.2f%% -> %.2f  \tTime: %.2f\"\n", SEQNUM, depth, global.results->LowerBound, 100*global.results->final_gap, global.results->UpperBound, elapsed(global.var->init_time)); }
	}	

	// Tell CPLEX if we had found new cuts:
	if (cuts_found == 0) { *useraction_p = CPX_CALLBACK_DEFAULT; }
	else                 { *useraction_p = CPX_CALLBACK_SET;     }

	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int CPXPUBLIC rounding_callback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, double *objval_p, double *x, int *checkfeas_p, int *useraction_p) {
   
	int		status = 0;
    int     i,j,e,k;
	double  aux, obj_value;
	int		depth, SEQNUM;
	
	GLOBAL_INFO global = *((GLOBAL_INFO_ptr) cbhandle);
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int		**index_e  = global.data->index_e;
	
	double	*lb    = NULL;
	double	*ub    = NULL;
	int		**info = NULL;

	
	// The first time that we enter here we must give to CPLEX the best heuristic solution found so far:
	if(global.var->flag_initial_sol == NO) {

		// OJO: podria ser que no fuera compatible con la 'info' que tenemos...

		for(e=0; e<E+K; e++) { x[e] = global.results->best_solution[e]; }
		(*objval_p) = global.results->UpperBound;
		global.var->flag_initial_sol = YES;		// Recuerda que YA hemos entrado aquí antes.

	} else {

		// Allocate Memory
		lb   = create_double_vector(E);
		ub   = create_double_vector(E);
		info = create_int_matrix(N, N);
		
		// Get callback depth
		status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
		assert(status == 0);
		
		// Get the Node Id
		status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &SEQNUM);
		assert(status == 0);

		// Get local LowerBound for each variable: 
		status = CPXgetcallbacknodelb(env, cbdata, wherefrom, lb, 0, E-1);
		assert(status == 0);
	
		// Get local UpperBound for each variable: 
		status = CPXgetcallbacknodeub(env, cbdata, wherefrom, ub, 0, E-1);
		assert(status == 0);

		// Build local info matrix:
		for (i=0; i<N; i++) {
			info[i][i] = 0;
			for (j=i+1; j<N; j++) {
				if      (ub[index_e[i][j]] < EPSILON)     { info[i][j] = info[j][i] =  0; }	// Variable fixed to 0
				else if (lb[index_e[i][j]] > 1.0-EPSILON) { info[i][j] = info[j][i] =  1; }	// Variable fixed to 1
				else									  { info[i][j] = info[j][i] = -1; }	// Variable with unknown value
			}
		}

		// Put the Master Variables in a Matrix:
		for (i=0; i<N; i++) {
			global.var->capacity[i][i] = 0.0;
			for (j=i+1; j<N; j++) {
				global.var->capacity[i][j] = global.var->capacity[j][i] = x[global.data->index_e[i][j]];
			}
		}

		// Round Fractional Solution:
		rounding_heuristic(N, global.var->capacity, info);
		tree_dist(N, global.data->C, global.var->capacity, global.var->dist);
	
		// Try to improve the integer solution:
		if (global.param->APPLY_AM_LOCAL_SEARCH == YES) { aux = AM_LS(N, global.data->C, global.data->W, global.var->capacity, global.var->dist, info); 	}

		// Return the new solution to CPLEX:
		for (e=0; e<E; e++) { x[e] = (global.var->capacity[global.data->index_i[e]][global.data->index_j[e]] > 0.5) ? (1.0) : (0.0); }

		// Evaluate the new integer solution:
		obj_value = 0.0;
		for (k=0; k<K; k++) {
			x[E+k]     = global.var->dist[global.data->Com[k].o][global.data->Com[k].d];		// Z_k
			obj_value += x[E+k] * global.data->Com[k].w;										// \Sum_k w^k Z_k
		}
		(*objval_p) = obj_value;

		// For every integer subproblem we must also check if the current solution is the best one and update the gap
		if (global.results->UpperBound > obj_value) {
			global.results->UpperBound = obj_value;
			for (e=0; e<E+K; e++) { global.results->best_solution[e] = x[e]; }
			// Update GAP
			global.results->final_gap = (global.results->UpperBound-global.results->LowerBound) / global.results->UpperBound;
			if (global.param->SCREEN_OUTPUT >= 1) { 	printf("Heur #%d depth:%d \t%.2f <- %.2f%% -> %.2f  \tTime: %.2f\"\n", SEQNUM, depth, global.results->LowerBound, 100*global.results->final_gap, global.results->UpperBound, elapsed(global.var->init_time)); }
		}	

		// Clean:
		free_int_matrix(&info, N);
		free_double_vector(&lb);
		free_double_vector(&ub);

	}

	// Tell CPLEX that the solution is already integer:
	*checkfeas_p = 0;

	// Tell CPLEX that a solution is being returned:
	*useraction_p = CPX_CALLBACK_SET;

	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int CPXPUBLIC set_branch_callback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt, int bdcnt, const int *nodebeg, const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p) {
	
	char	pause = 0;
	int		B_DEBUG = NO;

	int		i,ii,jj,j,k,kk,e,ee,d,c,status = 0;

	GLOBAL_INFO global = *((GLOBAL_INFO_ptr) cbhandle);
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	double  **C = global.data->C;
	double  **dist = global.var->dist;
	int		**index_e = global.data->index_e;
	int		**index_k = global.data->index_k;
	int		*index_i = global.data->index_i;
	int		*index_j = global.data->index_j;

	double	objval;
	
	double	*lb	   = create_double_vector(E+K);
	double	*ub	   = create_double_vector(E+K);
	int     *feas  = create_int_vector(E);
	int     *prior = create_int_vector(E);
	
	int     size, n_cc, n_edges, max_size, big_cc;
	int		best_e, best_size, best_prior;
	int		old_fixed_vars, new_fixed_vars;
	int     *sizes = create_int_vector(N);
    int     *stack = create_int_vector(N);
	int     *cc    = create_int_vector(N);
	int		**neighbors = create_int_matrix(N,N);
	int		*degree     = create_int_vector(N);
 	
	char    *sense = create_char_vector(E+K);
	double  *bound = create_double_vector(E+K);
	int		*index = create_int_vector(E+K);

	int     seqnum1, seqnum2, NODES_LEFT, SEQNUM, depth;	
	
	////////////////////////////////////////////////////////////////////////////

	if (B_DEBUG == YES) {
		printf("\n\nBRANCHING CALLBACK:\n");
		printf("\tGetting Information:\n");
	}

	// Get information: ////////////////////////////////////////////////////////
 	
	status = CPXgetcallbacknodeintfeas(env, cbdata, wherefrom, feas, 0, E-1);
	assert(status == 0);

	e=0;
	for (i=0; i<E; i++) { if (feas[i] == CPX_INTEGER_INFEASIBLE) { e++; } }
	if (B_DEBUG == YES) {
		status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
		assert(status == 0);
		status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &SEQNUM);
		assert(status == 0);
		printf("\t\tNode:      %d\n\t\tDepth:     %d\n", SEQNUM, depth);
		if (e==0) { printf("\t\tThis is an integer node!\n"); }
		else      { printf("\t\t%d integer variables have fractional values.\n", e); }
	}

	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	assert(status == 0);

	status = CPXgetcallbacknodelb(env, cbdata, wherefrom, lb, 0, K+E-1);
	assert(status == 0);

	status = CPXgetcallbacknodeub(env, cbdata, wherefrom, ub, 0, K+E-1);
	assert(status == 0); 
	
	status = CPXgetcallbackorder(env, cbdata, wherefrom, prior, NULL, 0, E-1);
	assert(status == 0);

	if (B_DEBUG == YES) { //Counting how many 
		jj = 0;
		for (e=0; e<E; e++) {
			if ((ub[e] < 1.0-2*EPSILON && ub[e]> 2*EPSILON) || (lb[e] < 1.0-2*EPSILON && lb[e]> 2*EPSILON)) { printf("\t\tERROR: Var(%d) -> %.3f <= %.3f\n", e, lb[e], ub[e]); }
			if ((ub[e] < 2*EPSILON) || (lb[e]>1.0-2*EPSILON)) { 	jj++; }
		}
		printf("\t\tThere are %d edges with LB[e]=1 or UB[e]=0\n", jj);
	}
	////////////////////////////////////////////////////////////////////////////
	
	
	// PRE-PROCESS /////////////////////////////////////////////////////////////
    
	// Get linked list version of the graph induced by the variables fixed to 1:
	for (i=0; i<N; i++) { degree[i] = 0; }
	for (i=0; i<N-1; i++) {
	    for (j=i+1; j<N; j++) {
            if (lb[index_e[i][j]] >= 0.5 && lb[index_e[j][i]] >= 0.5) {
				neighbors[i][degree[i]++] = j;
				neighbors[j][degree[j]++] = i;
            }
        } 
    }

    // Find Connected Components:
    for (i=0; i<N; i++) { cc[i] = NONE; }
	for (i=0; i<N; i++) { sizes[i] = 0; }
    n_cc = 0;
	max_size = 0;
    for (i=0; i<N; i++) {
		if (cc[i] == NONE) {
			cc[i] = n_cc;
			sizes[n_cc]++;
			size = 0;
			stack[size++] = i;
			while (size > 0) {
				j = stack[--size];
				for (d=0; d<degree[j]; d++) {
					k = neighbors[j][d];
					if (cc[k] == NONE) {
						cc[k] = n_cc;
						sizes[n_cc]++;
						stack[size++] = k;
					} else if (cc[k] != n_cc) { printf("ERROR: vertex %d already visited!\n", k); }
				}
			} 
			max_size = MAX(max_size, sizes[n_cc]);
			n_cc++;
		}
    }

	// Find Distances:
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { dist[i][j] = -1.0; } }
    for (i=0; i<N; i++) {
		size = 0;
		stack[size++] = i;
		dist[i][i] = 0.0;
		while (size > 0) {
			j = stack[--size];
			for (d=0; d<degree[j]; d++) {
				k = neighbors[j][d];
				if (dist[i][k] < -0.5) {
					dist[i][k] = dist[i][j] + C[j][k];
					stack[size++] = k;
				}
			}
		} 
    }

	if (B_DEBUG == YES) {
		printf("\tProcessing Information:\n");
		printf("\t\tAt least 1 Connected Component with %d vertices\n", max_size);
		if (max_size > 1) { printf("\t\tConnected Components with 2 or more vertices:\n"); }
		for (i=0; i<n_cc; i++) { 
			ii = 0;
			if (sizes[i] > 1) { printf("\t\t\t{ "); }
			for (j=0; j<N; j++) { if(cc[j]==i) { if (sizes[i] > 1){ printf("%d ", j); } ii++; } }
			if (sizes[i] > 1){ printf("}\n"); }
			if (sizes[i] != ii) { printf("ERROR: Size shoud be %d and is %d.\n", sizes[i], ii); }
		}
	}

	// Choose the branching variable:
	best_e = NONE;
	best_size = NONE;
	best_prior = NONE;
	for (e=0; e<E; e++) {
		if(feas[e] == CPX_INTEGER_INFEASIBLE) {
			ii = cc[index_i[e]];
			jj = cc[index_j[e]];
			size = 0;
			for (i=0; i<N-1; i++) {
				for (j=i+1; j<N; j++) {
					if ((cc[i]==ii && cc[j]==jj) || (cc[j]==ii && cc[i]==jj)){
						ee = index_e[i][j];
						kk = index_k[i][j];
						if (lb[ee]<0.5 && ub[ee]>0.5) { size++; }
						if (kk != NONE && ub[kk]-lb[kk] > 2*EPSILON) { size++; }
					}
				}
			}

			// Determine if this is the best variable-branching found so far:
			if (best_e == NONE) {
				best_e     = e;
				best_size  = size;
				best_prior = prior[e];
				if (B_DEBUG == YES) { printf("\tFinding best candidate:\n\t\tVAR(%d) joining components %d and %d into a component of size %d\n",e,ii,jj,size); }
			} else if (best_size < size) {
				best_e     = e;
				best_size  = size;
				best_prior = prior[e];
				if (B_DEBUG == YES) { printf("\t\tVAR(%d) joining components %d and %d into a component of size %d\n",e,ii,jj,size); }
			}  else if (best_size == size && prior[e] > best_prior) {
				best_e     = e;
				best_size  = size;
				best_prior = prior[e];
				if (B_DEBUG == YES) { printf("\t\tVAR(%d) joining components %d and %d into a component of size %d\n",e,ii,jj,size); }
			}
		}
	}
 
	////////////////////////////////////////////////////////////////////////////

	// Tell CPLEX

	if (best_e == NONE) { 

		// Set useraction to indicate the default branch
		*useraction_p = CPX_CALLBACK_DEFAULT;
		//printf("\tNO BRANCHEABLE VARIABLE FOUND\n");
		//getchar();

		if (B_DEBUG == YES) { printf("\tNO BRANCHEABLE VARIABLE FOUND\n"); }

	} else {

		// Set useraction to indicate a user-specified branch
		*useraction_p = CPX_CALLBACK_SET;
		/* Up node */
		size = 0;
		ii = cc[index_i[best_e]];
		jj = cc[index_j[best_e]];
		// Tree:
		for (i=0; i<N-1; i++) {
			for (j=i+1; j<N; j++) {
				if ((cc[i]==ii && cc[j]==jj) || (cc[j]==ii && cc[i]==jj)){
					ee = index_e[i][j];
					if (ee == best_e) {
						sense[size] = 'L';
						index[size] = ee;
						bound[size] = 1.0;
						size++; 
					} else if (lb[ee]<0.5 && ub[ee]>0.5) {
						sense[size] = 'U';
						index[size] = ee;
						bound[size] = 0.0;
						size++; 
		}	}	}	}
		// Dist
		for (i=0; i<N-1; i++) {
			for (j=i+1; j<N; j++) {
				if ((cc[i]==ii && cc[j]==jj) || (cc[j]==ii && cc[i]==jj)){
					kk = index_k[i][j];
					if (kk != NONE && ub[kk]-lb[kk] > 2*EPSILON) {
						sense[size] = 'B';
						index[size] = kk;
						bound[size] = dist[i][j];
						size++; 
		}	}	}	}
		status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, size, index, sense, bound, objval, NULL, &seqnum1);
		assert(status == 0);
		
		/* Down node */ 
		sense[0] = 'U';
		index[0] = best_e;
		bound[0] = 0.0;
		status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, index, sense, bound, objval, NULL, &seqnum2); 
		assert(status == 0);
	}

	////////////////////////////////////////////////////////////////////////////

	// Clean:
	free_char_vector(&sense);
	free_double_vector(&bound);
	free_int_vector(&index); 

	free_double_vector(&lb);
	free_double_vector(&ub);
	free_int_vector(&feas); 
	free_int_vector(&prior); 

	free_int_vector(&sizes); 
	free_int_vector(&stack);
	free_int_vector(&degree);
	free_int_matrix(&neighbors, N);
	free_int_vector(&cc);
	
	if (B_DEBUG == YES) {
		printf("\nPress enter to continue\n");
		pause = 0;
		while (pause != '\r' && pause != '\n') { pause = getchar(); }
	}

	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////	
int add_local_cuts(GLOBAL_INFO global, CPXCENVptr env, void *cbdata, int wherefrom, double *lb, double *ub) {
    
    int     i,j,k,ii,jj,c,d,e;
    int		N = global.data->N;
	int		E = global.data->E;
    double  **C = global.data->C;
    int     **index_e = global.data->index_e;
    int     **index_k = global.data->index_k;
	double  **dist = global.var->dist;
	int		status  = 0;
    int		added_cuts = 0;    
	int		dist_cuts, tree_cuts;

	int		numnz;
	double  rhs;
	char	sense;
	int		*cutind = create_int_vector(1);
	double	*cutval = create_double_vector(1);
    
	int     size, n_cc, n_edges;
    int     *stack = create_int_vector(N);
	int     *cc = create_int_vector(N);
    int     *comp = create_int_vector(N);
	int		**neighbors = create_int_matrix(N,N);
	int		*degree = create_int_vector(N);
	
	int		**local_edges = create_int_matrix(N, N);
	int		*tree_analysis_stack = create_int_vector(N);
	int		*is_component_tree = create_int_vector(N);
	int		count_edges;

    
    /////////////////
    // PRE-PROCESS //
    /////////////////

	// Get linked list version of the graph:
	n_edges = 0;
	for (i=0; i<N; i++) { 
		for (j = 0; j < N; j++) {
			local_edges[i][j] = 0;
		}
		degree[i] = 0; 
	}

	for (i=0; i<N-1; i++) {
	    for (j=i+1; j<N; j++) {
			if( (lb[index_e[i][j]] + EPSILON >= 1.0) && (lb[index_e[j][i]] + EPSILON >= 1.0)) {
				printf("Lower Value of arc(%d,%d)=%.4lf\n", i, j, lb[index_e[i][j]]);
				neighbors[i][degree[i]++] = j;
				neighbors[j][degree[j]++] = i;
				local_edges[i][j] = 1;
				n_edges++;
            }
        } 
    }
	printf("Number of edges %d\n", n_edges);

    // Find Connected Components:
    for (i=0; i<N; i++) { cc[i] = NONE; }
    n_cc = 0;
    for (i=0; i<N; i++) {
		if (cc[i] == NONE) {
			cc[i] = n_cc;
			size = 0;
			stack[size++] = i;
			while (size > 0) {
				j = stack[--size];
				for (d=0; d<degree[j]; d++) {
					k = neighbors[j][d];
					if (cc[k] == NONE) {
						cc[k] = n_cc;
						stack[size++] = k;
					}
				}
			} 
			n_cc++;
		}
    }
	printf("Number of connected components %d\n", n_cc);
	for (i = 0; i < N; i++) { printf("node %d is in component %d\n", i, cc[i]); }

	//Trying to detect whether the components are all trees
	for (i = 0; i < n_cc; i++) {
		size = 0;
		// Collecting the nodes in connected component i
		for (j = 0; j < N; j++) {
			if (cc[j] == i) {
				tree_analysis_stack[size++] = j;
			}
		}
		// Count how many edges are in the connected component i
		count_edges = 0;
		for (d = 0; d < size; d++) {
			for (e = 0; e < size; e++) {
				count_edges += local_edges[tree_analysis_stack[d]][tree_analysis_stack[e]];
			}
		}
		// Mark as tree
		if (count_edges == size - 1) { is_component_tree[i] = 1; }
		else { is_component_tree[i] = 0; printf("component %d is not a tree\n"); }
	}
	

	if (n_cc == 1 && n_edges==N-1) {
		goto Terminate;
	}

	// Find Distances:
    for (i=0; i<N; i++) { for (j=0; j<N; j++) { dist[i][j] = -1.0; } }
	for (i=0; i<N; i++) {
        dist[i][i] = 0.0;
        size = 0;
        stack[size++] = i;
        while (size > 0) {
			j = stack[--size];
			for (d=0; d<degree[j]; d++) {
				k = neighbors[j][d];
                if (dist[i][k] < 0.0) {
					dist[i][k] = dist[i][j] + C[j][k];
                    stack[size++] = k;
				}
            }
        } 
    }

    //////////////
    // ADD CUTS //
    //////////////
    
    // Now generate cuts:   
	added_cuts	= 0;
	dist_cuts	= 0;
	tree_cuts	= 0;
    numnz       = 1;        // This is always the same
    sense       = 'E';      // This is always the same
    cutval[0]   = 1.0;      // This is always the same    
    for (c=0; c<n_cc; c++) {
		if (is_component_tree[c] == 1) {
			printf("Writing cut for component %d\n", c);
			size = 0;
			for (i = 0; i < N; i++) {
				if (cc[i] == c) {
					stack[size++] = i;
				}
			}

			for (ii = 0; ii < size - 1; ii++) {
				i = stack[ii];
				for (jj = ii + 1; jj < size; jj++) {
					j = stack[jj];
					k = index_k[i][j];
					e = index_e[i][j];

					//printf("distance commodity (%d,%d)=%.4lf\n", i, j, dist[i][j]);
					// Add row: Z^k = dist[i,j]
					if (k != NONE && ub[E + k] - lb[E + k] > EPSILON) {
						cutind[0] = E + k;
						rhs = dist[i][j];
						status = CPXcutcallbackaddlocal(env, cbdata, wherefrom, numnz, rhs, sense, cutind, cutval);
						assert(status == 0);
						added_cuts++;
						dist_cuts++;

					}

					// Add row: Y_e = 0
					if (lb[e] + EPSILON < 1.0 && ub[e] - EPSILON > 0.0) {
						cutind[0] = e;
						rhs = 0.0;
						status = CPXcutcallbackaddlocal(env, cbdata, wherefrom, numnz, rhs, sense, cutind, cutval);
						assert(status == 0);
						added_cuts++;
						tree_cuts++;
					}
				}
			}
		}
    }    

	if (global.param->SCREEN_OUTPUT >= 2) { printf("\t\t\t %d local_dist and %d local_tree cuts found\n", dist_cuts, tree_cuts); }

    Terminate:
	// Remember
	global.results->n_local_cuts += added_cuts;
	
	// Clean:
	free_int_vector(&cutind);
	free_double_vector(&cutval);
    free_int_vector(&cc);
    free_int_vector(&comp);
    free_int_vector(&stack);
	free_int_matrix(&neighbors, N);
	free_int_vector(&degree);

	//if (added_cuts>0) { printf("\t%d+%d Local cuts found this iteration\n", tree_cuts, dist_cuts); }

	return added_cuts;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int add_OPT_root_cuts(GLOBAL_INFO global, CPXCENVptr env, CPXLPptr lp, int INTEGER) {

	int		e, k;
	int		E = global.data->E;
	int		K = global.data->K;
	int		status  = 0;
	int		numnz;
	double  *rhs = create_double_vector(1);
	char	*sense = create_char_vector(1);
	int		*matbeg = create_int_vector(1);
	int		*matind = create_int_vector(E+1);
	double	*matval = create_double_vector(E+1);
	double  *violation = create_double_vector(K);
	double  threshold;
	int     n_rows;
	int		added_cuts = 0;

	// Look for the most_violated_cut:
	for (k=0; k<K; k++) { violation[k] = OPT_violation(global.data, global.var->subproblem[k], k, global.var->master); }

	// Set the threshold:
	if (global.param->ADD_NON_VIOLATED_CUTS  == YES) { threshold = list_min(K, violation) - 1.0; }
	else if (global.var->flag_fist_master == NO) { threshold = 0.0; }
	else if (INTEGER == YES) { threshold = 0.01; }
	else if (INTEGER == NO)  { threshold = 0.1; }
	
	// Now add all optimality cuts that are over the threshold:
	for (k=0; k<K; k++) { 	
		if (violation[k] > threshold) {				
			numnz     = 0;
			matbeg[0] = 0;
			sense[0]  = 'G';
			rhs[0]	  = global.var->subproblem[k][E+global.data->Com[k].d] - global.var->subproblem[k][E+global.data->Com[k].o];

			// Y_e
			for (e=0; e<E; e++) {
				if (global.var->subproblem[k][e] > EPSILON) {
					matval[numnz] = global.var->subproblem[k][e];
					matind[numnz] = e;
					numnz++;
				}
			}

			// Z_k
			matval[numnz] = 1.0;
			matind[numnz] = E+k;
			numnz++;

			// Finally add a row to the master problem...
			n_rows = CPXgetnumrows(env, lp);
			status = CPXaddrows(env, lp, 0, 1, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
			assert(status == 0);
			assert(n_rows + 1 == CPXgetnumrows(env, lp));
			added_cuts++;
		}			
	}

 	// Remember
	global.results->n_opt_cuts += added_cuts;
	if (INTEGER == YES) { global.results->n_int_opt_cuts  += added_cuts; }
	else                { global.results->n_frac_opt_cuts += added_cuts; }
	
	// Clean:
	free_int_vector(&matind);
	free_int_vector(&matbeg);
	free_char_vector(&sense);
	free_double_vector(&rhs);
	free_double_vector(&matval);
	free_double_vector(&violation);
	
	return added_cuts;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int add_OPT_callback_cuts(GLOBAL_INFO global, CPXCENVptr env, void *cbdata, int wherefrom, int INTEGER) {

	int		e, k;
	int		E = global.data->E;
	int		K = global.data->K;
	int		status  = 0;
	int		numnz;
	double  *rhs = create_double_vector(1);
	char	*sense = create_char_vector(1);
	int		*matbeg = create_int_vector(1);
	int		*matind = create_int_vector(E+1);
	double	*matval = create_double_vector(E+1);
	double  *violation = create_double_vector(K);
	double  threshold;
	int		added_cuts = 0;

	// Look for the most_violated_cut:
	for (k=0; k<K; k++) { violation[k] = OPT_violation(global.data, global.var->subproblem[k], k, global.var->master); }
	
	// Set the threshold:
	if (global.param->ADD_NON_VIOLATED_CUTS  == YES) { threshold = list_min(K, violation) - 1.0; }
	else if (global.var->flag_fist_master == NO) { threshold = 0.0; }
	else if (INTEGER == YES) { threshold = 0.01; }
	else if (INTEGER == NO)  { threshold = 0.1; }

	// Now add all optimality cuts that are over the threshold:
	for (k=0; k<K; k++) { 	
		if (violation[k] > threshold) {				
			numnz     = 0;
			matbeg[0] = 0;
			sense[0]  = 'G';
			rhs[0]	  = global.var->subproblem[k][E+global.data->Com[k].d] - global.var->subproblem[k][E+global.data->Com[k].o];
					
			// Y_e
			for (e=0; e<E; e++) {
				if (global.var->subproblem[k][e] > EPSILON) {
					matval[numnz] = global.var->subproblem[k][e];
					matind[numnz] = e;
					numnz++;
				}
			}

			// Z_k
			matval[numnz] = 1.0;
			matind[numnz] = E+k;
			numnz++;

			// Finally add a row to the callback node!
			// Possible CPX_USECUT codes:
			// CPX_USECUT_FORCE		The cut is added to the relaxation and stays there
			// CPX_USECUT_PURGE		The cut is added to the relaxation but can be purged later on if CPLEX deems the cut ineffective.
			// CPX_USECUT_FILTER	The cut is treated exactly as cuts generated by CPLEX; that is, CPLEX applies its filtering process and may not even add the cut to the relaxation, for example, if CPLEX deems other cuts more effective, or if the cut is too dense.
			status = CPXcutcallbackadd(env, cbdata, wherefrom, numnz, *rhs, *sense, matind, matval, CPX_USECUT_FORCE);
			added_cuts++;
		}			
	}

	// Remember
	global.results->n_opt_cuts += added_cuts;
	if (INTEGER == YES) { global.results->n_int_opt_cuts  += added_cuts; }
	else                { global.results->n_frac_opt_cuts += added_cuts; }
	
	// Clean:
	free_int_vector(&matind);
	free_int_vector(&matbeg);
	free_char_vector(&sense);
	free_double_vector(&rhs);
	free_double_vector(&matval);
	free_double_vector(&violation);

	return added_cuts;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int add_FEAS_cuts(GLOBAL_INFO global, CPXCENVptr env, CPXLPptr lp, int INTEGER) {

	int		c, e;
	int		*cut_set;
	int		*index_i = global.data->index_i;
	int		*index_j = global.data->index_j;
	int		N = global.data->N;
	int		E = global.data->E;
	int		status  = 0;
	int		numnz;
	double  *rhs = create_double_vector(1);
	char	*sense = create_char_vector(1);
	int		*matbeg = create_int_vector(1);
	int		*matind = create_int_vector(E);
	double	*matval = create_double_vector(E);
	double  *violation = create_double_vector(N);
	int		n_rows;
	int		added_cuts = 0;
	double  threshold;

	// Look for the most_violated_cut:
	if (global.param->FEAS_CUT_TYPE == 0) {
		for (c=0; c<N; c++) { violation[c] = CAPACITY_violation(global.data, global.var->cut_sets[c], global.var->master); }
	} else if (global.param->FEAS_CUT_TYPE == 1) {
		for (c=0; c<N; c++) { violation[c] = SUBTOUR_violation(global.data, global.var->cut_sets[c], global.var->master); }
	} else {
		printf("Error: Unknown Feasibility Cut Type %d\n", global.param->FEAS_CUT_TYPE);
	}
	
	// Set the threshold:
	if (global.param->ADD_NON_VIOLATED_CUTS  == YES) { threshold = list_min(N, violation) - 1.0; }
	else if (global.var->flag_fist_master == NO) { threshold = 0.0; }
	else if (INTEGER == YES) { threshold = 0.01; }
	else if (INTEGER == NO)  { threshold = 0.1; }

	// Now add some of the most violated cuts:
	if (list_max(N, violation) > threshold) {
		for (c=0; c<N; c++) { 
			if (violation[c] >= threshold) {
				
				cut_set = global.var->cut_sets[c];
				numnz  = 0;
				matbeg[0] = 0;
				
				if (global.param->FEAS_CUT_TYPE == 0) {		// Sum y_e for e in delta(S) >= 1
					sense[0]  = 'G';
					rhs[0]	  =  1.0;
					// Y_e
					for (e=0; e<E; e++) {
						if (cut_set[index_i[e]] != cut_set[index_j[e]]) {
							matval[numnz] = 1.0;
							matind[numnz] = e;
							numnz++;
						}
					}

				} else if (global.param->FEAS_CUT_TYPE == 1) { 		// Sum y_e for e in E(S) <= |S|-1
					sense[0]  = 'L';
					rhs[0]    =  -1.0;
					for (e=0; e<N; e++) { if (cut_set[e] == 1) { rhs[0] += 1.0; } }
					// Y_e
					for (e=0; e<E; e++) {
						if (cut_set[index_i[e]] == 1 && cut_set[index_j[e]] == 1) {
							matval[numnz] = 1.0;
							matind[numnz] = e;
							numnz++;
						}
					}
				} else {
					printf("ERROR: Unkown Feasibility Cut Type %d", global.param->FEAS_CUT_TYPE);
				}

				// Finally add a row to the master problem...
				n_rows = CPXgetnumrows(env, lp);
				status = CPXaddrows(env, lp, 0, 1, numnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
				assert(status == 0);
				assert(n_rows + 1 == CPXgetnumrows(env, lp));
				added_cuts++;
			}			
		}
	}
	
	// Remember
	global.results->n_feas_cuts += added_cuts;
	if (INTEGER == YES) { global.results->n_int_feas_cuts  += added_cuts; }
	else                { global.results->n_frac_feas_cuts += added_cuts; }

	// Clean:
	free_int_vector(&matind);
	free_int_vector(&matbeg);
	free_char_vector(&sense);
	free_double_vector(&rhs);
	free_double_vector(&matval);
	free_double_vector(&violation);

	return added_cuts;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int add_FEAS_callback_cuts(GLOBAL_INFO global, CPXCENVptr env, void *cbdata, int wherefrom, int INTEGER) {

	int		c, e;
	int		*cut_set;
	int		*index_i = global.data->index_i;
	int		*index_j = global.data->index_j;
	int		N = global.data->N;
	int		E = global.data->E;
	int		status  = 0;
	int		numnz;
	double  *rhs = create_double_vector(1);
	char	*sense = create_char_vector(1);
	int		*matbeg = create_int_vector(1);
	int		*matind = create_int_vector(E);
	double	*matval = create_double_vector(E);
	double  *violation = create_double_vector(N);
	int		added_cuts = 0;
	double  threshold;

	// Look for the most_violated_cut:
	if (global.param->FEAS_CUT_TYPE == 0) {
		for (c=0; c<N; c++) { violation[c] = CAPACITY_violation(global.data, global.var->cut_sets[c], global.var->master); }
	} else if (global.param->FEAS_CUT_TYPE == 1) {
		for (c=0; c<N; c++) { violation[c] = SUBTOUR_violation(global.data, global.var->cut_sets[c], global.var->master); }
	} else {
		printf("Error: Unknown Feasibility Cut Type %d\n", global.param->FEAS_CUT_TYPE);
	}

	// Set the threshold:
	if (global.param->ADD_NON_VIOLATED_CUTS  == YES) { threshold = list_min(N, violation) - 1.0; }
	else if (global.var->flag_fist_master == NO) { threshold = 0.0; }
	else if (INTEGER == YES) { threshold = 0.01; }
	else if (INTEGER == NO)  { threshold = 0.01; }

	// Now add some of the most violated cuts:
	if (list_max(N, violation) > threshold) {
		for (c=0; c<N; c++) { 
			if (violation[c] >= threshold) {
				
				cut_set = global.var->cut_sets[c];
				numnz  = 0;
				matbeg[0] = 0;
				
				if (global.param->FEAS_CUT_TYPE == 0) {		// Sum y_e for e in delta(S) >= 1
					sense[0]  = 'G';
					rhs[0]	  =  1.0;
					// Y_e
					for (e=0; e<E; e++) {
						if (cut_set[index_i[e]] != cut_set[index_j[e]]) {
							matval[numnz] = 1.0;
							matind[numnz] = e;
							numnz++;
						}
					}

				} else if (global.param->FEAS_CUT_TYPE == 1) { 		// Sum y_e for e in E(S) <= |S|-1
					sense[0]  = 'L';
					rhs[0]    =  -1.0;
					for (e=0; e<N; e++) { if (cut_set[e] == 1) { rhs[0] += 1.0; } }
					// Y_e
					for (e=0; e<E; e++) {
						if (cut_set[index_i[e]] == 1 && cut_set[index_j[e]] == 1) {
							matval[numnz] = 1.0;
							matind[numnz] = e;
							numnz++;
						}
					}
				}

				// Finally add a row to the callback node!
				// Possible CPX_USECUT codes:
				// CPX_USECUT_FORCE		The cut is added to the relaxation and stays there
				// CPX_USECUT_PURGE		The cut is added to the relaxation but can be purged later on if CPLEX deems the cut ineffective.
				// CPX_USECUT_FILTER	The cut is treated exactly as cuts generated by CPLEX; that is, CPLEX applies its filtering process and may not even add the cut to the relaxation, for example, if CPLEX deems other cuts more effective, or if the cut is too dense.
				status = CPXcutcallbackadd(env, cbdata, wherefrom, numnz, *rhs, *sense, matind, matval, CPX_USECUT_FORCE);
				assert(status == 0);
				added_cuts++;
			}			
		}
	}

	// Remember
	global.results->n_feas_cuts += added_cuts;
	if (INTEGER == YES) { global.results->n_int_feas_cuts  += added_cuts; }
	else                { global.results->n_frac_feas_cuts += added_cuts; }

	// Clean:
	free_int_vector(&matind);
	free_int_vector(&matbeg);
	free_char_vector(&sense);
	free_double_vector(&rhs);
	free_double_vector(&matval);
	free_double_vector(&violation);

	return added_cuts;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void print_OPT_CUT(INSTANCE *data, double *sub_var, int k) {
	int		e;
	int		E = data->E;
	printf("\tZ[%d] ", k);
	for (e=0; e<E; e++) { if (sub_var[e] > EPSILON) { printf("+ %.2f Y[%d] ", sub_var[e], e); } }
	printf(">= %.2f\n", sub_var[E + data->Com[k].d] - sub_var[E + data->Com[k].o]);
}
double OPT_violation(INSTANCE *data, double *sub_var, int k, double *master_var) {
	int		e;
	int		E = data->E;
	double	lhs = master_var[E + k];											// z_k +
	for (e=0; e<E; e++) { lhs += sub_var[e]*master_var[e]; }					// SUM gamma_e*y_e
	return (sub_var[E + data->Com[k].d] - sub_var[E + data->Com[k].o]) - lhs; 	// >= beta_d - beta_o 
}
double SUBTOUR_violation(INSTANCE *data, int *cut_set, double *master_variables) {
	int		i, e;
	int		N = data->N;
	int		E = data->E;
	int		size  = 0;
	double	lhs_0 = 0.0;
	double	lhs_1 = 0.0;
	double	rhs_0;
	double	rhs_1;

	// first check if there is a cut at all:
	for (i=0; i<N; i++) { size += cut_set[i]; }		// size == |S_1|	
	if (size == 0 || size == N) { return -1.0; }		// Meaning: not violated

	// if this is the case we must compute and return the violation (of either side):
	for (e=0; e<E; e++) { 
		if (cut_set[data->index_i[e]] == 0 && cut_set[data->index_j[e]] == 0) {
			lhs_0 += master_variables[e];
		} else if (cut_set[data->index_i[e]] == 1 && cut_set[data->index_j[e]] == 1) {
			lhs_1 += master_variables[e];
		}
	}
	rhs_1 = ((double)(  size)) - 1.0;	//   size == |S_1|	
	rhs_0 = ((double)(N-size)) - 1.0;	// N-size == |S| - |S_1| == |S_0|

	// We want the "1" side to be the most violated one:
	if (lhs_0-rhs_0 > lhs_1-rhs_1) {
		for (i=0; i<N; i++) { cut_set[i] = 1-cut_set[i]; }
	}

	// return the biggest violation:
	return MAX(lhs_0-rhs_0, lhs_1-rhs_1); // SUM y_e \forall e \in \E(S) <= |S| - 1
}
double CAPACITY_violation(INSTANCE *data, int *cut_set, double *master_variables) {
	int		i, e;
	int		N = data->N;
	int		E = data->E;
	int		size = 0;
	double	lhs  = 0.0;
	double	rhs  = 1.0;

	// firs check if there is a cut at all:
	for (i=0; i<N; i++) { size += cut_set[i]; }		// size == |S_1|	
	if (size == 0 || size == N) { return -1.0; }		// Meaning: non violated
	
	// if this is the case we must compute and return the violation:
	for (e=0; e<E; e++) { 
		if (cut_set[data->index_i[e]] != cut_set[data->index_j[e]]) {
			lhs += master_variables[e];
		}	
	}		
	
	// SUM y_e \forall e \in \delta(S) >= 1
	return rhs - lhs;
}
////////////////////////////////////////////////////////////////////////////////


