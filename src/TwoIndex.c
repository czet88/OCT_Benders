#include "headers.h"


/// INFO ///////////////////////////////////////////////////////////////////////
//
// TwoIndex.c: Two Index Monolitic Formulation
//
// Author:  Carlos Luna-Mota
// Created: 2015 Jun 26
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
int TwoIndex(CONFIGURATION param, INSTANCE data, RESULTS *results) {

	int i, j, k, a, e, numvar;
	int status = 0;
	int N = data.N;
	int E = data.E;
	int K = data.K;
	int A = E+E;
	int IJK = (N*(N-1)*(N-2)) / 6; // = Combinations of N elements taken in groups of 3
	double **C = data.C;
	double **W = data.W;
	double **Y = create_double_matrix(N,N);
	double **dist = create_double_matrix(N,N);
	int	**index_e = data.index_e;
	int	*index_i = data.index_i;
	int	*index_j = data.index_j;
	COMMODITY *Com = data.Com;
	double Big_M;
	clock_t	start;
	
	/// CPLEX MASTER POBLEM ////////////////////////////////////////////////////
	CPXENVptr env = NULL; // pointer to Cplex enviroment .......................
	CPXLPptr  lp  = NULL; // pointer to Cplex Linear Program data structure ....
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
	char	  *ctype;  // tipo de variable ('C' continua o 'B' binaria) ........ 
    ////////////////////////////////////////////////////////////////////////////
	
	
	
	// Create enviroment ///////////////////////////////////////////////////////
	env = CPXopenCPLEX(&status);     // pointer to Cplex Enviroment
	assert(env != NULL);
	if (param.SCREEN_OUTPUT >= 2) { 
		status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);					// Output B&B info
		assert(status == 0);
	}
	status = CPXsetintparam(env, CPX_PARAM_THREADS,    1);							// Threads usados
	assert(status == 0);
    status = CPXsetdblparam(env, CPX_PARAM_TILIM,		param.MAX_CPU_TIME);		// Time Limit
	assert(status == 0);
	status = CPXsetdblparam(env, CPX_PARAM_EPGAP,		1.0e-4);						// Gap de Epsilon Optimalidad
	assert(status == 0);
    status = CPXsetintparam(env, CPX_PARAM_MIPSEARCH,  CPX_MIPSEARCH_TRADITIONAL);	// Turn on traditional search for use with control callbacks 
	assert(status == 0);
	status = CPXsetdblparam(env, CPX_PARAM_CUTSFACTOR,  1);	// no cuts
	assert(status == 0);
	status = CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, 0);							//	0	CPX_MIPEMPHASIS_BALANCED		Balance optimality and feasibility; default
	assert(status == 0);															//	1	CPX_MIPEMPHASIS_FEASIBILITY		Emphasize feasibility over optimality
																					//	2	CPX_MIPEMPHASIS_OPTIMALITY		Emphasize optimality over feasibility
																					//	3	CPX_MIPEMPHASIS_BESTBOUND		Emphasize moving best bound
																					//	4	CPX_MIPEMPHASIS_HIDDENFEAS		Emphasize finding hidden feasible solutions	
	////////////////////////////////////////////////////////////////////////////



	/// DEFINE LP //////////////////////////////////////////////////////////////
	objsen  = CPX_MIN;			// Dirección de la optimización (minimizar!)
	numcols = A+A+E;			// Number of variables of the 2-index formulation (X's, P's & D's)
	numrows = N+A+A+IJK+IJK;	// Number of constrains of the 2-index formulation
	numnz   = ((N-1)*N)+(2*A)+(2*A)+(4*IJK)+(5*IJK);		// Number of Non-Zero coefs of the 2-index formulation
	
	if (param.SCREEN_OUTPUT >= 1) { 
		printf("\n#############################################################\n\n");
		printf("2-index Formulation applied to %s\n\n", data.filename); 
		printf("Generating %10d rows ((N^3 + 9N^2 - 4N) / 6)\n", numrows); 
		printf("           %10d columns (5*E)\n", numcols); 
		printf("           %10d non-zero coeffitiens\n\n", numnz); 
	}
	
	obj		= create_double_vector(numcols);
    rhs		= create_double_vector(numrows);
    sense	= create_char_vector(numrows);
    matbeg	= create_int_vector(numcols);
    matcnt	= create_int_vector(numcols); 
    matind	= create_int_vector(numnz); 
    matval	= create_double_vector(numnz);
    lb		= create_double_vector(numcols);
    ub		= create_double_vector(numcols);
	ctype	= create_char_vector(numcols);

	// Compute Big M ///////////////////
	Big_M = 0.0;
	for (i=0; i<N; i++) { dist[0][i] = list_max(N, C[i]); }	// We just borrow the first row of 'dist'
	Big_M = list_sum(N, dist[0]) - list_min(N, dist[0]);

	/// CONSTRAINTS ////////////////////
	// Sum X_ij = 1-delta(j,0)		(first N rows)
	for (j=0; j<N; j++) {
		sense[j] = 'E';					// Sense of constraints (Equal)
		if (j == 0) { rhs[j] = 0.0; }	// Right Hand Side (0)
		else        { rhs[j] = 1.0; }	// Right Hand Side (1)
	}	
	// X_ij - P_ij <= 0 				(next A rows)
	for (a=0; a<A; a++) {
		sense[N + a] = 'L';				// Sense of constraints (Less or Equal)
		rhs[N + a]   = 0.0;				// Right Hand Side      (0)
	}
	// P_ij + P_ji <= 1 				(next A rows)
	for (a=0; a<A; a++) {
		sense[N+A + a] = 'L';			// Sense of constraints (Less or Equal)
		rhs[N+A + a]   = 1.0;			// Right Hand Side      (1)
	}
	// P_ij + P_jk + X_kj - P_ik <= 1 	(next IJK rows)
	for (i=0; i<IJK; i++) {
		sense[N+A+A + i] = 'L';			// Sense of constraints (Less or Equal)
		rhs[N+A+A + i]   = 1.0;			// Right Hand Side      (1)
	}
	//  	(next IJK rows)
	for (i=0; i<IJK; i++) {
		sense[N+A+A+IJK + i] = 'L';		// Sense of constraints (Less or Equal)
		rhs[N+A+A+IJK + i]   = Big_M;	// Right Hand Side      (1)
	}
	
	// ^^^   DONE   ^^^
	// |||          ||| 
	// VVV Not DONE VVV

	/// VARIABLES //////////////////////
	numnz = 0;	// index of the current nonzero element
	for (e=0; e<E; e++) {			// Y_e variables...
        lb[e]     = 0.0;				// lower bound
        ub[e]     = 1.0;				// upper bound
		obj[e]    = 0.0;				// Objective Function Coef.
        matcnt[e] = K + 1;				// Number of constraints it appears
		matbeg[e] = numnz;				// first NZ coeff related to this variable
        for (k=0; k<K; k++) {
			matind[numnz] = N*K + k*E + e;	// Capacity constraint (row:  y_e^k - x_e^k - x_-e^k >= 0)
			matval[numnz] = 1.0;			// Coef
			numnz++;
        }
		matind[numnz] = numrows-1;		// Cardinality constraint (last row)
        matval[numnz] = 1.0;			// Coef
        numnz++;
	}
	assert(numnz == E*K+E);
	numvar = 0;
	for(k=0; k<K; k++) {			// X_ij^k variables
		for (i=0; i<N; i++) {
			for (j=0; j<N; j++) {
				if (i!=j) {
					lb[E+numvar]     = 0.0;				// Lower Bound
					ub[E+numvar]     = 1.0;				// Upper Bound
					obj[E+numvar]    = (Com[k].w)*(data.C[i][j]);		// Obj. Fun. coeff
					matcnt[E+numvar] = 3;				// Number of constraints in which appear
					matbeg[E+numvar] = numnz;			// first NZ coeff related to this variable
					if (i<j) {
						matind[numnz] = N*k + i;		// Network Constraint (row:  Sum_a X_ai - Sum_b X_ib = supply(i))
						matval[numnz] = -1.0;			// Coef
						numnz++;
						matind[numnz] = N*k + j;		// Network Constraint (row:  Sum_a X_aj - Sum_b X_jb = supply(j))
						matval[numnz] =  1.0;			// Coef
						numnz++;
					} else {
						matind[numnz] = N*k + j;		// Network Constraint (row:  Sum_a X_aj - Sum_b X_jb = supply(j))
						matval[numnz] =  1.0;			// Coef
						numnz++;
						matind[numnz] = N*k + i;		// Network Constraint (row:  Sum_a X_ai - Sum_b X_ib = supply(i))
						matval[numnz] = -1.0;			// Coef
						numnz++;
					}
					e = index_e[i][j];
					matind[numnz] = N*K + k*E + e;	// Capacity constraint (row:  y_e^k - x_e^k - x_-e^k >= 0)
					matval[numnz] = -1.0;			// Coef
					numnz++;
					numvar++;
				}
			}
		}				
	}
	assert(numvar == K*E*2);
	assert(numnz  == K*E*7 + E);
	
	// Load the linear program to Cplex enviroment /////////////////////////////
	lp = CPXcreateprob(env, &status, "TwoIndex");
	assert(lp != NULL);
    status = CPXcopylp(env, lp, numcols, numrows, objsen, obj, rhs, sense, matbeg, matcnt, matind, matval, lb, ub, NULL);
	assert(status == 0);

	if (param.TO_INTEGRALITY == YES) {
		// Now tell CPLEX what kind of variables do we have:
		for (a=0; a<A; a++) { ctype[a]     = 'B'; }		// X_a (Binary)
		for (a=0; a<A; a++) { ctype[A+a]   = 'C'; }		// P_a (Continous)
		for (e=0; e<E; e++) { ctype[A+A+e] = 'C'; }		// D_e (Continous)
		status = CPXcopyctype(env, lp, ctype);
		assert(status == 0);
	}
	////////////////////////////////////////////////////////////////////////////


	
	/// SOLVE //////////////////////////////////////////////////////////////////
	if (param.SCREEN_OUTPUT >= 1) { printf("\nSolving to integrality:\n\n"); }
	start = clock();
	status = CPXmipopt(env, lp);
	assert(status == 0);
	results->master_time = elapsed(start);
	results->total_time  = results->master_time;
	
	if (param.SCREEN_OUTPUT >= 1) { printf("\nObtaining Solution Info:\n\n"); }

	// Obtain a valid Lower Bound for the integer problem
	status = CPXgetbestobjval(env, lp, &(results->LowerBound));
	assert(status == 0);
	
	// Obtain a valid Upper Bound for the integer problem
	status = CPXgetcutoff(env, lp, &(results->UpperBound));
	assert(status == 0);
	
	// Compute the gap
	results->final_gap = (results->UpperBound-results->LowerBound) / results->UpperBound;

	// Get the number of explored nodes
	results->n_explored_nodes = CPXgetnodecnt(env, lp);

	// Get solution:							// This may not work if time runs out!!!
	status = CPXgetx(env, lp, lb, 0, A+A+E-1);
	assert(status == 0);

	// Copy solution to Y
	a = 0;
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { Y[j][i] = 0.0; } }
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { if (i!=j) {
		Y[i][j] += lb[a];
		Y[j][i] += lb[a];
		a++;
	} } }
	tree_dist(N, C, Y, dist);

	// Copy solution to best_solution:
	if (results->best_solution == NULL) { results->best_solution = create_double_vector(E+K); }
	for (e=0; e<E; e++) { results->best_solution[e]   = Y[index_i[e]][index_j[e]]; }
	for (k=0; k<K; k++) { results->best_solution[E+k] = dist[Com[k].o][Com[k].d];  }
	////////////////////////////////////////////////////////////////////////////

	
	// Clean:
	if (param.SCREEN_OUTPUT >= 1) { printf("\nFreeing memory:\n\n"); }
	free_double_matrix(&Y, N);
	free_double_matrix(&dist, N);	
	free_double_vector(&obj);
	free_double_vector(&rhs);
	free_char_vector(&sense);
	free_int_vector(&matbeg);
	free_int_vector(&matcnt); 
	free_int_vector(&matind); 
	free_double_vector(&matval);
	free_double_vector(&lb);
	free_double_vector(&ub);
	free_char_vector(&ctype);
	
	return status;
}
////////////////////////////////////////////////////////////////////////////////
