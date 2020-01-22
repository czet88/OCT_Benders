#include "headers.h"


/// INFO ///////////////////////////////////////////////////////////////////////
//
// ThreeIndex.c: Three Index Monolitic Formulation
//
// Author:  Carlos Luna-Mota
// Created: 2015 Jun 26
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
int ThreeIndex(CONFIGURATION param, INSTANCE data, RESULTS *results) {

	int i, j, k, e, r, numvar;
	int status = 0;
	int N = data.N;
	int E = data.E;
	int K = data.K;
	double **C = data.C;
	double **W = data.W;
	double **Y = create_double_matrix(N,N);
	double **dist = create_double_matrix(N,N);
	int	**index_e = data.index_e;
	int	*index_i = data.index_i;
	int	*index_j = data.index_j;
	COMMODITY *Com = data.Com;
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
	status = CPXsetintparam(env, CPX_PARAM_THREADS,		1);							// Threads usados
	assert(status == 0);
    status = CPXsetdblparam(env, CPX_PARAM_TILIM,		param.MAX_CPU_TIME);		// Time Limit
	assert(status == 0);
	status = CPXsetdblparam(env, CPX_PARAM_EPGAP,		1.0e-4);						// Gap de Epsilon Optimalidad
	assert(status == 0);
    status = CPXsetintparam(env, CPX_PARAM_MIPSEARCH,	CPX_MIPSEARCH_TRADITIONAL);	// Turn on traditional search for use with control callbacks 
	assert(status == 0);
	status = CPXsetdblparam(env, CPX_PARAM_CUTSFACTOR,  1);							// no cuts
	assert(status == 0);
	status = CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, 0);							//	0	CPX_MIPEMPHASIS_BALANCED		Balance optimality and feasibility; default
	assert(status == 0);															//	1	CPX_MIPEMPHASIS_FEASIBILITY		Emphasize feasibility over optimality
																					//	2	CPX_MIPEMPHASIS_OPTIMALITY		Emphasize optimality over feasibility
																					//	3	CPX_MIPEMPHASIS_BESTBOUND		Emphasize moving best bound
																					//	4	CPX_MIPEMPHASIS_HIDDENFEAS		Emphasize finding hidden feasible solutions	
	////////////////////////////////////////////////////////////////////////////



	/// DEFINE LP //////////////////////////////////////////////////////////////
	objsen  = CPX_MIN;			// Dirección de la optimización (minimizar!)
	numcols = E + K*E*2;		// Number of variables of the 4-index formulation (Y's & X's)
	numrows = (N+E)*K + 1;		// Number of constrains of the 4-index formulation
	numnz   = K*E*7 + E;		// Number of Non-Zero coefs of the 4-index formulation
	
	if (param.SCREEN_OUTPUT >= 1) { 
		printf("\n#############################################################\n\n");
		printf("4-index Formulation applied to %s\n\n", data.filename); 
		printf("Generating %10d rows (K(V+E) + 1)\n", numrows); 
		printf("           %10d columns (2KE + E)\n", numcols); 
		printf("           %10d non-zero coeffitiens (7EK + E)\n\n", numnz); 
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

	
	/// CONSTRAINTS ////////////////////
	// Network Constraints (first N*K rows)
	for (k=0; k<K; k++) {
		for (j=0; j<N; j++) {
			sense[k*N+j] = 'E';								// Sense of constraints (Equal)
			if      (j == Com[k].d) { rhs[k*N+j] =  1.0; }	// Right Hand Side      ( 1)
			else if (j == Com[k].o) { rhs[k*N+j] = -1.0;	 }	// Right Hand Side      (-1)
			else                    { rhs[k*N+j] =  0.0; }	// Right Hand Side      ( 0)
		}
	}	
	// Capacity Constraints (next E*K rows)
	for (r=N*K; r<numrows-1; r++) {
	    sense[r] = 'G';					// Sense of constraints (Greater or Equal)
        rhs[r]   = 0.0;					// Right Hand Side      (0)
	}
	// Cardinality constraint (last row)
	sense[numrows-1] = 'E';				// Sense of constraints (Equality)
    rhs[numrows-1]   = ((double)(N-1));	// Right Hand Side      (N-1)
	
	
	
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
	lp = CPXcreateprob(env, &status, "ThreeIndex");
	assert(lp != NULL);
    status = CPXcopylp(env, lp, numcols, numrows, objsen, obj, rhs, sense, matbeg, matcnt, matind, matval, lb, ub, NULL);
	assert(status == 0);

	if (param.TO_INTEGRALITY == YES) {
		// Now tell CPLEX what kind of variables do we have:
		for (e=0; e<E;     e++) { ctype[e] = 'B'; 	}		// Y_e (Binary)
		for (e=0; e<K*E*2; e++) { ctype[E+e] = 'C'; }		// X_a^k (Continous)
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

	// Remember solution:
	if (results->best_solution == NULL) { results->best_solution = create_double_vector(E+K); }
	status = CPXgetx(env, lp, results->best_solution, 0, E-1);
	assert(status == 0);
	for (i=0; i<N; i++) { Y[i][i] = 0.0; }
	for (e=0; e<E; e++) { Y[index_i[e]][index_j[e]] = Y[index_j[e]][index_i[e]] = results->best_solution[e]; }
	tree_dist(N, C, Y, dist);
	for (k=0; k<K; k++) { results->best_solution[E+k] = dist[Com[k].o][Com[k].d]; }
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
