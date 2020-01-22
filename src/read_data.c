#include "headers.h"

/// INFO ///////////////////////////////////////////////////////////////////////
//
// read_data.c: Input & Memory handling functions
//
// Author:  Carlos Luna-Mota
// Created: 2015 Jan 21
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// safe file open function
FILE *open_file(const char *filename, const char *mode) {
    FILE *file;
    if ((file = fopen(filename, mode)) == NULL) { 
		fprintf(stderr, "ERROR: Unable to open file.\n");
    }
    return file;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Read the parameters and stores them in "param". 
int read_CONFIGURATION(char *filename, CONFIGURATION *param) {
	int status = 0;

	// default PARAMETERS: ////////////////////////////////////////////////////
	
	// Select output configurations
	param->SCREEN_OUTPUT			= ALL;	// NO: No; YES: Minimal; ALL: Verbose
	param->FILE_OUTPUT				= ALL; 	// NO: No; YES: Minimal; ALL: Verbose
	
	// Select solving algorithm
	param->TO_INTEGRALITY			= YES;	// NO: Solve only the linear relaxation. YES: Solve the integer problem
	param->ALGORITHM				= 1;	// -1: None, apply only Heuristic procedures. 
											//  0: Iterative Benders (for the 4-index formulation)
											//  1: Branch&Cut Benders (for the 4-index formulation)
											//  2: Use the 2-index formulation
											//  3: Use the 3-index formulation 
											//  4: Use the 4-index formulation 

	// Select stopping criterions
	param->MIN_ABS_GAP				= 0.01;	// Algorithm stops when the absolute optimality gap is smaller than this.
	param->MAX_CPU_TIME				= 24*3600;	// In seconds
	param->UPPER_CUTOFF				= NO;	// NO: do not use Upper CutOff; YES: use Upper CutOff

	// Heuristic Phase:
	param->DO_MST_HEUR				= NO;
	param->DO_AMT_HEUR				= YES;
	param->DO_GHT_HEUR				= NO;
	param->DO_STAR_HEUR				= NO;	
	param->DO_MCMCT_HEUR			= NO;
	param->DO_BUDaC_HEUR			= NO;
	param->ADD_LOCAL_SEARCH_HEUR	= NO;
	
	// General Benders Settings:
	param->ROOT_CLEANUP				= 0.3;	// 0.0: Erase all non-binding root-node constraint
											// 0.X: Erase all non-binding root-node constraints with slack > (0.X*Numrows)'th smallest slack
											// 1.0: Keep all root-node constraints (binding or not)
	
	param->USE_ROUNDING_CALLBACK	= YES;	// YES: Use the CPLEX Heuristic Callback with a rounding heuristic for frac solutions
	param->APPLY_AM_LOCAL_SEARCH	= YES;	// YES:	Apply Ahuja-Murty Local search after rounding
	param->ADD_LOCAL_CUTS			= YES;	// YES: Add local cuts in the callbacks!
	param->USE_BRANCHING_CALLBACK   = YES;	// YES: Use the CPLEX Branching Callback
	param->ADD_NON_VIOLATED_CUTS	= NO;	// NO: Add only violated cuts. YES: Add them all.
	
	// Select initial formulation
	param->GLOBAL_LOWER_BOUND		= YES;	// Add a single "global lower bound" cut
	param->COMMODITY_LOWER_BOUND	= YES;	// Add a "commodity lower bound" for each commodity
	param->TREE_3_LOWER_BOUND		= NO;	// Add "3-tree lower bounds"
	param->TREE_4_LOWER_BOUND		= NO;	// Add "4-tree lower bounds"
	param->MST_OPT_CUTS				= YES;	// Add the optimality cuts associated to the MST
	param->AMT_OPT_CUTS				= YES;	// Add the optimality cuts associated to the Ahuja-Murty Tree
	param->GHT_OPT_CUTS				= YES;	// Add the optimality cuts associated to the Gomory-Hu Tree
	param->MCMCT_OPT_CUTS			= NO;	// Add the optimality cuts associated to the MinCutMinCost Tree
	param->STARS_OPT_CUTS			= NO;	// Add the optimality cuts associated to each Star-tree

	// Feasibility Cuts:
	param->FEAS_CUT_TYPE			= 1;	// 0: Cutset Inequalities, 1: Packing Inequalities
	 
	// Optimality Cuts:
	param->OPT_CUT_TYPE				= 1;	// 0: Solve Dual Subproblem using CPLEX.
											// 1: Solve Magnanti-Wong Pareto Dual Subproblem using CPLEX 
											// 2: Solve Papadakos     Pareto Dual Subproblem using CPLEX 
	param->STAR_LIKE_CORE_POINT		= YES;	// YES: Reset core point before each update!
											// NO:	Update from previous core point, never reset it!

	// Integer Cuts:
	param->HEUR_LIFTING				= NO;	// YES: Try to lift all optimality cuts.
	
	// Fractional Cuts:
	param->FRAC_CUTS_FREQ			= 5;	// User cut callback is called only in levels at depths that are multiples of this number.
	param->ADD_FRAC_FEAS_CUTS		= NO;   // YES: Add fractional cuts.
	param->ADD_FRAC_OPT_CUTS		= YES;  // YES: Add fractional cuts.
	////////////////////////////////////////////////////////////////////////////

	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Read and pre-process and instance from a file.
int read_INSTANCE(INSTANCE *instance_ptr) {
	
	int i,j,k,e;
	int N, E, K;
	int status = 0;
	long long int aux;
    FILE *input;
	
	input = open_file(instance_ptr->filename, "r");

    // Read the instance size (=N) /////////////////////////////////////////////
	fscanf(input, "%d ", &N);
    instance_ptr->N = N;

	// Make enought room to just read the rest of the instance /////////////////
    instance_ptr->C  = create_double_matrix(N, N);
	instance_ptr->W  = create_double_matrix(N, N);
	
	// Read the distace matrix (=C) ////////////////////////////////////////////
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (fscanf(input, "%lld", &aux) != YES) {
		        fprintf(stderr, "ERROR: Unable to read dist matrix (%d,%d)\n",i,j);
	        } else {
				instance_ptr->C[i][j] = ((double)(aux));
			}

		}
	}

    // Read the communication requests matrix (=W) /////////////////////////////
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (fscanf(input, "%lld", &aux) != YES) {
	        fprintf(stderr, "ERROR: Unable to read request matrix (%d,%d)\n",i,j);
	        } else {
				instance_ptr->W[i][j] = ((double)(aux));
			}
		}
	}

	// We don't need more information:
	fclose(input);

	// PRE PROCESS INSTANCE: ///////////////////////////////////////////////////

    // Assert that the distace matrix C symetric:
    for (i=0; i<N; i++) {
        assert(instance_ptr->C[i][i] == 0.0);									// There are 0's in the main diagonal				
		for (j=i+1; j<N; j++) {
            assert(instance_ptr->C[j][i] == instance_ptr->C[i][j]);         	// It's symetric
        }
	}

	// Make the communication requests matrix W upper-triangular:
    for (i=0; i<N; i++) {
        instance_ptr->W[i][i] = 0.0;					// Put 0's in the main diagonal
		for (j=i+1; j<N; j++) {		
			instance_ptr->W[i][j] += instance_ptr->W[j][i];			// Make upper-triangular
			instance_ptr->W[j][i] = 0.0;
		}
	}

	// Compute the characteristic sizes of this instance:
	E =  (N*(N-1))/2;		// Number of edges
	instance_ptr->E = E;

	K = 0;	                // Number of NonZero comunication requests
    for (i=0; i<N; i++) { 
		for (j=0; j<N; j++) {
			if (instance_ptr->W[i][j] > 0) { K++; }
		} 
	}
	instance_ptr->K = K;

	instance_ptr->index_i = create_int_vector(E);
	instance_ptr->index_j = create_int_vector(E);
	instance_ptr->index_e = create_int_matrix(N, N);

	for (i=0; i<N; i++) { instance_ptr->index_e[i][i] = NONE; }
	e=0;
	for (i=0; i<N-1; i++) {
		for (j=i+1; j<N; j++) {
			instance_ptr->index_i[e] = i;
			instance_ptr->index_j[e] = j;
			instance_ptr->index_e[i][j] = e;
			instance_ptr->index_e[j][i] = e;
			e++;
		}
	}
	assert(e == E);

    // Compute the C1 matrix:
    instance_ptr->C1 = create_double_matrix(N, N);
    for (i=0; i<N; i++) { for (j=0; j<N; j++) {
        instance_ptr->C1[i][j] = instance_ptr->C[i][j];
    } }
	Floyd_Warshall(N, instance_ptr->C1);

	// Compute the C2 matrix:
	instance_ptr->C2 = create_double_matrix(N, N);
	for (i=0; i<N; i++) {
        instance_ptr->C2[i][i] = 0.0;									// Put 0's in the main diagonal				
		for (j=i+1; j<N; j++) {
			instance_ptr->C2[i][j] = -1.0;
			for (k=0; k<N; k++) {
				if(k!=i && k!=j) {
					if (instance_ptr->C2[i][j] == -1.0) {
						instance_ptr->C2[i][j] = instance_ptr->C1[i][k] + instance_ptr->C1[k][j];
					} else {
						instance_ptr->C2[i][j] = MIN(instance_ptr->C2[i][j], instance_ptr->C1[i][k]+instance_ptr->C1[k][j]);
					}
				}
			}
			instance_ptr->C2[j][i] = instance_ptr->C2[i][j];		// Make it symmetric
		}
	}

    // Determine if the instance is euclidean:
    instance_ptr->IS_EUCLIDEAN = YES;
    for (i=0; i<N; i++) { for (j=0; j<N; j++) {
        if (instance_ptr->C1[i][j] != instance_ptr->C[i][j]) {
            instance_ptr->IS_EUCLIDEAN = NO;        
        }
    } }

	// Compute MIN_CUT and MIN_CUT_MIN_COST matrices:
	instance_ptr->MIN_CUT = create_double_matrix(N, N);
	instance_ptr->MIN_CUT_MIN_COST = create_double_matrix(N, N);
	compute_min_cut_data(N, instance_ptr->C, instance_ptr->W, instance_ptr->MIN_CUT, instance_ptr->MIN_CUT_MIN_COST);
	
	// Now store all Comunication request by decreasing order of MIN_CUT[i][j]*MIN_CUT_MIN_COST[i][j]:
	// I will use here an O(NN^4) selection sort for simplicity (this will be done just once!)
	instance_ptr->Com = create_COMMODITY_vector(K);
	k = 0;
	while (k<K) {
		instance_ptr->Com[k].w = -1.0;
		instance_ptr->Com[k].o = NONE;
		instance_ptr->Com[k].d = NONE;
		for (i=0; i<N-1; i++) {
		    for (j=i+1; j<N; j++) {
			    if (instance_ptr->W[i][j] > 0.0) {
					if (instance_ptr->Com[k].w == -1.0) {
						instance_ptr->Com[k].w = instance_ptr->W[i][j];	// Value
						instance_ptr->Com[k].o = i;					// Origin
						instance_ptr->Com[k].d = j;					// Destination
						instance_ptr->W[i][j] = -1.0;
					// ALTERNATIVE: Use just the MIN_CUT data.
					// } else if (instance_ptr->MIN_CUT[i][j] > instance_ptr->MIN_CUT[instance_ptr->Com[k].o][instance_ptr->Com[k].d]) {
					} else if (instance_ptr->MIN_CUT_MIN_COST[i][j] * instance_ptr->MIN_CUT[i][j] > instance_ptr->MIN_CUT_MIN_COST[instance_ptr->Com[k].o][instance_ptr->Com[k].d] * instance_ptr->MIN_CUT[instance_ptr->Com[k].o][instance_ptr->Com[k].d]) {
						instance_ptr->W[instance_ptr->Com[k].o][instance_ptr->Com[k].d] = instance_ptr->Com[k].w;	// Restore previous
						instance_ptr->Com[k].w = instance_ptr->W[i][j];	// Value
						instance_ptr->Com[k].o = i;					// Origin
						instance_ptr->Com[k].d = j;					// Destination
						instance_ptr->W[i][j] = NONE;
					}
				}
			}
		}
		k++;
	}
	assert(k==K);

	// Now restore the W matrix...
	for (k=0; k<K; k++) { instance_ptr->W[instance_ptr->Com[k].o][instance_ptr->Com[k].d] = instance_ptr->Com[k].w; }
	
	// ...and build the index_k:
	instance_ptr->index_k = create_int_matrix(N, N);
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { instance_ptr->index_k[i][j] = NONE; } }
	for (k=0; k<K; k++) { 
		instance_ptr->index_k[instance_ptr->Com[k].o][instance_ptr->Com[k].d] = k;
		instance_ptr->index_k[instance_ptr->Com[k].d][instance_ptr->Com[k].o] = k;
	}
	////////////////////////////////////////////////////////////////////////////
	
	return status;
}
////////////////////////////////////////////////////////////////////////////////



// VECTORS /////////////////////////////////////////////////////////////////////
char *create_char_vector(int cells) {
    char *ptr = (char *) calloc(cells, sizeof(char));
    if (ptr == NULL) {
        fprintf(stderr, "ERROR: Unable to allocate memory for char_vector\n");
    }
    return ptr;
}
void free_char_vector(char **ptr) {
	if (*ptr == NULL) {
		fprintf(stderr, "ERROR: Unable to free memory from char_vector\n");
	} else {
		free(*ptr);
		*ptr = NULL;
	}
}

double *create_double_vector(int cells) {
    double *ptr = (double *) calloc(cells, sizeof(double));
    if (ptr == NULL) {
        fprintf(stderr, "ERROR: Unable to allocate memory for double_vector\n");
    }
    return ptr;
}
void free_double_vector(double **ptr) {
	if (*ptr == NULL) {
		fprintf(stderr, "ERROR: Unable to free memory from double_vector\n");
	} else {
		free(*ptr);
		*ptr = NULL;
	}
}

int *create_int_vector(int cells) {
    int *ptr = (int *) calloc(cells, sizeof(int));
    if (ptr == NULL) {
        fprintf(stderr, "ERROR: Unable to allocate memory for int_vector\n");
    }
    return ptr;
}
void free_int_vector(int **ptr) {
	if (*ptr == NULL) {
		fprintf(stderr, "ERROR: Unable to free memory from int_vector\n");
	} else {
		free(*ptr);
		*ptr = NULL;
	}
}

INSTANCE *create_INSTANCE_vector(int *cells_ptr, const char *filename) {
	int		i;
	FILE	*instance_list;      
	INSTANCE *ptr;
	const int MAX_FILEPATH_SIZE = 256;


	// Open file in read mode:
	instance_list = open_file(filename, "r");

	// Read the number of instances to solve.
	fscanf(instance_list, "%d", cells_ptr);	
	
	// Allocate vector of INSTANCES and store filenames:
	ptr = (INSTANCE *) calloc((*cells_ptr), sizeof(INSTANCE));
	if (ptr == NULL) {
        fprintf(stderr, "ERROR: Unable to allocate memory for INSTANCE_vector\n");
	} else {		
		for (i=0; i<(*cells_ptr); i++) { 			
			ptr[i].filename = create_char_vector(MAX_FILEPATH_SIZE);
			fscanf(instance_list, "%s", ptr[i].filename);
		}
    }
    
	// Close the input file
	fclose(instance_list);

	return ptr;
}
void free_INSTANCE(INSTANCE *instance_ptr) {
	free(instance_ptr->filename);
	free_double_matrix(&(instance_ptr->C), instance_ptr->N);
	free_double_matrix(&(instance_ptr->C1), instance_ptr->N);
    free_double_matrix(&(instance_ptr->C2), instance_ptr->N);
	free_double_matrix(&(instance_ptr->W), instance_ptr->N);
	free_double_matrix(&(instance_ptr->MIN_CUT), instance_ptr->N);
	free_double_matrix(&(instance_ptr->MIN_CUT_MIN_COST), instance_ptr->N);
	free_COMMODITY_vector(&(instance_ptr->Com));
	free_int_vector(&(instance_ptr->index_i));
	free_int_vector(&(instance_ptr->index_j));
	free_int_matrix(&(instance_ptr->index_e), instance_ptr->N);
	free_int_matrix(&(instance_ptr->index_k), instance_ptr->N);
}
void free_INSTANCE_vector(INSTANCE **ptr) {
	if (*ptr == NULL) {
		fprintf(stderr, "ERROR: Unable to free memory from INSTANCE_vector\n");
	} else {
		free(*ptr);
		*ptr = NULL;
	}
}

COMMODITY *create_COMMODITY_vector(int cells) {
	COMMODITY *ptr = (COMMODITY *) calloc(cells, sizeof(COMMODITY));
    if (ptr == NULL) {
        fprintf(stderr, "ERROR: Unable to allocate memory for COMMODITY_vector\n");
    }
    return ptr;
}
void free_COMMODITY_vector(COMMODITY **ptr) {
	if (*ptr == NULL) {
		fprintf(stderr, "ERROR: Unable to free memory from COMMODITY_vector\n");
	} else {
		free(*ptr);
		*ptr = NULL;
	}
}

RESULTS *create_RESULTS_vector(int cells) {
	int i;
	RESULTS *ptr = (RESULTS *) calloc(cells, sizeof(RESULTS));
    if (ptr == NULL) {
        fprintf(stderr, "ERROR: Unable to allocate memory for RESULTS_vector\n");
    } else {
		for (i=0; i<cells; i++) {

			ptr[i].total_time = 0.0;
			ptr[i].root_time = 0.0;

			ptr[i].master_time = 0.0;
			ptr[i].sub_time = 0.0;

			ptr[i].feas_time = 0.0;
			ptr[i].opt_time = 0.0;

			ptr[i].initial_cons = 0;
			ptr[i].root_cons = 0;
			ptr[i].final_cons = 0;

			ptr[i].LowerBound = 0.0;
			ptr[i].UpperBound = 0.0;

			ptr[i].initial_gap = 0.0;
			ptr[i].root_gap    = 0.0;
			ptr[i].final_gap   = 0.0;

			ptr[i].n_opt_cuts = 0;
			ptr[i].n_feas_cuts = 0;
			ptr[i].n_local_cuts = 0;
			
			ptr[i].n_int_opt_cuts = 0;
			ptr[i].n_int_feas_cuts = 0;
			ptr[i].n_frac_opt_cuts = 0;
			ptr[i].n_frac_feas_cuts = 0;
			
			ptr[i].n_opt_subproblems = 0;
			ptr[i].n_feas_subproblems = 0;
			
			ptr[i].n_int_opt_subproblems = 0;
			ptr[i].n_int_feas_subproblems = 0;
			ptr[i].n_frac_opt_subproblems = 0;
			ptr[i].n_frac_feas_subproblems = 0;

			ptr[i].best_solution = NULL;
			ptr[i].n_explored_nodes = 0;
		}
	}
    return ptr;
}
void free_RESULTS_vector(RESULTS **ptr, int cells) {
	int i;
	if (*ptr == NULL) {
		fprintf(stderr, "ERROR: Unable to free memory from RESULTS_vector\n");
	} else {
		for (i=0; i<cells; i++) {
			if ((*ptr)[i].best_solution == NULL) {
				fprintf(stderr, "ERROR: Unable to free memory from RESULTS.double_vector\n");
			} else {
				free_double_vector(&((*ptr)[i].best_solution));	
			}
		}
		free(*ptr);
		*ptr = NULL;
	}
}
////////////////////////////////////////////////////////////////////////////////


// MATRICES ////////////////////////////////////////////////////////////////////
double **create_double_matrix(int rows, int columns) {
    int i;
    double **ptr = (double **) calloc(rows, sizeof(double *));
    if (ptr == NULL) {
        fprintf(stderr, "ERROR: Unable to allocate memory for double_matrix\n");
    } else {
		for (i=0; i<rows; i++) { ptr[i] = create_double_vector(columns); }
	}
    return ptr;
}
void free_double_matrix(double ***ptr, int rows) {
	int i;
	if (*ptr == NULL) {
		fprintf(stderr, "ERROR: Unable to free memory from double_matrix\n");
	} else {
		for (i=0; i<rows; i++) { free_double_vector( &(*ptr)[i] ); }
		*ptr = NULL;
	}
}

int **create_int_matrix(int rows, int columns) {
    int i;
    int **ptr = (int **) calloc(rows, sizeof(int *));
    if (ptr == NULL) {
        fprintf(stderr, "ERROR: Unable to allocate memory for int_matrix\n");
    } else {
		for (i=0; i<rows; i++) { ptr[i] = create_int_vector(columns); }
	}
    return ptr;
}
void free_int_matrix(int ***ptr, int rows) {
	int i;
	if (*ptr == NULL) {
		fprintf(stderr, "ERROR: Unable to free memory from int_matrix\n");
	} else {
		for (i=0; i<rows; i++) { free_int_vector( &(*ptr)[i] ); }
		*ptr = NULL;
	}
}
////////////////////////////////////////////////////////////////////////////////
