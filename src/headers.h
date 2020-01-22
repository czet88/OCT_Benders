#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include "cplex.h"



/// INFO ///////////////////////////////////////////////////////////////////////
//
// headers.h: Call for external libraries, pre-processor directives and
//            function declarations for everything.
//
// Author:  Carlos Luna-Mota
// Created: 2015 Jan 21
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////



/// C PREPROCESSOR CODE ////////////////////////////////////////////////////////

// Boolean values simulated with signed integers
#define ALL		2
#define YES     1	// True
#define NO      0	// False
#define NONE   -1	// Unknown

// General precision for dealing with non integer comparisons
#define EPSILON 0.0001

// Some basic functions:
#define ABS(x) (((x) > 0 ) ? (x) : -(x))	
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) > (b)) ? (b) : (a))
#define SIZE(x) (sizeof (x) / sizeof (*(x)))
////////////////////////////////////////////////////////////////////////////////



/// Typedefs & structs /////////////////////////////////////////////////////////

// Generic commodity data structure
typedef struct COMMODITY { 
    int		o;				// Origin
    int		d;				// Destination
    double	w;				// Communication Request
} COMMODITY;

// OCSTP Instance data structure
typedef struct INSTANCE {

	// Raw Data:
	char	*filename;	// Instance (file)name
	int		N;			// Number of Nodes
	double	**C;		// Edge Cost (symmetric)
	double	**W;		// Communication Requirement (upper triangular)
	
	// Pre-processed data:
	int			E;		// Number of Edges
	int			K;		// Number of Commodities
	COMMODITY	*Com;	// List of Commodities (for the sake of sparsenest)

	int		*index_i;				// Y[e] = Y[index_i[e]][index_j[e]]
	int		*index_j;				// Y[e] = Y[index_i[e]][index_j[e]]
	int		**index_e;				// e = index_e[i][j]
	int		**index_k;				// k = index_k[i][j] <=> Com[k].o==i && Com[k].d==j

    int     IS_EUCLIDEAN;           // Tells if C is euclidean.
	double  **C1;					// Shortest path distance matrix
	double  **C2;					// Second shortest path distance matrix

	double	**MIN_CUT;				// Value of the (i,j)-MinCut
	double	**MIN_CUT_MIN_COST;		// Minimum Cost of any edge crossing the (i,j)-MinCut

} INSTANCE;

// OCSTP Experimental settings data structure
typedef struct CONFIGURATION {
	
	// Select output configurations
	int		SCREEN_OUTPUT;			// NO: No; YES: Minimal; ALL: Verbose
	int		FILE_OUTPUT;			// NO: No; YEs: Minimal; ALL: Verbose
	
	// Select solving algorithm
	int		TO_INTEGRALITY;			// NO: Solve only the linear relaxation. YES: Solve the integer problem
	int		ALGORITHM;				// -1: None, apply only Heuristic procedures. 
									//  0: Iterative Benders (for the 4-index formulation)
									//  1: Branch&Cut Benders (for the 4-index formulation)
									//  2: Use the 2-index formulation
									//  3: Use the 3-index formulation 
									//  4: Use the 4-index formulation 
		
	// Select stopping criterions
	double	MIN_ABS_GAP;			// Algorithm stops when optimality gap is smaller than this
	double	MAX_CPU_TIME;			// In seconds
	int		UPPER_CUTOFF;			// NO: do not use Upper CutOff; YES: use Upper CutOff

	// Select Heuristics for initial UpperBound
	int		DO_MST_HEUR;
	int		DO_AMT_HEUR;
	int		DO_GHT_HEUR;
	int		DO_STAR_HEUR;
	int     DO_MCMCT_HEUR;
	int		DO_BUDaC_HEUR;
	int		ADD_LOCAL_SEARCH_HEUR;

	// General Benders Settings:
	double  ROOT_CLEANUP;			// 0.0  Erase all non-binding root-node constraint
									// 0.X  Erase all non-binding root-node constraints with slack > (0.X*Numrows)'th smallest slack
									// 1.0  Keep all root-node constraints (binding or not)

	int		USE_ROUNDING_CALLBACK;	// YES: Use the CPLEX Heuristic Callback with a rounding heuristic
	int		APPLY_AM_LOCAL_SEARCH;	// YES:	Apply Ahuja-Murty Local search after rounding.
	int		ADD_LOCAL_CUTS;			// YES: Add local cuts in the callbacks!
	int		USE_BRANCHING_CALLBACK; // YES: Use the CPLEX Branching Callback
	int		ADD_NON_VIOLATED_CUTS;	// NO: Add only violated cuts. YES: Add them all.

	// Select initial formulation
	int		GLOBAL_LOWER_BOUND;		// Add a single "global lower bound" cut
	int		COMMODITY_LOWER_BOUND;	// Add a "commodity lower bound" for each commodity
	int		TREE_3_LOWER_BOUND;		// Add "3-tree lower bounds"
	int		TREE_4_LOWER_BOUND;		// Add "4-tree lower bounds"
	int		MST_OPT_CUTS;			// Add the optimality cuts associated to the MST
	int		AMT_OPT_CUTS;			// Add the optimality cuts associated to the Ahuja-Murty Tree
	int		GHT_OPT_CUTS;			// Add the optimality cuts associated to the Gomory-Hu Tree
	int		MCMCT_OPT_CUTS;			// Add the optimality cuts associated to the MinCutMinCost Tree
	int		STARS_OPT_CUTS;			// Add the optimality cuts associated to each Star-tree

	// Feasibility Cuts:
	int		FEAS_CUT_TYPE;			// 0: Cutset Inequalities, 1: Subtour Elimination Inequalities
	
	// Optimality Cuts:
	int		OPT_CUT_TYPE;			// 0: Solve Dual Subproblem using CPLEX.
									// 1: Solve Magnanti-Wong Pareto Dual Subproblem using CPLEX.
									// 2: Solve Papadakos Pareto Dual Subproblem using CPLEX.
	int		STAR_LIKE_CORE_POINT;	// YES: Star-Like update strategy. NO: Path-like update strategy.

	// Integer Cuts:
	int		HEUR_LIFTING;			// YES: Try to lift all integer optimality cuts.
	
	// Fractional Cuts:	
	int		FRAC_CUTS_FREQ;			// From then on, cuts are included every FRAC_CUTS_FREQ levels of depth.
	int		ADD_FRAC_FEAS_CUTS;		// YES: Add fractional cuts.
	int		ADD_FRAC_OPT_CUTS;		// YES: Add fractional cuts.
	
} CONFIGURATION;

// Generic Benders Decomposition results
typedef struct RESULTS {
		
	double	heuristic_times[30];
	double  heuristic_values[30];

	double	total_time;		// total cputime (in seconds)
	double  root_time;		// cputime (in seconds) spent solving the root node

	double  master_time;	// cputime (in seconds) spent solving the master problem
	double	sub_time;		// cputime (in seconds) spent solving the subproblem problems

	double  feas_time;		// cputime (in seconds) spent solving feasibility subproblems
	double	opt_time;		// cputime (in seconds) spent solving optimality subproblems

	int		initial_cons;
	int		root_cons;
	int		final_cons;

	double	UpperBound;
	double	LowerBound;

	double  initial_gap;
	double  root_gap;
	double  final_gap;

	int		n_opt_cuts;			// Total
	int		n_feas_cuts;		// Total
	int		n_local_cuts;

	int		n_int_opt_cuts;
	int		n_int_feas_cuts;
	int		n_frac_opt_cuts;
	int		n_frac_feas_cuts;
	
	int		n_opt_subproblems;	// Total
	int		n_feas_subproblems;	// Total
	
	int		n_int_opt_subproblems;
	int		n_int_feas_subproblems;
	int		n_frac_opt_subproblems;
	int		n_frac_feas_subproblems;
	
	double  *best_solution;

	int		n_explored_nodes;

} RESULTS;

// Variables for Benders Algorithm
typedef struct BENDERS_VARIABLES {
	double		*master;		// Master Problem variables
	double		**subproblem;	// Subproblem variables
	int			**cut_sets;		// Cut-Sets
	double		**capacity;		// Capacity auxiliary matrix
	double		**dist;			// Distances auxiliary matrix

	clock_t		init_time;		// Initial time

	int			last_id;		// Id of last branch&bound node visited
	int			last_id_visits;	// # of times that we have visited last B&B node
	int			flag_initial_sol;
	int			flag_generate_cuts;
	int			flag_fist_master;

	double		*core_point;	// Core Point 
	double		lambda;			// Core Point Update Factor
} BENDERS_VARIABLES;

// Variables Globales
typedef struct GLOBAL_INFO {
	// General Global Variables:
	CONFIGURATION		*param;
	INSTANCE			*data;
	RESULTS				*results;
	BENDERS_VARIABLES	*var;
} GLOBAL_INFO, *GLOBAL_INFO_ptr;
////////////////////////////////////////////////////////////////////////////////



/// read_data.c ////////////////////////////////////////////////////////////////
FILE	*open_file(const char *filename, const char *mode);
int		read_CONFIGURATION(char *filename, CONFIGURATION *param);
int		read_INSTANCE(INSTANCE *instance_ptr);
char	*create_char_vector(int cells);
void	free_char_vector(char **ptr);
double	*create_double_vector(int cells);
void	free_double_vector(double **ptr);
int		*create_int_vector(int cells);
void	free_int_vector(int **ptr);
INSTANCE *create_INSTANCE_vector(int *cells, const char *filename);
void	free_INSTANCE(INSTANCE *instance);
void	free_INSTANCE_vector(INSTANCE **ptr);
COMMODITY *create_COMMODITY_vector(int cells);
void	free_COMMODITY_vector(COMMODITY **ptr);
RESULTS *create_RESULTS_vector(int cells);
void	free_RESULTS_vector(RESULTS **ptr, int cells);
double	**create_double_matrix(int rows, int columns);
void	free_double_matrix(double ***ptr, int rows);
int		**create_int_matrix(int rows, int columns);
void	free_int_matrix(int ***ptr, int rows);
////////////////////////////////////////////////////////////////////////////////



// output.c ////////////////////////////////////////////////////////////////////
void show_help(char *program_name);
void output_RESULTS(char *filename, CONFIGURATION param, INSTANCE instance, RESULTS results);
void output_averaged_RESULTS(char *filename, CONFIGURATION param, int n_of_instances, RESULTS *results);
void output_settings(char *filename, CONFIGURATION param, int n_of_instances);
////////////////////////////////////////////////////////////////////////////////



/// Heuristics.c ///////////////////////////////////////////////////////////////
int		Heuristics(CONFIGURATION param, INSTANCE data, RESULTS *results);
////////////////////////////////////////////////////////////////////////////////



/// TwoIndex.c /////////////////////////////////////////////////////////////////
int		TwoIndex(CONFIGURATION param, INSTANCE data, RESULTS *results);
////////////////////////////////////////////////////////////////////////////////



/// ThreeIndex.c ///////////////////////////////////////////////////////////////
int		ThreeIndex(CONFIGURATION param, INSTANCE data, RESULTS *results);
////////////////////////////////////////////////////////////////////////////////



/// FourIndex.c ////////////////////////////////////////////////////////////////
int		FourIndex(CONFIGURATION param, INSTANCE data, RESULTS *results);
////////////////////////////////////////////////////////////////////////////////



// Benders.c ///////////////////////////////////////////////////////////////////
int		Benders(CONFIGURATION param, INSTANCE instance, RESULTS *results);
////////////////////////////////////////////////////////////////////////////////



// Benders_master.c ////////////////////////////////////////////////////////////
int		create_CPLEX_master_enviroment(CPXENVptr *env, CONFIGURATION param, GLOBAL_INFO info);
int		create_CPLEX_master_lp(CPXENVptr env, CPXLPptr *lp, GLOBAL_INFO info);
int		solve_root_node(CPXENVptr env, CPXLPptr lp, GLOBAL_INFO info);
int		clean_root_node(CPXENVptr env, CPXLPptr lp, GLOBAL_INFO global);
int		create_CPLEX_master_milp(CPXENVptr env, CPXLPptr lp, GLOBAL_INFO info);
int		solve_to_integrality(CPXENVptr env, CPXLPptr lp, GLOBAL_INFO info);
int		CPXPUBLIC lazy_cut_callback(CPXCENVptr env, void *cbdata, int wherefrom, void *info, int *useraction_p);
int		CPXPUBLIC user_cut_callback(CPXCENVptr env, void *cbdata, int wherefrom, void *info, int *useraction_p);
int		CPXPUBLIC rounding_callback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, double *objval_p, double *x, int *checkfeas_p, int *useraction_p);
int CPXPUBLIC set_branch_callback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int brtype, int sos, int nodecnt, int bdcnt, const int *nodebeg, const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p);
int		add_local_cuts(GLOBAL_INFO global, CPXCENVptr env, void *cbdata, int wherefrom, double *lb, double *ub);
int		add_OPT_root_cuts(GLOBAL_INFO global, CPXCENVptr env, CPXLPptr lp, int INTEGER);
int		add_OPT_callback_cuts(GLOBAL_INFO global, CPXCENVptr env, void *cbdata, int wherefrom, int INTEGER);
int		add_FEAS_cuts(GLOBAL_INFO global, CPXCENVptr env, CPXLPptr lp, int INTEGER);
int		add_FEAS_callback_cuts(GLOBAL_INFO global, CPXCENVptr env, void *cbdata, int wherefrom, int INTEGER);
double	OPT_violation(INSTANCE *data, double *subproblem_variables, int k, double *master_variables);
void	print_OPT_CUT(INSTANCE *data, double *subproblem_variables, int k);
double	SUBTOUR_violation(INSTANCE *data, int *cut_set, double *master_variables);
double	CAPACITY_violation(INSTANCE *data, int *cut_set, double *master_variables);
////////////////////////////////////////////////////////////////////////////////



// Benders_subproblem.c ////////////////////////////////////////////////////////
int		solve_feasibility_subproblem(GLOBAL_INFO global, int INTEGER);
int		solve_optimality_subproblems(GLOBAL_INFO global, int INTEGER);
int		get_INTEGER_obj_values(GLOBAL_INFO global, double *obj_values);
int		solve_benders_subproblems(GLOBAL_INFO global, double *obj_values);
int		solve_INTEGER_magnanti_wong_subproblems(GLOBAL_INFO global, double *obj_values);
int		solve_magnanti_wong_subproblems(GLOBAL_INFO global, double *obj_values, int FLag);
int		solve_papadakos_subproblems(GLOBAL_INFO global);
void	update_core_point(GLOBAL_INFO global, double **master);
double	lift_coeff(GLOBAL_INFO global, int e, int k, double **master_variables);
void	print_tree(GLOBAL_INFO global);
void	print_encoded_tree(GLOBAL_INFO global);
////////////////////////////////////////////////////////////////////////////////



// CLM_algorithm.c /////////////////////////////////////////////////////////////
int		randint(int ini, int end);
double	list_max(int N, double *list);
double	list_min(int N, double *list);
double	list_sum(int N, double *list);
double	list_prod(int N, double *list);
double	find_smallest(int N, double *list, int k);
int		ARE_EQUAL(int N, int *arr1, int *arr2);
double	elapsed(clock_t init_time);
////////////////////////////////////////////////////////////////////////////////



// CLM_graph.c /////////////////////////////////////////////////////////////////
void	Floyd_Warshall(int N, double **D);
double	LongestHamiltonianPathUB(int N, double **C);
void	tree_dist(int N, double **C, double **y, double **dist);
void	compute_min_cut_data(int N, double **C, double **D, double **min_cut, double **min_cut_min_cost);
void	GH_tree(int N, double **Capacity, int *Tree, double *Flow);
double	min_cap_cuts(int N, double **capacity, int **cut_sets);
double	connected_components(int N, double **capacity, int **cut_sets);
double	maxFlowAndCut(int N, int s, int t, double **Capacity, int *visited);
void	random_tree(int N, double **Y);
void	tree_to_code(int N, int **neighbours, int *degree, int *code);
void	matrix_to_code(int N, double **Y, int *code);
void	code_to_tree(int N, int *code, int **neighbours, int *degree);
void	code_to_matrix(int N, int *code, double **Y);
////////////////////////////////////////////////////////////////////////////////



/// CLM_OCSTP.c ////////////////////////////////////////////////////////////////
double	evaluate_tree(int N, double **C, double **D, double **Y);
void	rounding_heuristic(int N, double **Y, int **info);
double	Reinelt_lower_bound(int N, double **W, double **C, double **C2);
double	Elena_lower_bound(int N, double **W, double **C);
//void	compute_local_C_C2(INSTANCE data, double **local_C1, double **local_C2, int **INFO);
void	MST(int N, double **C, double **Y, int **info);
void	GHT(int N, double **W, double **arc);
void	STARS(int N, double **C, double **W, double **Y);
void	AMT(int N, double **C, double **W, double **Y, int **info);
void	MCMCT(int N, double **CUT, double **COST, double **Y, int **info);
double	OCSTP_DaC(int N, double **C, double **W, double **capacity, int SPLIT_TYPE, int IMPROVE, int FINAL_IMPROVE);
double  OCSTP_BUDaC(int N, double **C, double **W, double **capacity, int IMPROVE, int FINAL_IMPROVE);
void	rec_OCSTP_DaC(int N, double **C, double **W, int **T, double **Dist, int *V, int SPLIT_TYPE, int IMPROVE);
void	OCSTP_DaC_SPLIT_0(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_SPLIT_1(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_SPLIT_2(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_SPLIT_3(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_SPLIT_4(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_SPLIT_5(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_SPLIT_6(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_SPLIT_7(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_SPLIT_8(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_SPLIT_9(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_SPLIT_10(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_SPLIT_11(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_SPLIT_12(int N, double **C, double **W, int *S, int *S1, int *S2);
void	OCSTP_DaC_MERGE(int N, double **C, double **W, int **T, double **Dist, int *S1, int *S2);
void	OCSTP_DaC_IMPROVE(int N, double **C, double **W, int **T, double **Dist, int *S);
double	AM_LS(int N, double **C, double **W, double **Y, double **Dist, int **info);
double  AM_LS2(GLOBAL_INFO global, int **info, int BEST_IMPROVEMENT);
double	AM_LS3(int N, double **C, double **W, double **Y, int **info, int BEST_IMPROVEMENT);
double  DLS(int N, double **C, double **W, double **Y);
double  DLS2(int N, double **C, double **W, double **Y);
////////////////////////////////////////////////////////////////////////////////

double FinalLP;
int countfails;
