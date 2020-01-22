#include "headers.h"

/// INFO ///////////////////////////////////////////////////////////////////////
//
// benders.c: Main Benders Procedure
//
// Author:  Carlos Luna-Mota
// Created: 2015 Jan 21
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
int Benders(CONFIGURATION param, INSTANCE data, RESULTS *results) {

	int status = 0;
	int N = data.N;
	int E = data.E;
	int K = data.K;

	CPXENVptr env = NULL;     // pointer to Cplex Enviroment
	CPXLPptr  lp  = NULL;     // pointer to Cplex Linear Program data structure
    
	/////////////////////////////////////////////////////////////////////////////
	// Structs of Global Variables:
	BENDERS_VARIABLES var;
	GLOBAL_INFO	global;
	var.master			= create_double_vector(E+K);
	var.subproblem		= create_double_matrix(K, E+N);
	var.cut_sets		= create_int_matrix(N, N);
	var.capacity		= create_double_matrix(N, N);
	var.dist			= create_double_matrix(N, N);	
	var.init_time		= clock();
	var.last_id				= NONE;
	var.last_id_visits		= 0;
	var.flag_initial_sol	= NO;
	var.flag_generate_cuts	= YES;
	var.flag_fist_master	= YES;
	var.core_point		= create_double_vector(E);
	var.lambda			= 0.0;
	global.data			= &data;
	global.param		= &param;
	global.results		= results;		
	global.var			= &var;
	update_core_point(global, NULL);
	/////////////////////////////////////////////////////////////////////////////
    	
	// Create the CPLEX Enviroment
	if (elapsed(var.init_time) < param.MAX_CPU_TIME) {
		status = create_CPLEX_master_enviroment(&env, param, global);
		assert(status == 0);
	}

	// Generate the initial LP:
	if (elapsed(var.init_time) < param.MAX_CPU_TIME) {
		if (param.SCREEN_OUTPUT >= 1) { printf("\nGENERATING BASE MODEL:\n"); }
		status = create_CPLEX_master_lp(env, &lp, global);
		assert(status == 0);	
	}
	
	// Now solve the root node:
	if (elapsed(var.init_time) < param.MAX_CPU_TIME) {
		if (param.SCREEN_OUTPUT >= 1) { printf("\nSOLVING ROOT NODE:\n"); }
		status = solve_root_node(env, lp, global);
		assert(status == 0);

		if (param.SCREEN_OUTPUT >= 1) { printf("\nCLEANING ROOT NODE:\n"); }
		status = clean_root_node(env, lp, global);
		assert(status == 0);

		results->root_time = elapsed(var.init_time);
	}

	// Solve to integrality
	if (elapsed(var.init_time) < param.MAX_CPU_TIME							// Unless we run out of time...
		&& param.TO_INTEGRALITY == YES										// or we are told to just solve the LP...
		&& results->UpperBound - results->LowerBound > param.MIN_ABS_GAP) {	// or we have already solved the problem...

		// Transform the LP into a MILP and set the callbacks
		status = create_CPLEX_master_milp(env, lp, global);
		assert(status == 0);

		if (param.SCREEN_OUTPUT >= 1) { printf("\nSOLVING INTEGER PROBLEM:\n"); }

		// And then solve the problem to integrality
		status = solve_to_integrality(env, lp, global);
		assert(status == 0);
	}

	results->total_time = elapsed(var.init_time);
		
	// Clean:
	free_double_vector(&(var.master)); 
	free_double_matrix(&(var.subproblem), K);
	free_int_matrix(&(var.cut_sets), N);
	free_double_matrix(&(var.capacity), N);
	free_double_matrix(&(var.dist), N);
	free_double_vector(&(var.core_point));

	status = CPXfreeprob(env, &lp);
	assert(status == 0);

	status = CPXcloseCPLEX(&env);
    assert(status == 0);
    
	return status;
}
////////////////////////////////////////////////////////////////////////////////
