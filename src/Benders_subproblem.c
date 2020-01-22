#include "headers.h"


/// INFO ///////////////////////////////////////////////////////////////////////
//
// subproblem.c: SubProblem Functions
//
// Author:  Carlos Luna-Mota
// Created: 2015 Jan 21
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int solve_feasibility_subproblem(GLOBAL_INFO global, int INTEGER) {

	int	i, j;
	int status = 0;
	int N = global.data->N;
	double NN = ((double)(N));
	int **index_e = global.data->index_e;

	double **capacity = create_double_matrix(N, N);
	double min_flow;

	// Build the capacity matrix:
	for (i=0; i<N; i++) {
		capacity[i][i] = 0.0; 
		for (j=i+1; j<N; j++) {
			capacity[i][j] = capacity[j][i] = global.var->master[index_e[i][j]];
		}
	}

	// Now find the candidates to min_cuts:
	if (INTEGER == YES) { min_flow = connected_components(N, capacity, global.var->cut_sets); }
	else                { min_flow = min_cap_cuts(N, capacity, global.var->cut_sets);   }

	// Remember: 
	global.results->n_feas_subproblems++;
	if (INTEGER == YES) { global.results->n_int_feas_subproblems++;  }
	else                { global.results->n_frac_feas_subproblems++; }

	// clean
	free_double_matrix(&capacity, N);

	return status;
} 
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int solve_optimality_subproblems(GLOBAL_INFO global, int INTEGER) {

	double	obj_value = 0.0;
	int		i, j, e, k, d;
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	double	*obj_values = create_double_vector(K);
	int		status = 0;
	int		RESET_CORE = NO;

	// Put the Master Variables in a Matrix:
	for (i=0; i<N; i++) {
		global.var->capacity[i][i] = 0.0;
		for (j=i+1; j<N; j++) {
			global.var->capacity[i][j] = global.var->master[global.data->index_e[i][j]];
			global.var->capacity[j][i] = global.var->capacity[i][j]; 
		}
	}	
	
	// Now solve the optimality subproblem:
	if (global.param->OPT_CUT_TYPE == 0) {			// Regular Primal Subproblem
				
		status = solve_benders_subproblems(global, obj_values);
	
	} else if(global.param->OPT_CUT_TYPE == 1) {		// Magnanti-Wong Pareto Dual Subproblem
		
		// Update Core-Point:
		if (global.param->STAR_LIKE_CORE_POINT == YES) { update_core_point(global, NULL); }		
		update_core_point(global, global.var->capacity); 

		// Compute Obj-Function:
		if (INTEGER == YES) { get_INTEGER_obj_values(global, obj_values);             }
		//else                { status = solve_benders_subproblems(global, obj_values); }

		// Generate the cuts:
		if (global.var->flag_generate_cuts == YES) {			
			//if(INTEGER == YES) status = solve_benders_subproblems(global, obj_values);
			status = solve_magnanti_wong_subproblems(global, obj_values, INTEGER);
		}

		if (INTEGER == YES && global.var->flag_generate_cuts == NO) {
			global.var->flag_generate_cuts = YES;
			status = solve_benders_subproblems(global, obj_values);
		}
		
	} else if(global.param->OPT_CUT_TYPE == 2) {		// Papadakos Pareto Primal Subproblem
				
		// Update Core-Point:
		if (global.param->STAR_LIKE_CORE_POINT == YES) { update_core_point(global, NULL); }		
		update_core_point(global, global.var->capacity); 

		// Compute Obj-Function:
		if (INTEGER == YES) { get_INTEGER_obj_values(global, obj_values); } 

		// Compute Cuts
		status = solve_papadakos_subproblems(global);

		// For integer subproblems, if Papadakos is infeasible we should generate regular benders cuts:
		if (INTEGER == YES && global.var->flag_generate_cuts == NO) {
			global.var->flag_generate_cuts = YES;
			status = solve_benders_subproblems(global, obj_values);
		}

	} else {

		// Non Implemented:
		printf("ERROR: Unknown Optimality Cut Type: %d ", global.param->OPT_CUT_TYPE);
		status = 1;
	}

	// Try to lift some coeff
	if (global.param->HEUR_LIFTING == YES && INTEGER == YES && global.data->IS_EUCLIDEAN == YES) {
		for (k=0; k<K; k++) {
			d = global.data->Com[k].d;
			for (e=0; e<E; e++) {
				if (global.var->subproblem[k][e] < EPSILON) {
					global.var->subproblem[k][e]    = lift_coeff(global, e, k, global.var->capacity);	// Solve lifting problem
					global.var->subproblem[k][E+d] += global.var->subproblem[k][e];	// Update RHS
				}
			}
		}
	}

	// For every integer subproblem we must also check if the current solution is the best one and update the gap
	if (INTEGER == YES) {
		obj_value = 0.0;
		for (k=0; k<K; k++) { obj_value += global.data->Com[k].w * obj_values[k]; }
		if (global.results->UpperBound >= obj_value) {
			global.results->UpperBound = obj_value;
			for (e=0; e<E; e++) { global.results->best_solution[e] = global.var->capacity[global.data->index_i[e]][global.data->index_j[e]]; }
			for (k=0; k<K; k++) { global.results->best_solution[E+k] = obj_values[k]; }

			// Update GAP
			global.results->final_gap = (global.results->UpperBound-global.results->LowerBound) / global.results->UpperBound;
			if (global.param->SCREEN_OUTPUT >= 1) { 	printf("SubPro\t\t\t%.2f <- %.2f%% -> %.2f  \tTime: %.2f\"\n", global.results->LowerBound, 100*global.results->final_gap, global.results->UpperBound, elapsed(global.var->init_time)); }
		}
	}


	// Remember:
	global.results->n_opt_subproblems++;      
	if (INTEGER == YES) { global.results->n_int_opt_subproblems++;  }
	else                { global.results->n_frac_opt_subproblems++; }

	// clean
	free_double_vector(&obj_values);
	
	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int get_INTEGER_obj_values(GLOBAL_INFO global, double *obj_values) {
	
	int		status = 0;
	int		k;
	int		N = global.data->N;
	int		K = global.data->K;
	COMMODITY *Com = global.data->Com;

	// Compute Distances:
	tree_dist(N, global.data->C, global.var->capacity, global.var->dist);	// Compute distances by DFS assuming that capacity holds a tree.
		
	// Compute Obj Values:
	for (k=0; k<K; k++) { obj_values[k] = global.var->dist[Com[k].o][Com[k].d]; }

	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int solve_benders_subproblems(GLOBAL_INFO global, double *obj_values){
	
	int		solstat, status = 0;
	int		i, k, e;

	// Problem Data:
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int		**index_e = global.data->index_e;
	int		*index_i  = global.data->index_i;
	int		*index_j  = global.data->index_j;
	COMMODITY *Com = global.data->Com;

	// Allocate memory for CPLEX PROBLEM DEFINITION
	CPXENVptr env = NULL;
	CPXNETptr net = NULL;
	double	 *supply = create_double_vector(N);
	int		 *tail	 = create_int_vector(2*E);
	int	     *head	 = create_int_vector(2*E);
	double	 *obj	 = create_double_vector(2*E);
	double	 *ub	 = create_double_vector(2*E);
	double	 *lb	 = create_double_vector(2*E);
	// AND OUTPUT:
	double   *beta    = create_double_vector(N);
	
	double	available_time = MAX(0.0, global.param->MAX_CPU_TIME - elapsed(global.var->init_time)) / ((double)(K));
	///////////////////////////////////////////////////////////////////////////
	
	
	/// OPEN CPLEX ////////////////////////////////////////////////////////////
	env = CPXopenCPLEX(&status);
	assert(env != NULL);

	net = CPXNETcreateprob(env, &status, "PrimalSubproblem");
	assert(net != NULL);
	///////////////////////////////////////////////////////////////////////////
	
	
	// Define network structure: //////////////////////////////////////////////
	for (e=0; e<E; e++) {
		tail[e] = index_i[e];
		tail[e+E] = index_j[e];
		head[e] = index_j[e];
		head[e+E] = index_i[e];
		obj[e] = obj[e+E] = global.data->C[index_i[e]][index_j[e]];
		lb[e] = lb[e+E] = 0.0;
		ub[e] = ub[e+E] = global.var->capacity[index_i[e]][index_j[e]];
	}
	///////////////////////////////////////////////////////////////////////////


	/// SOLVE K Subproblems ///////////////////////////////////////////////////
	for (k=0; k<K; k++) {

		// Define the Supply:
		for (i=0; i<N; i++) { supply[i] = 0.0; }
		supply[Com[k].o] =  1.0;
		supply[Com[k].d] = -1.0;

		// Delete existing network (if any). 
		if (CPXNETgetnumnodes(env, net) > 0) {
			status = CPXNETdelnodes (env, net, 0, CPXNETgetnumnodes (env, net)-1);
			assert(status==0);
		}

		// Set optimization sense:
		status = CPXNETchgobjsen(env, net, CPX_MIN);
		assert(status==0);

		// Add nodes to network along with their supply values:
		status = CPXNETaddnodes (env, net, N, supply, NULL);
		assert(status==0);

		// Add arcs to network along with their objective values and bounds:
		status = CPXNETaddarcs(env, net, 2*E, tail, head, lb, ub, obj, NULL);
		assert(status==0);
			
		// Optimize the problem and obtain solution:
		status = CPXNETprimopt(env, net);
		assert(status == 0);

		// Check network dimensions
		assert(2*E == CPXNETgetnumarcs(env, net));
		assert( N  == CPXNETgetnumnodes(env, net));

		// Get solution data:
		status = CPXNETsolution(env, net, &solstat, &(obj_values[k]), NULL, beta, NULL, NULL);
		assert(status == 0);	
		
		// Check it!
		if (solstat != CPX_STAT_OPTIMAL) { 
			printf("Error: Benders cut[%d] failed.\n", k);
			global.var->flag_generate_cuts = NO;
			break;
		} else {
			// store it:
			for (i=0; i<N; i++) { global.var->subproblem[k][E+i] = -beta[i]; }
			for (e=0; e<E; e++) {
				global.var->subproblem[k][e] = MAX(0, ABS(beta[index_i[e]] - beta[index_j[e]])-global.data->C[index_i[e]][index_j[e]]);
				if(global.var->subproblem[k][e] < 0.01) { global.var->subproblem[k][e] = 0; }
			}
		}
	}

	// Clean:	   
	free_double_vector(&supply);
	free_double_vector(&obj);
	free_double_vector(&ub);
	free_double_vector(&lb);
	free_int_vector(&tail);
	free_int_vector(&head);
	
	free_double_vector(&beta);
	
	status = CPXNETfreeprob(env, &net);
	assert(status == 0);

	status = CPXcloseCPLEX(&env);
	assert(status == 0);

	return status;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int solve_magnanti_wong_subproblems(GLOBAL_INFO global, double *obj_values, int flag) {

    int		solstat, status = 0;
	int		i, k, e;
	double	obj_primal_sub;
	double x0;

	// Problem Data:
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int		**index_e = global.data->index_e;
	int		*index_i  = global.data->index_i;
	int		*index_j  = global.data->index_j;
	COMMODITY *Com = global.data->Com;

	// Allocate memory for CPLEX PROBLEM DEFINITION
	CPXENVptr env = NULL;
	CPXNETptr net = NULL;
	double	 *supply = create_double_vector(N);
	int		 *tail	 = create_int_vector(2*E);
	int	     *head	 = create_int_vector(2*E);
	double	 *obj	 = create_double_vector(2*E);
	double	 *ub	 = create_double_vector(2*E);
	double	 *lb	 = create_double_vector(2*E);
	// AND OUTPUT:
	double   *beta    = create_double_vector(N);
	
	double	available_time = MAX(0.0, global.param->MAX_CPU_TIME - elapsed(global.var->init_time)) / ((double)(K));
	///////////////////////////////////////////////////////////////////////////
	
	
	/// OPEN CPLEX ////////////////////////////////////////////////////////////
	env = CPXopenCPLEX(&status);
	assert(env != NULL);

	net = CPXNETcreateprob(env, &status, "PrimalSubproblem");
	assert(net != NULL);
	///////////////////////////////////////////////////////////////////////////
	x0=0;
	//Calculating x0////////////////////////////////////////////////////////////
	for(e=0;e<E;e++){
		x0=x0+ global.var->core_point[e];
	}
	///////////////////////////////////////////////////////////////////////////
	// Define network structure: //////////////////////////////////////////////
	for (e=0; e<E; e++) {
		tail[e] = index_i[e];
		tail[e+E] = index_j[e];
		head[e] = index_j[e];
		head[e+E] = index_i[e];
		obj[e] = obj[e+E] = global.data->C[index_i[e]][index_j[e]];
		lb[e] = lb[e+E] = 0.0;
		ub[e] = ub[e+E] = global.var->core_point[e]+x0*global.var->master[e];
		//if(ub[e]>1)printf("ub[%d]=%lf\n", e, ub[e]);
	}
	//getchar();
	///////////////////////////////////////////////////////////////////////////


	/// SOLVE K Subproblems ///////////////////////////////////////////////////
	for (k=0; k<K; k++) {

		// Define the Supply:
		for (i=0; i<N; i++) { supply[i] = 0.0; }
		supply[Com[k].o] =  1.0+x0;
		supply[Com[k].d] = -1.0-x0;

		// Delete existing network (if any). 
		if (CPXNETgetnumnodes(env, net) > 0) {
			status = CPXNETdelnodes (env, net, 0, CPXNETgetnumnodes (env, net)-1);
			assert(status==0);
		}

		// Set optimization sense:
		status = CPXNETchgobjsen(env, net, CPX_MIN);
		assert(status==0);

		// Add nodes to network along with their supply values:
		status = CPXNETaddnodes (env, net, N, supply, NULL);
		assert(status==0);

		// Add arcs to network along with their objective values and bounds:
		status = CPXNETaddarcs(env, net, 2*E, tail, head, lb, ub, obj, NULL);
		assert(status==0);
			
		// Optimize the problem and obtain solution:
		status = CPXNETprimopt(env, net);
		assert(status == 0);

		// Check network dimensions
		assert(2*E == CPXNETgetnumarcs(env, net));
		assert( N  == CPXNETgetnumnodes(env, net));

		// Get solution data:
		status = CPXNETsolution(env, net, &solstat, &obj_primal_sub, NULL, beta, NULL, NULL);
		assert(status == 0);
		if (solstat != CPX_STAT_OPTIMAL) { 
			printf("Warning: Papadakos[%d] infeasible for !\n", k);
			
			countfails++;
			global.var->flag_generate_cuts = NO;
			break;
		} else {
			for (i=0; i<N; i++) { global.var->subproblem[k][E+i] = -beta[i]; }
			for (e=0; e<E; e++) {
				global.var->subproblem[k][e] = MAX(0, ABS(beta[index_i[e]] - beta[index_j[e]])-global.data->C[index_i[e]][index_j[e]]);
				if(global.var->subproblem[k][e] < 0.01) { global.var->subproblem[k][e] = 0; }
			}
		}
	}

	// Clean:	   
	free_double_vector(&supply);
	free_double_vector(&obj);
	free_double_vector(&ub);
	free_double_vector(&lb);
	free_int_vector(&tail);
	free_int_vector(&head);
	
	free_double_vector(&beta);
	
	status = CPXNETfreeprob(env, &net);
	assert(status == 0);

	status = CPXcloseCPLEX(&env);
	assert(status == 0);

	return status;

}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
int solve_papadakos_subproblems(GLOBAL_INFO global) {
		
	int		solstat, status = 0;
	int		i, k, e;
	double	obj_primal_sub;

	// Problem Data:
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int		**index_e = global.data->index_e;
	int		*index_i  = global.data->index_i;
	int		*index_j  = global.data->index_j;
	COMMODITY *Com = global.data->Com;

	// Allocate memory for CPLEX PROBLEM DEFINITION
	CPXENVptr env = NULL;
	CPXNETptr net = NULL;
	double	 *supply = create_double_vector(N);
	int		 *tail	 = create_int_vector(2*E);
	int	     *head	 = create_int_vector(2*E);
	double	 *obj	 = create_double_vector(2*E);
	double	 *ub	 = create_double_vector(2*E);
	double	 *lb	 = create_double_vector(2*E);
	// AND OUTPUT:
	double   *beta    = create_double_vector(N);
	
	double	available_time = MAX(0.0, global.param->MAX_CPU_TIME - elapsed(global.var->init_time)) / ((double)(K));
	///////////////////////////////////////////////////////////////////////////
	
	
	/// OPEN CPLEX ////////////////////////////////////////////////////////////
	env = CPXopenCPLEX(&status);
	assert(env != NULL);

	net = CPXNETcreateprob(env, &status, "PrimalSubproblem");
	assert(net != NULL);
	///////////////////////////////////////////////////////////////////////////
	
	// Define network structure: //////////////////////////////////////////////
	for (e=0; e<E; e++) {
		tail[e] = index_i[e];
		tail[e+E] = index_j[e];
		head[e] = index_j[e];
		head[e+E] = index_i[e];
		obj[e] = obj[e+E] = global.data->C[index_i[e]][index_j[e]];
		lb[e] = lb[e+E] = 0.0;
		ub[e] = ub[e+E] = global.var->core_point[e];
	}
	///////////////////////////////////////////////////////////////////////////


	/// SOLVE K Subproblems ///////////////////////////////////////////////////
	for (k=0; k<K; k++) {

		// Define the Supply:
		for (i=0; i<N; i++) { supply[i] = 0.0; }
		supply[Com[k].o] =  1.0;
		supply[Com[k].d] = -1.0;

		// Delete existing network (if any). 
		if (CPXNETgetnumnodes(env, net) > 0) {
			status = CPXNETdelnodes (env, net, 0, CPXNETgetnumnodes (env, net)-1);
			assert(status==0);
		}

		// Set optimization sense:
		status = CPXNETchgobjsen(env, net, CPX_MIN);
		assert(status==0);

		// Add nodes to network along with their supply values:
		status = CPXNETaddnodes (env, net, N, supply, NULL);
		assert(status==0);

		// Add arcs to network along with their objective values and bounds:
		status = CPXNETaddarcs(env, net, 2*E, tail, head, lb, ub, obj, NULL);
		assert(status==0);
			
		// Optimize the problem and obtain solution:
		status = CPXNETprimopt(env, net);
		assert(status == 0);

		// Check network dimensions
		assert(2*E == CPXNETgetnumarcs(env, net));
		assert( N  == CPXNETgetnumnodes(env, net));

		// Get solution data:
		status = CPXNETsolution(env, net, &solstat, &obj_primal_sub, NULL, beta, NULL, NULL);
		assert(status == 0);
		if (solstat != CPX_STAT_OPTIMAL) { 
			printf("Warning: Papadakos[%d] infeasible!\n", k); 
			global.var->flag_generate_cuts = NO;
			break;
		} else {
			for (i=0; i<N; i++) { global.var->subproblem[k][E+i] = -beta[i]; }
			for (e=0; e<E; e++) {
				global.var->subproblem[k][e] = MAX(0, ABS(beta[index_i[e]] - beta[index_j[e]])-global.data->C[index_i[e]][index_j[e]]);
				if(global.var->subproblem[k][e] < 0.01) { global.var->subproblem[k][e] = 0; }
			}
		}
	}

	// Clean:	   
	free_double_vector(&supply);
	free_double_vector(&obj);
	free_double_vector(&ub);
	free_double_vector(&lb);
	free_int_vector(&tail);
	free_int_vector(&head);
	
	free_double_vector(&beta);
	
	status = CPXNETfreeprob(env, &net);
	assert(status == 0);

	status = CPXcloseCPLEX(&env);
	assert(status == 0);

	return status;

}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void update_core_point(GLOBAL_INFO global, double **master) {
	int e;
	double valrisk=1;
	int conta=0;

	double NN = ((double) (global.data->N));

	if (master == NULL) {
		for (e=0; e<global.data->E; e++) { 
			//global.var->core_point[e] = 2/NN;
			if(global.results->best_solution[e]>0.5){
				global.var->core_point[e] = valrisk;
				//conta++;
				//printf("arc %d=%lf\n", conta,global.var->core_point[e]);
			}
			else{
				global.var->core_point[e] = .3/*((1-valrisk)/**(NN-1))/(global.data->E-NN+1)*/;
			}
		}
		
	} else {
		for (e=0; e<global.data->E; e++) { 
			global.var->core_point[e] *= (1.0 - global.var->lambda);
			global.var->core_point[e] += global.var->lambda * master[global.data->index_i[e]][global.data->index_j[e]];
		}
	}
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
double lift_coeff(GLOBAL_INFO global, int e, int k, double **master_variables){
	
	int		N = global.data->N;
	int		E = global.data->E;
	int		K = global.data->K;
	int		**index_e = global.data->index_e;
	int		o = global.data->Com[k].o;
	int		d = global.data->Com[k].d;
	int		i = global.data->index_i[e];
	int		j = global.data->index_j[e];
	int		index_min;

	double	*dist	 = create_double_vector(N);
	int		*visited = create_int_vector(N);

	double gamma_e = 0.0;
	
	if (global.var->subproblem[k][e] <= EPSILON && master_variables[i][j] >= 1.0-EPSILON) {
		
		// Gamma_e = C2[o][d] - (Beta[d]-Beta[o]) + ...   <=>  e == {o,d}
		//         = C1[o][d] - (Beta[d]-Beta[o]) + ...   <=>  e != {o,d}
		if (index_e[o][d] == e) { gamma_e += global.data->C2[o][d]; }
		else                    { gamma_e += global.data->C1[o][d];  }
		gamma_e -= global.var->subproblem[k][E+d] - global.var->subproblem[k][E+o];
	
		// Gamma_e = ... + MST(Gamma) such that e not in the Tree
		for (i=0; i<N; i++) { visited[i] = NO; }	// Nobody is in the tree
		visited[0] = YES;							// 0 is in the tree
		for (i=1; i<N; i++) { dist[i] = global.var->subproblem[k][index_e[0][i]]; }
		dist[0]    = 0.0;							// Distances to a tree node		
		for (i=1; i<N; i++) {			// now add the other N-1 nodes to the tree
			index_min = NONE;			// best candidate (closer)
			for (j=1; j<N; j++) {
				if (visited[j] == NO) {
					if (index_min == NONE || dist[index_min] > dist[j]) {
						index_min = j; 
					}
				}
			}
			assert(index_min != NONE);
			// Now add this node and update distances
			visited[index_min] = YES;
			for (j=1; j<N; j++) {
    			if (visited[j] == NO) {
    				dist[j] = MIN(dist[j], global.var->subproblem[k][index_e[index_min][j]]);
    			}
			}
		}
		for (i=1; i<N; i++) { gamma_e += dist[i]; } // Compute Spanning tree Cost
	}

	// clean:
	free_double_vector(&dist);
	free_int_vector(&visited);
	
	return MAX(0.0, gamma_e);
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
void print_encoded_tree(GLOBAL_INFO global) {
	
	int i,j, min;
	int N = global.data->N;
	printf("\t\t0");
	for (i=1; i<N; i++) {
		min = N;
		for (j=N-1; j>=0; j--) {
			if (global.var->master[global.data->index_e[i][j]] > EPSILON) {
				min = j;
			}
		}
		printf(" %d",min); 
	}
	printf("\n");	
}
void print_tree(GLOBAL_INFO global) {
	
	int e, ee, new_line;

	ee = 0;
	new_line = YES;
	printf("\tTree:");
	for (e=0; e<global.data->E; e++) {
		if (global.var->master[e] > EPSILON) {
			printf("\t(%d,%d)",global.data->index_i[e],global.data->index_j[e]); 
			ee++;
			new_line = NO;
		}	
		if (ee%7 == 0 && new_line == NO && ee+1 <= global.data->N) { 
			new_line = YES;
			printf("\n\t");
		}
	} printf("\n");	
}
////////////////////////////////////////////////////////////////////////////////
