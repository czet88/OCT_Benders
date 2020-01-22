#include "headers.h"



/// INFO ///////////////////////////////////////////////////////////////////////
//
// CLM_OCSTP.c: Optimum Communication Spanning Tree Problem related functions
//
// Author:  Carlos Luna-Mota
// Created: 2015 Jan 21
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////


/////////////
// GENERAL //
/////////////


////////////////////////////////////////////////////////////////////////////////
// Evaluates a tree stored in the symmetric binary matrix Y using the edge-cost
// matrix C and the communication requests matrix D.
double evaluate_tree(int N, double **C, double **D, double **Y) {
	int		i,j;
	double	**dist = create_double_matrix(N,N);
	double	result = 0.0;
	
	// Compute the distances for every pair of nodes:
	tree_dist(N, C, Y, dist);

	// Compute the Communication Cost:
	for (i=0; i<N; i++) { for(j=0; j<N; j++) { result += dist[i][j] * D[i][j]; } }

	// clean
	free_double_matrix(&dist, N);

	return result;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Uses Prim's algorithm to find the MaxSpanningTree in O(n^2) time
void rounding_heuristic(int N, double **Y, int **info) {
    int		i, j, iteration, index_max;
    int		*pred = create_int_vector(N);
	int		*in_T = create_int_vector(N);
	double	*dist = create_double_vector(N);

	// Initialice everything:
	for (i=0; i<N; i++) { in_T[i] = NO; }			// all the others haven't been visited
	for (i=0; i<N; i++) { pred[i] = 0; }			// and have 0 as a predecesor
	for (i=0; i<N; i++) { dist[i] = Y[0][i]; }	// current distances to the tree
    
	in_T[0] = YES;	// 0 is the only visited node so far
    
	// during N-1 iterations add a new edge to the MST
	for (iteration=1; iteration<N; iteration++) {

		index_max = NONE;

		// look for a mandatory edges first:		
		if (info != NULL) {
			for (j=1; j<N; j++) {
				if (in_T[j] == NO && info[pred[j]][j] == YES) {
					index_max = j;
					break;
				}
			}
		}

		// If there aren't any mandatory edges look for the furthest one:
		if (index_max == NONE) {
			for (j=1; j<N; j++) {
				// Check if feasible:
				if (in_T[j] == NO && (info==NULL || info[pred[j]][j]!=NO)) {	
					// Check if optimal:
					if (index_max == NONE || dist[index_max] < dist[j]) {
						index_max = j; 
					}
				}
			}
		}
		assert(index_max != NONE);

		// We add it to the MST:
		in_T[index_max] = YES;

		// and update the distances:
    	for (j=1; j<N; j++) {
			// Easy case: we have no restrictions, just check for the best one:
			if (info == NULL) {
				if (in_T[j] == NO && dist[j] < Y[index_max][j]) {
    				dist[j] = Y[index_max][j];
    				pred[j] = index_max;
    			}
			// Difficult case: we have some restrictions:
			} else {
				// Check for feasibility:
    			if (in_T[j] == NO && info[index_max][j]!=NO) {
					// Check for optimality:
					if (info[index_max][j]==YES || info[j][pred[j]]==NO || (info[j][pred[j]]==NONE && dist[j] < Y[index_max][j])) {
    					dist[j] = Y[index_max][j];
    					pred[j] = index_max;
					}
    			}
			}
		}
	}

	// output the result in the "y" matrix
	for (i=0; i<N; i++) {
		Y[i][i] = 0.0;		// Main diagonal
		for (j=i+1; j<N; j++) {
			if (pred[i] == j || pred[j] == i) { Y[i][j] = Y[j][i] = 1.0; } // Edge in the tree 
			else                              { Y[i][j] = Y[j][i] = 0.0; } // Edge not in the tree
		}
	}

	free_int_vector(&pred);
	free_double_vector(&dist);
	free_int_vector(&in_T);
}
////////////////////////////////////////////////////////////////////////////////



//////////////////
// LOWER BOUNDS //
//////////////////


////////////////////////////////////////////////////////////////////////////////
// Computes Reinelt's Lower Bound for OCSTP:
double Reinelt_lower_bound(int N, double **W, double **C, double **C2) {
	int i,j;
	double lwr = 0.0;

	double	**D	= create_double_matrix(N,N);
	double  **Y	= create_double_matrix(N,N);
	
	/// REINELT'S LwrBound /////////////////////////////////////////////////////
	
	// Initialize
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { Y[i][j] = Y[j][i] = 0.0; } }
	for (i=0; i<N; i++) {
		D[i][i] = 0.0;
		for (j=i+1; j<N; j++) { D[i][j] = D[j][i] = (W[i][j] + W[j][i])*(C[i][j] - C2[i][j]); }
	}

	// Compute
	MST(N, D, Y, NULL);
	for (i=0; i<N; i++) {
		for (j=i+1; j<N; j++) {			
			lwr += (W[i][j] + W[j][i])* C2[i][j];
			lwr += Y[i][j]*D[i][j];
		}
	}
	
	// clean:
	free_double_matrix(&D, N);
	free_double_matrix(&Y, N);
	
	// Return the Reinelt's lower bound:
	return lwr;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Computes Elena's Lower Bound for the OCSTP
double Elena_lower_bound(int N, double **W, double **C) {
	int i,j,e,k;
	double lwr;

	double	**Capacity	= create_double_matrix(N,N);
	double  **Y			= create_double_matrix(N,N);
	double	*Dist		= create_double_vector(N);
	double	*Flow		= create_double_vector(N);
	int		*Tree		= create_int_vector(N);	
	int		*VisitF		= create_int_vector(N);
	int		*VisitD		= create_int_vector(N);
	int		best_f, best_d;

	// ELENA'S LwrBound ////////////////////////////////////////////////////////
	
	// Initialize
	for (i=0; i<N; i++) { Tree[i] = NONE; }
	for (i=0; i<N; i++) { Flow[i] =  0.0; }
	for (i=0; i<N; i++) { Dist[i] = -1.0; }
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { Y[i][j] = 0.0; } }
	for (i=0; i<N; i++) {
		Capacity[i][i] = 0.0;
		for (j=i+1; j<N; j++) { Capacity[i][j] = Capacity[j][i] = (W[i][j] + W[j][i]); 	}
	}
	
	// Compute Gomory-Hu Tree and MST
	GH_tree(N, Capacity, Tree, Flow);
	MST(N, C, Y, NULL);
	e = 0;
	for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { if (Y[i][j] >= 0.5) { Dist[e++] = C[i][j]; } } }
	assert(e == N-1);
	
	// Compute Lower Bound
	lwr = 0.0;
	for (i=0; i<N; i++) { VisitF[i] = NO; }
	for (i=0; i<N; i++) { VisitD[i] = NO; }
	VisitF[0]   = YES;
	VisitD[N-1] = YES;
	for (k=0; k<N-1; k++) {
		best_f = NONE;
		best_d = NONE;
		for (i=0; i<N; i++) { 
			if (VisitF[i] == NO && (best_f == NONE || Flow[i] < Flow[best_f])) { best_f = i; }
			if (VisitD[i] == NO && (best_d == NONE || Dist[i] > Dist[best_d])) { best_d = i; }
		}
		assert(best_f != NONE && best_d != NONE);
		VisitD[best_d] = YES;
		VisitF[best_f] = YES;
		lwr += Flow[best_f] * Dist[best_d];
	}
	
	////////////////////////////////////////////////////////////////////////////

	// Clean:
	free_double_matrix(&Capacity, N);
	free_double_matrix(&Y, N);
	free_double_vector(&Dist);
	free_double_vector(&Flow);
	free_int_vector(&Tree);
	free_int_vector(&VisitF);
	free_int_vector(&VisitD);

	// Return the Elena's lower bound:
	return lwr;

}
////////////////////////////////////////////////////////////////////////////////



/////////////////////////////
// CONSTRUCTIVE HEURISTICS //
/////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Good old O(N^2) Prim's algorithm for finding the MST.
void MST(int N, double **cost, double **arc, int **info) {
    int i, j, iteration, index_min;
    int		*pred = create_int_vector(N);
	int		*in_T = create_int_vector(N);
	double	*dist = create_double_vector(N);

	// Initialice everything:
	for (i=0; i<N; i++) { in_T[i] = NO; }			// all the others haven't been visited
	for (i=0; i<N; i++) { pred[i] = 0; }			// and have 0 as a predecesor
	for (i=0; i<N; i++) { dist[i] = cost[0][i]; }	// current distances to the tree
    
	in_T[0] = YES;	// 0 is the only visited node so far
    
	// during N-1 iterations add a new edge to the MST
	for (iteration=1; iteration<N; iteration++) {

		index_min = NONE;

		// look for a mandatory edges first:		
		if (info != NULL) {
			for (j=1; j<N; j++) {
				if (in_T[j] == NO && info[pred[j]][j] == YES) {
					index_min = j;
					break;
				}
			}
		}

		// If there aren't any mandatory edges look for the furthest one:
		if (index_min == NONE) {
			for (j=1; j<N; j++) {
				// Check if feasible:
				if (in_T[j] == NO && (info==NULL || info[pred[j]][j]!=NO)) {	
					// Check if optimal:
					if (index_min == NONE || dist[index_min] > dist[j]) {
						index_min = j; 
					}
				}
			}
		}
		assert(index_min != NONE);

		// We add it to the MST:
		in_T[index_min] = YES;

		// and update the distances:
    	for (j=1; j<N; j++) {
			// Easy case: we have no restrictions, just check for the best one:
			if (info == NULL) {
				if (in_T[j] == NO && dist[j] > cost[index_min][j]) {
    				dist[j] = cost[index_min][j];
    				pred[j] = index_min;
    			}
			// Difficult case: we have some restrictions:
			} else {
				// Check for feasibility:
    			if (in_T[j] == NO && info[index_min][j]!=NO) {
					// Check for optimality:
					if (info[index_min][j]==YES || info[j][pred[j]]==NO || (info[j][pred[j]]==NONE && dist[j] > cost[index_min][j])) {
    					dist[j] = cost[index_min][j];
    					pred[j] = index_min;
					}
    			}
			}
		}
	}

	// output the result in the "y" matrix
	for (i=0; i<N; i++) {
		arc[i][i] = 0.0;		// Main diagonal
		for (j=i+1; j<N; j++) {
			if (pred[i] == j || pred[j] == i) { arc[i][j] = arc[j][i] = 1.0; } // Edge in the tree 
			else                              { arc[i][j] = arc[j][i] = 0.0; } // Edge not in the tree
		}
	}

	free_int_vector(&pred);
	free_double_vector(&dist);
	free_int_vector(&in_T);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Computes the Gomory-Hu tree (with respect to R) and stores it in Y
void GHT(int N, double **W, double **Y) {
    int		i, j;
	double	**capacity = create_double_matrix(N, N);
    double	*flow	   = create_double_vector(N);
	int		*tree      = create_int_vector(N);

	assert(N > 0);
	assert(W != NULL);
	assert(Y != NULL);

	// Init everything:
	for (i=0; i<N; i++) { tree[i] = NONE; }
	for (i=0; i<N; i++) { flow[i] = 0.0;  }
	for (i=0; i<N; i++) {		// capacity must be symmetric :-/
		capacity[i][i] = 0.0;
		for (j=i+1; j<N; j++) {
			capacity[i][j] = capacity[j][i] = W[i][j]+W[j][i];
		}
	}
	
	// compute the Gomory-Hu Tree
	GH_tree(N, capacity, tree, flow);
		
	// output the result in the "y" matrix
	for (i=0; i<N; i++) {
        Y[i][i] = 0.0;			// Main diagonal
		for (j=i+1; j<N; j++) {
			if (tree[i] == j || tree[j] == i) { Y[i][j] = Y[j][i] = 1.0; } // Edge in the tree 
			else                              { Y[i][j] = Y[j][i] = 0.0; } // Edge not in the tree
		}
	}
	
	// Clean:
	free_double_matrix(&capacity, N); 
	free_double_vector(&flow);
	free_int_vector(&tree);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Finds the best star-shapped tree
void STARS(int N, double **C, double **W, double **Y) {
    int i,j, c, best_center;
	double value, best_value;

	assert(N > 0);
	assert(C != NULL);
	assert(W != NULL);
	assert(Y != NULL);

	// Evaluate 0-star:
	c = best_center = 0;
	value = 0.0;
	for (i=0; i<N; i++) {
		for (j=i+1; j<N; j++) {
			if (i==c || j==c) { value += (W[i][j]+W[j][i]) * (C[i][j]); }
			else      { value += (W[i][j]+W[j][i]) * (C[i][c]+C[c][j]); }
		}
	}
	best_value = value;

	// Evaluate all other stars:
	for (c=1; c<N; c++) {
		value = 0.0;
		for (i=0; i<N; i++) {
			for (j=i+1; j<N; j++) {
				if (i==c || j==c) { value += (W[i][j]+W[j][i]) * (C[i][j]); }
				else      { value += (W[i][j]+W[j][i]) * (C[i][c]+C[c][j]); }
			}
		}
		if (best_value > value) {
			best_value = value;
			best_center = c;
		}
	}

	// Remember the best:
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			if (i==j) { Y[i][j] = 0.0; }
			else if(i==best_center || j==best_center) { Y[i][j] = 1.0; }
			else { Y[i][j] = 0.0; }
		}
	}
}
////////////////////////////////////////////////////////////////////////////////
 

////////////////////////////////////////////////////////////////////////////////
// Constructive Heuristic for the OCSTP from Ahuja & Murty in O(N^3)
// Stores tree in Y
void AMT(int N, double **C, double **W, double **Y, int **info) {
	int		i, j, k, best_i, best_j, n_in_T;
	double	delta, best_delta, sum_w;
	double	**dist_in_G = create_double_matrix(N,N);
	double	**dist_in_T = create_double_matrix(N,N);
	double	*w = create_double_vector(N);
	double	*cc = create_double_vector(N);
	int		*in_T = create_int_vector(N);
	
	// Starting Point...
	int	SEED = 0;				// Deterministic version
	//int	SEED = randint(0, N-1); // Randomized version

	assert(N > 0);
	assert(C != NULL);
	assert(W != NULL);
	assert(Y != NULL);

	// Compute the minimum distances over G:
	for (i=0; i<N; i++) { for(j=0; j<N; j++) { dist_in_G[i][j] = C[i][j]; } }

	// Build the null graph:
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { Y[i][j] = 0.0; } }
	
	// Step 1: Initialization
	for (i=0; i<N; i++) { in_T[i] = NO; }	// Nobody is in the Tree
	in_T[SEED] = YES;						// SEED is in the Tree
	n_in_T = 1;								// Nodes in the Tree = 1
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { dist_in_T[i][j] = 0.0; } }
	
	while (n_in_T < N) {	// While someone outside the tree:
		// Compute w[i] and sum_w:
		sum_w = 0.0;			// Communication between tree and the outside
		for (i=0; i<N; i++) {
			w[i] = 0.0;			// Comminucation between 'i' and the other side
			for (j=0; j<N; j++) {
				if (in_T[i] != in_T[j]) { w[i] += W[i][j]+W[j][i]; }
			}
			if (in_T[i] == YES) { sum_w += w[i]; }
		}
		// Compute cc[i]
		for (i=0; i<N; i++) {
			cc[i] = 0.0;			// How much it cost to move the information of a whole side to 'i'
			for (j=0; j<N; j++) {
				if (in_T[i] == YES && in_T[j] == YES) {   cc[i] += w[j]*dist_in_T[i][j]; }
				else if (in_T[i] == NO && in_T[j] == NO){ cc[i] += w[j]*dist_in_G[i][j]; }
			}
		}

		// Find the best edge to enlarge the tree:
		best_i = NONE;
		best_j = NONE;
		best_delta = 0.0;
		for (i=0; i<N; i++) {
			for (j=i+1; j<N; j++) {
				if (in_T[i] != in_T[j] && (info == NULL || info[i][j] != NO)) {	// if the edge is allowed.
					delta = cc[i] + C[i][j]*sum_w + cc[j];
					if ((info != NULL && info[i][j] == YES) || best_i == NONE || delta < best_delta) { // If the edge is mandatory or if the edge is the best one
						best_i = i;
						best_j = j;
						best_delta = delta;
					}
					if (info != NULL && info[i][j] == YES) { break; }			// Stop searching!
				}
			}
		}
		assert(best_i!=NONE && best_j!=NONE);

		// Add it to the tree:
		Y[best_i][best_j] = Y[best_j][best_i] = 1.0;

		if (in_T[best_i] == NO) {		// add best_i to the tree
			for (k=0; k<N; k++) {
				if (in_T[k] == YES) {
					dist_in_T[k][best_i] = dist_in_T[k][best_j] + C[best_j][best_i];
					dist_in_T[best_i][k] = dist_in_T[k][best_i];
				}
			}
			in_T[best_i] = YES; 
		} else {						// otherwise: add best_j to the tree
			for (k=0; k<N; k++) {
				if (in_T[k] == YES) {
					dist_in_T[k][best_j] = dist_in_T[k][best_i] + C[best_i][best_j];
					dist_in_T[best_j][k] = dist_in_T[k][best_j];
				}
			}
			in_T[best_j] = YES; 
		}
		n_in_T++;
	}

	// clean:
	free_double_matrix(&dist_in_G, N);
	free_double_matrix(&dist_in_T, N);
	free_int_vector(&in_T);
	free_double_vector(&cc);
	free_double_vector(&w);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////	// DONE but may improve?
// Constructive Heur. for the OCSTP based on the MinCut & MinCutMinCost matrices
void MCMCT(int N, double **CUT, double **COST, double **Y, int **info) {
	int i, j;
	double **MCMC = create_double_matrix(N, N);

	// Compute the MinCutMinCost matrix:
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { MCMC[i][j] = CUT[i][j]*COST[i][j]; } }

	// Compute the MST of MCMC:
	MST(N, MCMC, Y, info);

	free_double_matrix(&MCMC, N);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Divide and Conquer Heuristic for the OCSTP:
//
// Total time = max{ (n-1)*time(SPLIT), O(n^3), (improvements)*O(n^4) }  (if IMPROVE==YES)
//            = max{ (n-1)*time(SPLIT), O(n^3) }                         (if IMPROVE==NO ) 
//
// Total time += (improvements)*O(n^3)  if FINAL_IMPROVE==YES
//
// The basic structure of this heuristic is:
//
//		* Split the set of vertices into 2 subsets
//		* Recurse to solve each subset independently and get 2 sub-trees
//		* Merge both sub-trees into a single tree
//
// The Merging phase is always performed in a greedy fashion using a variation 
// of the Ahuja-Murty Local Search for the OCSTP. Therefore, the real work should
// be done during the Split phase.
// 
// Different Splitting heuristics will lead to different overall behaviour.
//
double OCSTP_DaC(int N, double **C, double **W, double **capacity, int SPLIT_TYPE, int IMPROVE, int FINAL_IMPROVE) {
	int i, j;
	int *V  = create_int_vector(N);
	int *S1 = create_int_vector(N);
	int *S2 = create_int_vector(N);
	int **T = create_int_matrix(N,N);
	double **Dist = create_double_matrix(N,N);
	double value;

	// Initialize:
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { T[i][j]    = NO; } }
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { Dist[i][j] = 0.0; } }
	for (i=0; i<N; i++) { V[i] = YES; }
	
	// Solve recursively:
	rec_OCSTP_DaC(N, C, W, T, Dist, V, SPLIT_TYPE, IMPROVE);

	// Try to improve as a post-processing step if it hasn't been already done:
	if (IMPROVE==NO && FINAL_IMPROVE == YES) { OCSTP_DaC_IMPROVE(N, C, W, T, Dist, V); }

	// Prepare Output:
	for (i=0; i<N; i++) { for (j=0; j<N; j++) {
		if (T[i][j] == YES) { capacity[i][j] = 1.0; }
		else                { capacity[i][j] = 0.0; }
	} }

	// Compute value:
	value = 0.0;
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { value += Dist[i][j] * W[i][j]; } }

	// Clean:
	free_int_vector(&V);
	free_int_vector(&S1);
	free_int_vector(&S2);
	free_int_matrix(&T, N);
	free_double_matrix(&Dist, N);

	// Return value:
	return value;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Recursive Call to the Divide and Conquer heuristic for the OCSTP
void rec_OCSTP_DaC(int N, double **C, double **W, int **T, double **Dist, int *V, int SPLIT_TYPE, int IMPROVE) {

	int i;
	int size = 0;
	int *S1;
	int *S2;
	int v[2];	

	// Compute the size of V:
	for (i=0; i<N; i++) { if (V[i] == YES) { size++; } }

	if (size == 2) {
		// This is a trivial case and it is not necessary to recourse:
		size = 0;
		for (i=0; i<N; i++) { if (V[i] == YES) { v[size++] = i; } }
		T[v[0]][v[1]] = YES;
		T[v[1]][v[0]] = YES;
		Dist[v[0]][v[1]] = C[v[0]][v[1]];
		Dist[v[1]][v[0]] = C[v[1]][v[0]];

	} else if (size > 2) {
		// This is the general case
		S1 = create_int_vector(N);
		S2 = create_int_vector(N);

		// SPLIT 
		if      (SPLIT_TYPE == 0)  { OCSTP_DaC_SPLIT_0(N, C, W, V, S1, S2); } // Even/Odd splitting
		else if (SPLIT_TYPE == 1)  { OCSTP_DaC_SPLIT_1(N, C, W, V, S1, S2); } // Single node splitting argmax(sum_{j in S} W[i][j])
		else if (SPLIT_TYPE == 2)  { OCSTP_DaC_SPLIT_2(N, C, W, V, S1, S2); } // Single node splitting argmax(sum_{j in V} W[i][j])
		else if (SPLIT_TYPE == 3)  { OCSTP_DaC_SPLIT_3(N, C, W, V, S1, S2); } // Single node splitting argmin(sum_{j in S} W[i][j])
		else if (SPLIT_TYPE == 4)  { OCSTP_DaC_SPLIT_4(N, C, W, V, S1, S2); } // Single node splitting argmin(sum_{j in V} W[i][j])
		else if (SPLIT_TYPE == 5)  { OCSTP_DaC_SPLIT_5(N, C, W, V, S1, S2); } // Heuristic Max Cut
		else if (SPLIT_TYPE == 6)  { OCSTP_DaC_SPLIT_6(N, C, W, V, S1, S2); } // Heuristic Clustering MaxEdge of MinimumSpanningTree(W)
		else if (SPLIT_TYPE == 7)  { OCSTP_DaC_SPLIT_7(N, C, W, V, S1, S2); } // Heuristic Clustering MaxEdge of MinimumSpanningTree(W)
		else if (SPLIT_TYPE == 8)  { OCSTP_DaC_SPLIT_8(N, C, W, V, S1, S2); } // Heuristic Clustering MaxEdge of MinimumSpanningTree(W)
		else if (SPLIT_TYPE == 9)  { OCSTP_DaC_SPLIT_9(N, C, W, V, S1, S2); } // Heuristic Clustering MaxEdge of MinimumSpanningTree(W)
		else if (SPLIT_TYPE == 10) { OCSTP_DaC_SPLIT_10(N, C, W, V, S1, S2); } // Heuristic Clustering MaxEdge of MinimumSpanningTree(W)
		else if (SPLIT_TYPE == 11) { OCSTP_DaC_SPLIT_11(N, C, W, V, S1, S2); } // Heuristic Clustering MaxEdge of MinimumSpanningTree(W)
		else                       { OCSTP_DaC_SPLIT_0(N, C, W, V, S1, S2); } // DEFAULT
		
		// RECOURSE
		rec_OCSTP_DaC(N, C, W, T, Dist, S1, SPLIT_TYPE, IMPROVE);
		rec_OCSTP_DaC(N, C, W, T, Dist, S2, SPLIT_TYPE, IMPROVE);

		// MERGE
		OCSTP_DaC_MERGE(N, C, W, T, Dist, S1, S2);

		// TRY TO IMPROVE
		if (IMPROVE == YES) { OCSTP_DaC_IMPROVE(N, C, W, T, Dist, V); }
	
		// Clean:
		free_int_vector(&S1);
		free_int_vector(&S2);
	}
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// - As long as the set S contains at least 2 elements this function
//   should return 2 non-empty sets S1 and S2.
// - It should take no more than O(n^2) time to do so.
//
// Divide in halves. Even/odd split.
void OCSTP_DaC_SPLIT_0(int N, double **C, double **W, int *S, int *S1, int *S2) {
	int i, count;
	for (i=0; i<N; i++) { S1[i] = NO; }
	for (i=0; i<N; i++) { S2[i] = NO; }
	count = 0;
	
	// Odd/Even Split:
	for (i=0; i<N; i++) {
		if (S[i] == YES) {
			if (count%2 == 0) { S1[i] = YES; }
			else              { S2[i] = YES; }
			count++;
		}
	}
}
// Split one node from the rest: The i = argmax(sum_{j in S} W[i][j])
void OCSTP_DaC_SPLIT_1(int N, double **C, double **W, int *S, int *S1, int *S2) {
	int i, j, best_i;
	double sum_W, best_sum_W;
	for (i=0; i<N; i++) { S1[i] = NO; }
	for (i=0; i<N; i++) { S2[i] = S[i]; }
	best_i = NONE;
	for (i=0; i<N; i++) {
		if (S[i] == YES) {
			sum_W = 0.0;
			for (j=0; j<N; j++) { if (S[j] == YES) { sum_W += W[i][j]; } }
			if (best_i == NONE || best_sum_W < sum_W) {
				best_i = i;
				best_sum_W = sum_W;
			}
		}
	}
	S2[best_i] = NO;
	S1[best_i] = YES;
}
// Split one node from the rest: The i = argmax(sum_{j in V} W[i][j])
void OCSTP_DaC_SPLIT_2(int N, double **C, double **W, int *S, int *S1, int *S2) {
	int i, j, best_i;
	double sum_W, best_sum_W;
	for (i=0; i<N; i++) { S1[i] = NO; }
	for (i=0; i<N; i++) { S2[i] = S[i]; }
	best_i = NONE;
	for (i=0; i<N; i++) {
		if (S[i] == YES) {
			sum_W = 0.0;
			for (j=0; j<N; j++) { sum_W += W[i][j]; }
			if (best_i == NONE || best_sum_W < sum_W) {
				best_i = i;
				best_sum_W = sum_W;
			}
		}
	}
	S2[best_i] = NO;
	S1[best_i] = YES;
}
// Split one node from the rest: The i = argmin(sum_{j in S} W[i][j])
void OCSTP_DaC_SPLIT_3(int N, double **C, double **W, int *S, int *S1, int *S2) {
	int i, j, best_i;
	double sum_W, best_sum_W;
	for (i=0; i<N; i++) { S1[i] = NO; }
	for (i=0; i<N; i++) { S2[i] = S[i]; }
	best_i = NONE;
	for (i=0; i<N; i++) {
		if (S[i] == YES) {
			sum_W = 0.0;
			for (j=0; j<N; j++) { if (S[j] == YES) { sum_W += W[i][j]; } }
			if (best_i == NONE || best_sum_W > sum_W) {
				best_i = i;
				best_sum_W = sum_W;
			}
		}
	}
	S2[best_i] = NO;
	S1[best_i] = YES;
}
// Split one node from the rest: The i = argmin(sum_{j in V} W[i][j])
void OCSTP_DaC_SPLIT_4(int N, double **C, double **W, int *S, int *S1, int *S2) {
	int i, j, best_i;
	double sum_W, best_sum_W;
	for (i=0; i<N; i++) { S1[i] = NO; }
	for (i=0; i<N; i++) { S2[i] = S[i]; }
	best_i = NONE;
	for (i=0; i<N; i++) {
		if (S[i] == YES) {
			sum_W = 0.0;
			for (j=0; j<N; j++) { sum_W += W[i][j]; }
			if (best_i == NONE || best_sum_W > sum_W) {
				best_i = i;
				best_sum_W = sum_W;
			}
		}
	}
	S2[best_i] = NO;
	S1[best_i] = YES;
}
// Heuristic max cut
void OCSTP_DaC_SPLIT_5(int N, double **C, double **W, int *S, int *S1, int *S2) {
	int i, j, k, best_i, best_j;
	double best_W;
	double *W1 = create_double_vector(N);
	double *W2 = create_double_vector(N);
	int size, size1, size2;

	// Compute set-size:
	size = 0;
	for (i=0; i<N; i++) { if (S[i] == YES) { size++; } }

	// Find the pair of nodes with the biggest communication request:
	best_W = -1.0;
	for (i=0; i<N; i++) {
		if (S[i] == YES) {
			for (j=0; j<N; j++) {
				if (S[j] == YES  && W[i][j] > best_W) {
					best_W = W[i][j];
					best_i = i;
					best_j = j;
				}
			}
		}
	}

	// Initialize S1 and S2:
	for (i=0; i<N; i++) { S1[i] = NO; }
	for (i=0; i<N; i++) { S2[i] = NO; }
	S1[best_i] = YES;
	S2[best_j] = YES;
	size1 = 1;
	size2 = 1;

	// Initialize W1 and W2:
	for (k=0; k<N; k++) { 
		if (S[k] == YES && S1[k] == NO && S2[k] == NO) { 
			W1[k] = W[best_i][k]; 
			W2[k] = W[best_j][k];
		} else {
			W1[k] = -1.0;
			W2[k] = -1.0;
		}
	}

	// Main LOOP:
	while (size1+size2 < size) {

		// Find the next candidate
		best_i = NONE;
		for (i=0; i<N; i++) {
			if (S[i] == YES && S1[i] == NO && S2[i] == NO) {
				if (best_i == NONE || best_W < MAX(W1[i], W2[i])) {
					best_i = i;
					best_W = MAX(W1[i], W2[i]);
				}
			}
		}

		// Include it to the corresponding subset:
		if (W1[best_i] > W2[best_i] || (W1[best_i]==W2[best_i] && size2>size1)) {
			S1[best_i] = YES;
			size1++;
			for (i=0; i<N; i++) { if (S[i]==YES && S1[i]==NO && S2[i]==NO) { W1[i] += W[best_i][i]; } }
		} else { 
			S2[best_i] = YES,
			size2++;
			for (i=0; i<N; i++) { if (S[i]==YES && S1[i]==NO && S2[i]==NO) { W2[i] += W[best_i][i]; } }
		}
		W1[best_i] = -1.0;
		W2[best_i] = -1.0;
	}

	// Clean:
	free_double_vector(&W1);
	free_double_vector(&W2);
}
// Heuristic Clustering: MaxEdge of MinSpanningTree(W)
void OCSTP_DaC_SPLIT_6(int N, double **C, double **W, int *S, int *S1, int *S2) {
	int		i, j, best_i;
	double	best_W;
	int		iteration, index_min;
	int		*stack = create_int_vector(N);
	int		*pred  = create_int_vector(N);
	double	*dist  = create_double_vector(N);
	int		stack_size;

	int size = 0;
	for (i=0; i<N; i++) { if (S[i] == YES) { size++; best_i=i; } }

	// Initialice everything:
	for (i=0; i<N; i++) { pred[i] = best_i; }
	for (i=0; i<N; i++) { dist[i] = W[best_i][i]; }
    S1[best_i] = YES;	// best_i is the only visited node so far
    
	// during size-1 iterations add a new edge to the MST
	for (iteration=1; iteration<size; iteration++) {
		// look for the closest node to the current tree:
		index_min = NONE;
		for (j=0; j<N; j++) {
			if (S[j]==YES && S1[j] == NO) {
				if (index_min == NONE || dist[index_min] > dist[j]) {
					index_min = j; 
				}
			}
		}
		assert(index_min != NONE);

		// We add it to the MST:
		S1[index_min] = YES;

		// and update the distances:
    	for (j=0; j<N; j++) {
    		if (S[j]==YES && S1[j] == NO ) {
				if (dist[j] > W[index_min][j]) {
    				dist[j] = W[index_min][j];
    				pred[j] = index_min;
				}
    		}
		}
	}

	// Now select the biggest edge and split the tree:
	best_i = NONE;
	for (i=0; i<N; i++) { 
		if (S1[i] == YES && (best_i==NONE || best_W < W[i][pred[i]])) {
			best_i = i;
			best_W = W[i][pred[i]];
		}
	}

	// Transfer to S2 the half containing best_i:
	S1[best_i] = NO;
	S2[best_i] = YES;
	stack_size = 0;
	stack[stack_size++] = best_i;
	while (stack_size>0) {
		i = stack[--stack_size];
		for (j=0; j<N; j++) {
			if (S1[j]==YES && (i==pred[j] || j==pred[i]) && j!= pred[best_i]) {
				stack[stack_size++] = j;
				S1[j] = NO;
				S2[j] = YES;
			}
		}
	}

	// Clean:
	free_int_vector(&pred);
	free_double_vector(&dist);
	free_int_vector(&stack);
}
// Heuristic Clustering: MaxEdge of MaxSpanningTree(W)
void OCSTP_DaC_SPLIT_7(int N, double **C, double **W, int *S, int *S1, int *S2) {
	int		i, j, best_i;
	double	best_W;
	int		iteration, index_min;
	int		*stack = create_int_vector(N);
	int		*pred  = create_int_vector(N);
	double	*dist  = create_double_vector(N);
	int		stack_size;

	int size = 0;
	for (i=0; i<N; i++) { if (S[i] == YES) { size++; best_i=i; } }

	// Initialice everything:
	for (i=0; i<N; i++) { pred[i] = best_i; }
	for (i=0; i<N; i++) { dist[i] = W[best_i][i]; }
    S1[best_i] = YES;	// best_i is the only visited node so far
    
	// during size-1 iterations add a new edge to the MST
	for (iteration=1; iteration<size; iteration++) {
		// look for the farthest node to the current tree:
		index_min = NONE;
		for (j=0; j<N; j++) {
			if (S[j]==YES && S1[j] == NO) {
				if (index_min == NONE || dist[index_min] < dist[j]) {
					index_min = j; 
				}
			}
		}
		assert(index_min != NONE);

		// We add it to the MST:
		S1[index_min] = YES;

		// and update the distances:
    	for (j=0; j<N; j++) {
    		if (S[j]==YES && S1[j] == NO ) {
				if (dist[j] < W[index_min][j]) {
    				dist[j] = W[index_min][j];
    				pred[j] = index_min;
				}
    		}
		}
	}

	// Now select the biggest edge and split the tree:
	best_i = NONE;
	for (i=0; i<N; i++) { 
		if (S1[i] == YES && (best_i==NONE || best_W < W[i][pred[i]])) {
			best_i = i;
			best_W = W[i][pred[i]];
		}
	}

	// Transfer to S2 the half containing best_i:
	S1[best_i] = NO;
	S2[best_i] = YES;
	stack_size = 0;
	stack[stack_size++] = best_i;
	while (stack_size>0) {
		i = stack[--stack_size];
		for (j=0; j<N; j++) {
			if (S1[j]==YES && (i==pred[j] || j==pred[i]) && j!= pred[best_i]) {
				stack[stack_size++] = j;
				S1[j] = NO;
				S2[j] = YES;
			}
		}
	}

	// Clean:
	free_int_vector(&pred);
	free_double_vector(&dist);
	free_int_vector(&stack);
}
// Heuristic Clustering: MinEdge of MinSpanningTree(W)
void OCSTP_DaC_SPLIT_8(int N, double **C, double **W, int *S, int *S1, int *S2) {
	int		i, j, best_i;
	double	best_W;
	int		iteration, index_min;
	int		*stack = create_int_vector(N);
	int		*pred  = create_int_vector(N);
	double	*dist  = create_double_vector(N);
	int		stack_size;

	int size = 0;
	for (i=0; i<N; i++) { if (S[i] == YES) { size++; best_i=i; } }

	// Initialice everything:
	for (i=0; i<N; i++) { pred[i] = best_i; }
	for (i=0; i<N; i++) { dist[i] = W[best_i][i]; }
    S1[best_i] = YES;	// best_i is the only visited node so far
    
	// during size-1 iterations add a new edge to the MST
	for (iteration=1; iteration<size; iteration++) {
		// look for the closest node to the current tree:
		index_min = NONE;
		for (j=0; j<N; j++) {
			if (S[j]==YES && S1[j] == NO) {
				if (index_min == NONE || dist[index_min] > dist[j]) {
					index_min = j; 
				}
			}
		}
		assert(index_min != NONE);

		// We add it to the MST:
		S1[index_min] = YES;

		// and update the distances:
    	for (j=0; j<N; j++) {
    		if (S[j]==YES && S1[j] == NO ) {
				if (dist[j] > W[index_min][j]) {
    				dist[j] = W[index_min][j];
    				pred[j] = index_min;
				}
    		}
		}
	}

	// Now select the biggest edge and split the tree:
	best_i = NONE;
	best_W = -1.0;
	for (i=0; i<N; i++) { 
		if (S1[i] == YES && (best_i==NONE || best_W > dist[i])) {
			best_i = i;
			best_W = dist[i];
		}
	}
	assert(best_i != NONE);

	// Transfer to S2 the half containing best_i:
	S1[best_i] = NO;
	S2[best_i] = YES;
	stack_size = 0;
	stack[stack_size++] = best_i;
	while (stack_size>0) {
		i = stack[--stack_size];
		for (j=0; j<N; j++) {
			if (S1[j]==YES && (i==pred[j] || j==pred[i]) && j!= pred[best_i]) {
				stack[stack_size++] = j;
				S1[j] = NO;
				S2[j] = YES;
			}
		}
	}

	// Clean:
	free_int_vector(&pred);
	free_double_vector(&dist);
	free_int_vector(&stack);
}
// Heuristic Clustering: MinEdge of MaxSpanningTree(W)
void OCSTP_DaC_SPLIT_9(int N, double **C, double **W, int *S, int *S1, int *S2) {
	int		i, j, best_i;
	double	best_W;
	int		iteration, index_min;
	int		*stack = create_int_vector(N);
	int		*pred  = create_int_vector(N);
	double	*dist  = create_double_vector(N);
	int		stack_size;

	int size = 0;
	for (i=0; i<N; i++) { if (S[i] == YES) { size++; best_i=i; } }

	// Initialice everything:
	for (i=0; i<N; i++) { pred[i] = best_i; }
	for (i=0; i<N; i++) { dist[i] = W[best_i][i]; }
    S1[best_i] = YES;	// best_i is the only visited node so far
    
	// during size-1 iterations add a new edge to the MST
	for (iteration=1; iteration<size; iteration++) {
		// look for the farthest node to the current tree:
		index_min = NONE;
		for (j=0; j<N; j++) {
			if (S[j]==YES && S1[j] == NO) {
				if (index_min == NONE || dist[index_min] < dist[j]) {
					index_min = j; 
				}
			}
		}
		assert(index_min != NONE);

		// We add it to the MST:
		S1[index_min] = YES;

		// and update the distances:
    	for (j=0; j<N; j++) {
    		if (S[j]==YES && S1[j] == NO ) {
				if (dist[j] < W[index_min][j]) {
    				dist[j] = W[index_min][j];
    				pred[j] = index_min;
				}
    		}
		}
	}

	// Now select the smallest edge and split the tree:
	best_i = NONE;
	for (i=0; i<N; i++) { 
		if (S1[i] == YES && (best_i==NONE || best_W > W[i][pred[i]])) {
			best_i = i;
			best_W = W[i][pred[i]];
		}
	}

	// Transfer to S2 the half containing best_i:
	S1[best_i] = NO;
	S2[best_i] = YES;
	stack_size = 0;
	stack[stack_size++] = best_i;
	while (stack_size>0) {
		i = stack[--stack_size];
		for (j=0; j<N; j++) {
			if (S1[j]==YES && (i==pred[j] || j==pred[i]) && j!= pred[best_i]) {
				stack[stack_size++] = j;
				S1[j] = NO;
				S2[j] = YES;
			}
		}
	}

	// Clean:
	free_int_vector(&pred);
	free_double_vector(&dist);
	free_int_vector(&stack);
}
// Divide in aproximate halves randomly.
void OCSTP_DaC_SPLIT_10(int N, double **C, double **W, int *S, int *S1, int *S2) {
	int i;
	int size1 = 0;
	int size2 = 0;
	for (i=0; i<N; i++) { S1[i] = NO; }
	for (i=0; i<N; i++) { S2[i] = NO; }
	
	// aproximate random split:
	for (i=0; i<N; i++) {
		if (S[i] == YES) {
			if      (size1 == 0)                      { S1[i] = YES; size1++; }
			else if (size2 == 0)                      { S2[i] = YES; size2++; }
			else if (randint(0,size1+size2) <= size2) { S1[i] = YES; size1++; }
			else                                      { S2[i] = YES; size2++; }
		}
	}
}
// Divide in exact halves randomly.
void OCSTP_DaC_SPLIT_11(int N, double **C, double **W, int *S, int *S1, int *S2) {
	int i, j, aux;
	int size = 0;
	int *perm = create_int_vector(N);
	
	// Init:
	for (i=0; i<N; i++) { S1[i] = NO; }
	for (i=0; i<N; i++) { S2[i] = NO; }
	for (i=0; i<N; i++) { if (S[i] == YES) { perm[size++] = i; } }

	// Fisher Yates Shuffle	
	for (i=0; i<size-1; i++) {
		j = randint(i,size-1);
		aux = perm[i]; perm[i] = perm[j]; perm[j] = aux;
	}

	// Store it properly
	for (i=0;      i<size/2; i++) { S1[perm[i]] = YES; }
	for (i=size/2; i<size;   i++) { S2[perm[i]] = YES; }
	
	// Clean
	free_int_vector(&perm);

}
// Divide By the whole Gomory Hu Tree:															// Must pre-compute all splits and pass them to the recursive call via C and W matrices
void OCSTP_DaC_SPLIT_12(int N, double **C, double **W, int *S, int *S1, int *S2) {

	int		i,j,k;						// index variables for loops
    int		s,t;						// source and sink nodes
    double	auxFlow;					// aux variable to make a switch
	int		*tree    = create_int_vector(N);
	int		*stack   = create_int_vector(N);
	int		stack_size, size_S, size_S1;
	int		e, best_e, best_dif;

	// compute size of S
	size_S = 0;
	for(i=0; i<N; i++) { if (S[i] == YES) { size_S++; } }
	assert(size_S >= 2);

    // COMPUTE Gomory-Hu Tree:
    for (i=0; i<N; i++) { tree[i] = 0; }
    for (s=1; s<N; s++) {
        t = tree[s];                                 
		auxFlow = maxFlowAndCut(N, s, t, W, stack);   
        for (i=0; i<N; i++) { if (stack[i] == YES  &&  i != s  &&  tree[i] == t) { tree[i] = s; } }
		if (stack[tree[t]] == YES) { tree[s] = tree[t]; tree[t] = s; }
    }

	// Try to remove each edge of the tree and rememeber the best one:
	best_e   = NONE;
	best_dif = NONE;
	for (e=1; e<N; e++) {
		if (S[e] == YES && S[tree[e]] == YES) {
			// Put the 'i' side into S1:
			for (i=0; i<N; i++) { S1[i] = NO; }
			size_S1	   = 1;
			S1[e]	   = YES;
			stack_size = 0;
			stack[stack_size++] = e;
			while (stack_size > 0) {
				i = stack[--stack_size];
				j = tree[i];
				// Case 1:
				if (S[j] == YES && S1[j] == NO) { 
					stack[stack_size++] = j;
					S1[j] = YES;
					size_S1++;
				}			
				// Case 2:
				for (k=1; k<N; k++) {
					if (tree[k] == i && S[k] == YES && S1[k] == NO) {
						stack[stack_size++] = k;
						S1[k] = YES;
						size_S1++;	
					}
				}
			}
			if (best_e = NONE || best_dif > ABS(2*size_S1 - size_S)) {
				best_e = e;
				best_dif = ABS(2*size_S1 - size_S);
			}
		}
	}
	assert(best_e != NONE);
	assert(best_dif <= (size_S-2));

	// Compute the final S1:
	for (i=0; i<N; i++) { S1[i] = NO; }
	S1[best_e] = YES;
	stack_size = 0;
	stack[stack_size++] = best_e;
	while (stack_size > 0) {
		i = stack[--stack_size];
		j = tree[i];
		// Case 1:
		if (S[j] == YES && S1[j] == NO) { 
			stack[stack_size++] = j;
			S1[j] = YES;
		}			
		// Case 2:
		for (k=1; k<N; k++) {
			if (tree[k] == i && S[k] == YES && S1[k] == NO) {
				stack[stack_size++] = k;
				S1[k] = YES;
			}
		}
	}

	// Compute the final S2:
	for (i=0; i<N; i++) { 
		if (S[i]==YES && S1[i] == NO) { S2[i] = YES; }
		else                          { S2[i] = NO;  } 
	}

	// clean:
	free_int_vector(&stack);
	free_int_vector(&tree);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////	// Must be improved
// merge S1 and S2 and update the T and D matrices
void OCSTP_DaC_MERGE(int N, double **C, double **W, int **T, double **Dist, int *S1, int *S2) {

	int i, j;
	int *X = create_int_vector(N);
	int *Y = create_int_vector(N);
	int size1 = 0;
	int size2 = 0;
	int x, y, best_x, best_y;
	double value, best_value;
	double *w = create_double_vector(N);
	double *cc = create_double_vector(N);
	double sum_w;
		
	// Get the list of nodes of each subset and their sizes:
	for (i=0; i<N; i++) {
		if (S1[i] == YES) { X[size1++] = i; }
		if (S2[i] == YES) { Y[size2++] = i; }
	}
	assert(size1>0);
	assert(size2>0);

	// PRECOMPUTE //////////////////////////////////////////////////////////////
	// Now compute the amount of communication 'w[x]' from each node 'x'
	// to the other side of the tree:
	for (x=0; x<N; x++) { w[x] = 0.0; }

	for (i=0; i<size1; i++) {
		x = X[i];
		for (j=0; j<size2; j++) {
			y = Y[j];
			w[x] += W[x][y] + W[y][x];
			w[y] += W[x][y] + W[y][x];
		}
	}
	// And the total exchange of communication between both sides:
	sum_w = 0.0;
	for (i=0; i<size1; i++) { sum_w += w[X[i]]; }
	// And the communication cost of moving all the info to the exit point of each side:
	for (x=0; x<N; x++) { cc[x] = 0.0; }
	for (i=0; i<size1-1; i++) {	for (j=i+1; j<size1; j++) {
			cc[X[i]] += w[X[j]]*Dist[X[i]][X[j]];
			cc[X[j]] += w[X[i]]*Dist[X[j]][X[i]];
	}	}
	for (i=0; i<size2-1; i++) { for (j=i+1; j<size2; j++) {
			cc[Y[i]] += w[Y[j]]*Dist[Y[i]][Y[j]];
			cc[Y[j]] += w[Y[i]]*Dist[Y[j]][Y[i]];
	}	}
	////////////////////////////////////////////////////////////////////////////

	best_x = NONE;
	best_y = NONE;
	// Look for the best edge joining S1 and S2:
	for (i=0; i<size1; i++) {
		x = X[i];
		for (j=0; j<size2; j++) {
			y = Y[j];
			// Get value and remember the best one!
			value = cc[x] + C[x][y]*sum_w + cc[y];
			if (best_x == NONE || value < best_value) {
				best_value = value;
				best_x = x;
				best_y = y;
			}
		}
	}
	assert(best_x != NONE);
	assert(best_y != NONE);

	// Add the edge and update everything:
	T[best_x][best_y] = YES;
	T[best_y][best_x] = YES;
	for (i=0; i<size1; i++) {
		x = X[i];
		for (j=0; j<size2; j++) {
			y = Y[j];
			Dist[x][y] = Dist[x][best_x] + C[best_x][best_y] + Dist[best_y][y];
			Dist[y][x] = Dist[x][y];
		}
	}

	// Clean
	free_int_vector(&X);
	free_int_vector(&Y);
	free_double_vector(&w);
	free_double_vector(&cc);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Applies a Ahuja-Murty local search to the sub-tree formed by the nodes in 
// in the subset 'S' and then updates 'T' and 'D'.
// Assumes that the nodes in 'S' form a tree in 'T'
void OCSTP_DaC_IMPROVE(int N, double **C, double **W, int **T, double **Dist, int *S) {

	int		i, j, x, y, edge, best_x, best_y;
	double	sum_w;
	int		*index_x = create_int_vector(N-1);
	int		*index_y = create_int_vector(N-1);
	int		*X = create_int_vector(N);
	int		*Y = create_int_vector(N);
	double	*w  = create_double_vector(N);
	double	*cc = create_double_vector(N);
	int		*stack = create_int_vector(N);
	int		stack_size, size, size_x, size_y;
	double	improvement, best_improvement, cost, best_cost;
	int		old_edge, new_x, new_y;
	
	// Compute size:
	size = 0;
	for (x=0; x<N; x++) { if (S[x] == YES) { size++; } }

	// Store the edges in a list
	edge = 0;
	for (x=0; x<N; x++) {
		if (S[x] == YES) {
			for (y=x+1; y<N; y++) {
				if (S[y]==YES && T[x][y] == YES) {
					index_x[edge] = x;
					index_y[edge] = y;
					edge++;
				}
			}
		}
	}
	assert(edge == size-1);
	
	/// MAIN LOOP //////////////////////////////////////////////////////////////
	best_improvement = 1.0;			// This is just for entering in the main loop
	while (best_improvement > EPSILON) {


		for (edge=0; edge<size-1; edge++) { assert(T[index_x[edge]][index_y[edge]] == YES); }

		// Examine all the edges to find the best one
		old_edge = NONE;
		new_x    = NONE;
		new_y    = NONE;
		best_improvement = -1.0;
		for (edge=0; edge<size-1; edge++) {
			// we are going to "remove" the edge (x,y):
			x = index_x[edge];
			y = index_y[edge];

			// Set X contains the edges in the 'x' side of the T\{(x,y)}
			for (i=0; i<N; i++) { X[i] = NO; }
			// We retrieve X by means of a DFS:
			X[x]       = YES;
			size_x     = 1;
			stack_size = 0;
			stack[stack_size++] = x;
			while (stack_size>0) {
				i = stack[--stack_size];
				for (j=0; j<N; j++) {
					if (j!=y && T[i][j]==YES && X[j]==NO) {
						stack[stack_size++] = j;
						X[j] = YES;
						size_x++;
					}
				}
			}

			// Set Y contains the edges in the 'y' side of the T\{(x,y)}
			size_y = 0;
			for (i=0; i<N; i++) { 
				if (S[i] == YES && X[i] == NO) { Y[i] = YES; size_y++; }
				else                           { Y[i] = NO;            }
			}			
			assert(size_x+size_y == size);

			// Now make them compact!
			size_x = 0;
			for (i=0; i<N; i++) { if (X[i] == YES) { X[size_x++] = i; } }
			size_y = 0;
			for (i=0; i<N; i++) { if (Y[i] == YES) { Y[size_y++] = i; } }
			assert(size_x+size_y == size);

			// Now compute the amount of communication 'w[x]' from each node 'x'
			// to the other side of the tree:
			for (i=0; i<N; i++) { w[i] = 0.0; }
			for (i=0; i<size_x; i++) {
				for (j=0; j<size_y; j++) {
					w[X[i]] += W[X[i]][Y[j]] + W[Y[j]][X[i]];
					w[Y[j]] += W[X[i]][Y[j]] + W[Y[j]][X[i]];
				}
			}
			// And the total exchange of communication between both sides:
			sum_w = 0.0;
			for (i=0; i<size_x; i++) { sum_w += w[X[i]]; }
			// And the communication cost of moving all the info to the exit point of each side:
			for (i=0; i<N; i++) { cc[i] = 0.0; }
			for (i=0; i<size_x-1; i++) {	for (j=i+1; j<size_x; j++) {
					cc[X[i]] += w[X[j]]*Dist[X[i]][X[j]];
					cc[X[j]] += w[X[i]]*Dist[X[j]][X[i]];
			}	}
			for (i=0; i<size_y-1; i++) { for (j=i+1; j<size_y; j++) {
					cc[Y[i]] += w[Y[j]]*Dist[Y[i]][Y[j]];
					cc[Y[j]] += w[Y[i]]*Dist[Y[j]][Y[i]];
			}	}


			// now find the edge that minimices the CC:
			best_cost = cc[x] + C[x][y]*sum_w + cc[y];
			best_x = x;
			best_y = y;
			for (i=0; i<size_x; i++) {
				for (j=0; j<size_y; j++) {
					cost = cc[X[i]] + C[X[i]][Y[j]]*sum_w + cc[Y[j]];
					if (cost < best_cost) {
						best_cost = cost;
						best_x = X[i];
						best_y = Y[j];
					}
				}
			}

			// If we have improved, remember how!
			if (best_x != x || best_y != y) {
				improvement = (cc[x] + C[x][y]*sum_w + cc[y]) - best_cost;
				if (improvement > best_improvement) {
					best_improvement = improvement;
					old_edge = edge;
					new_x = best_x;
					new_y = best_y;
				}
			} 
		}

		if (best_improvement > EPSILON) {

			x = index_x[old_edge];
			y = index_y[old_edge];

			// Get X and Y:

			// Set X contains the edges in the 'x' side of the T\{(x,y)}
			for (i=0; i<N; i++) { X[i] = NO; }
			// We retrieve X by means of a DFS:
			X[x]       = YES;
			size_x     = 1;
			stack_size = 0;
			stack[stack_size++] = x;
			while (stack_size>0) {
				i = stack[--stack_size];
				for (j=0; j<N; j++) {
					if (j!=y && T[i][j]==YES && X[j]==NO) {
						stack[stack_size++] = j;
						X[j] = YES;
						size_x++;
					}
				}
			}

			// Set Y contains the edges in the 'y' side of the T\{(x,y)}
			size_y = 0;
			for (i=0; i<N; i++) { 
				if (S[i] == YES && X[i] == NO) { Y[i] = YES; size_y++; }
				else                           { Y[i] = NO;            }
			}			
			assert(size_x+size_y == size);

			// Now make them compact!
			size_x = 0;
			for (i=0; i<N; i++) { if (X[i] == YES) { X[size_x++] = i; } }
			size_y = 0;
			for (i=0; i<N; i++) { if (Y[i] == YES) { Y[size_y++] = i; } }
			assert(size_x+size_y == size);

			// Update Distance matrix!
			for (i=0; i<size_x; i++) {
				for (j=0; j<size_y; j++) {
					Dist[X[i]][Y[j]] = Dist[X[i]][new_x] + C[new_x][new_y] + Dist[new_y][Y[j]];
					Dist[Y[j]][X[i]] = Dist[X[i]][Y[j]];
				}
			}

			// remove edge (x,y) from T:
			T[x][y] = NO;
			T[y][x] = NO;

			// add edge (new_x,new_y) to T:
			T[new_x][new_y] = YES;
			T[new_y][new_x] = YES;

			// Update edge list:
			index_x[old_edge] = new_x;
			index_y[old_edge] = new_y;
		}	
	}

	// Clean:
	free_int_vector(&index_x);
	free_int_vector(&index_y);
	free_int_vector(&X);
	free_int_vector(&Y);
	free_double_vector(&w);
	free_double_vector(&cc);
	free_int_vector(&stack);
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Bottom Up version of the Divide & Conquer Heuristic:
double OCSTP_BUDaC(int N, double **C, double **W, double **capacity, int IMPROVE, int FINAL_IMPROVE) {
	int i, j, n_cc, c_1, c_2;
	int *cc = create_int_vector(N);
	int *S1 = create_int_vector(N);
	int *S2 = create_int_vector(N);
	int *next = create_int_vector(N);

	double	**comm = create_double_matrix(N,N);
	int		**T = create_int_matrix(N,N);
	double	**Dist = create_double_matrix(N,N);
	double	value, max_comm;

	// Initialize the tree:
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { T[i][j]    = NO; } }
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { Dist[i][j] = 0.0; } }
	
	// Each vertex starts in a different connected component
	for (i=0; i<N; i++) { cc[i] = i; }
	// This holds the communication between components
	for (i=0; i<N; i++) { for (j=i+1; j<N; j++) { comm[i][j] = comm[j][i] = W[i][j]; } }
	// This holds the closer component:
	for (i=0; i<N; i++) {
		next[i] = NONE;
		for (j=0; j<N; j++) { if (next[i] == NONE || comm[i][next[i]] < comm[i][j]) { next[i] = j; } }
	}
	
	for (n_cc=N; n_cc>1; n_cc--) {

		// find the next 2 components to get merged in O(n)
		max_comm = -1.0;
		c_1 = c_2 = -1;
		for (i=0; i<N; i++) {
			if (cc[i] == i) {
				if (max_comm < comm[i][next[i]]) {
					max_comm = comm[i][next[i]];
					c_1 = i; 
					c_2 = next[i];
				}
			}
		}
		assert(c_1 >= 0);
		
		// build the explicit S_1 and S_2 in O(n)
		for (i=0; i<N; i++) {
			if (cc[i] == c_1) { S1[i] = YES; }
			else              { S1[i] = NO; }
			if (cc[i] == c_2) { S2[i] = YES; }
			else              { S2[i] = NO; }
		}

		// Merge them:
		OCSTP_DaC_MERGE(N, C, W, T, Dist, S1, S2);

		// TRY TO IMPROVE
		if (IMPROVE == YES) {
			for (i=0; i<N; i++) { S1[i] += S2[i]; }
			OCSTP_DaC_IMPROVE(N, C, W, T, Dist, S1);
		}	

		// Update CC:
		for (i=0; i<N; i++) { if (cc[i] == c_2) { cc[i] = c_1; } }
		// Update COMM:
		for (i=0; i<N; i++) {
			if (cc[i] == i) {
				comm[i][c_1] += comm[i][c_2];
				comm[c_1][i] += comm[c_2][i];
				comm[i][c_2] = 0.0;
				comm[c_2][i] = 0.0;
			}
			comm[i][i] = 0.0;
		}
		// Update NEXT:
		for (i=0; i<N; i++) {
			if (cc[i] == i) {
				if (comm[i][c_1] > comm[i][next[i]]){
					next[i] = c_1;
				}
			}
		}
	}
	
	// Try to improve as a post-processing step if it hasn't been already done:
	if (IMPROVE==NO && FINAL_IMPROVE == YES) { 
		for (i=0; i<N; i++) { S1[i] = YES; }
		OCSTP_DaC_IMPROVE(N, C, W, T, Dist, S1);
	}

	// Prepare Output:
	for (i=0; i<N; i++) { for (j=0; j<N; j++) {
		if (T[i][j] == YES) { capacity[i][j] = 1.0; }
		else                { capacity[i][j] = 0.0; }
	} }

	// Compute value:
	value = 0.0;
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { value += Dist[i][j] * W[i][j]; } }

	// Clean:
	free_int_vector(&cc);
	free_int_vector(&next);
	free_int_vector(&S1);
	free_int_vector(&S2);
	free_int_matrix(&T, N);
	free_double_matrix(&Dist, N);
	free_double_matrix(&comm, N);

	// Return value:
	return value;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////
// IMPROVEMENT HEURISTICS //
////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// - Given a cost matrix C, a communication requests matrix W, the adjacency 
// matrix of a tree T and its corresponding distance matrix Dist, this function
// returns the result of applying Ahuja-Murty Local Search to the tree.
// - If the 'info' matrix is provided (!=NULL), the tree will be enforced to
// follow the constraints of that matrix regarding mandatory edges (=YES),
// forbiden edges (=NO) or otherwise optional edges (=NONE).
// - If no info matrix is provided, the result will always be non-negative,
// with improvement > 0 if some improvement has been made and improvement = 0
// if the tree is left unchanged.
// - If the info matrix is provided, negative improvement can be returned as a
// result of the sanitation process that removes forbiden edges and adds 
// mandatory edges without enforcing and improvement in the objective function.
double AM_LS(int N, double **C, double **W, double **T, double **Dist, int **info) {

	int		i, j, x, y, edge, best_x, best_y;
	double	sum_w;
	int		*index_x = create_int_vector(N-1);
	int		*index_y = create_int_vector(N-1);
	int		*X = create_int_vector(N);
	int		*Y = create_int_vector(N);
	double	*w  = create_double_vector(N);
	double	*cc = create_double_vector(N);
	int		*stack = create_int_vector(N);
	int		stack_size, size_x, size_y;
	double	total_improvement, improvement, best_improvement, cost, best_cost;
	int		old_edge, new_x, new_y;
	
	// Store the edges in a list
	edge = 0;
	for (x=0; x<N; x++) {
		for (y=x+1; y<N; y++) {
			if (T[x][y] > 0.5) {
				index_x[edge] = x;
				index_y[edge] = y;
				edge++;
			}
		}
	}
	assert(edge == N-1);
	
	/// MAIN LOOP //////////////////////////////////////////////////////////////
	total_improvement = 0.0;
	best_improvement = 1.0;			// This is just for entering in the main loop
	while (best_improvement > EPSILON) {

		// Examine all the edges to find the best one
		old_edge = NONE;
		new_x    = NONE;
		new_y    = NONE;
		best_improvement = -1.0;
		for (edge=0; edge<N-1; edge++) {
			// we are going to "remove" the edge (x,y):
			x = index_x[edge];
			y = index_y[edge];

			if (info != NULL && info[x][y] == YES) { continue; } // Mandatory edge, examine the next one!

			// Set X contains the edges in the 'x' side of the T\{(x,y)}
			for (i=0; i<N; i++) { X[i] = NO; }
			// We retrieve X by means of a DFS:
			X[x]       = YES;
			size_x     = 1;
			stack_size = 0;
			stack[stack_size++] = x;
			while (stack_size>0) {
				i = stack[--stack_size];
				for (j=0; j<N; j++) {
					if (j!=y && T[i][j]>0.5 && X[j]==NO) {
						stack[stack_size++] = j;
						X[j] = YES;
						size_x++;
					}
				}
			}

			// Set Y contains the edges in the 'y' side of the T\{(x,y)}
			size_y = 0;
			for (i=0; i<N; i++) { 
				if (X[i] == NO) { Y[i] = YES; size_y++; }
				else            { Y[i] = NO;            }
			}			

			// Now make them compact!
			size_x = 0;
			for (i=0; i<N; i++) { if (X[i] == YES) { X[size_x++] = i; } }
			size_y = 0;
			for (i=0; i<N; i++) { if (Y[i] == YES) { Y[size_y++] = i; } }
			assert(size_x+size_y == N);

			// Now compute the amount of communication 'w[x]' from each node 'x'
			// to the other side of the tree:
			for (i=0; i<N; i++) { w[i] = 0.0; }
			for (i=0; i<size_x; i++) {
				for (j=0; j<size_y; j++) {
					w[X[i]] += W[X[i]][Y[j]] + W[Y[j]][X[i]];
					w[Y[j]] += W[X[i]][Y[j]] + W[Y[j]][X[i]];
				}
			}
			// And the total exchange of communication between both sides:
			sum_w = 0.0;
			for (i=0; i<size_x; i++) { sum_w += w[X[i]]; }
			// And the communication cost of moving all the info to the exit point of each side:
			for (i=0; i<N; i++) { cc[i] = 0.0; }
			for (i=0; i<size_x-1; i++) {	for (j=i+1; j<size_x; j++) {
					cc[X[i]] += w[X[j]]*Dist[X[i]][X[j]];
					cc[X[j]] += w[X[i]]*Dist[X[j]][X[i]];
			}	}
			for (i=0; i<size_y-1; i++) { for (j=i+1; j<size_y; j++) {
					cc[Y[i]] += w[Y[j]]*Dist[Y[i]][Y[j]];
					cc[Y[j]] += w[Y[i]]*Dist[Y[j]][Y[i]];
			}	}


			// now find the feasible edge that minimices the CC:
			best_cost = cc[x] + C[x][y]*sum_w + cc[y];
			best_x = x;
			best_y = y;
			for (i=0; i<size_x; i++) {
				for (j=0; j<size_y; j++) {
					cost = cc[X[i]] + C[X[i]][Y[j]]*sum_w + cc[Y[j]];
					// Easy case, without info:
					if (info == NULL) {
						if (cost < best_cost) {
							best_cost = cost;
							best_x = X[i];
							best_y = Y[j];
						}
					// Difficult case, taking info into account:
					} else if (info[best_x][best_y] == NO || info[X[i]][Y[j]] == YES || (info[X[i]][Y[j]] != NO && cost < best_cost)) {
						best_cost = cost;
						best_x = X[i];
						best_y = Y[j];
					} 
					if (info != NULL && info[best_x][best_y] == YES) { break; }
				}
				if (info != NULL && info[best_x][best_y] == YES) { break; }
			}

			// If we have improved, remember how!
			if (best_x != x || best_y != y) {
				improvement = (cc[x] + C[x][y]*sum_w + cc[y]) - best_cost;
				// Easy case, no info:
				if (info == NULL) {
					if (improvement > best_improvement) {
						best_improvement = improvement;
						old_edge = edge;
						new_x = best_x;
						new_y = best_y;
					}
				// Difficult case, taking info into account:
				} else if (improvement > best_improvement) {
					best_improvement = improvement;
					old_edge = edge;
					new_x = best_x;
					new_y = best_y;
				}
			} 
		}

		if (best_improvement > EPSILON) {

			x = index_x[old_edge];
			y = index_y[old_edge];

			// Get X and Y:

			// Set X contains the edges in the 'x' side of the T\{(x,y)}
			for (i=0; i<N; i++) { X[i] = NO; }
			// We retrieve X by means of a DFS:
			X[x]       = YES;
			size_x     = 1;
			stack_size = 0;
			stack[stack_size++] = x;
			while (stack_size>0) {
				i = stack[--stack_size];
				for (j=0; j<N; j++) {
					if (j!=y && T[i][j]>0.5 && X[j]==NO) {
						stack[stack_size++] = j;
						X[j] = YES;
						size_x++;
					}
				}
			}

			// Set Y contains the edges in the 'y' side of the T\{(x,y)}
			size_y = 0;
			for (i=0; i<N; i++) { 
				if (X[i] == NO) { Y[i] = YES; size_y++; }
				else            { Y[i] = NO;            }
			}			

			// Now make them compact!
			size_x = 0;
			for (i=0; i<N; i++) { if (X[i] == YES) { X[size_x++] = i; } }
			size_y = 0;
			for (i=0; i<N; i++) { if (Y[i] == YES) { Y[size_y++] = i; } }
			assert(size_x+size_y == N);

			// Update Distance matrix!
			for (i=0; i<size_x; i++) {
				for (j=0; j<size_y; j++) {
					Dist[X[i]][Y[j]] = Dist[X[i]][new_x] + C[new_x][new_y] + Dist[new_y][Y[j]];
					Dist[Y[j]][X[i]] = Dist[X[i]][Y[j]];
				}
			}

			// remove edge (x,y) from T:
			T[x][y] = 0.0;
			T[y][x] = 0.0;

			// add edge (new_x,new_y) to T:
			T[new_x][new_y] = 1.0;
			T[new_y][new_x] = 1.0;

			// Update edge list:
			index_x[old_edge] = new_x;
			index_y[old_edge] = new_y;
		}	
	}

	// Clean:
	free_int_vector(&index_x);
	free_int_vector(&index_y);
	free_int_vector(&X);
	free_int_vector(&Y);
	free_double_vector(&w);
	free_double_vector(&cc);
	free_int_vector(&stack);

	return total_improvement;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Ahuja Murty Local Search: returns improvement >0.0 if found, 0.0 otherwise.
double AM_LS2(GLOBAL_INFO global, int **info, int BEST_IMPROVEMENT) {
	
	int		i, j, e, u, v, d, edge, found, min_i, min_j;
	int		N = global.data->N;
	double  **C = global.data->C;
	double  **W = global.data->W;
	int		best_edge, best_new_i, best_new_j;
	double	delta = 0.0;
	double	best_delta, sum_w, current_cost, best_cost;

	// Efficient ways to work with trees:
	int		**neighbours = create_int_matrix(N,N-1);
	int		*degree = create_int_vector(N);
	int		*index_i = create_int_vector(N-1);
	int		*index_j = create_int_vector(N-1);
	
	int		*S = create_int_vector(N);
	int		*stack = create_int_vector(N);
	int		stack_size;
	
	double	**dist = create_double_matrix(N,N);
	double	*w  = create_double_vector(N);
	double	*cc = create_double_vector(N);
	
	// Copy Y to neighbours/degree data structure:
	for (i=0; i<N; i++) { degree[i] = 0; }
	for (i=0; i<N-1; i++) { 
		for (j=i+1; j<N; j++) {
			if (global.var->master[global.data->index_e[i][j]] > 0.5) {
				neighbours[i][degree[i]++] = j;
				neighbours[j][degree[j]++] = i;
			}
		}
	}
	for (i=0; i<N; i++) { for (d=degree[i]; d<N-1; d++) { neighbours[i][d] = NONE; } }
	
	best_delta = 1.0;				// Just to enter the loop;
	while (best_delta > EPSILON) {	// While we are improving...

		// Update current edge list:
		edge = 0;
		for (i=0; i<N-1; i++) {
			for (d=0; d<degree[i]; d++) {
				if (i < neighbours[i][d]) {
					index_i[edge] = i;
					index_j[edge] = neighbours[i][d];
					edge++;
				}
			}
		}				
		assert(edge == N-1);
		
		// Compute the distances over the current tree
		for (i=0; i<N; i++) {
			for (j=0; j<N; j++) { dist[i][j] = -1.0; }
			dist[i][i] = 0.0;
			// use a DFS:
			stack_size = 0;
			stack[stack_size++] = i;
			while (stack_size > 0) {
				u = stack[--stack_size];
				for (d=0; d<degree[u]; d++) {
					j = neighbours[u][d];
					if (dist[i][j] < -0.5) {
						dist[i][j] = dist[i][u] + C[u][j];
						stack[stack_size++] = j;
					}
				}
			}
		}

		// Examine all the edges to find the best one
		best_edge  = NONE;
		best_new_i = NONE;
		best_new_j = NONE;
		best_delta = 0.0;
		for (edge=0; edge<N-1; edge++) {
			
			// we are going to "remove" the edge (i,j):
			i = index_i[edge];
			j = index_j[edge];

			if (info != NULL && info[i][j]==YES) { continue; } // This edge can not be removed

			// Set S contains the edges in the 'i' side of the T\{(i,j)}
			for (u=0; u<N; u++) { S[u] = NO; }
			// We retrieve S by means of a DFS:
			S[i] = YES;
			stack_size = 0;
			stack[stack_size++] = i;
			while (stack_size>0){
				u = stack[--stack_size];
				for (d=0; d<degree[u]; d++) {
					if (neighbours[u][d]!=j && S[neighbours[u][d]]==NO) {
						stack[stack_size++] = neighbours[u][d];
						S[neighbours[u][d]] = YES;
					}
				}
			}

			// Now compute the amount of communication 'w[u]' from each node 'u'
			// to the other side of the tree:
			for (u=0; u<N; u++) {
				w[u] = 0;
				for (v=0; v<N; v++) {
					if (S[u] != S[v]) {
						w[u] += W[u][v] + W[v][u];
					}
				}
			}
			sum_w = 0;
			for (u=0; u<N; u++) {
				if (S[u] == YES) { sum_w += w[u]; }
			}

			// And the communication cost from nodes of the same side:
			for (u=0; u<N; u++) {
				cc[u] = 0;
				for (v=0; v<N; v++) {
					if (u!=v && S[u] == S[v]) {
						cc[u] += w[v] * dist[u][v];
					}
				}
			}

			// now find the edge that minimices the CC:
			best_cost = cc[i] + C[i][j]*sum_w + cc[j];
			min_i = i;
			min_j = j;
			for (u=0; u<N; u++) {
				for (v=u+1; v<N; v++) {
					if (S[u] != S[v] && (info == NULL || info[u][v] != NO)) {
						current_cost = cc[u] + C[u][v]*sum_w + cc[v];
						if ((info != NULL && info[i][j] == YES) || current_cost < best_cost) {
							current_cost = current_cost;
							min_i = u;
							min_j = v;
						}
						if (info != NULL && info[u][v] == YES) { break; }	//
					}
				}
			}

			// If we have improved, remember how!
			if (min_i != i || min_j != j) {
				if (((cc[i]+C[i][j]*sum_w+cc[j])-best_cost) > best_delta) {
					best_delta = (cc[i] + C[i][j]*sum_w + cc[j]) - best_cost;
					best_edge = edge;
					best_new_i = min_i;
					best_new_j = min_j;
				}
				if (BEST_IMPROVEMENT == NO) { break; }
			} 

		}

		if (best_delta > 0.0) {

			i = index_i[best_edge];
			j = index_j[best_edge];

			// remove j from new_neighbours[i]:
			found = NO;
			degree[i]--;
			for (d=0; d<degree[i]; d++) {
				if (found == NO && neighbours[i][d] == j) { found = YES; }
				if (found == YES) { neighbours[i][d] = neighbours[i][d+1]; }
			}
			neighbours[i][degree[i]] = NONE;

			// remove i from new_neighbours[j]:
			found = NO;
			degree[j]--;
			for (d=0; d<degree[j]; d++) {
				if (found == NO && neighbours[j][d] == i) { found = YES; }
				if (found == YES) { neighbours[j][d] = neighbours[j][d+1]; }
			}
			neighbours[j][degree[j]] = NONE;

			// add best_new_j to new_neighbours[best_new_i]:
			neighbours[best_new_i][degree[best_new_i]++] = best_new_j;

			// add best_new_i to new_neighbours[best_new_j]:
			neighbours[best_new_j][degree[best_new_j]++] = best_new_i;

			// add current delta to the global delta: 
			delta += best_delta;
		}	
	}

	// Copy best_solution to Y:
	for (e=0; e<global.data->E; e++) { global.var->master[e] = 0.0; }
	for (i=0; i<N; i++) {
		for (d=0; d<degree[i]; d++) {
			global.var->master[global.data->index_e[i][neighbours[i][d]]] = 1.0;
		}
	}

	// clean:
	free_int_matrix(&neighbours, N);
	free_int_vector(&degree);
	free_int_vector(&index_i);
	free_int_vector(&index_j);
	free_int_vector(&S);
	free_int_vector(&stack);
	free_double_matrix(&dist, N);
	free_double_vector(&w);
	free_double_vector(&cc);

	return delta;
}
double AM_LS3(int N, double **C, double **W, double **Y, int **info, int BEST_IMPROVEMENT) {
	
	int		i, j, u, v, d, edge, found, min_i, min_j;
	int		best_edge, best_new_i, best_new_j;
	double	delta = 0.0;
	double	best_delta, sum_w, current_cost, best_cost;

	// Efficient ways to work with trees:
	int		**neighbours = create_int_matrix(N,N-1);
	int		*degree = create_int_vector(N);
	int		*index_i = create_int_vector(N-1);
	int		*index_j = create_int_vector(N-1);
	
	int		*S = create_int_vector(N);
	int		*stack = create_int_vector(N);
	int		stack_size;
	
	double	**dist = create_double_matrix(N,N);
	double	*w  = create_double_vector(N);
	double	*cc = create_double_vector(N);
	
	// Copy Y to neighbours/degree data structure:
	for (i=0; i<N; i++) { degree[i] = 0; }
	for (i=0; i<N-1; i++) { 
		for (j=i+1; j<N; j++) {
			if (Y[i][j]+Y[j][i] >= 1.0-EPSILON) {
				neighbours[i][degree[i]++] = j;
				neighbours[j][degree[j]++] = i;
			}
		}
	}
	for (i=0; i<N; i++) { for (d=degree[i]; d<N-1; d++) { neighbours[i][d] = NONE; } }
	
	best_delta = 1.0;				// Just to enter the loop;
	while (best_delta > EPSILON) {	// While we are improving...

		// Update current edge list:
		edge = 0;
		for (i=0; i<N-1; i++) {
			for (d=0; d<degree[i]; d++) {
				if (i < neighbours[i][d]) {
					index_i[edge] = i;
					index_j[edge] = neighbours[i][d];
					edge++;
				}
			}
		}				
		assert(edge == N-1);
		
		// Compute the distances over the current tree
		for (i=0; i<N; i++) {
			for (j=0; j<N; j++) { dist[i][j] = -1.0; }
			dist[i][i] = 0.0;
			// use a DFS:
			stack_size = 0;
			stack[stack_size++] = i;
			while (stack_size > 0) {
				u = stack[--stack_size];
				for (d=0; d<degree[u]; d++) {
					j = neighbours[u][d];
					if (dist[i][j] < -0.5) {
						dist[i][j] = dist[i][u] + C[u][j];
						stack[stack_size++] = j;
					}
				}
			}
		}

		// Examine all the edges to find the best one
		best_edge  = NONE;
		best_new_i = NONE;
		best_new_j = NONE;
		best_delta = -1.0;
		for (edge=0; edge<N-1; edge++) {
			
			// we are going to "remove" the edge (i,j):
			i = index_i[edge];
			j = index_j[edge];

			if (info != NULL && info[i][j]==YES) { continue; } // This edge can not be removed

			// Set S contains the edges in the 'i' side of the T\{(i,j)}
			for (u=0; u<N; u++) { S[u] = NO; }
			// We retrieve S by means of a DFS:
			S[i] = YES;
			stack_size = 0;
			stack[stack_size++] = i;
			while (stack_size>0){
				u = stack[--stack_size];
				for (d=0; d<degree[u]; d++) {
					if (neighbours[u][d]!=j && S[neighbours[u][d]]==NO) {
						stack[stack_size++] = neighbours[u][d];
						S[neighbours[u][d]] = YES;
					}
				}
			}

			// Now compute the amount of communication 'w[u]' from each node 'u'
			// to the other side of the tree:
			for (u=0; u<N; u++) {
				w[u] = 0;
				for (v=0; v<N; v++) {
					if (S[u] != S[v]) {
						w[u] += W[u][v] + W[v][u];
					}
				}
			}
			sum_w = 0;
			for (u=0; u<N; u++) {
				if (S[u] == YES) { sum_w += w[u]; }
			}

			// And the communication cost from nodes of the same side:
			for (u=0; u<N; u++) {
				cc[u] = 0;
				for (v=0; v<N; v++) {
					if (u!=v && S[u] == S[v]) {
						cc[u] += w[v] * dist[u][v];
					}
				}
			}

			// now find the edge that minimices the CC:
			best_cost = cc[i] + C[i][j]*sum_w + cc[j];
			min_i = i;
			min_j = j;
			for (u=0; u<N; u++) {
				for (v=u+1; v<N; v++) {
					if (S[u] != S[v] && (info == NULL || info[u][v] != NO)) {
						current_cost = cc[u] + C[u][v]*sum_w + cc[v];
						if ((info != NULL && info[i][j] == YES) || current_cost < best_cost) {
							current_cost = current_cost;
							min_i = u;
							min_j = v;
						}
						if (info != NULL && info[u][v] == YES) { break; }	//
					}
				}
			}

			// If we have improved, remember how!
			if (min_i != i || min_j != j) {
				if (((cc[i]+C[i][j]*sum_w+cc[j])-best_cost) > best_delta) {
					best_delta = (cc[i] + C[i][j]*sum_w + cc[j]) - best_cost;
					best_edge = edge;
					best_new_i = min_i;
					best_new_j = min_j;
				}
				if (BEST_IMPROVEMENT == NO) { break; }
			} 

		}

		if (best_delta > 0.0) {

			i = index_i[best_edge];
			j = index_j[best_edge];

			// remove j from new_neighbours[i]:
			found = NO;
			degree[i]--;
			for (d=0; d<degree[i]; d++) {
				if (found == NO && neighbours[i][d] == j) { found = YES; }
				if (found == YES) { neighbours[i][d] = neighbours[i][d+1]; }
			}
			neighbours[i][degree[i]] = NONE;

			// remove i from new_neighbours[j]:
			found = NO;
			degree[j]--;
			for (d=0; d<degree[j]; d++) {
				if (found == NO && neighbours[j][d] == i) { found = YES; }
				if (found == YES) { neighbours[j][d] = neighbours[j][d+1]; }
			}
			neighbours[j][degree[j]] = NONE;

			// add best_new_j to new_neighbours[best_new_i]:
			neighbours[best_new_i][degree[best_new_i]++] = best_new_j;

			// add best_new_i to new_neighbours[best_new_j]:
			neighbours[best_new_j][degree[best_new_j]++] = best_new_i;

			// add current delta to the global delta: 
			delta += best_delta;
		}	
	}

	// Copy best_solution to Y:
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { Y[i][j] = 0.0; } }  
	for (i=0; i<N; i++) {
		for (d=0; d<degree[i]; d++) {
			Y[i][neighbours[i][d]] = 1.0;
		}
	}

	// clean:
	free_int_matrix(&neighbours, N);
	free_int_vector(&degree);
	free_int_vector(&index_i);
	free_int_vector(&index_j);
	free_int_vector(&S);
	free_int_vector(&stack);
	free_double_matrix(&dist, N);
	free_double_vector(&w);
	free_double_vector(&cc);

	return delta;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Best-improve local seach for the 1-character exchange neighborhood:
double DLS(int N, double **C, double **W, double **Y) {
	double **dist = create_double_matrix(N,N);
	int **neighbors = create_int_matrix(N,N);
	int *degree = create_int_vector(N);
	int	*code = create_int_vector(N);
	int *best = create_int_vector(N);
	int *stack = create_int_vector(N);

	int i,j,k,l,m,d, size, old;
	double val, best_val, old_val;
	int STOP = NO;

	for (i=0; i<N; i++) { code[i] = NONE;}
	for (i=0; i<N; i++) { best[i] = NONE;}

	// INPUT:
	matrix_to_code(N, Y, best);
	tree_dist(N, C, Y, dist);
	old_val = 0.0;
	for (i=0; i<N; i++) { for(j=0; j<N; j++) { old_val += dist[i][j] * W[i][j]; } }
	best_val = old_val;

	// DLS
	while (STOP == NO) {
		
		STOP = YES;
		for (i=0; i<N; i++) { code[i] = best[i]; }

		for (i=0; i<N-2; i++) {

			// Remember:
			old = code[i];	

			for (j=0; j<N; j++) {
				if (j!= old) {

					// Decode tree:
					code[i] = j;
					code_to_tree(N, code, neighbors, degree);
					
					// Compute Distances:
					for (k=0; k<N; k++) {
						for (l=0; l<N; l++) { dist[k][l] = -1.0; }
						dist[k][k] = 0.0;
						size = 0;
						stack[size++] = k;
						while (size > 0) {
							l = stack[--size];
							for (d=0; d<degree[l]; d++) {
								m = neighbors[l][d];
								if (dist[k][m] == -1.0) {
									dist[k][m] = dist[k][l] + C[l][m];
									stack[size++] = m;
								}
							}
						}
					}

					// Compute obj-function
					val = 0.0;
					for (k=0; k<N; k++) { for(l=0; l<N; l++) { val += dist[k][l] * W[k][l]; } }

					// Update best code:
					if (val < best_val) {
						STOP = NO;
						for (k=0; k<N; k++) { best[k] = code[k]; }
						best_val = val;
					}
				}
			}

			// Restore:
			code[i] = old;
		}
	}
	
	// OUTPUT:
	code_to_matrix(N, best, Y);
	
	// Clean:
	free_int_vector(&stack);
	free_int_vector(&best);
	free_int_vector(&code);
	free_int_vector(&degree);
	free_int_matrix(&neighbors, N);
	free_double_matrix(&dist, N);

	return old_val - best_val;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// First-improve local search for the 2-character-exchange neighborhood 
double DLS2(int N, double **C, double **W, double **Y) {
	double **dist = create_double_matrix(N,N);
	int **neighbors = create_int_matrix(N,N);
	int *degree = create_int_vector(N);
	int	*code = create_int_vector(N);
	int *best = create_int_vector(N);
	int *stack = create_int_vector(N);

	int i,j,k,l,m,d, size, old, oldold, ii, jj;
	double val, best_val, old_val;
	int STOP = NO;
	int STOP2 = NO;

	for (i=0; i<N; i++) { code[i] = NONE;}
	for (i=0; i<N; i++) { best[i] = NONE;}

	// INPUT:
	matrix_to_code(N, Y, best);
	tree_dist(N, C, Y, dist);
	old_val = 0.0;
	for (i=0; i<N; i++) { for(j=0; j<N; j++) { old_val += dist[i][j] * W[i][j]; } }
	best_val = old_val;

	// DLS + DLS2
	while (STOP2 == NO) {
		STOP2 = YES;

		// DLS
		while (STOP == NO) {
			STOP = YES;
			for (i=0; i<N; i++) { code[i] = best[i]; }
			for (i=0; i<N-2; i++) {
				old = code[i];	
				for (j=0; j<N; j++) {
					if (j!= old) {
						code[i] = j;
						code_to_tree(N, code, neighbors, degree);
						for (k=0; k<N; k++) {
							for (l=0; l<N; l++) { dist[k][l] = -1.0; }
							dist[k][k] = 0.0;
							size = 0;
							stack[size++] = k;
							while (size > 0) {
								l = stack[--size];
								for (d=0; d<degree[l]; d++) {
									m = neighbors[l][d];
									if (dist[k][m] == -1.0) {
										dist[k][m] = dist[k][l] + C[l][m];
										stack[size++] = m;
									}
								}
							}
						}
						val = 0.0;
						for (k=0; k<N; k++) { for(l=0; l<N; l++) { val += dist[k][l] * W[k][l]; } }
						if (val < best_val) {
							STOP = NO;
							for (k=0; k<N; k++) { best[k] = code[k]; }
							best_val = val;
						}
					}
				}
				code[i] = old;
			}
		}

		// If we arrive here we have reached a local optimum of the 1-character-exchange neighborhood
		// Try to scape!

		// DLS 2  (first improve!)
		for (i=0; i<N; i++) { code[i] = best[i]; }
		for (i=0; i<N-2; i++) {
			if (STOP2 == NO) { break; }
			old = code[i];	
			for (ii=0; ii<N-2; ii++) {
				if (STOP2 == NO) { break; }
				if (ii!=i) {
					oldold = code[ii];
					for (j=0; j<N; j++) {
						if (STOP2 == NO) { break; }
						if (j!= old) {
							code[i] = j;
							for (jj=0; jj<N; jj++) {
								if (STOP2 == NO) { break; }
								if (jj!=oldold) {
									code[ii] = jj;
									code_to_tree(N, code, neighbors, degree);
									for (k=0; k<N; k++) {
										for (l=0; l<N; l++) { dist[k][l] = -1.0; }
										dist[k][k] = 0.0;
										size = 0;
										stack[size++] = k;
										while (size > 0) {
											l = stack[--size];
											for (d=0; d<degree[l]; d++) {
												m = neighbors[l][d];
												if (dist[k][m] == -1.0) {
													dist[k][m] = dist[k][l] + C[l][m];
													stack[size++] = m;
												}
											}
										}
									}
									val = 0.0;
									for (k=0; k<N; k++) { for(l=0; l<N; l++) { val += dist[k][l] * W[k][l]; } }
									if (val < best_val) {
										STOP2 = NO;
										for (k=0; k<N; k++) { best[k] = code[k]; }
										best_val = val;
									}
								} 
							}
						}
					}
					code[ii] = oldold;
				}
			}
			code[i] = old;
		}
	}
	
	// OUTPUT:
	code_to_matrix(N, best, Y);
	
	// Clean:
	free_int_vector(&stack);
	free_int_vector(&best);
	free_int_vector(&code);
	free_int_vector(&degree);
	free_int_matrix(&neighbors, N);
	free_double_matrix(&dist, N);

	return old_val - best_val;
}
////////////////////////////////////////////////////////////////////////////////
