#include "headers.h"

/// INFO ///////////////////////////////////////////////////////////////////////
//
// CLM_graph.c: Graph Theory related functions
//
// Author:  Carlos Luna-Mota
// Created: 2015 Jan 21
// Updated: 2017 Jun 23
//
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Applyes the Floyd-Warshall algorithm (in place) to a distance matrix D.
// Takes O(N^3) time (optimal for dense graphs)
void Floyd_Warshall(int N, double **D) {
	int	i,j,k;			// matrix indices
	
	// Floyd Warshall:
	for (k=0; k<N; k++) {
		for (i=0; i<N; i++) {
			for (j=0; j<N; j++) {
				D[i][j] = MIN(D[i][j], D[i][k]+D[k][j]);
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Computes an UpperBound over the value of the longest Hamiltonian Path 
// in a graph of N nodes with a cost matrix 'C'
double LongestHamiltonianPathUB(int N, double **C) {
	double UB, row_max, min_row_max;
	int row,col;

	// Compute for first row to initialize everything:
	row_max = C[0][0];
	for (col=1; col<N; col++) { row_max = MAX(row_max, C[0][col]); }
	min_row_max = row_max;
	UB = row_max;

	// Compute for the rest of rows:
	for (row=1; row<N; row++) {
		row_max = C[row][0];
		for (col=1; col<N; col++) { row_max = MAX(row_max, C[row][col]); }
		min_row_max = MIN(min_row_max, row_max);
		UB += row_max;
	}

	// Return the sum of the N-1 biggest edges of the graph, 
	// provided that each row contains at most one of them:
	return UB - min_row_max;
}
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
// Computes the min-cut matrix and the min-cut-min-cost matrix for all 0<=i!=j<N
void compute_min_cut_data(int N, double **C, double **D, double **min_cut, double **min_cut_min_cost) {
	int		i, j;
	double **capacity = create_double_matrix(N, N);
    double *flow	  = create_double_vector(N);
	int    *tree      = create_int_vector(N);

    int		tail,head;
	int		queueSize, current;
	int		*visited = create_int_vector(N);
	int		*queue	 = create_int_vector(N);
	
	// Init everything:
	for (i=0; i<N; i++) { tree[i] = NONE; }
	for (i=0; i<N; i++) { flow[i] = 0.0;  }
	for (i=0; i<N; i++) {
		capacity[i][i] = 0.0;
		for (j=i+1; j<N; j++) {
			capacity[i][j] = capacity[j][i] = D[i][j]+D[j][i];
		}
	}
	for (i=0; i<N; i++) {
		min_cut[i][i] = 0.0;
		for (j=i+1; j<N; j++) {
			min_cut[i][j] = min_cut[j][i] = -1.0;
		}
	}
	for (i=0; i<N; i++) {
		min_cut_min_cost[i][i] = 0.0;
		for (j=i+1; j<N; j++) {
			min_cut_min_cost[i][j] = min_cut_min_cost[j][i] = -1.0;
		}
	}

	// compute the Gomory-Hu Tree
	GH_tree(N, capacity, tree, flow);
	
	for (head=1; head<N; head++) {

		tail = tree[head];
		tree[head] = head;

		queue[0] = tail;				// we are going to recover the subset S that contains tail
		queueSize = 1;                  // so far it only contains one element
		current = 0;                    // I use the cutNodes array and this int as a queue of size nCutNodes
		for (i=0; i<N; i++) { visited[i] = NO; }	// nobody is visited yet...
		visited[tail] = YES;						// nobody but tail, of course.
		
		// now iterate following the edges of the tree defined by the GHtree array until find S
		while (current < queueSize) {
			for (i=0; i<N; i++) {
				// case 1: queue[current] -> GHtree[queue[current]]
				if (i == queue[current] && visited[tree[i]] == NO) {
					visited[tree[i]] = YES;       // mark as visited
					queue[queueSize] = tree[i];   // put it in the next empty slot of the queue
					queueSize++;                    // update queue size
				}
				// case 2: i -> queue[current]
				if( tree[i] == queue[current] && visited[i] == NO ) {
					visited[i] = YES;               // mark as visited
					queue[queueSize] = i;			// put it in the next empty slot of the queue
					queueSize++;                    // update queue size    
				}    
			}        
			current++;
		}

		// Retore the edge!!!
		tree[head] = tail;			
			
		// Now, for every edge crossing the cut: assign the value of the flow (if minimal)
		for (i=0; i<N-1; i++) {
			for (j=i+1; j<N; j++) {
				if (visited[i] != visited[j]) {
					if (min_cut[i][j] == -1.0) {
						min_cut[i][j] = min_cut[j][i] = flow[head];
					} else {
						min_cut[i][j] = min_cut[j][i] = MIN(min_cut[i][j], flow[head]);
					}
					if (min_cut_min_cost[i][j] == -1.0) {
						min_cut_min_cost[i][j] = min_cut_min_cost[j][i] = ((double)(C[i][j]));
					} else {
						min_cut_min_cost[i][j] = min_cut_min_cost[j][i] = MIN(min_cut_min_cost[i][j], ((double)(C[i][j])));
					}
				}
			}
		}
	}

	// Clean:
	free_double_matrix(&capacity, N);
	free_double_vector(&flow);
	free_int_vector(&tree);
	free_int_vector(&visited);
	free_int_vector(&queue);
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Take as argument the capacity matrix 'Capacity' and return the arrays
// 'Tree' and 'Flow'. It cost as much as O( n*f(n) )  where f(n) is how much 
// it cost to make a max_flow_and_cut query.
void GH_tree(int N, double **Capacity, int *Tree, double *Flow){  
    int		i;							// index variables for loops
    int		s,t;						// source and sink nodes
    double	auxFlow;					// aux variable to make a switch
    int		*in_s_side = create_int_vector(N);
	
    // initialize everything
    for (i=0; i<N; i++) {                // the starting tree is star-like
        Tree[i] = 0;                     //  with the '0' as central node
        Flow[i] = 0.0;                   //  and null flow in each edge
    }
    
    // visit the nodes and update 'Tree' and 'Flow'
    for (s=1; s<N; s++) {                // the '0' node is already 'visited'

		// in this iteration we are going to update the edge 's-Tree[s]'
        t = Tree[s];                                 
		Flow[s] = maxFlowAndCut(N, s, t, Capacity, in_s_side);   
		
        // we have to update all the nodes connected with 't' that are on the 's' side
        for (i=0; i<N; i++) {
            if (in_s_side[i] == YES  &&  i != s  &&  Tree[i] == t) { Tree[i] = s; }
		}

        // if the predecessor of 't' is in the connected component of 's' (='X')
        // we must swich the direction of the 's'-'t' link and the 's'-'t' flow!
        if (in_s_side[Tree[t]] == YES) {
            Tree[s] = Tree[t];
            Tree[t] = s;
            auxFlow = Flow[s];
            Flow[s] = Flow[t];
            Flow[t] = auxFlow;
        }
    }

	// clean:
	free_int_vector(&in_s_side);
}
////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////// // OJO!!!! Aun debugando!
// Returns up to N-1 cuts represented in the binary NxN-matrix cut_sets
// Uses a simplified variant of the Gusfield's Gomory-Hu-Tree routine.
double min_cap_cuts(int N, double **capacity, int **cut_sets){  
    int		i,j;							// index variables for loops
    int		s,t;							// source and sink nodes
	int		*tree = cut_sets[N-1];			// using the last row as a tree structure
	double	min_flow;

	// Compute dummy min_flow:
	min_flow = 0.0;
	for (i=1; i<N; i++) { min_flow += capacity[0][i]; }

	// RESET cut_sets!
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			cut_sets[i][j] = 0;
		}
	}

	for (i=0; i<N; i++) { tree[i] = 0; }    // we start with the 0-star tree
    
    for (s=1; s<N; s++) {					// visit the other nodes and update 'tree'
    
		// in this iteration we are going to update the edge 's-tree[s]'
        t = tree[s];                                 
		min_flow = MIN(min_flow, maxFlowAndCut(N, s, t, capacity, cut_sets[s-1]));   
		
        // we have to update all the nodes connected with 't' that are on the 's' side
        for (i=0; i<N; i++) {
            if (cut_sets[s-1][i] == YES  &&  i != s  &&  tree[i] == t) { tree[i] = s; }
		}

        // if the predecessor of 't' is in the connected component of 's' (='X')
        // we must swich the direction of the 's'-'t' link!
        if (cut_sets[s-1][tree[t]] == YES) {
            tree[s] = tree[t];
            tree[t] = s;
        }
    }

	// clean:
	for (i=0; i<N; i++) { tree[i] = 0; }

	return min_flow;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Returns up to N-1 cuts represented in the binary NxN-matrix cut_sets
double connected_components(int N, double **capacity, int **cut_sets) {
	int	i,j;
    int	queueSize, current;
	int	starting_point, n_conn_comp;
	int	*visited = create_int_vector(N);
	int	*queue	 = create_int_vector(N);
	
	// RESET cut_sets!
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			cut_sets[i][j] = 0;
		}
	}

	// Init everything:
	for (i=0; i<N; i++) { visited[i] = NO;  }
	queueSize = 0;
	starting_point = 0;
	n_conn_comp = 0;

	// While not all nodes visited:
	while (queueSize<N) {
		// look for an unvisited node and add it to the queue:
		for (i=starting_point; i<N; i++){
			if (visited[i] == NO) {
				starting_point = i+1;
				visited[i] = YES;
				current = queueSize;
				queue[queueSize++] = i;
				break;
			}
		}
		// We want to rememenber the connected component of this element:
		for (i=0; i<N; i++) { cut_sets[n_conn_comp][i] = NO; }
		// Lets visit a whole connected component using a BFS:
		while (queueSize>current) {
			i = queue[current++];
			cut_sets[n_conn_comp][i] = YES;
			for (j=0; j<N; j++) {
				if (visited[j] == NO && capacity[i][j] > 0.5) {	// we are assuming integer capacities {0.0, 1.0}
					visited[j] = YES;
					queue[queueSize++] = j;
				}
			}
		}
		n_conn_comp++;
	}
	
	// We can erase the last one: its the complementary of the other one!
	if (n_conn_comp==2) for (i=0; i<N; i++) { cut_sets[n_conn_comp-1][i] = NO; }

	// clean
	free_int_vector(&visited);
	free_int_vector(&queue);

	// Now compute the min flow:
	if (n_conn_comp == 1) { return 1.0; }
	else                  { return 0.0; }
		
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Performs the Edmonds-Karp Algorithm and returns the maximum flow between
// 's' and 't'. The algorithm also remembers the connected component that
// contains 's' in the binary array 'visited'
double maxFlowAndCut(int N, int s, int t, double **Capacity, int *visited){
    int		i,j;							// index variables for loops
    int		u,v;                            // index variables for nodes
	int		queueIndex, queueSize;          // just for use an array as a queue
    double	flow;							// the augmenting flow
    double	maxFlow = 0.0;                  // the total amount of flow

    double  **currentFlow = create_double_matrix(N,N);
	int		*previous = create_int_vector(N);
	int		*queue = create_int_vector(N);	// basic queue

	int DEBUG = s<0 ? YES : NO;
	s = ABS(s);

    // initialize the 'currentFlow'    
    for (i=0; i<N; i++) { for (j=0; j<N; j++) { currentFlow[i][j] = 0.0; } }
    
    // perform as many flow augmentations as you can!
    flow = 1.0;                                     // this is just to start the loop
    while (flow > 0.0) {
        flow = 0.0;
        // initialize everything else!
		for (i=0; i<N; i++) { visited[i] = NO; }    // nobody has been visited yet
        visited[s] = YES;                           //  ...nobody but 's'
        queueSize = 0;                              // the queue contains 0 elements
		queue[queueSize++] = s;                     // now it contains s
        queueIndex = 0;                             // next element to leave the queue
      
        // perform a BFS to find a augmenting path between 's' and 't'
        while (queueIndex < queueSize){
            u = queue[queueIndex++];                // get the next element
            // now look if you can extend the augmenting path from 'u' to any other 'v'
            for (v=0; v<N; v++) {
                if (visited[v]==NO && Capacity[u][v]-currentFlow[u][v] > 0.0) {
                    visited[v] = YES;               // if you can, visit the new node
                    queue[queueSize++] = v;           // put it in the queue
                    previous[v] = u;                // remeber where you come from
                    if (v==t) break;                // and avoid do extra work! 
                }
			}

            // now, before continue, check if we have finished the BFS!
            if (visited[t] == YES) {
                // we are going to trace back our steps looking for how much flow we can send
                v = t;
                u = previous[t];
                flow = Capacity[u][v]-currentFlow[u][v];
                while (u!=s) {
                    v = u;
                    u = previous[v];
                    if (flow > Capacity[u][v]-currentFlow[u][v])
                        flow = Capacity[u][v]-currentFlow[u][v];
                }
                // now we can update the currentFlowMatrix with that amount of flow
                v = t;
                u = previous[t];
                currentFlow[u][v] += flow; // we add the flow in one direction
                currentFlow[v][u] -= flow; // and remove it from the other one!
                while (u!=s) {
                    v = u;
                    u = previous[v];
                    currentFlow[u][v] += flow; // we add the flow in one direction
                    currentFlow[v][u] -= flow; // and remove it from the other one!
                }        
                maxFlow += flow;                    // we have improved so we update
                break;                              // avoid do extra work!
            }
        }  
    }

	if (DEBUG == YES) {
		printf("Last Flow = %.12f\n",flow);
		printf("Max Flow  = %.12f\n",maxFlow);
		flow = 0.0;
		for (i=0; i<N; i++) {
			for (j=0; j<N; j++) {
				if (currentFlow[i][j]>0.0) { 
					printf("(%d,%d)\t-> %.12f / %.12f\n", i,j,currentFlow[i][j], Capacity[i][j]); 
					if (j==t) { flow += currentFlow[i][j]; }
				}
			}
		}
		printf("Flow Sum = %.12f\n",flow);
	}

	// clean:
	free_double_matrix(&currentFlow, N);
	free_int_vector(&previous);
	free_int_vector(&queue);

    return maxFlow;
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// we assume that the structure stored in the adjacency matrix 'y' is a tree and 
// compute its distance matrix 'dist' based on the edge cost matrix 'C' in O(N^2)
void tree_dist(int N, double **C, double **y, double **dist) {
	int i,j,k,d;
	int **neighbours = create_int_matrix(N, N);
	int *degree		 = create_int_vector(N);
	int *stack		 = create_int_vector(N);
	int size;

	// Build a neighbour-list representation of the tree for speed-up
	for (i=0; i<N; i++) { degree[i] = 0; }	
	for (i=0; i<N; i++) {	
		for (j=i+1; j<N; j++) {
			if (y[i][j] + y[j][i] > 1.0-EPSILON) {
				neighbours[i][degree[i]++] = j;
				neighbours[j][degree[j]++] = i;
			}
		}
	}

	// perform N DFS to find the distance from each node to the rest:
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) { dist[i][j] = -1.0; }
		dist[i][i] = 0.0;
		size = 0;
		stack[size++] = i;
		while (size > 0) {
			j = stack[--size];
			for (d=0; d<degree[j]; d++) {
				k = neighbours[j][d];
				if (dist[i][k] == -1.0) {
					dist[i][k] = dist[i][j] + C[j][k];
					stack[size++] = k;
				}
			}
		}
	}

	// clean:
	free_int_matrix(&neighbours, N);
	free_int_vector(&degree);
	free_int_vector(&stack);	
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Puts an uniformly random tree in Y
void random_tree(int N, double **Y) {

	int i;
	int *code = create_int_vector(N-2);
	
	// Generate a Dandelion Code uniformly at random
	for (i=0; i<N-2; i++) { 	code[i] = randint(0, N-1); }

	// Translate the code to an adjacency matrix:
	code_to_matrix(N, code, Y);

	// Clean
	free_int_vector(&code);
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Takes the tree stored in the adjacency list (neighbours, degree) and
// stores its corresponding dandelion code in 'code'. 
// Time: O(N)
void tree_to_code(int N, int **neighbours, int *degree, int *code) {
	int i,queue_size, queue_index, current, previous, path_len;
	int current_min;
	int *succ  = create_int_vector(N);
	int *visit  = create_int_vector(N);
	int *queue = create_int_vector(N);
	int *path  = create_int_vector(N);
	const int FIRST = 0;
	const int LAST = N-1;
	
	for (i=0; i<N; i++) { visit[i] = NO; }

	// BFS from LAST to fill the Successor array:
	succ[LAST] = LAST;
	visit[LAST] = YES;
	
	queue[0] = LAST;
	queue_size = 1;
	queue_index = 0;

	// This is O(E) (= O(N) for a Tree!)
	while (queue_index < N) {
		current = queue[queue_index];
		for (i=0; i<degree[current]; i++) {
			previous = neighbours[current][i];
			if (visit[previous] == NO) {
				queue[queue_size] = previous;
				queue_size++;
				visit[previous] = YES;
				succ[previous] = current;
			}
		}
		queue_index++;
	}
    
	// Recover the FIRST-LAST-path
	path_len = 0;
	current = FIRST;
	while (current != LAST) {
		path[path_len] = current;
		current = succ[current];
		path_len++;
	}
	path[path_len] = LAST;
	
	// initialice the Dandelion Code:
	for (i=0; i<N-2; i++) { code[i] = NONE; }

	// break the path in cicles by finding left-to-right-mins
	current = path_len-1;
	previous = path_len;
	current_min = path[current];
	
	while (current>0) {
		if (current==1 || path[current-1] < current_min) {  // we have the cycle [current, previous)
	
			// every element in a cycle must point to the next one:
			for (i=current; i<previous-1; i++) { code[path[i]-1] = path[i+1]; }
			code[path[previous-1]-1] = path[current]; // the last elem of the cycle must point to the first!
	
			previous = current;
			current_min = path[current-1];
		}
		current--;
	}
	
	// the rest of elements must point to its successor:
	for (i=0; i<N-2; i++) { if (code[i] == NONE) { code[i] = succ[i+1]; } }

	// Clean
	free_int_vector(&succ);
	free_int_vector(&visit);
	free_int_vector(&queue);
	free_int_vector(&path);
}
// Time O(N^2)
void matrix_to_code(int N, double **Y, int *code) {
	int i,j,e,queue_size, queue_index, current, previous, path_len;
	int current_min;
	int *degree = create_int_vector(N);
	int **neighbours = create_int_matrix(N,N);
	int *succ  = create_int_vector(N);
	int *visit  = create_int_vector(N);
	int *queue = create_int_vector(N);
	int *path  = create_int_vector(N);
	const int FIRST = 0;
	const int LAST = N-1;
	
	// Copy the tree structure at Y to the neighbours-degree structure	// THIS IS THE ONLY O(N^2) part!
	e = 0;
	for (i=0; i<N; i++) { degree[i] = 0; }
	for (i=0; i<N-1; i++) { for (j=i+1; j<N; j++) {
		if (Y[i][j] + Y[j][i] > 1.0-EPSILON) {
			neighbours[i][degree[i]++] = j;
			neighbours[j][degree[j]++] = i;
			e++;
		}
	} }
	assert(e == N-1);
	//////////////////////////////////////////////////////////////////

	for (i=0; i<N; i++) { visit[i] = NO; }

	// BFS from LAST to fill the Successor array:
	succ[LAST] = LAST;
	visit[LAST] = YES;
	
	queue[0] = LAST;
	queue_size = 1;
	queue_index = 0;

	// This is O(E) (= O(N) for a Tree!)
	while (queue_index < N) {
		current = queue[queue_index];
		for (i=0; i<degree[current]; i++) {
			previous = neighbours[current][i];
			if (visit[previous] == NO) {
				queue[queue_size] = previous;
				queue_size++;
				visit[previous] = YES;
				succ[previous] = current;
			}
		}
		queue_index++;
	}
    
	// Recover the FIRST-LAST-path
	path_len = 0;
	current = FIRST;
	while (current != LAST) {
		path[path_len] = current;
		current = succ[current];
		path_len++;
	}
	path[path_len] = LAST;
	
	// initialice the Dandelion Code:
	for (i=0; i<N-2; i++) { code[i] = NONE; }

	// break the path in cicles by finding left-to-right-mins
	current = path_len-1;
	previous = path_len;
	current_min = path[current];
	
	while (current>0) {
		if (current==1 || path[current-1] < current_min) {  // we have the cycle [current, previous)
	
			// every element in a cycle must point to the next one:
			for (i=current; i<previous-1; i++) { code[path[i]-1] = path[i+1]; }
			code[path[previous-1]-1] = path[current]; // the last elem of the cycle must point to the first!
	
			previous = current;
			current_min = path[current-1];
		}
		current--;
	}
	
	// the rest of elements must point to its successor:
	for (i=0; i<N-2; i++) { if (code[i] == NONE) { code[i] = succ[i+1]; } }

	// Clean
	free_int_vector(&succ);
	free_int_vector(&visit);
	free_int_vector(&queue);
	free_int_vector(&path);
	free_int_vector(&degree);
	free_int_matrix(&neighbours, N);
}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// Takes the dandelion code 'code' and stores its corresponding tree in the
// adjacency list (neighbours, degree). 
// Time: O(N)
void code_to_tree(int N, int *code, int **neighbours, int *degree) {
	int i,j, current, current_min;
	int *orbit  = create_int_vector(N);
	int *isMin  = create_int_vector(N);
	int *long_code = create_int_vector(N);
	const int FIRST = 0;
	const int LAST = N-1;

	// Extend the code: (0, (c0, c1, c2, ... cN-3), N-1)
	long_code[FIRST] = FIRST;
	long_code[LAST] = LAST;
	for (i=0; i<N-2; i++) { long_code[i+1] = code[i]; }

	// Now compute the isMin array in linear time:
	for (i=0; i<N; i++) { isMin[i] = NO; }
	for (i=0; i<N; i++) { orbit[i] = NONE; }
	for (i=0; i<N; i++) {
		j = i;
		while (orbit[j] == NONE) {
			orbit[j] = i;
			j = long_code[j];
		}
		if (orbit[j] == i) {
			current = j;
			current_min = j;
			while (long_code[current] != j) {
				current = long_code[current];
				current_min = MIN(current_min, current);
			}
			isMin[current_min] = YES;
		}
	}

	// finally store the tree:
	for (i=0; i<N; i++) { degree[i] = 0; }
	current_min = FIRST;
	for (current=1; current<N; current++) {
		if (isMin[current] == YES) {
			i = long_code[current];
			j = current_min;
			neighbours[i][degree[i]] = j;
			neighbours[j][degree[j]] = i;
			degree[i]++;
			degree[j]++;
			current_min = current;
		} else {
			i = long_code[current];
			j = current;
			neighbours[i][degree[i]] = j;
			neighbours[j][degree[j]] = i;
			degree[i]++;
			degree[j]++;
		}
	}
	
	// Clean
	free_int_vector(&long_code);
	free_int_vector(&orbit);
	free_int_vector(&isMin);
}
// Time: O(N^2)
void code_to_matrix(int N, int *code, double **Y) {
	int i,j, current, current_min;
	int *orbit  = create_int_vector(N);
	int *isMin  = create_int_vector(N);
	int *long_code = create_int_vector(N);
	const int FIRST = 0;
	const int LAST = N-1;

	// Extend the code: (0, (c0, c1, c2, ... cN-3), N-1)
	long_code[FIRST] = FIRST;
	long_code[LAST] = LAST;
	for (i=0; i<N-2; i++) { long_code[i+1] = code[i]; }

	// Now compute the isMin array in linear time:
	for (i=0; i<N; i++) { isMin[i] = NO; }
	for (i=0; i<N; i++) { orbit[i] = NONE; }
	for (i=0; i<N; i++) {
		j = i;
		while (orbit[j] == NONE) {
			orbit[j] = i;
			j = long_code[j];
		}
		if (orbit[j] == i) {
			current = j;
			current_min = j;
			while (long_code[current] != j) {
				current = long_code[current];
				current_min = MIN(current_min, current);
			}
			isMin[current_min] = YES;
		}
	}

	// finally store the tree:
	for (i=0; i<N; i++) { for (j=0; j<N; j++) { Y[i][j] = 0.0; } }
	current_min = FIRST;
	for (current=1; current<N; current++) {
		if (isMin[current] == YES) {
			i = long_code[current];
			j = current_min;
			Y[i][j] = Y[j][i] = 1.0;
			current_min = current;
		} else {
			i = long_code[current];
			j = current;
			Y[i][j] = Y[j][i] = 1.0;
		}
	}
	
	// Clean
	free_int_vector(&long_code);
	free_int_vector(&orbit);
	free_int_vector(&isMin);
}
////////////////////////////////////////////////////////////////////////////////


