/*
    Copyright Vladimir Kolmogorov vnk@ist.ac.at 2013

    This file is part of POTTS.

    POTTS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    POTTS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef JGASFJASGNASJAIUHADGD
#define JGASFJASGNASJAIUHADGD

#define __VERBOSE__H_ // enable this for verbose output in this package

#ifdef __VERBOSE__H_
#define exec_r(r) r;
#else
#define exec_r(r)
#endif

#define REUSE_GRAPH_STRUCTURE // if this option is on, the alpha expansion with graph reusing uses a bit less memory

///////////////////////////////////////////////////////////////////////////////////////////////

#include "maxflow/graph.h"
#ifdef REUSE_GRAPH_STRUCTURE 
#include "maxflow/graphX.h"
#endif

class Potts
{
public:
	typedef int REAL;
	typedef int LabelId;
	typedef int NodeId;
	
	// 'unaries_fn': pointer to the function that computes unary cost for a given node i and label k.
	// This method of specifying unaries saves memory (since no internal array for unaries is allocated),
	// but could be a bit slower (??).
	//


	// It is possible to add more edges than edge_num_max 
	// (and edge_num_max can be zero). However, if the count is exceeded, then 
	// the internal memory is reallocated (increased by 50%) which is expensive. 
	// Also, temporarily the amount of allocated memory would be more than twice than needed.
	Potts(int node_num, int edge_num_max, int label_num, REAL (*unaries_fn)(NodeId i, LabelId k) = NULL); 
	~Potts();

	//////////////////////////////////////////////////////////////////////////////
	////////////////////////////// CREATING INSTANCE /////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	// Returns a pointer to an array of size label_num (which the user can then modify).
	// Initially all entries are zeros.
	// Cannot be called if unaries_fn != NULL in the constructor.
	REAL* GetUnariesPtr(NodeId i) { return unaries + i*K; }


	void AddEdge(NodeId i, NodeId j, REAL lambda);

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////     ALGORITHMS     ////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

// Options. 
#define MAXIMIZE_PERSISTENCY      (0x01) // if true, run a postprocessing to label as many nodes as possible in a single Kovtun reduction.
                                         // Matters if e.g. there are many ties in the unary costs.

#define INIT_WITH_KOVTUN_SOLUTION (0x02) // ways to initialize alpha-expansion.
#define INIT_WITH_KOVTUN_FLOW     (0x04) // if off, init with zero solution/flow

#define POTTS_DEBUG               (0x08) // if on, test that Kovtun's phase produces correct solution and flow

	// First run Kovtun's reduction using parametric maxflow,
	// then run 'iter_num_max' iterations of alpha-expansion (iter_num_max >= 0).
	// Returns the cost of the solution.
	// Second parameter:
	//   reuse_graphs=false: naive alpha-expansion
	//   reuse_graphs=true: starting with the second iteration, reuse graphs & flow, as in [Alahari, Kohli, Torr, PAMI'10]
	//                    (should be faster, but needs much more memory)
	// Third parameter: the default is recommended
	REAL Solve(int iter_num_max, bool reuse_graphs, unsigned options = INIT_WITH_KOVTUN_SOLUTION);


	// After calling Solve(), you may also call expansions manually for individual labels
	REAL DoExpansion(LabelId k);


	// Read solution
	LabelId GetSolution(NodeId i) { return solution[i] % K; } 
	bool isPersistent(NodeId i) { return (solution[i] >= K); } // true if labeled by Kovtun's reduction

	int GetUnlabeledNum() { return unlabeled_num; } // returns total number of non-persistent nodes









	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////


private:

	//////////// binary tree whose nodes are intervals (= subsets of labels) ////////
	typedef int IntervalId;
	struct Interval
	{
		LabelId i, j, split;
		IntervalId parent;
		IntervalId children[2]; // left and right
		int depth;
		bool is_leaf;
	};
	void ConstructSubtree(IntervalId I);
	IntervalId GetLCA(LabelId i, LabelId j);
	int GetSign(LabelId i, IntervalId j);
	int GetSign(LabelId a, LabelId i, LabelId j) { return GetSign(a, GetLCA(i, j)); }

	Interval* intervals; // 0,...,K-1 correspond to leaves, K is the root
	int interval_num;
	int** LCA; // LCA[i][j] = GetLCA(i, j)   i<K, j<K
	int** SIGN;    // SIGN[i][j] = GetSign(i, j)  i<K, j<interval_num
	///////////////////////////////////////////

	typedef Graph<REAL,REAL,REAL> G;
#ifdef REUSE_GRAPH_STRUCTURE 
	typedef GraphX<REAL,REAL,REAL> GX;
#else
	typedef Graph<REAL,REAL,REAL> GX;
#endif
	struct Node
	{
		NodeId	next;
		int		k_min;
		REAL	f_delta; // second smallest - smallest ( >= 0 )
		REAL	lambda; // = -g0
		union
		{
			REAL	current_unary;
			NodeId	flag; // used in the last phase of Solve (in max_persistent)
		};
	};

	int node_num, K;
	Node* nodes;
	LabelId* solution;
	NodeId* label_lists; // of size K
	G* graph;

	REAL* unaries;
	REAL (*unaries_fn)(NodeId i, LabelId k);

	void InitKovtun();
	void SolveKovtun();
	void SolveKovtun(Interval I, NodeId list);
	void MaximizePersistency();
	void TestKovtun(); // for debugging

	int unlabeled_num;

	///////////////////////////////////
	// second stage - after Kovtun
	struct Node2
	{
		NodeId orig_id; // \in [0,node_num-1]
		LabelId kovtun_label; // needed for extracting flow
		LabelId solution;
	};
	struct Edge2
	{
		NodeId i, j;
		REAL lambda;
		//REAL message0; // before final maxflow in SolveKovtun()
		//REAL message1; // after final maxflow in SolveKovtun()
	};
	struct LabelInfo
	{
		GX* graph;
		LabelId* old_labels; // of size unlabeled_num;
		REAL source_cost;
	};

	Node2* nodes2; // of size unlabeled_num
	Edge2* edges2;
	REAL* messages;
	int edge_num2;
	LabelInfo* label_info; // of size K
#ifdef REUSE_GRAPH_STRUCTURE
	GX* graph_ptr; // pointer, don't delete it!
#endif

	REAL cost0, current_cost;
	bool last_expansion_successful;

	bool reuse_graphs;
	bool maximize_persistency;
	bool init_with_kovtun_solution;
	bool init_with_kovtun_flow;
	bool potts_debug;

	void InitExpansion();
	void SolveExpansion(int iter_num_max);
	void DoExpansionNaive(LabelId k);
	void DoExpansionReuse_first(LabelId k);
	void DoExpansionReuse_second(LabelId k);

	void ComputeCurrentCost();

	////////////////////////////////////
	// if unaries_fn is used, then we need to store explicitly unaries that arose after Kovtun's reduction
	struct Unary2
	{
		LabelId k;
		REAL lambda;
	};
	struct Unary2Ptr
	{
		union
		{
			Unary2* ptr;
			int num; // temporary, used during construction
		};
		int k_cashed;
		REAL u_cashed;
	};
	Unary2Ptr* unaries2; // of size unlabeled_num+1
	REAL unaries2_fn(int i, LabelId k, bool is_frequent=false)
	{
		if (k == unaries2[i].k_cashed) return unaries2[i].u_cashed;

		REAL u = (*unaries_fn)(nodes2[i].orig_id, k);
		Unary2* ptr;
		for (ptr=unaries2[i].ptr; ptr<unaries2[i+1].ptr; ptr++)
		{
			if (ptr->k == k) u -= ptr->lambda;
		}

		if (is_frequent)
		{
			unaries2[i].k_cashed = k;
			unaries2[i].u_cashed = u;
		}

		return u;
	}
};

inline void Potts::AddEdge(NodeId i, NodeId j, REAL lambda)
{
	graph->add_edge(i, j, lambda, lambda);
}


#endif
