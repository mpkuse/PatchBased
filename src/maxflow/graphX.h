/* graph.h */
/*
	version 3.02

	This software library implements the maxflow algorithm
	described in

		"An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision."
		Yuri Boykov and Vladimir Kolmogorov.
		In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), 
		September 2004

	This algorithm was developed by Yuri Boykov and Vladimir Kolmogorov
	at Siemens Corporate Research. To make it available for public use,
	it was later reimplemented by Vladimir Kolmogorov based on open publications.

	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication.

	----------------------------------------------------------------------

	REUSING TREES:

	Starting with version 3.0, there is a also an option of reusing search
	trees from one maxflow computation to the next, as described in

		"Efficiently Solving Dynamic Markov Random Fields Using GraphX Cuts."
		Pushmeet Kohli and Philip H.S. Torr
		International Conference on Computer Vision (ICCV), 2005

	If you use this option, you should cite
	the aforementioned paper in any resulting publication.
*/
	


/*
	For description, license, example usage see README.TXT.
*/

#ifndef __GRAPHX_H__
#define __GRAPHX_H__

#include <string.h>
#include "block.h"

#include <assert.h>
// NOTE: in UNIX you need to use -DNDEBUG preprocessor option to supress assert's!!!



// captype: type of edge capacities (excluding t-links)
// tcaptype: type of t-links (edges between nodes and terminals)
// flowtype: type of total flow
//
// Current instantiations are in instances.inc
template <typename captype, typename tcaptype, typename flowtype> class GraphX
{
public:
	typedef enum
	{
		SOURCE	= 0,
		SINK	= 1
	} termtype; // terminals 
	typedef int node_id;
	typedef int arc_id;
#define NULL_ARC_ID ((arc_id)-1)
#define SISTER(a) ((arc_id) ((a) + ((a&1) ? -1 : +1)))

	/////////////////////////////////////////////////////////////////////////
	//                     BASIC INTERFACE FUNCTIONS                       //
	//              (should be enough for most applications)               //
	/////////////////////////////////////////////////////////////////////////

	// Constructor. 
	// The first argument gives an estimate of the maximum number of nodes that can be added
	// to the graph, and the second argument is an estimate of the maximum number of edges.
	// The last (optional) argument is the pointer to the function which will be called 
	// if an error occurs; an error message is passed to this function. 
	// If this argument is omitted, exit(1) will be called.
	//
	// IMPORTANT: It is possible to add more nodes to the graph than node_num_max 
	// (and node_num_max can be zero). However, if the count is exceeded, then 
	// the internal memory is reallocated (increased by 50%) which is expensive. 
	// Also, temporarily the amount of allocated memory would be more than twice than needed.
	// Similarly for edges.
	// If you wish to avoid this overhead, you can download version 2.2, where nodes and edges are stored in blocks.
	GraphX(int node_num_max, int edge_num_max, void (*err_function)(const char *) = NULL);

	// copy graph from g (except for residual capacities)
	GraphX(GraphX* g);

	// Destructor
	~GraphX();

	// Adds node(s) to the graph. By default, one node is added (num=1); then first call returns 0, second call returns 1, and so on. 
	// If num>1, then several nodes are added, and node_id of the first one is returned.
	// IMPORTANT: see note about the constructor 
	node_id add_node(int num = 1);

	// Adds a bidirectional edge between 'i' and 'j' with the weights 'cap' and 'rev_cap'.
	// IMPORTANT: see note about the constructor 
	void add_edge(node_id i, node_id j, captype cap, captype rev_cap);

	// Adds new edges 'SOURCE->i' and 'i->SINK' with corresponding weights.
	// Can be called multiple times for each node.
	// Weights can be negative.
	// NOTE: the number of such edges is not counted in edge_num_max.
	//       No internal memory is allocated by this call.
	void add_tweights(node_id i, tcaptype cap_source, tcaptype cap_sink);


	// Computes the maxflow. Can be called several times.
	// FOR DESCRIPTION OF reuse_trees, SEE mark_node().
	// FOR DESCRIPTION OF changed_list, SEE remove_from_changed_list().
	flowtype maxflow(bool reuse_trees = false, Block<node_id>* changed_list = NULL);

	// After the maxflow is computed, this function returns to which
	// segment the node 'i' belongs (GraphX<captype,tcaptype,flowtype>::SOURCE or GraphX<captype,tcaptype,flowtype>::SINK).
	//
	// Occasionally there may be several minimum cuts. If a node can be assigned
	// to both the source and the sink, then default_segm is returned.
	termtype what_segment(node_id i, termtype default_segm = SOURCE);



	//////////////////////////////////////////////
	//       ADVANCED INTERFACE FUNCTIONS       //
	//      (provide access to the graph)       //
	//////////////////////////////////////////////

//private:
	struct node;
	struct arc;

public:

	////////////////////////////
	// 1. Reallocating graph. //
	////////////////////////////

	// Removes all nodes and edges. 
	// After that functions add_node() and add_edge() must be called again. 
	//
	// Advantage compared to deleting GraphX and allocating it again:
	// no calls to delete/new (which could be quite slow).
	//
	// If the graph structure stays the same, then an alternative
	// is to go through all nodes/edges and set new residual capacities
	// (see functions below).
	void reset();

	////////////////////////////////////////////////////////////////////////////////
	// 2. Functions for getting pointers to arcs and for reading graph structure. //
	//    NOTE: adding new arcs may invalidate these pointers (if reallocation    //
	//    happens). So it's best not to add arcs while reading graph structure.   //
	////////////////////////////////////////////////////////////////////////////////

	// The following two functions return arcs in the same order that they
	// were added to the graph. NOTE: for each call add_edge(i,j,cap,cap_rev)
	// the first arc returned will be i->j, and the second j->i.
	// If there are no more arcs, then the function can still be called, but
	// the returned arc_id is undetermined.
	arc_id get_first_arc() { return 0; }
	arc_id get_next_arc(arc_id a) { return a+1; }

	// we can also query arcs insident to a specific node.
	arc_id get_first_incident_arc(node_id i); // may return NULL_ARC_ID
	arc_id get_next_incident_arc(arc_id a); // may return NULL_ARC_ID
	arc_id get_reverse_arc(arc_id a);

	// other functions for reading graph structure
	int get_node_num() { return node_num; }
	int get_arc_num() { return arc_num; }
	void get_arc_ends(arc_id a, node_id& i, node_id& j); // returns i,j so that a = i->j
	void get_arc_end(arc_id a, node_id& j);

	///////////////////////////////////////////////////
	// 3. Functions for reading residual capacities. //
	///////////////////////////////////////////////////

	// returns residual capacity of SOURCE->i minus residual capacity of i->SINK
	tcaptype get_trcap(node_id i); 
	// returns residual capacity of arc a
	captype get_rcap(arc_id a);

	/////////////////////////////////////////////////////////////////
	// 4. Functions for setting residual capacities.               //
	//    NOTE: If these functions are used, the value of the flow //
	//    returned by maxflow() will not be valid!                 //
	/////////////////////////////////////////////////////////////////

	void set_trcap(node_id i, tcaptype trcap); 
	void add_trcap(node_id i, tcaptype trcap); // similar to add_tweights(i, trcap, 0) but doesn't care about the value of total flow
	tcaptype add_trcap_X(node_id i, tcaptype trcap); // see implementation
	void set_rcap(arc_id a, captype rcap);

	void reset_flow(); // set flow to zero

	////////////////////////////////////////////////////////////////////
	// 5. Functions related to reusing trees & list of changed nodes. //
	////////////////////////////////////////////////////////////////////

	// If flag reuse_trees is true while calling maxflow(), then search trees
	// are reused from previous maxflow computation. 
	// In this case before calling maxflow() the user must
	// specify which parts of the graph have changed by calling mark_node():
	//   add_tweights(i),set_trcap(i)    => call mark_node(i)
	//   add_edge(i,j),set_rcap(a)       => call mark_node(i); mark_node(j)
	//
	// This option makes sense only if a small part of the graph is changed.
	// The initialization procedure goes only through marked nodes then.
	// 
	// mark_node(i) can either be called before or after graph modification.
	// Can be called more than once per node, but calls after the first one
	// do not have any effect.
	// 
	// NOTE: 
	//   - This option cannot be used in the first call to maxflow().
	//   - It is not necessary to call mark_node() if the change is ``not essential'',
	//     i.e. sign(trcap) is preserved for a node and zero/nonzero status is preserved for an arc.
	//   - To check that you marked all necessary nodes, you can call maxflow(false) after calling maxflow(true).
	//     If everything is correct, the two calls must return the same value of flow. (Useful for debugging).
	void mark_node(node_id i);

	// If changed_list is not NULL while calling maxflow(), then the algorithm
	// keeps a list of nodes which could potentially have changed their segmentation label.
	// Nodes which are not in the list are guaranteed to keep their old segmentation label (SOURCE or SINK).
	// Example usage:
	//
	//		typedef GraphX<int,int,int> G;
	//		G* g = new GraphX(nodeNum, edgeNum);
	//		Block<G::node_id>* changed_list = new Block<G::node_id>(128);
	//
	//		... // add nodes and edges
	//
	//		g->maxflow(); // first call should be without arguments
	//		for (int iter=0; iter<10; iter++)
	//		{
	//			... // change graph, call mark_node() accordingly
	//
	//			g->maxflow(true, changed_list);
	//			G::node_id* ptr;
	//			for (ptr=changed_list->ScanFirst(); ptr; ptr=changed_list->ScanNext())
	//			{
	//				G::node_id i = *ptr; assert(i>=0 && i<nodeNum);
	//				g->remove_from_changed_list(i);
	//				// do something with node i...
	//				if (g->what_segment(i) == G::SOURCE) { ... }
	//			}
	//			changed_list->Reset();
	//		}
	//		delete changed_list;
	//		
	// NOTE:
	//  - If changed_list option is used, then reuse_trees must be used as well.
	//  - In the example above, the user may omit calls g->remove_from_changed_list(i) and changed_list->Reset() in a given iteration.
	//    Then during the next call to maxflow(true, &changed_list) new nodes will be added to changed_list.
	//  - If the next call to maxflow() does not use option reuse_trees, then calling remove_from_changed_list()
	//    is not necessary. ("changed_list->Reset()" or "delete changed_list" should still be called, though).
	void remove_from_changed_list(node_id i) 
	{ 
		assert(i>=0 && i<node_num && nodes[i].is_in_changed_list); 
		nodes[i].is_in_changed_list = 0;
	}

	// Removes arcs that go from the source to the sink component, and returns the sum of their residual costs (>=0).
	// If there is such an arc (i,j), then this function must be called both for i and j.
	// Should be called after maxflow().
	captype remove_crossing_arcs(node_id i);

	void restore_crossing_arcs();

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
	
//private:
	// internal variables and functions

	struct node
	{
		arc_id		first;		// first outcoming arc

		arc_id		parent;	// node's parent
		node		*next;		// pointer to the next active node
								//   (or to itself if it is the last node in the list)
		int			TS;			// timestamp showing when DIST was computed
		int			DIST;		// distance to the terminal
		int			is_sink : 1;	// flag showing whether the node is in the source or in the sink tree (if parent!=NULL)
		int			is_marked : 1;	// set by mark_node()
		int			is_in_changed_list : 1; // set by maxflow if 

		tcaptype	tr_cap;		// if tr_cap > 0 then tr_cap is residual capacity of the arc SOURCE->node
								// otherwise         -tr_cap is residual capacity of the arc node->SINK 

	};

	struct arc
	{
		node_id		head;		// node the arc points to
		arc_id		next;		// next arc with the same originating node
	};

	struct nodeptr
	{
		node    	*ptr;
		nodeptr		*next;
	};
	static const int NODEPTR_BLOCK_SIZE = 128;

	node				*nodes;
	arc					*arcs;
	bool				external_arcs; // don't delete them
	captype*			r_caps;

	int					node_num, node_num_max, arc_num, arc_num_max;

	DBlock<nodeptr>		*nodeptr_block;

	void	(*error_function)(const char *);	// this function is called if a error occurs,
										// with a corresponding error message
										// (or exit(1) is called if it's NULL)

	flowtype			flow;		// total flow

	// reusing trees & list of changed pixels
	int					maxflow_iteration; // counter
	Block<node_id>		*changed_list;

	/////////////////////////////////////////////////////////////////////////

	node				*queue_first[2], *queue_last[2];	// list of active nodes
	nodeptr				*orphan_first, *orphan_last;		// list of pointers to orphans
	int					TIME;								// monotonically increasing global counter

	/////////////////////////////////////////////////////////////////////////

	// functions for processing active list
	void set_active(node *i);
	node *next_active();

	// functions for processing orphans list
	void set_orphan_front(node* i); // add to the beginning of the list
	void set_orphan_rear(node* i);  // add to the end of the list

	void add_to_changed_list(node* i);

	void maxflow_init();             // called if reuse_trees == false
	void maxflow_reuse_trees_init(); // called if reuse_trees == true
	void augment(arc_id middle_arc);
	void process_source_orphan(node *i);
	void process_sink_orphan(node *i);

	termtype what_segment(node* i, termtype default_segm = SOURCE);

	void test_consistency(node* current_node=NULL); // debug function
};











///////////////////////////////////////
// Implementation - inline functions //
///////////////////////////////////////



template <typename captype, typename tcaptype, typename flowtype> 
	inline typename GraphX<captype,tcaptype,flowtype>::node_id GraphX<captype,tcaptype,flowtype>::add_node(int num)
{
	assert(num > 0);
	assert(node_num+num <= node_num_max);


	int _i = node_num;
	node_num += num;
	node* i = nodes + _i;
	memset(i, 0, num*sizeof(node));
	for ( ; i<nodes+node_num; i++)
	{
		i->first = NULL_ARC_ID;
		i->parent = NULL_ARC_ID;
	}

	return _i;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::add_tweights(node_id i, tcaptype cap_source, tcaptype cap_sink)
{
	assert(i >= 0 && i < node_num);

	tcaptype delta = nodes[i].tr_cap;
	if (delta > 0) cap_source += delta;
	else           cap_sink   -= delta;
	flow += (cap_source < cap_sink) ? cap_source : cap_sink;
	nodes[i].tr_cap = cap_source - cap_sink;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::add_edge(node_id _i, node_id _j, captype cap, captype rev_cap)
{
	assert(_i >= 0 && _i < node_num);
	assert(_j >= 0 && _j < node_num);
	assert(_i != _j);
	assert(cap >= 0);
	assert(rev_cap >= 0);

	assert(arc_num < arc_num_max);

	arc_id a = arc_num ++;
	arc_id a_rev = arc_num ++;

	node* i = nodes + _i;
	node* j = nodes + _j;

	arcs[a].next = i -> first;
	i -> first = a;
	arcs[a_rev].next = j -> first;
	j -> first = a_rev;
	arcs[a].head = _j;
	arcs[a_rev].head = _i;
	r_caps[a] = cap;
	r_caps[a_rev] = rev_cap;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline typename GraphX<captype,tcaptype,flowtype>::arc_id GraphX<captype,tcaptype,flowtype>::get_first_incident_arc(node_id i)
{
	return nodes[i].first;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline typename GraphX<captype,tcaptype,flowtype>::arc_id GraphX<captype,tcaptype,flowtype>::get_next_incident_arc(arc_id a) 
{
	return arcs[a].next; 
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline typename GraphX<captype,tcaptype,flowtype>::arc_id GraphX<captype,tcaptype,flowtype>::get_reverse_arc(arc_id a) 
{
	return SISTER(a);
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::get_arc_ends(arc_id a, node_id& i, node_id& j)
{
	assert(a >= 0 && a < arc_num);
	i = arcs[SISTER(a)].head;
	j = arcs[a].head;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::get_arc_end(arc_id a, node_id& j)
{
	assert(a >= 0 && a < arc_num);
	j = arcs[a].head;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline tcaptype GraphX<captype,tcaptype,flowtype>::get_trcap(node_id i)
{
	assert(i>=0 && i<node_num);
	return nodes[i].tr_cap;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline captype GraphX<captype,tcaptype,flowtype>::get_rcap(arc_id a)
{
	assert(a >= 0 && a < arc_num);
	return r_caps[a];
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::set_trcap(node_id i, tcaptype trcap)
{
	assert(i>=0 && i<node_num); 
	nodes[i].tr_cap = trcap;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::add_trcap(node_id i, tcaptype trcap)
{
	assert(i>=0 && i<node_num);
	nodes[i].tr_cap += trcap;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline tcaptype GraphX<captype,tcaptype,flowtype>::add_trcap_X(node_id i, tcaptype trcap)
{
	assert(i>=0 && i<node_num);
	tcaptype v = (nodes[i].tr_cap < 0) ? nodes[i].tr_cap : 0;
	nodes[i].tr_cap += trcap;
	if (nodes[i].tr_cap < 0) v -= nodes[i].tr_cap;
	return v;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::set_rcap(arc_id a, captype rcap)
{
	assert(a >= 0 && a < arc_num);
	r_caps[a] = rcap;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::reset_flow()
{
	flow = 0;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline typename GraphX<captype,tcaptype,flowtype>::termtype GraphX<captype,tcaptype,flowtype>::what_segment(node_id i, termtype default_segm)
{
	if (nodes[i].parent != NULL_ARC_ID)
	{
		return (nodes[i].is_sink) ? SINK : SOURCE;
	}
	else
	{
		return default_segm;
	}
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline typename GraphX<captype,tcaptype,flowtype>::termtype GraphX<captype,tcaptype,flowtype>::what_segment(node* i, termtype default_segm)
{
	if (i->parent != NULL_ARC_ID)
	{
		return (i->is_sink) ? SINK : SOURCE;
	}
	else
	{
		return default_segm;
	}
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::mark_node(node_id _i)
{
	node* i = nodes + _i;
	if (!i->next)
	{
		/* it's not in the list yet */
		if (queue_last[1]) queue_last[1] -> next = i;
		else               queue_first[1]        = i;
		queue_last[1] = i;
		i -> next = i;
	}
	i->is_marked = 1;
}


template <typename captype, typename tcaptype, typename flowtype> 
	inline captype GraphX<captype,tcaptype,flowtype>::remove_crossing_arcs(node_id _i)
{
	node* i = nodes + _i;
	arc_id* a_ptr;
	captype sum = 0;
	if (what_segment(i) == SOURCE)
	{
		for (a_ptr=&i->first; *a_ptr!=NULL_ARC_ID; )
		{
			arc_id a = *a_ptr;
			if (what_segment(arcs[a].head) == SOURCE) a_ptr = &arcs[a].next;
			else
			{
				sum += r_caps[SISTER(a)];
				*a_ptr = arcs[a].next;
				arcs[a].next = a;
			}
		}
	}
	else
	{
		for (a_ptr=&i->first; *a_ptr!=NULL_ARC_ID; )
		{
			arc_id a = *a_ptr;
			if (what_segment(arcs[a].head) != SOURCE) a_ptr = &arcs[a].next;
			else
			{
				sum += r_caps[a];
				*a_ptr = arcs[a].next;
				arcs[a].next = a;
			}
		}
	}
	return sum;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::restore_crossing_arcs()
{
	arc_id a;
	for (a=0; a<arc_num; a+=2)
	{
		if (arcs[a].next != a && arcs[a+1].next != a+1) continue;
		node* i = nodes + arcs[a+1].head;
		node* j = nodes + arcs[a].head;
		arcs[a].next = i->first;
		i->first = a;
		arcs[a+1].next = j->first;
		j->first = a+1;
	}
}

#endif
