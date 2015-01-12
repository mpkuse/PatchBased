/* graph.cpp */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graphX.h"

/*
	special constants for node->parent. Duplicated in maxflow.cpp, both should match!
*/
#define TERMINAL ( (arc_id) -2 )		/* to terminal */
#define ORPHAN   ( (arc_id) -3 )		/* orphan */

template <typename captype, typename tcaptype, typename flowtype> 
	GraphX<captype, tcaptype, flowtype>::GraphX(int _node_num_max, int _edge_num_max, void (*err_function)(const char *))
	: node_num(0), node_num_max(_node_num_max), arc_num(0), arc_num_max(2*_edge_num_max),
	  external_arcs(false),
	  nodeptr_block(NULL),
	  error_function(err_function)
{
	nodes = (node*) malloc(node_num_max*sizeof(node));
	arcs = (arc*) malloc(arc_num_max*sizeof(arc));
	r_caps = (captype*) malloc(arc_num_max*sizeof(captype));
	if (!nodes || !arcs || !r_caps) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	maxflow_iteration = 0;
	flow = 0;
}

template <typename captype, typename tcaptype, typename flowtype> 
	GraphX<captype, tcaptype, flowtype>::GraphX(GraphX<captype, tcaptype, flowtype>* g)
	: node_num(g->node_num), node_num_max(g->node_num), arc_num(g->arc_num), arc_num_max(g->arc_num),
	  external_arcs(true),
	  nodeptr_block(NULL),
	  error_function(g->error_function)
{
	nodes = (node*) malloc(node_num_max*sizeof(node));
	arcs = g->arcs;
	r_caps = (captype*) malloc(arc_num_max*sizeof(captype));
	if (!nodes || !arcs || !r_caps) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	node_id i;

	memset(nodes, 0, node_num*sizeof(node));
	for (i=0; i<node_num; i++)
	{
		nodes[i].first = g->nodes[i].first;
		nodes[i].parent = NULL_ARC_ID;
	}

	maxflow_iteration = 0;
	flow = 0;
}

template <typename captype, typename tcaptype, typename flowtype> 
	GraphX<captype,tcaptype,flowtype>::~GraphX()
{
	if (nodeptr_block) 
	{ 
		delete nodeptr_block; 
		nodeptr_block = NULL; 
	}
	free(nodes);
	if (!external_arcs) free(arcs);
	free(r_caps);
}

template <typename captype, typename tcaptype, typename flowtype> 
	void GraphX<captype,tcaptype,flowtype>::reset()
{
	node_num = arc_num = 0;

	if (nodeptr_block) 
	{ 
		delete nodeptr_block; 
		nodeptr_block = NULL; 
	}

	maxflow_iteration = 0;
	flow = 0;
}


#include "instancesX.inc"
