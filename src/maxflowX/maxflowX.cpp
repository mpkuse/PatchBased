/* maxflow.cpp */


#include <stdio.h>
#include "graphX.h"


/*
	special constants for node->parent. Duplicated in graph.cpp, both should match!
*/
#define TERMINAL ( (arc_id) -2 )		/* to terminal */
#define ORPHAN   ( (arc_id) -3 )		/* orphan */


#define INFINITE_D ((int)(((unsigned)-1)/2))		/* infinite distance to the terminal */

/***********************************************************************/

/*
	Functions for processing active list.
	i->next points to the next node in the list
	(or to i, if i is the last node in the list).
	If i->next is NULL iff i is not in the list.

	There are two queues. Active nodes are added
	to the end of the second queue and read from
	the front of the first queue. If the first queue
	is empty, it is replaced by the second queue
	(and the second queue becomes empty).
*/


template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::set_active(node *i)
{
	if (!i->next)
	{
		/* it's not in the list yet */
		if (queue_last[1]) queue_last[1] -> next = i;
		else               queue_first[1]        = i;
		queue_last[1] = i;
		i -> next = i;
	}
}

/*
	Returns the next active node.
	If it is connected to the sink, it stays in the list,
	otherwise it is removed from the list
*/
template <typename captype, typename tcaptype, typename flowtype> 
	inline typename GraphX<captype,tcaptype,flowtype>::node* GraphX<captype,tcaptype,flowtype>::next_active()
{
	node *i;

	while ( 1 )
	{
		if (!(i=queue_first[0]))
		{
			queue_first[0] = i = queue_first[1];
			queue_last[0]  = queue_last[1];
			queue_first[1] = NULL;
			queue_last[1]  = NULL;
			if (!i) return NULL;
		}

		/* remove it from the active list */
		if (i->next == i) queue_first[0] = queue_last[0] = NULL;
		else              queue_first[0] = i -> next;
		i -> next = NULL;

		/* a node in the list is active iff it has a parent */
		if (i->parent != NULL_ARC_ID) return i;
	}
}

/***********************************************************************/

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::set_orphan_front(node *i)
{
	nodeptr *np;
	i -> parent = ORPHAN;
	np = nodeptr_block -> New();
	np -> ptr = i;
	np -> next = orphan_first;
	orphan_first = np;
}

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::set_orphan_rear(node *i)
{
	nodeptr *np;
	i -> parent = ORPHAN;
	np = nodeptr_block -> New();
	np -> ptr = i;
	if (orphan_last) orphan_last -> next = np;
	else             orphan_first        = np;
	orphan_last = np;
	np -> next = NULL;
}

/***********************************************************************/

template <typename captype, typename tcaptype, typename flowtype> 
	inline void GraphX<captype,tcaptype,flowtype>::add_to_changed_list(node *i)
{
	if (changed_list && !i->is_in_changed_list)
	{
		node_id* ptr = changed_list->New();
		*ptr = (node_id)(i - nodes);
		i->is_in_changed_list = true;
	}
}

/***********************************************************************/

template <typename captype, typename tcaptype, typename flowtype> 
	void GraphX<captype,tcaptype,flowtype>::maxflow_init()
{
	node *i;

	queue_first[0] = queue_last[0] = NULL;
	queue_first[1] = queue_last[1] = NULL;
	orphan_first = NULL;

	TIME = 0;

	for (i=nodes; i<nodes+node_num; i++)
	{
		i -> next = NULL;
		i -> is_marked = 0;
		i -> is_in_changed_list = 0;
		i -> TS = TIME;
		if (i->tr_cap > 0)
		{
			/* i is connected to the source */
			i -> is_sink = 0;
			i -> parent = TERMINAL;
			set_active(i);
			i -> DIST = 1;
		}
		else if (i->tr_cap < 0)
		{
			/* i is connected to the sink */
			i -> is_sink = 1;
			i -> parent = TERMINAL;
			set_active(i);
			i -> DIST = 1;
		}
		else
		{
			i -> parent = NULL_ARC_ID;
		}
	}
}

template <typename captype, typename tcaptype, typename flowtype> 
	void GraphX<captype,tcaptype,flowtype>::maxflow_reuse_trees_init()
{
	node* i;
	node* j;
	node* queue = queue_first[1];
	arc_id a;
	nodeptr* np;

	queue_first[0] = queue_last[0] = NULL;
	queue_first[1] = queue_last[1] = NULL;
	orphan_first = orphan_last = NULL;

	TIME ++;

	while ((i=queue))
	{
		queue = i->next;
		if (queue == i) queue = NULL;
		i->next = NULL;
		i->is_marked = 0;
		set_active(i);

		if (i->tr_cap == 0)
		{
			if (i->parent != NULL_ARC_ID) set_orphan_rear(i);
			continue;
		}

		if (i->tr_cap > 0)
		{
			if (i->parent == NULL_ARC_ID || i->is_sink)
			{
				i->is_sink = 0;
				for (a=i->first; a!=NULL_ARC_ID; a=arcs[a].next)
				{
					j = nodes + arcs[a].head;
					if (!j->is_marked)
					{
						if (j->parent == SISTER(a)) set_orphan_rear(j);
						if (j->parent!=NULL_ARC_ID && j->is_sink && r_caps[a] > 0) set_active(j);
					}
				}
				add_to_changed_list(i);
			}
		}
		else
		{
			if (i->parent==NULL_ARC_ID || !i->is_sink)
			{
				i->is_sink = 1;
				for (a=i->first; a!=NULL_ARC_ID; a=arcs[a].next)
				{
					j = nodes + arcs[a].head;
					if (!j->is_marked)
					{
						if (j->parent == SISTER(a)) set_orphan_rear(j);
						if (j->parent!=NULL_ARC_ID && !j->is_sink && r_caps[SISTER(a)] > 0) set_active(j);
					}
				}
				add_to_changed_list(i);
			}
		}
		i->parent = TERMINAL;
		i -> TS = TIME;
		i -> DIST = 1;
	}

	//test_consistency();

	/* adoption */
	while ((np=orphan_first))
	{
		orphan_first = np -> next;
		i = np -> ptr;
		nodeptr_block -> Delete(np);
		if (!orphan_first) orphan_last = NULL;
		if (i->is_sink) process_sink_orphan(i);
		else            process_source_orphan(i);
	}
	/* adoption end */

	//test_consistency();
}

template <typename captype, typename tcaptype, typename flowtype> 
	void GraphX<captype,tcaptype,flowtype>::augment(arc_id middle_arc)
{
	node *i;
	arc_id a;
	tcaptype bottleneck;


	/* 1. Finding bottleneck capacity */
	/* 1a - the source tree */
	bottleneck = r_caps[middle_arc];
	for (i=nodes+arcs[SISTER(middle_arc)].head; ; i=nodes+arcs[a].head)
	{
		a = i -> parent;
		if (a == TERMINAL) break;
		if (bottleneck > r_caps[SISTER(a)]) bottleneck = r_caps[SISTER(a)];
	}
	if (bottleneck > i->tr_cap) bottleneck = i -> tr_cap;
	/* 1b - the sink tree */
	for (i=nodes+arcs[middle_arc].head; ; i=nodes+arcs[a].head)
	{
		a = i -> parent;
		if (a == TERMINAL) break;
		if (bottleneck > r_caps[a]) bottleneck = r_caps[a];
	}
	if (bottleneck > - i->tr_cap) bottleneck = - i -> tr_cap;


	/* 2. Augmenting */
	/* 2a - the source tree */
	r_caps[SISTER(middle_arc)] += bottleneck;
	r_caps[middle_arc] -= bottleneck;
	for (i=nodes+arcs[SISTER(middle_arc)].head; ; i=nodes+arcs[a].head)
	{
		a = i -> parent;
		if (a == TERMINAL) break;
		r_caps[a] += bottleneck;
		r_caps[SISTER(a)] -= bottleneck;
		if (r_caps[SISTER(a)] <= 0)
		{
			set_orphan_front(i); // add i to the beginning of the adoption list
		}
	}
	i -> tr_cap -= bottleneck;
	if (!i->tr_cap)
	{
		set_orphan_front(i); // add i to the beginning of the adoption list
	}
	/* 2b - the sink tree */
	for (i=nodes+arcs[middle_arc].head; ; i=nodes+arcs[a].head)
	{
		a = i -> parent;
		if (a == TERMINAL) break;
		r_caps[SISTER(a)] += bottleneck;
		r_caps[a] -= bottleneck;
		if (r_caps[a] <= 0)
		{
			set_orphan_front(i); // add i to the beginning of the adoption list
		}
	}
	i -> tr_cap += bottleneck;
	if (!i->tr_cap)
	{
		set_orphan_front(i); // add i to the beginning of the adoption list
	}


	flow += bottleneck;
}

/***********************************************************************/

template <typename captype, typename tcaptype, typename flowtype> 
	void GraphX<captype,tcaptype,flowtype>::process_source_orphan(node *i)
{
	node *j;
	arc_id a0, a0_min = NULL_ARC_ID, a;
	int d, d_min = INFINITE_D;

	/* trying to find a new parent */
	for (a0=i->first; a0!=NULL_ARC_ID; a0=arcs[a0].next)
	if (r_caps[SISTER(a0)] > 0)
	{
		j = nodes + arcs[a0].head;
		if (!j->is_sink && ((a=j->parent)!=NULL_ARC_ID))
		{
			/* checking the origin of j */
			d = 0;
			while ( 1 )
			{
				if (j->TS == TIME)
				{
					d += j -> DIST;
					break;
				}
				a = j -> parent;
				d ++;
				if (a==TERMINAL)
				{
					j -> TS = TIME;
					j -> DIST = 1;
					break;
				}
				if (a==ORPHAN) { d = INFINITE_D; break; }
				j = nodes + arcs[a].head;
			}
			if (d<INFINITE_D) /* j originates from the source - done */
			{
				if (d<d_min)
				{
					a0_min = a0;
					d_min = d;
				}
				/* set marks along the path */
				for (j=nodes+arcs[a0].head; j->TS!=TIME; j=nodes+arcs[j->parent].head)
				{
					j -> TS = TIME;
					j -> DIST = d --;
				}
			}
		}
	}

	if ((i->parent = a0_min) != NULL_ARC_ID)
	{
		i -> TS = TIME;
		i -> DIST = d_min + 1;
	}
	else
	{
		/* no parent is found */
		add_to_changed_list(i);

		/* process neighbors */
		for (a0=i->first; a0!=NULL_ARC_ID; a0=arcs[a0].next)
		{
			j = nodes+arcs[a0].head;
			if (!j->is_sink && ((a=j->parent)!=NULL_ARC_ID))
			{
				if (r_caps[SISTER(a0)] > 0) set_active(j);
				if (a!=TERMINAL && a!=ORPHAN && nodes+arcs[a].head==i)
				{
					set_orphan_rear(j); // add j to the end of the adoption list
				}
			}
		}
	}
}

template <typename captype, typename tcaptype, typename flowtype> 
	void GraphX<captype,tcaptype,flowtype>::process_sink_orphan(node *i)
{
	node *j;
	arc_id a0, a0_min = NULL_ARC_ID, a;
	int d, d_min = INFINITE_D;

	/* trying to find a new parent */
	for (a0=i->first; a0!=NULL_ARC_ID; a0=arcs[a0].next)
	if (r_caps[a0] > 0)
	{
		j = nodes + arcs[a0].head;
		if (j->is_sink && ((a=j->parent)!=NULL_ARC_ID))
		{
			/* checking the origin of j */
			d = 0;
			while ( 1 )
			{
				if (j->TS == TIME)
				{
					d += j -> DIST;
					break;
				}
				a = j -> parent;
				d ++;
				if (a==TERMINAL)
				{
					j -> TS = TIME;
					j -> DIST = 1;
					break;
				}
				if (a==ORPHAN) { d = INFINITE_D; break; }
				j = nodes + arcs[a].head;
			}
			if (d<INFINITE_D) /* j originates from the sink - done */
			{
				if (d<d_min)
				{
					a0_min = a0;
					d_min = d;
				}
				/* set marks along the path */
				for (j=nodes+arcs[a0].head; j->TS!=TIME; j=nodes+arcs[j->parent].head)
				{
					j -> TS = TIME;
					j -> DIST = d --;
				}
			}
		}
	}

	if ((i->parent = a0_min) != NULL_ARC_ID)
	{
		i -> TS = TIME;
		i -> DIST = d_min + 1;
	}
	else
	{
		/* no parent is found */
		add_to_changed_list(i);

		/* process neighbors */
		for (a0=i->first; a0!=NULL_ARC_ID; a0=arcs[a0].next)
		{
			j = nodes + arcs[a0].head;
			if (j->is_sink && ((a=j->parent)!=NULL_ARC_ID))
			{
				if (r_caps[a0] > 0) set_active(j);
				if (a!=TERMINAL && a!=ORPHAN && nodes+arcs[a].head==i)
				{
					set_orphan_rear(j); // add j to the end of the adoption list
				}
			}
		}
	}
}

/***********************************************************************/

template <typename captype, typename tcaptype, typename flowtype> 
	flowtype GraphX<captype,tcaptype,flowtype>::maxflow(bool reuse_trees, Block<node_id>* _changed_list)
{
	node *i, *j, *current_node = NULL;
	arc_id a;
	nodeptr *np, *np_next;

	if (!nodeptr_block)
	{
		nodeptr_block = new DBlock<nodeptr>(NODEPTR_BLOCK_SIZE, error_function);
	}

	changed_list = _changed_list;
	if (maxflow_iteration == 0 && reuse_trees) { if (error_function) (*error_function)("reuse_trees cannot be used in the first call to maxflow()!"); exit(1); }
	if (changed_list && !reuse_trees) { if (error_function) (*error_function)("changed_list cannot be used without reuse_trees!"); exit(1); }

	if (reuse_trees) maxflow_reuse_trees_init();
	else             maxflow_init();

	// main loop
	while ( 1 )
	{
		// test_consistency(current_node);

		if ((i=current_node))
		{
			i -> next = NULL; /* remove active flag */
			if (i->parent == NULL_ARC_ID) i = NULL;
		}
		if (!i)
		{
			if (!(i = next_active())) break;
		}

		/* growth */
		if (!i->is_sink)
		{
			/* grow source tree */
			for (a=i->first; a!=NULL_ARC_ID; a=arcs[a].next)
			if (r_caps[a] > 0)
			{
				j = nodes + arcs[a].head;
				if (j->parent == NULL_ARC_ID)
				{
					j -> is_sink = 0;
					j -> parent = SISTER(a);
					j -> TS = i -> TS;
					j -> DIST = i -> DIST + 1;
					set_active(j);
					add_to_changed_list(j);
				}
				else if (j->is_sink) break;
				else if (j->TS <= i->TS &&
				         j->DIST > i->DIST)
				{
					/* heuristic - trying to make the distance from j to the source shorter */
					j -> parent = SISTER(a);
					j -> TS = i -> TS;
					j -> DIST = i -> DIST + 1;
				}
			}
		}
		else
		{
			/* grow sink tree */
			for (a=i->first; a!=NULL_ARC_ID; a=arcs[a].next)
			if (r_caps[SISTER(a)] > 0)
			{
				j = nodes + arcs[a].head;
				if (j->parent == NULL_ARC_ID)
				{
					j -> is_sink = 1;
					j -> parent = SISTER(a);
					j -> TS = i -> TS;
					j -> DIST = i -> DIST + 1;
					set_active(j);
					add_to_changed_list(j);
				}
				else if (!j->is_sink) { a = SISTER(a); break; }
				else if (j->TS <= i->TS &&
				         j->DIST > i->DIST)
				{
					/* heuristic - trying to make the distance from j to the sink shorter */
					j -> parent = SISTER(a);
					j -> TS = i -> TS;
					j -> DIST = i -> DIST + 1;
				}
			}
		}

		TIME ++;

		if (a != NULL_ARC_ID)
		{
			i -> next = i; /* set active flag */
			current_node = i;

			/* augmentation */
			augment(a);
			/* augmentation end */

			/* adoption */
			while ((np=orphan_first))
			{
				np_next = np -> next;
				np -> next = NULL;

				while ((np=orphan_first))
				{
					orphan_first = np -> next;
					i = np -> ptr;
					nodeptr_block -> Delete(np);
					if (!orphan_first) orphan_last = NULL;
					if (i->is_sink) process_sink_orphan(i);
					else            process_source_orphan(i);
				}

				orphan_first = np_next;
			}
			/* adoption end */
		}
		else current_node = NULL;
	}
	// test_consistency();

	if (!reuse_trees || (maxflow_iteration % 64) == 0)
	{
		delete nodeptr_block; 
		nodeptr_block = NULL; 
	}

	maxflow_iteration ++;
	return flow;
}

/***********************************************************************/


template <typename captype, typename tcaptype, typename flowtype> 
	void GraphX<captype,tcaptype,flowtype>::test_consistency(node* current_node)
{
	node *i;
	arc_id a;
	int r;
	int num1 = 0, num2 = 0;

	// test whether all nodes i with i->next!=NULL are indeed in the queue
	for (i=nodes; i<nodes+node_num; i++)
	{
		if (i->next || i==current_node) num1 ++;
	}
	for (r=0; r<3; r++)
	{
		i = (r == 2) ? current_node : queue_first[r];
		if (i)
		for ( ; ; i=i->next)
		{
			num2 ++;
			if (i->next == i)
			{
				if (r<2) assert(i == queue_last[r]);
				else     assert(i == current_node);
				break;
			}
		}
	}
	assert(num1 == num2);

	for (i=nodes; i<nodes+node_num; i++)
	{
		// test whether all edges in seach trees are non-saturated
		if (i->parent == NULL_ARC_ID) {}
		else if (i->parent == ORPHAN) {}
		else if (i->parent == TERMINAL)
		{
			if (!i->is_sink) assert(i->tr_cap > 0);
			else             assert(i->tr_cap < 0);
		}
		else
		{
			if (!i->is_sink) assert (r_caps[SISTER(i->parent)] > 0);
			else             assert (r_caps[i->parent] > 0);
		}
		// test whether passive nodes in search trees have neighbors in
		// a different tree through non-saturated edges
		if (i->parent!=NULL_ARC_ID && !i->next)
		{
			if (!i->is_sink)
			{
				assert(i->tr_cap >= 0);
				for (a=i->first; a!=NULL_ARC_ID; a=arcs[a].next)
				{
					if (r_caps[a] > 0) assert(nodes[arcs[a].head].parent!=NULL_ARC_ID && !nodes[arcs[a].head].is_sink);
				}
			}
			else
			{
				assert(i->tr_cap <= 0);
				for (a=i->first; a!=NULL_ARC_ID; a=arcs[a].next)
				{
					if (r_caps[SISTER(a)] > 0) assert(nodes[arcs[a].head].parent!=NULL_ARC_ID && nodes[arcs[a].head].is_sink);
				}
			}
		}
		// test marking invariants
		if (i->parent!=NULL_ARC_ID && i->parent!=ORPHAN && i->parent!=TERMINAL)
		{
			assert(i->TS <= nodes[arcs[i->parent].head].TS);
			if (i->TS == nodes[arcs[i->parent].head].TS) assert(i->DIST > nodes[arcs[i->parent].head].DIST);
		}
	}
}

#include "instancesX.inc"
