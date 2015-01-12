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

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "Potts.h"

Potts::REAL POTTS_TOTAL_FLOW = 0; // for debugging purposes


Potts::Potts(int _node_num, int edge_num_max, int label_num, REAL (*_unaries_fn)(NodeId i, LabelId k))
	: node_num(_node_num), K(label_num),
	  unaries_fn(_unaries_fn), unlabeled_num(_node_num),
	  nodes2(NULL), edges2(NULL), messages(NULL), unaries2(NULL), label_info(NULL)
#ifdef REUSE_GRAPH_STRUCTURE
	  , graph_ptr(NULL)
#endif
{
	int i, j;

	graph = new G(node_num, edge_num_max);
	graph->add_node(node_num);
	nodes = new Node[node_num];
	solution = new LabelId[node_num];
	unaries = NULL;
	if (unaries_fn) unaries = NULL;
	else
	{
		unaries = new REAL[node_num*K];
		memset(unaries, 0, node_num*K*sizeof(REAL));
	}
	label_lists = new NodeId[K];

	// construct tree
	intervals = new Interval[2*K];
	if (K == 1) return;
	intervals[K].i = 0;
	intervals[K].j = K-1;
	intervals[K].parent = -1;
	intervals[K].depth = 0;
	interval_num = K+1;
	ConstructSubtree(K);

	LCA = new int*[2*K];
	LCA[0] = new int[(K+interval_num)*K];
	SIGN = LCA + K;
	SIGN[0] = LCA[0] + K*K;
	for (i=1; i<K; i++)
	{
		LCA[i] = LCA[i-1] + K;
		SIGN[i] = SIGN[i-1] + interval_num;
	}
	for (i=0; i<K; i++)
	{
		for (j=0; j<K; j++) LCA[i][j] = GetLCA(i, j);
		for (j=0; j<interval_num; j++) SIGN[i][j] = (j == i) ? 0 : GetSign(i, j);
	}
}

Potts::~Potts()
{
	delete [] LCA[0];
	delete [] LCA;
	if (graph) delete graph;
	if (label_info)
	{
		int k;
		for (k=0; k<K; k++)
		{
			if (label_info[k].graph) delete label_info[k].graph;
			if (label_info[k].old_labels) delete [] label_info[k].old_labels;
		}
		delete [] label_info;
	}
	if (nodes) delete [] nodes;
	if (nodes2) delete [] nodes2;
	if (edges2) delete [] edges2;
	if (messages) delete messages;
	delete [] solution;
	if (unaries) delete [] unaries;
	if (unaries2)
	{
		delete [] unaries2[0].ptr;
		delete [] unaries2;
	}
	delete [] intervals;
	delete [] label_lists;
}

void Potts::ConstructSubtree(IntervalId I)
{
	if (I < K)
	{
		intervals[I].children[0] = intervals[I].children[1] = -1;
		intervals[I].is_leaf = true;
		return;
	}
	intervals[I].is_leaf = false;

	int i = intervals[I].i, j = intervals[I].j, split = (i + j + 1) / 2;
	intervals[I].split = split;
	intervals[I].children[0]  = (i == split-1) ? i : (interval_num ++);
	intervals[I].children[1] = (j == split) ? j : (interval_num ++);
	Interval* L = &intervals[intervals[I].children[0]];
	Interval* R = &intervals[intervals[I].children[1]];
	L->i = i; L->j = split-1;
	R->i = split; R->j = j;
	L->depth = R->depth = intervals[I].depth + 1;
	L->parent = R->parent = I;
	ConstructSubtree(intervals[I].children[0]);
	ConstructSubtree(intervals[I].children[1]);
}

Potts::IntervalId Potts::GetLCA(LabelId i, LabelId j)
{
	while (intervals[i].depth > intervals[j].depth) i = intervals[i].parent;
	while (intervals[i].depth < intervals[j].depth) j = intervals[j].parent;
	while (i != j)
	{
		i = intervals[i].parent;
		j = intervals[j].parent;
	}
	return i;
}

int Potts::GetSign(LabelId i, IntervalId j)
{
	IntervalId i_prev = i;
	i = intervals[i].parent;
	while (intervals[i].depth > intervals[j].depth) { i_prev = i; i = intervals[i].parent; }
	while (intervals[i].depth < intervals[j].depth) j = intervals[j].parent;
	while (i != j)
	{
		i_prev = i;
		i = intervals[i].parent;
		j = intervals[j].parent;
	}
	return (intervals[i].children[0] == i_prev) ? +1 : -1;
}

Potts::REAL Potts::Solve(int iter_num_max, bool _reuse_graphs, unsigned options)
{
	reuse_graphs = _reuse_graphs;
	init_with_kovtun_solution = (options & INIT_WITH_KOVTUN_SOLUTION) ? true : false;
	init_with_kovtun_flow     = (options & INIT_WITH_KOVTUN_FLOW) ? true : false;
	maximize_persistency      = (options & MAXIMIZE_PERSISTENCY) ? true : false;
	potts_debug               = (options & POTTS_DEBUG) ? true : false;

	SolveKovtun();
	exec_r(printf( "[Potts] Done with Kovtun's Initial Solution\n" ));
	if (unlabeled_num == 0)
	{
		// compute cost
		current_cost = 0;
		G::arc_id a = graph->get_first_arc();
		int i, j, e, E = graph->get_arc_num() / 2;
		for (i=0; i<node_num; i++)
		{
			if (unaries) current_cost += unaries[i*K + solution[i]-K];
			else         current_cost += (*unaries_fn)(i, solution[i]-K);
		}
		for (e=0; e<E; e++)
		{
			graph->get_arc_ends(a, i, j);
			REAL cij = graph->get_rcap(a);
			a = graph->get_next_arc(a);
			REAL cji = graph->get_rcap(a);
			a = graph->get_next_arc(a);
			if (solution[i] != solution[j]) current_cost += (cij + cji) / 2;
		}
	}
	else
	{
		exec_r(printf( "[Potts] Solving expansion\n" ));
		SolveExpansion(iter_num_max);
		exec_r(printf( "\n[Potts] Done Optimizing....\n"));
	}
	return current_cost;
}


void Potts::InitKovtun()
{
	NodeId i;
	LabelId k;

	if (K < 2)
	{
		for (i=0; i<node_num; i++) solution[i] = 0;
		unlabeled_num = 0;
		return;
	}

	for (i=0; i<node_num; i++)
	{
		LabelId k_min;
		REAL X, Y; // min and second min
		if (unaries)
		{
			const REAL* u = unaries + i*K;
			if (u[0] < u[1]) { k_min = 0; X = u[0]; Y = u[1]; }
			else             { k_min = 1; X = u[1]; Y = u[0]; }
			for (k=2; k<K; k++)
			{
				if (Y <= u[k]) continue;
				if (X <= u[k]) Y = u[k];
				else { k_min = k; Y = X; X = u[k]; }
			}
		}
		else
		{
			REAL u;
			X = (*unaries_fn)(i, 0);
			Y = (*unaries_fn)(i, 1);
			if (X <= Y) k_min = 0;
			else { k_min = 1; u = X; X = Y; Y = u; }
			for (k=2; k<K; k++)
			{
				u = (*unaries_fn)(i, k);
				if (Y <= u) continue;
				if (X <= u) Y = u;
				else { k_min = k; Y = X; X = u; }
			}
		}
		nodes[i].k_min = k_min;
		nodes[i].f_delta = Y-X;
		nodes[i].lambda = 0;
	}
}

//////////////////////////////////////////////////////////////

void Potts::SolveKovtun()
{
	NodeId i;
	LabelId k;

	InitKovtun();
	if (K < 2) return;

	for (i=0; i<node_num; i++)
	{
		nodes[i].current_unary =  (nodes[i].k_min < intervals[K].split) ? nodes[i].f_delta : -nodes[i].f_delta;
		graph->add_trcap(i, nodes[i].current_unary);

		nodes[i].next = i+1;
	}
	nodes[node_num-1].next = -1;

	graph->maxflow();
	SolveKovtun(intervals[K], 0);

	if (init_with_kovtun_flow || maximize_persistency || potts_debug)
	{
		// read flows to 'messages'
		int e, E = graph->get_arc_num() / 2;
		messages = new REAL[E];
		G::arc_id a = graph->get_first_arc();
		for (e=0; e<E; e++)
		{
			REAL cij = graph->get_rcap(a);
			a = graph->get_next_arc(a);
			REAL cji = graph->get_rcap(a);
			a = graph->get_next_arc(a);
			messages[e] = (cij - cji) / 2;
		}
	}

	for (k=0; k<K; k++)
	{
		for (i=label_lists[k]; i>=0; i=nodes[i].next)
		{
			REAL u;
			if (k == nodes[i].k_min) u = - nodes[i].lambda + nodes[i].f_delta;
			else
			{
				if (unaries) u = - nodes[i].lambda + unaries[i*K+nodes[i].k_min] - unaries[i*K+k];
				else         u = - nodes[i].lambda + (*unaries_fn)(i, nodes[i].k_min) - (*unaries_fn)(i, k);
			}
			u -= nodes[i].current_unary;
			if (u != 0)
			{
				graph->add_trcap(i, u);
				graph->mark_node(i);
			}
		}

		graph->maxflow(true);

		for (i=label_lists[k]; i>=0; i=nodes[i].next)
		{
			if (graph->what_segment(i) == G::SOURCE) { solution[i] = K+k; unlabeled_num --; }
			else                                       solution[i] = k;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////

	if (maximize_persistency) MaximizePersistency();

	if (potts_debug)
	{
		if (!maximize_persistency) for (i=0; i<node_num; i++) nodes[i].flag = solution[i];
		TestKovtun();
	}

	if (!init_with_kovtun_flow && messages)
	{
		delete [] messages;
		messages = NULL;
	}

	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////


	delete [] nodes;
	nodes = NULL;
}

void Potts::MaximizePersistency()
{
	graph->restore_crossing_arcs();

	NodeId i, j;
	LabelId k;
	G::arc_id a, a_rev;
	NodeId first_unlabeled;
	NodeId* prev_ptr = &first_unlabeled;
	for (i=0; i<node_num; i++)
	{
		if (solution[i] >= K)
		{
			nodes[i].flag = solution[i];
			nodes[i].next = -solution[i];
		}
		else
		{
			*prev_ptr = i;
			prev_ptr = &nodes[i].next;
			if (unaries) nodes[i].lambda = unaries[i*K + nodes[i].k_min];
			else         nodes[i].lambda = (*unaries_fn)(i, + nodes[i].k_min);
			nodes[i].flag = -1;
		}
	}
	*prev_ptr = -1;

	for (k=0; k<K; k++)
	{
		for (i=first_unlabeled; i>=0; i=nodes[i].next)
		{
			if (nodes[i].flag >= 0) continue;

			REAL u;
			if (k == nodes[i].k_min) u = -nodes[i].f_delta;
			else
			{
				if (unaries) u = unaries[i*K + k] - nodes[i].lambda;
				else         u = (*unaries_fn)(i, k) - nodes[i].lambda;
			}

			LabelId _xi = solution[i] % K;
			for (a=graph->get_first_incident_arc(i); a; a=graph->get_next_incident_arc(a))
			{
				graph->get_arc_end(a, j);
				LabelId _xj = solution[j] % K;
				LabelId _xij = LCA[_xi][_xj];

				REAL cij = graph->get_rcap(a);
				a_rev = graph->get_reverse_arc(a);
				REAL cji = graph->get_rcap(a_rev);

				REAL message;
				if (k == _xij) message = (cij - cji) / 2;
				else
				{
					int e = graph->get_arc_cardinal_num(a);
					message = messages[e/2]*SIGN[k][_xij];
					if (e & 1) message = -message;
				}

				if (-nodes[j].next >= K && -nodes[j].next != k+K && 2*message > -(cij + cji))
				{
					u = 1;
					break;
				}

				u -= message;
			}
			if (u <= 0) continue;

			// traverse all nodes that can be reached from i
			NodeId i0 = i, first = i;
			nodes[i].flag = node_num;
			while ( first < node_num )
			{
				i = first;
				first = nodes[i].flag;

				LabelId _xi = solution[i] % K;
				for (a=graph->get_first_incident_arc(i); a; a=graph->get_next_incident_arc(a))
				{
					graph->get_arc_end(a, j);
					if (nodes[j].flag >= 0) continue;

					LabelId _xj = solution[j] % K;
					LabelId _xij = LCA[_xi][_xj];

					REAL cij = graph->get_rcap(a);
					a_rev = graph->get_reverse_arc(a);
					REAL cji = graph->get_rcap(a_rev);

					REAL message;
					if (k == _xij) message = (cij - cji) / 2;
					else
					{
						int e = graph->get_arc_cardinal_num(a);
						message = messages[e/2]*SIGN[k][_xij];
						if (e & 1) message = -message;
					}
					if (2*message < cij + cji)
					{
						nodes[j].flag = first;
						first = j;
					}
				}
			}
			i = i0;
		}

		// almost done. nodes that haven't been visited can be assigned label 'k'
		prev_ptr = &first_unlabeled;
		for (i=*prev_ptr; i>=0; i=*prev_ptr)
		{
			if (nodes[i].flag >= 0)
			{
				nodes[i].flag = -1;
				prev_ptr = &nodes[i].next;
			}
			else
			{
				nodes[i].flag = k + K;
				*prev_ptr = nodes[i].next;
				nodes[i].next = - (k+K);
			}
		}
	}

	for (i=0; i<node_num; i++)
	{
		if (solution[i] >= K) continue;
		int tmp = solution[i];
		if (nodes[i].flag >= 0) { solution[i] = nodes[i].flag; unlabeled_num --; }
		nodes[i].flag = tmp; // for TestKovtun()
	}
}


void Potts::SolveKovtun(Interval I, NodeId list_current)
{
	NodeId i;

	if (I.is_leaf)
	{
		label_lists[I.i] = list_current;
		return;
	}

	NodeId list[2];
	NodeId* i_ptr[2];

	i_ptr[0] = &list[0];
	i_ptr[1] = &list[1];
	for (i=list_current; i>=0; i=nodes[i].next)
	{
		int s = graph->what_segment(i); // 0 or 1
		REAL sum = graph->remove_crossing_arcs(i) / 2;
		nodes[i].lambda += sum;
		nodes[i].current_unary += (2*s-1)*sum;
		*i_ptr[s] = i;
		i_ptr[s] = &nodes[i].next;

		// set new unary term
		if ( !intervals[I.children[s]].is_leaf &&
		     I.i <= nodes[i].k_min &&
		     nodes[i].k_min <= I.j )
		{
			REAL new_unary_min, new_unary_max, new_unary;
			if (nodes[i].k_min < intervals[I.children[s]].split)
			{
				new_unary_min = - nodes[i].lambda + nodes[i].f_delta; new_unary_max = nodes[i].lambda + nodes[i].f_delta;
			}
			else
			{
				new_unary_min = - nodes[i].lambda - nodes[i].f_delta; new_unary_max = nodes[i].lambda - nodes[i].f_delta;
			}


			if      (new_unary_min > nodes[i].current_unary) new_unary = new_unary_min;
			else if (new_unary_max < nodes[i].current_unary) new_unary = new_unary_max;
			else continue;

			graph->add_trcap(i, new_unary - nodes[i].current_unary);
			nodes[i].current_unary = new_unary;
			graph->mark_node(i);
		}
	}
	*i_ptr[0] = *i_ptr[1] = -1;

	graph->maxflow(true);
	SolveKovtun(intervals[I.children[0]], list[0]);
	SolveKovtun(intervals[I.children[1]], list[1]);
}


void Potts::TestKovtun()
{
	NodeId i;
	LabelId k0, k;

	int e, E = graph->get_arc_num() / 2;
	Edge2* edges2 = new Edge2[E];
	REAL* M = new REAL[2*E];
	G::arc_id a = graph->get_first_arc();
	for (e=0; e<E; e++)
	{
		graph->get_arc_ends(a, edges2[e].i, edges2[e].j);
		REAL cij = graph->get_rcap(a);
		a = graph->get_next_arc(a);
		REAL cji = graph->get_rcap(a);
		a = graph->get_next_arc(a);
		edges2[e].lambda = (cij + cji) / 2;
		M[2*e] = messages[e];
		M[2*e+1] = (cij - cji) / 2;
	}

	REAL* rep = new REAL[node_num];
	for (k0=0; k0<K; k0++)
	{
		for (i=0; i<node_num; i++)
		{
			REAL u = 1000000;
			for (k=0; k<K; k++)
			{
				if (k == k0) continue;
				REAL uk;
				if (unaries) uk = unaries[i*K + k];
				else         uk = (*unaries_fn)(i, k);
				if (u > uk) u = uk;
			}
			if (unaries) rep[i] = unaries[i*K + k0] - u;
			else         rep[i] = (*unaries_fn)(i, k0) - u;
		}
		for (e=0; e<E; e++)
		{
			LabelId _xi = nodes[edges2[e].i].flag % K, _xj = nodes[edges2[e].j].flag % K;
			LabelId _xij = LCA[_xi][_xj];
			REAL message = (k0 == _xij) ? M[2*e+1] : M[2*e]*SIGN[k0][_xij];

			rep[edges2[e].i] -= message;
			rep[edges2[e].j] += message;
			if (  (solution[edges2[e].i]==k0+K && solution[edges2[e].j]!=k0+K && message != -edges2[e].lambda)
			   || (solution[edges2[e].i]!=k0+K && solution[edges2[e].j]==k0+K && message != +edges2[e].lambda) )
			{
				printf("!");
				getchar();
			}
		}
		for (i=0; i<node_num; i++)
		{
			if ((solution[i]!=k0+K && rep[i] < 0) || (solution[i]==k0+K && rep[i] > 0))
			{
				printf("?");
				getchar();
			}
		}
	}

	delete [] rep;
	delete [] edges2;
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


void Potts::SolveExpansion(int iter_num_max)
{
	int k, iter;

	InitExpansion();

	int unsuccessful_num = 0;
	for (iter=0; iter<iter_num_max; iter++)
	{exec_r(printf( "\n[Potts] Expansion Iteration (CurrCost=%d) %d of %d", current_cost, iter+1, iter_num_max ));
		for (k=0; k<K; k++)
		{exec_r( printf( "."); fflush(stdout); );
			DoExpansion(k);
			if (last_expansion_successful) unsuccessful_num = 0;
			else unsuccessful_num ++;
			if (unsuccessful_num >= K) return;
		}
	}
}

void Potts::InitExpansion()
{
	int i, j, k, e;

	cost0 = 0;
	int* mapping = new int[node_num];
	nodes2 = new Node2[unlabeled_num];
	j = 0;
	for (i=0; i<node_num; i++)
	{
		if (solution[i] < K)
		{
			nodes2[j].orig_id = i;
			nodes2[j].kovtun_label = solution[i];
			if (!init_with_kovtun_solution) solution[i] = 0;
			nodes2[j].solution = solution[i];
			mapping[i] = j;
			if (unaries) memcpy(unaries+j*K, unaries+i*K, K*sizeof(REAL));
			j ++;
		}
		else
		{
			mapping[i] = -1;
			if (unaries) cost0 += unaries[i*K + solution[i] - K];
			else         cost0 += (*unaries_fn)(i, solution[i] - K);
		}
	}

	edge_num2 = 0;
	int unaries_num2 = 0;
	if (!unaries)
	{
		unaries2 = new Unary2Ptr[unlabeled_num+1];
		for (i=0; i<=unlabeled_num; i++)
		{
			unaries2[i].num = 0;
			unaries2[i].k_cashed = -1;
		}
	}
	int E = graph->get_arc_num() / 2;
	G::arc_id a = graph->get_first_arc();
	for (e=0; e<E; e++)
	{
		graph->get_arc_ends(a, i, j);
		if (mapping[i]>=0 && mapping[j]>=0) edge_num2 ++;
		else if (unaries2)
		{
			if      (mapping[i]>=0) { unaries2[mapping[i]].num ++; unaries_num2 ++; }
			else if (mapping[j]>=0) { unaries2[mapping[j]].num ++; unaries_num2 ++; }
		}
		a = graph->get_next_arc(a);
		a = graph->get_next_arc(a);
	}
	edges2 = new Edge2[edge_num2];
	Unary2* unaries2_ptr;
	if (unaries2)
	{
		int prev_num = unaries2[0].num;
		unaries2_ptr = unaries2[0].ptr = new Unary2[unaries_num2];
		for (i=1; i<=unlabeled_num; i++)
		{
			Unary2* ptr = unaries2[i-1].ptr + prev_num;
			prev_num = unaries2[i].num;
			unaries2[i].ptr = ptr;
		}
	}
	REAL* messages_old = messages;
	if (messages)
	{
		messages = new REAL[2*edge_num2];
	}

	edge_num2 = 0;
	a = graph->get_first_arc();
	for (e=0; e<E; e++)
	{
		graph->get_arc_ends(a, i, j);
		REAL cij = graph->get_rcap(a);
		a = graph->get_next_arc(a);
		REAL cji = graph->get_rcap(a);
		a = graph->get_next_arc(a);
		REAL lambda = (cij + cji) / 2;

		if (mapping[i]>=0 && mapping[j]>=0)
		{
			edges2[edge_num2].i = mapping[i];
			edges2[edge_num2].j = mapping[j];
			edges2[edge_num2].lambda = lambda;
			if (messages_old)
			{
				messages[2*edge_num2]   = messages_old[e];
				messages[2*edge_num2+1] = cij - lambda;
			}
			edge_num2 ++;
		}
		else if (mapping[i]>=0)
		{
			cost0 += lambda;
			if (unaries) unaries[mapping[i]*K + solution[j] - K] -= lambda;
			else
			{
				unaries2[mapping[i]].ptr->k = solution[j] - K;
				unaries2[mapping[i]].ptr->lambda = lambda;
				unaries2[mapping[i]].ptr ++;
			}
		}
		else if (mapping[j]>=0)
		{
			cost0 += lambda;
			if (unaries) unaries[mapping[j]*K + solution[i] - K] -= lambda;
			else
			{
				unaries2[mapping[j]].ptr->k = solution[i] - K;
				unaries2[mapping[j]].ptr->lambda = lambda;
				unaries2[mapping[j]].ptr ++;
			}
		}
		else
		{
			if (solution[i] != solution[j]) cost0 += lambda;
		}
	}

	if (messages_old) delete [] messages_old;

	if (unaries2)
	{
		for (i=unlabeled_num-1; i>0; i--) unaries2[i].ptr = unaries2[i-1].ptr;
		unaries2[0].ptr = unaries2_ptr;
	}

	delete [] mapping;
	delete graph;
	graph = NULL;

	//////////////////////////////////

	if (!reuse_graphs)
	{
		graph = new G(unlabeled_num, edge_num2);
		graph->add_node(unlabeled_num);
		for (e=0; e<edge_num2; e++)
		{
			graph->add_edge(edges2[e].i, edges2[e].j, 0, 0);
		}
	}
	else
	{
		label_info = new LabelInfo[K];
		for (k=0; k<K; k++)
		{
			label_info[k].graph = NULL;
			label_info[k].old_labels = NULL;
		}
	}

	//////////////////////////////////

	ComputeCurrentCost();
}


void Potts::ComputeCurrentCost()
{
	int i, e;
	current_cost = cost0;
	if (unaries)
	{
		for (i=0; i<unlabeled_num; i++)
		{
			current_cost += unaries[i*K + nodes2[i].solution];
		}
	}
	else
	{
		for (i=0; i<unlabeled_num; i++)
		{
			current_cost += unaries2_fn(i, nodes2[i].solution);
		}
	}
	for (e=0; e<edge_num2; e++)
	{
		if (nodes2[edges2[e].i].solution != nodes2[edges2[e].j].solution) current_cost += edges2[e].lambda;
	}
}

Potts::REAL Potts::DoExpansion(LabelId k)
{
	if (unlabeled_num == 0) return current_cost;
	if (reuse_graphs)
	{
		if (!label_info[k].graph) DoExpansionReuse_first(k);
		else                      DoExpansionReuse_second(k);
	}
	else DoExpansionNaive(k);
	return current_cost;
}

void Potts::DoExpansionNaive(LabelId k)
{
	int i, j, e;
	LabelId xi, xj;
	LabelId _xi, _xj, _xij;
	G::arc_id a;

	// source = orig label
	// sink = label k

	for (i=0; i<unlabeled_num; i++)
	{
		xi = nodes2[i].solution;
		if (xi == k) graph->set_trcap(i, 0);
		else
		{
			if (unaries) graph->set_trcap(i, unaries[i*K + k] - unaries[i*K + xi]);
			else         graph->set_trcap(i, unaries2_fn(i, k) - unaries2_fn(i, xi, true));
		}
	}
	a = graph->get_first_arc();
	for (e=0; e<edge_num2; e++)
	{
		i = edges2[e].i; j = edges2[e].j;
		xi = nodes2[i].solution; xj = nodes2[j].solution;
		REAL lambda = edges2[e].lambda;

		if (xi == k || xj == k)
		{
			graph->set_rcap(a, 0);
			a = graph->get_next_arc(a);
			graph->set_rcap(a, 0);
			a = graph->get_next_arc(a);
			if (xi != k) graph->add_trcap(i, -lambda);
			else if (xj != k) graph->add_trcap(j, -lambda);
			continue;
		}

		REAL m;

		if (init_with_kovtun_flow)
		{
			_xi = nodes2[i].kovtun_label; _xj = nodes2[j].kovtun_label;
			_xij = LCA[_xi][_xj];

			m = (k == _xij) ? messages[2*e+1] : messages[2*e]*SIGN[k][_xij];

			if      (m < -lambda) m = -lambda;
			else if (m > +lambda) m = +lambda;
		}
		else
		{
			m = 0;
		}

		if (xi == xj)
		{
			graph->set_rcap(a, lambda-m);
			a = graph->get_next_arc(a);
			graph->set_rcap(a, lambda+m);
			a = graph->get_next_arc(a);
			graph->add_tweights(i, -m, 0);
			graph->add_tweights(j, +m, 0);
		}
		else
		{
			REAL cij = (REAL)((lambda - m) / 2);
			graph->set_rcap(a, cij);
			a = graph->get_next_arc(a);
			graph->set_rcap(a, lambda-cij);
			a = graph->get_next_arc(a);
			graph->add_tweights(i, cij-lambda, 0);
			graph->add_tweights(j, -cij, 0);
		}
	}

	REAL source_cost = 0;
	for (i=0; i<unlabeled_num; i++)
	{
		REAL ci = graph->get_trcap(i);
		if (ci < 0) source_cost -= ci;
	}
	graph->reset_flow();
	REAL flow = graph->maxflow();
	source_cost -= flow;

	POTTS_TOTAL_FLOW += flow; // for debugging

	if (source_cost > 1e-10)
	{
		for (i=0; i<unlabeled_num; i++)
		{
			if (graph->what_segment(i) == G::SINK) nodes2[i].solution = solution[nodes2[i].orig_id] = k;
		}
		current_cost -= source_cost;

		last_expansion_successful = true;
	}
	else last_expansion_successful = false;
}

void Potts::DoExpansionReuse_first(LabelId k)
{
	int i, j, e;
	LabelId xi, xj;

	// source = orig label
	// sink = label k

	GX* g;
#ifdef REUSE_GRAPH_STRUCTURE
	if (graph_ptr)
	{
		g = new GX(graph_ptr);
	}
	else
#endif
	{
		g = new GX(unlabeled_num, edge_num2);
		g->add_node(unlabeled_num);
	}
	label_info[k].graph = g;
	LabelId* old_labels = label_info[k].old_labels = new int[unlabeled_num];
	for (i=0; i<unlabeled_num; i++)
	{
		old_labels[i] = xi = nodes2[i].solution;
		if (xi == k) g->set_trcap(i, 0);
		else
		{
			if (unaries) g->set_trcap(i, unaries[i*K + k] - unaries[i*K + xi]);
			else         g->set_trcap(i, unaries2_fn(i, k) - unaries2_fn(i, xi, true));
		}
	}
	for (e=0; e<edge_num2; e++)
	{
		i = edges2[e].i; j = edges2[e].j;
		xi = nodes2[i].solution; xj = nodes2[j].solution;
		REAL lambda = edges2[e].lambda;

		if (xi == k || xj == k)
		{
			if (xi != k) g->add_trcap(i, -lambda);
			else if (xj != k) g->add_trcap(j, -lambda);

#ifdef REUSE_GRAPH_STRUCTURE
			if (graph_ptr) g->r_caps[2*e] = g->r_caps[2*e+1] = 0;
			else 
#endif
			g->add_edge(i, j, 0, 0);

			continue;
		}

		REAL m;

		if (init_with_kovtun_flow)
		{
			LabelId _xi = nodes2[i].kovtun_label, _xj = nodes2[j].kovtun_label;
			LabelId _xij = LCA[_xi][_xj];

			m = (k == _xij) ? messages[2*e+1] : messages[2*e]*SIGN[k][_xij];

			if      (m < -lambda) m = -lambda;
			else if (m > +lambda) m = +lambda;
		}
		else
		{
			m = 0;
		}

		if (xi == xj)
		{
			g->add_tweights(i, -m, 0);
			g->add_tweights(j, +m, 0);

#ifdef REUSE_GRAPH_STRUCTURE
			if (graph_ptr)
			{
				g->r_caps[2*e] = lambda-m;
				g->r_caps[2*e+1] = lambda+m;
			}
			else
#endif
			g->add_edge(i, j, lambda-m, lambda+m);
		}
		else
		{
			REAL cij = (REAL)((lambda - m) / 2);
			g->add_tweights(i, cij-lambda, 0);
			g->add_tweights(j, -cij, 0);

#ifdef REUSE_GRAPH_STRUCTURE
			if (graph_ptr)
			{
				g->r_caps[2*e] = cij;
				g->r_caps[2*e+1] = lambda-cij;
			}
			else
#endif
			g->add_edge(i, j, cij, lambda-cij);
		}
	}

#ifdef REUSE_GRAPH_STRUCTURE
	if (!graph_ptr) graph_ptr = g;
#endif

	label_info[k].source_cost = 0;
	for (i=0; i<unlabeled_num; i++)
	{
		REAL ci = g->get_trcap(i);
		if (ci < 0) label_info[k].source_cost -= ci;
	}
	g->reset_flow();
	REAL flow = g->maxflow();
	label_info[k].source_cost -= flow;

	POTTS_TOTAL_FLOW += flow; // for debugging

	if (label_info[k].source_cost > 1e-10)
	{
		for (i=0; i<unlabeled_num; i++)
		{
			if (g->what_segment(i) == GX::SINK) nodes2[i].solution = solution[nodes2[i].orig_id] = k;
		}
		current_cost -= label_info[k].source_cost;

		last_expansion_successful = true;
	}
	else last_expansion_successful = false;
}

void Potts::DoExpansionReuse_second(LabelId k)
{
	int i, j, e;
	LabelId xi, xj;

	// source = orig label
	// sink = label k

	GX* g = label_info[k].graph;
	LabelId* old_labels = label_info[k].old_labels;
	for (i=0; i<unlabeled_num; i++)
	{
		int xi_old = old_labels[i];
		xi = nodes2[i].solution;
		if (xi_old == xi) continue;

		REAL delta;
		if (xi_old == k)
		{
			if (unaries) delta = unaries[i*K + k] - unaries[i*K + xi];
			else         delta = unaries2_fn(i, k) - unaries2_fn(i, xi, true);
		}
		else if (xi == k)
		{
			if (unaries) delta = -unaries[i*K + k] + unaries[i*K + xi_old];
			else         delta = -unaries2_fn(i, k) + unaries2_fn(i, xi_old);
		}
		else
		{
			if (unaries) delta = unaries[i*K + xi_old] - unaries[i*K + xi];
			else         delta = unaries2_fn(i, xi_old) - unaries2_fn(i, xi, true);
		}
		label_info[k].source_cost += g->add_trcap_X(i, delta);
		g->mark_node(i);
	}
	for (e=0; e<edge_num2; e++)
	{
		i = edges2[e].i; j = edges2[e].j;
		int xi_old = old_labels[i], xj_old = old_labels[j];
		xi = nodes2[i].solution; xj = nodes2[j].solution;
		if (xi_old == xi && xj_old == xj) continue;

		REAL lambda = edges2[e].lambda;
		REAL delta00 = ( ((xi == xj) ? 0 : 1) - ((xi_old == xj_old) ? 0 : 1) ) * lambda;
		REAL delta01 = ( ((xi ==  k) ? 0 : 1) - ((xi_old ==      k) ? 0 : 1) ) * lambda;
		REAL delta10 = ( ((k  == xj) ? 0 : 1) - ((k      == xj_old) ? 0 : 1) ) * lambda;
		REAL ci = 0, cj = 0;

		if (delta00 != 0)
		{
			REAL c = delta00 / 2;
			delta01 -= c;
			delta10 -= delta00-c;
			ci = -c;
			cj = -(delta00-c);
		}
#ifdef REUSE_GRAPH_STRUCTURE
		REAL& rij = g->r_caps[2*e];
		REAL& rji = g->r_caps[2*e+1];
#else
		REAL& rij = g->arcs[2*e].r_cap;
		REAL& rji = g->arcs[2*e+1].r_cap;
#endif
		rij += delta01;
		rji += delta10;
		if (rij < 0)
		{
			rji += rij;
			ci -= rij;
			cj += rij;
			rij = 0;
		}
		else if (rji < 0)
		{
			rij += rji;
			ci += rji;
			cj -= rji;
			rji = 0;
		}
		label_info[k].source_cost += g->add_trcap_X(i, ci);
		label_info[k].source_cost += g->add_trcap_X(j, cj);
		g->mark_node(i);
		g->mark_node(j);
	}

	for (i=0; i<unlabeled_num; i++)
	{
		old_labels[i] = nodes2[i].solution;
	}

	g->reset_flow();
	REAL new_cost = g->maxflow(true);
	label_info[k].source_cost -= new_cost;

	POTTS_TOTAL_FLOW += new_cost; // for debugging

	if (label_info[k].source_cost > 1e-10)
	{
		for (i=0; i<unlabeled_num; i++)
		{
			if (g->what_segment(i) == GX::SINK) nodes2[i].solution = solution[nodes2[i].orig_id] = k;
		}
		current_cost -= label_info[k].source_cost;

		last_expansion_successful = true;
	}
	else last_expansion_successful = false;
}

