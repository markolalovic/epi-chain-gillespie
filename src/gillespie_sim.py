# src/gillespie_sim.py

import numpy as np
from copy import deepcopy

# main simulation function
def gillespie_sim(graph_g, model_params, time_max, rng, n_init_infected=1):
    """
    Simulates SIR dynamics on graph_g.

    - Adapted to track two-type dynamics
    - Assumes graph_g has node attribute "group" set to "A" or "B"
      - See `example_contact_network()` in `src/contact_networks.py`
    
    """
    graph = deepcopy(graph_g)
    node_order = get_node_order(graph)

    # initialize, aggregate counts and rates
    init_info = initialize_states(graph, rng, n_init_infected)
    counts = initialize_counts(graph)
    total_rate = initialize_rates(graph, model_params)

    # output
    t = 0.0
    times = [t]
    counts_t = [counts.copy()]
    events = []

    # Gillespie loop
    while t < time_max:
        # if no  more events possible
        if np.isclose(total_rate, 0.0):
            if t < time_max:
                times.append(time_max)
                counts_t.append(counts.copy())
            break
        
        v_selected, tau = draw_next_event(graph, total_rate, rng)

        # if next event occurs after horizon, then stop without state change
        if t + tau > time_max:
            times.append(time_max)
            counts_t.append(counts.copy())
            break

        t += tau
        event_type, counts, total_rate, parent = update_states(
            graph, v_selected, counts, total_rate, model_params, rng
        )

        events.append({
            't': t,
            'event_type': event_type,
            'node': v_selected,
            'parent': parent,
            'states': snapshot_states(graph, node_order),
        })

        times.append(t)
        counts_t.append(counts.copy())

    return {
        'times': times,
        'counts': counts_t,
        'events': events,
        'node_order': node_order,
        'init': init_info,
        'graph': graph,
    }


### core functions
def initialize_rates(graph, model_params):
    beta = model_params['beta']
    mu = model_params['mu']
    total_rate = 0.0

    for v in graph.nodes:
        s = graph.nodes[v]['state']

        if s == 'S':
            
            w_inf = 0.0 # sum of weights to infected neighbors
            for u in graph.neighbors(v):
                if graph.nodes[u]['state'] == 'I':
                    w = graph[v][u].get('weight', 1.0)
                    w_inf += w
            rate_v = beta * w_inf

        elif s == 'I':
            rate_v = mu

        else: # recovered
            rate_v = 0.0

        graph.nodes[v]['rate'] = rate_v
        total_rate += rate_v

    return total_rate


def draw_next_event(graph, total_rate, rng):
    u1, u2 = rng.random(2)

    # draw waiting time tau
    tau = -np.log(1.0 - u1) / total_rate

    # select v_selected by linear search
    target = u2 * total_rate
    running_sum = 0.0
    v_selected = None

    for i, attrs in graph.nodes(data=True):
        running_sum += attrs['rate']
        # if u2 == 0, then target = 0, we should not select the first vtx
        if running_sum > target: 
            v_selected = i
            break

    if v_selected is None: # pick last vtx
        v_selected = list(graph.nodes)[-1]

    return v_selected, tau


def update_states(graph, v_selected, counts, total_rate, model_params, rng):
    """ If infection of 'v': 
      - sample parent infector of 'v' from infected neighbors 'u' 
      with probability proportionally to weight ( = number of edges) of 'u'
    """
    beta = model_params['beta']
    mu = model_params['mu']

    # S -> I (infection)
    if graph.nodes[v_selected]['state'] == 'S':
        event_type = 'infection'
        parent = sample_parent(graph, v_selected, rng)

        counts['S'] -= 1
        counts['I'] += 1

        old_rate = graph.nodes[v_selected]['rate']
        graph.nodes[v_selected]['state'] = 'I'
        graph.nodes[v_selected]['rate'] = mu
        total_rate += mu - old_rate

        for j in graph.neighbors(v_selected):
            if graph.nodes[j]['state'] == 'S':
                w = graph[v_selected][j].get('weight', 1.0)
                dw = beta * w
                graph.nodes[j]['rate'] += dw
                total_rate += dw
    
    # NOTE: not 'S' we treat as 'I'
    # This holds if selected_i has positive rate and states are only S, I, R
    else:   # I -> R (recovery)
        event_type = 'recovery'
        parent = None # if recovery happens, there is no infector parent

        counts['I'] -= 1
        counts['R'] += 1

        old_rate = graph.nodes[v_selected]['rate']
        graph.nodes[v_selected]['state'] = 'R'
        graph.nodes[v_selected]['rate'] = 0.0
        total_rate -= old_rate

        for j in graph.neighbors(v_selected):
            if graph.nodes[j]['state'] == 'S':
                w = graph[v_selected][j].get('weight', 1.0)
                dw = beta * w
                graph.nodes[j]['rate'] -= dw
                total_rate -= dw

    return event_type, counts, total_rate, parent


def sample_parent(graph, v_selected, rng):
    """ Samples parent infector from infected neighbors of v_selected. """
    infected_neighbors = [u for u in graph.neighbors(v_selected) if graph.nodes[u]['state'] == 'I']

    # sample parent from infected neighbors with probability proportional to w(u_i, v)
    weights = np.array([graph[v_selected][u].get('weight', 1.0) for u in infected_neighbors], dtype=float)
    probs = weights / weights.sum()
    parent = rng.choice(infected_neighbors, p = probs)
    
    # for cleaner output, e.g., instead of np.str_('a1'), to get plain 'a1'
    parent = parent.item() if hasattr(parent, "item") else parent

    return parent


### additional helpers
def initialize_states(graph, rng, n_init_infected = 1):
    """ Sets n_init_infected chosen uniformly at random, 
    NOTE: only from group A.
    """
    nodes = list(graph.nodes())
    nodes_A = [i for i in nodes if graph.nodes[i]["group"] == "A"]

    for i in nodes:
        graph.nodes[i]["state"] = "S"
        graph.nodes[i]["rate"] = 0.0

    infected = []
    if n_init_infected > 0:
        infected = list(
            rng.choice(nodes_A, size=n_init_infected, replace=False)
        )
        for i in infected:
            graph.nodes[i]["state"] = "I"

    return {"initial_infected_nodes": infected}

def initialize_counts(graph):
    """ Counts of current compartments. """
    counts = {'S': 0, 'I': 0, 'R': 0}
    for i in graph.nodes:
        counts[graph.nodes[i]['state']] += 1
    return counts

def get_node_order(graph):
    """ Returns fixed node ordering. """
    return sorted(graph.nodes())

def print_states(graph):
    """ Prints index, state, intensities. """
    node_order = get_node_order(graph)
    for i in node_order:
        print(i, graph.nodes[i]['state'], graph.nodes[i]['rate'])

def snapshot_states(graph, node_order):
    """ Returns states in fixed node order e.g.: ['I', 'S', 'I', 'S', 'S']. """
    return [graph.nodes[i]['state'] for i in node_order]

