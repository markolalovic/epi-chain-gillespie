# src/gillespie_sim.py

import numpy as np
from copy import deepcopy

# main simulation function
def gillespie_sim(graph_g, model_params, time_max, rng, n_init_infected=1):
    r"""
    Simulates linear-chain SEIR on a weighted two-group contact network graph_g.
    UPDATE: I_pre block added

    - `model_params` is a dictionary with keys:
      - "beta_AA", "beta_BB", "beta_AB"

    - Now with global updates only:
      - after each event, we recompute all rates from scratch

    - Assumes graph_g has node attribute "group" set to "A" or "B"
    - See `example_contact_network()` in `src/contact_networks.py`
  
    - Node states are strings:
      "S", "R",
      "E:1", ..., "E:K1",
      "Ip:1", ..., "Ip:K_pre",  <- UPDATE
      "Ia:1", ..., "Ia:K2",
      "Is:1", ..., "Is:K3"

    - Infection intensity for susceptible v:
      lambda_v = \sum_{u~v} w_uv beta_{g(u)g(v)} 1_{u is infectious, i.e.: Ipre or Iasym or Isym} <- UPDATE
    
    - Exposed is a block of stages:
      - E:1, ..., E:K1
    
    - I_p:K_pre branches to Ia:1 or Is:1 with probability \alpha and 1 - \alpha
    
    - Detection is another bookkeeping process:
      - upon entry of `v` to "Is:1", do Bernoulli(p_detect)
      - and set `v` attribute to `detected = True` if success
    
    """
    graph = deepcopy(graph_g)
    node_order = get_node_order(graph)

    # initialize, aggregate counts and rates
    init_info = initialize_states(graph, rng, n_init_infected)
    
    counts = compute_counts(graph)
    total_rate = compute_all_rates(graph, model_params)

    # do checks
    do_checks(graph, counts, total_rate, prev_D=None)

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

        # mutate state, no local rate or count updates here
        event_type, parent, detected_inc = update_states(graph, v_selected, model_params, rng)

        # global recompute
        counts = compute_counts(graph)
        total_rate = compute_all_rates(graph, model_params)

        # do checks
        prev_D = counts_t[-1].get("D", 0) if counts_t else None
        do_checks(graph, counts, total_rate, prev_D = prev_D)

        events.append({
            "t": t,
            "event_type": event_type,
            "node": v_selected,
            "parent": parent,
            "detected_inc": detected_inc,
            "states": snapshot_states(graph, node_order),
        })

        times.append(t)
        counts_t.append(counts.copy())

    return {
        "times": times,
        "counts": counts_t,
        "events": events,
        "node_order": node_order,
        "init": init_info,
        "graph": graph
    }


## core functions 
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
        # inverse CDF selection
        if running_sum > target: 
            v_selected = i
            break

    if v_selected is None: # pick last vtx
        v_selected = list(graph.nodes)[-1]

    return v_selected, tau


def initialize_states(graph, rng, n_init_infected=1):
    """ 
    Samples n_init_infected chosen uniformly at random from group A.

    Before: set them to "Ia:1"
    UPDATE: With I_pre block: set them to "Ip:1"
    """
    nodes = list(graph.nodes())
    nodes_A = [i for i in nodes if graph.nodes[i]["group"] == "A"]

    for i in nodes:
        graph.nodes[i]["state"] = "S"
        graph.nodes[i]["rate"] = 0.0
        graph.nodes[i]["detected"] = False

    infected = []
    if n_init_infected > 0:
        # TODO: add guard: n_init_infected > len(nodes_A)
        infected = list(rng.choice(nodes_A, size=n_init_infected, replace=False))
        for i in infected:
            # graph.nodes[i]["state"] = "Ia:1"
            graph.nodes[i]["state"] = "Ip:1" # UPDATE
    
    # clean up output e.g. `np.str_('a2')` with .item() as in sample_parent
    infected = [inf.item() for inf in infected]

    return {"initial_infected_nodes": infected}

def update_states(graph, v_selected, model_params, rng):
    """
    This function only mutates node states, rates and counts are recomputed outside.

    UPDATE: adding node state mutation for staged I_pre block.

    Global update of state: no local neighbor updates to generalize more easily.

    Returns: event_type, parent, detected_inc
    """
    alpha = model_params["alpha"]
    p_detect = model_params["p_detect"]

    # TODO: perhaps rename also: K1 -> K_expo, K_2 -> K_asym, K3 -> K_sym
    K1 = int(model_params["K1"])
    K_pre = int(model_params["K_pre"])
    K2 = int(model_params["K2"])
    K3 = int(model_params["K3"])

    s = graph.nodes[v_selected]["state"]
    parent = None
    detected_inc = 0

    # S -> E:1 (infection)
    if s == "S":
        event_type = "infection"
        parent = sample_parent(graph, v_selected, model_params, rng)
        graph.nodes[v_selected]["state"] = "E:1"
        return event_type, parent, detected_inc
    
    # R has zero rate and should never be selected
    if s == "R":
        raise RuntimeError("Error: update_states selected a recovered node with zero rate!")

    # progress through blocks E, I_a, I_s
    # detection bookkeeping on progression to I_s:1
    block, k = parse_state(s)

    # E block: progression only, no branching here with added I_pre block
    if block == "E":
        if k < K1:
            event_type = "E_progress"
            graph.nodes[v_selected]["state"] = f"E:{k+1}"
        else:
            # k == K1
            event_type = "E_to_Ip"
            graph.nodes[v_selected]["state"] = "Ip:1"
        return event_type, parent, detected_inc

    # Ip block: progression and branching to Ia:1, Is:1
    if block == "Ip":
        if k < K_pre:
            event_type = "Ip_progress"
            graph.nodes[v_selected]["state"] = f"Ip:{k+1}"
        else:
            # k == K_pre            
            u = rng.random()
            if u < alpha:
                event_type = "Ip_to_Ia"
                graph.nodes[v_selected]["state"] = "Ia:1"
            else:
                event_type = "Ip_to_Is"
                graph.nodes[v_selected]["state"] = "Is:1"

                # detection bookkeeping upon entry to Is:1
                # if not yet detected run Bernoulli(p_detect)
                if (not graph.nodes[v_selected].get("detected", False)) and (rng.random() < p_detect):
                    graph.nodes[v_selected]["detected"] = True
                    detected_inc = 1

        return event_type, parent, detected_inc

    # I_a infectious progression -> R
    if block == "Ia":
        if k < K2:
            event_type = "Ia_progress"
            graph.nodes[v_selected]["state"] = f"Ia:{k+1}"
        else:
            event_type = "Ia_to_R"
            graph.nodes[v_selected]["state"] = "R"
        return event_type, parent, detected_inc

    # I_s infectious progression -> R
    if block == "Is":
        if k < K3:
            event_type = "Is_progress"
            graph.nodes[v_selected]["state"] = f"Is:{k+1}"
        else:
            event_type = "Is_to_R"
            graph.nodes[v_selected]["state"] = "R"
        return event_type, parent, detected_inc

    # else: throw unknown block state error
    raise ValueError(f"Error: update_states encountered unknown state block: state={s}, block={block}")


## global rate computation 
def compute_all_rates(graph, model_params):
    """ UPDATE: adding lambda_pre, K_pre to handle block I_pre = Ip """

    sigma = model_params["sigma"]
    lambda_pre = model_params["lambda_pre"]
    mu = model_params["mu"]

    K1 = int(model_params["K1"])
    K_pre = int(model_params["K_pre"])
    K2 = int(model_params["K2"])
    K3 = int(model_params["K3"])

    total_rate = 0.0

    for v in graph.nodes:
        s = graph.nodes[v]["state"]

        if s == "S":
            rate_v = infection_intensity(graph, v, model_params)

        elif s == "R":
            rate_v = 0.0

        else:
            block, k = parse_state(s)
            if block == "E":
                rate_v = sigma * K1
            elif block == "Ip":
                rate_v = lambda_pre * K_pre                
            elif block == "Ia":
                rate_v = mu * K2
            elif block == "Is":
                rate_v = mu * K3
            else:
                raise ValueError(f"Error: unknown state block in compute_all_rates: {s}")

        graph.nodes[v]["rate"] = rate_v
        total_rate += rate_v

    return total_rate


def infection_intensity(graph, v, model_params):
    lam = 0.0
    gv = graph.nodes[v]["group"]

    for u in graph.neighbors(v):
        if is_infectious(graph.nodes[u]["state"]):
            gu = graph.nodes[u]["group"]
            beta = beta_of_groups(gu, gv, model_params)
            w = graph[v][u].get("weight", 1.0)
            lam += w * beta

    return lam


def beta_of_groups(gu, gv, model_params):
    r""" Returns \beta_{g(u), g(v)}, where \beta_AB = \beta_BA. """
    if gu == "A" and gv == "A":
        return model_params["beta_AA"]
    if gu == "B" and gv == "B":
        return model_params["beta_BB"]
    return model_params["beta_AB"]  # = beta_BA (symmetric)


def sample_parent(graph, v_selected, model_params, rng):
    """
    Samples parent infector from infectious neighbors of v_selected
    with weights proportional to w_uv * beta_{g(u)g(v)}.    
    """
    gv = graph.nodes[v_selected]["group"]

    infectious_neighbors = []
    weights = []

    for u in graph.neighbors(v_selected):
        if is_infectious(graph.nodes[u]["state"]):
            gu = graph.nodes[u]["group"]
            beta = beta_of_groups(gu, gv, model_params)
            w = graph[v_selected][u].get("weight", 1.0)
            infectious_neighbors.append(u)
            weights.append(w * beta)

    if not infectious_neighbors:
        raise RuntimeError("Error infection event selected but no infectious neighbors.")

    weights = np.array(weights, dtype=float)
    s = weights.sum()
    if s <= 0.0:
        raise RuntimeError(f"Error: non-positive parent-weight sum {s}.")

    probs = weights / s
    parent = rng.choice(infectious_neighbors, p=probs)

    # for cleaner output, e.g., instead of np.str_('a1'), to get plain 'a1'
    parent = parent.item() if hasattr(parent, "item") else parent

    return parent


### additional helpers
def parse_state(state):
    # returns (block, k) given staged state "{state}:k"
    block, k = state.split(":")
    return block, int(k)

def is_infectious(state):
    # definition of infectious set \mathcal{I} = I_pre, I_asym, I_sym
    return state.startswith("Ip:") or state.startswith("Ia:") or state.startswith("Is:")


def compute_counts(graph):
    """
    Aggregate counts (substage-collapsed):
      S, E, I_pre, I_asym, I_sym, R plus bookkeeping D (cumulative detected).
    """
    counts = {"S": 0, "E": 0, "I_pre": 0, "I_asym": 0, "I_sym": 0, "R": 0, "D": 0}

    for v in graph.nodes:
        s = graph.nodes[v]["state"]

        if s == "S":
            counts["S"] += 1
        elif s == "R":
            counts["R"] += 1
        elif s.startswith("E:"):
            counts["E"] += 1
        elif s.startswith("Ip:"):
            counts["I_pre"] += 1            
        elif s.startswith("Ia:"):
            counts["I_asym"] += 1
        elif s.startswith("Is:"):
            counts["I_sym"] += 1
        else:
            raise ValueError(f"Error: unknown state block in counts: {s}")

        if graph.nodes[v].get("detected", False):
            counts["D"] += 1

    return counts


def do_checks(graph, counts, total_rate, prev_D=None, tol=1e-12):
    # total rate = sum of node rates
    s = 0.0
    for v in graph.nodes:
        s += graph.nodes[v].get("rate", 0.0)
    if abs(s - total_rate) > tol:
        raise RuntimeError(f"Error rate mismatch: sum(node.rate)={s}, total_rate={total_rate}")

    # compartments invariant
    N = len(graph.nodes)
    nsum = counts["S"] + counts["E"] + counts["I_pre"] + counts["I_asym"] + counts["I_sym"] + counts["R"]
    if nsum != N:
        raise RuntimeError(f"Error in block count invariant: S + E + Ip + Ia + Is + R = {nsum} neq N = {N}")

    # bounds on D and monotonicity of D
    D = counts["D"]
    if not (0 <= D <= N):
        raise RuntimeError(f"Error D out of bounds: D={D}, N={N}")
    if prev_D is not None and D < prev_D:
        raise RuntimeError(f"Error D decreased: prev_D={prev_D}, D={D}")


def get_node_order(graph):
    return sorted(graph.nodes())


def print_states(graph):
    node_order = get_node_order(graph)
    for i in node_order:
        print(i, graph.nodes[i]["state"], graph.nodes[i].get("rate", None))


def snapshot_states(graph, node_order):
    return [graph.nodes[i]["state"] for i in node_order]

