# src/gillespie_sim.py

import numpy as np
from copy import deepcopy

# main simulation function
def gillespie_sim(
        graph_g, 
        model_params, 
        time_max, 
        rng, 
        n_init_infected=1,
        wave=None,
        detection_params=None,        
    ):
    r"""
    Simulates linear-chain SEIR on a weighted two-group contact network graph_g.

    UPDATE: added detection layer for testing and observation process

    If wave is:
      - None: no detection layer
      - "wave2": regular testing using PCR
      - "omicron": reduced testing using LFD

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
    
    UPDATE: Returned dict contains out["detected"] with the simulated 
    cumulative detected resident cases D(t)
    """
    graph = deepcopy(graph_g)
    node_order = get_node_order(graph)

    # initialize disease states
    init_info = initialize_states(graph, rng, n_init_infected)
    
    # initialize disease counts and rates
    counts = compute_counts(graph)
    total_rate = compute_all_rates(graph, model_params)

    # initialize detection layer
    detection_params = build_detection_params(wave, detection_params)
    detected = initialize_detected_output(graph, detection_params)
    next_routine_test_time = next_initial_routine_test_time(detection_params, rng)

    outbreak_test_times = []
    outbreak_state = {
        "cooldown_until": -np.inf,
        "n_triggered": 0,
    }

    # checks
    do_checks(graph, counts, total_rate)
    do_detection_checks(graph, detected, detection_params)

    # output
    t = 0.0
    times = [t]
    counts_t = [counts.copy()]
    events = []

    # Gillespie loop and testing loop
    while t < time_max:
        # if no more events possible
        # no disease events remain, stop disease process
        # NOTE: We assume S and R never test positive
        # so further testing cannot create new detected cases
        if np.isclose(total_rate, 0.0):
            if t < time_max:
                times.append(time_max)
                counts_t.append(counts.copy())

                detected["times"].append(float(time_max))
                detected["counts"].append(compute_detected_counts(graph, detection_params))
            break
        
        # draw candidate time of next disease event
        tau = draw_waiting_time(total_rate, rng)
        next_disease_time = t + tau

        next_outbreak_time = next_outbreak_test_time(outbreak_test_times)
        next_test_time = min(next_routine_test_time, next_outbreak_time)

        # if the next routine testing time happens first,
        # process the testing event
        if (
            detection_enabled(detection_params)
            and next_test_time <= time_max
            and next_test_time <= next_disease_time
        ):
            t = next_test_time

            if next_outbreak_time <= next_routine_test_time:
                # outbreak testing round
                pop_next_outbreak_test_time(outbreak_test_times)

                test_event = process_outbreak_testing_event(
                    graph=graph,
                    t=t,
                    detection_params=detection_params,
                    rng=rng,
                )

            else:
                # routine testing round
                test_event = process_routine_testing_event(
                    graph=graph,
                    t=t,
                    detection_params=detection_params,
                    rng=rng,
                )

                period = float(detection_params["routine_testing"]["period"])
                next_routine_test_time += period

            if test_event is not None:
                append_detected_event(detected, graph, detection_params, test_event)
                do_detection_checks(graph, detected, detection_params)

                # routine or symptom test can trigger an outbreak protocol
                if test_event.get("n_counted_positive", 0) > 0:
                    maybe_trigger_outbreak_testing(
                        t=t,
                        detected=detected,
                        detection_params=detection_params,
                        rng=rng,
                        outbreak_state=outbreak_state,
                        outbreak_test_times=outbreak_test_times,
                    )

            # no disease state changed, and testing does not affect rates
            continue

        # if next disease event occurs after horizon, stop without state change
        if next_disease_time > time_max:
            times.append(time_max)
            counts_t.append(counts.copy())

            detected["times"].append(float(time_max))
            detected["counts"].append(compute_detected_counts(graph, detection_params))
            break

        # otherwise process disease event
        t = next_disease_time
        v_selected = draw_event_node(graph, total_rate, rng)

        # mutate disease state only
        event_type, parent = update_states(graph, v_selected, model_params, rng)

        # symptom-triggered testing, if a resident enters Is:1
        if detection_enabled(detection_params) and event_type == "Ip_to_Is":
            symptom_event = process_symptom_testing_event(
                graph=graph,
                t=t,
                v=v_selected,
                detection_params=detection_params,
                rng=rng,
            )

            if symptom_event is not None:
                append_detected_event(detected, graph, detection_params, symptom_event)
                do_detection_checks(graph, detected, detection_params)

                if symptom_event.get("n_counted_positive", 0) > 0:
                    maybe_trigger_outbreak_testing(
                        t=t,
                        detected=detected,
                        detection_params=detection_params,
                        rng=rng,
                        outbreak_state=outbreak_state,
                        outbreak_test_times=outbreak_test_times,
                    )
        
        # global recompute of disease rates and disease counts
        counts = compute_counts(graph)
        total_rate = compute_all_rates(graph, model_params)

        # checks
        do_checks(graph, counts, total_rate)

        events.append({
            "t": float(t),
            "event_type": event_type,
            "node": v_selected,
            "parent": parent,
            "states": snapshot_states(graph, node_order),
        })

        times.append(float(t))
        counts_t.append(counts.copy())

    return {
        "times": times,
        "counts": counts_t,
        "events": events,
        "node_order": node_order,
        "init": init_info,
        "graph": graph,
        "detected": detected,
    }


## core functions 

## UPDATE: split draw_next_event() into two functions
# - draw waiting time
# - draw event node
def draw_waiting_time(total_rate, rng):
    u = rng.random()
    return -np.log(1.0 - u) / total_rate

def draw_event_node(graph, total_rate, rng):
    u = rng.random()
    target = u * total_rate
    running_sum = 0.0
    v_selected = None

    for i, attrs in graph.nodes(data=True):
        running_sum += attrs["rate"]
        if running_sum > target:
            v_selected = i
            break

    if v_selected is None:
        v_selected = list(graph.nodes)[-1]

    return v_selected

def draw_next_event(graph, total_rate, rng):
    tau = draw_waiting_time(total_rate, rng)
    v_selected = draw_event_node(graph, total_rate, rng)
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

        # observation attributes
        # UPDATE: flag if already contributed to observed resident (group B) case count
        graph.nodes[i]["detected"] = False
        graph.nodes[i]["detected_time"] = None
        graph.nodes[i]["detected_by"] = None
        graph.nodes[i]["detected_test_type"] = None        

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
    This function mutates disease state only, rates and counts are recomputed outside.

    UPDATE: detection / testing is handled separatelly by observation layer

    UPDATE: adding node state mutation for staged I_pre block.

    Global update of state: no local neighbor updates to generalize more easily.

    Returns: event_type, parent
    """
    alpha = model_params["alpha"]

    K1 = int(model_params["K1"])
    K_pre = int(model_params["K_pre"])
    K2 = int(model_params["K2"])
    K3 = int(model_params["K3"])

    s = graph.nodes[v_selected]["state"]
    parent = None
    
    # S -> E:1 (infection)
    if s == "S":
        event_type = "infection"
        parent = sample_parent(graph, v_selected, model_params, rng)
        graph.nodes[v_selected]["state"] = "E:1"
        return event_type, parent
    
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
        return event_type, parent

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
        return event_type, parent

    # I_a infectious progression -> R
    if block == "Ia":
        if k < K2:
            event_type = "Ia_progress"
            graph.nodes[v_selected]["state"] = f"Ia:{k+1}"
        else:
            event_type = "Ia_to_R"
            graph.nodes[v_selected]["state"] = "R"
        return event_type, parent

    # I_s infectious progression -> R
    if block == "Is":
        if k < K3:
            event_type = "Is_progress"
            graph.nodes[v_selected]["state"] = f"Is:{k+1}"
        else:
            event_type = "Is_to_R"
            graph.nodes[v_selected]["state"] = "R"
        return event_type, parent

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
      S, E, I_pre, I_asym, I_sym, R

    UPDATE: Detected cases are now handled separately in out["detected"]
    """
    counts = {"S": 0, "E": 0, "I_pre": 0, "I_asym": 0, "I_sym": 0, "R": 0}

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

    return counts


def do_checks(graph, counts, total_rate, tol=1e-12):
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

def get_node_order(graph):
    return sorted(graph.nodes())


def print_states(graph):
    node_order = get_node_order(graph)
    for i in node_order:
        print(i, graph.nodes[i]["state"], graph.nodes[i].get("rate", None))


def snapshot_states(graph, node_order):
    return [graph.nodes[i]["state"] for i in node_order]



# UPDATE: wave specific detection defaults
# NOTE: placeholders
# TODO: this should not be calibrated, we have too many parameters
# TODO: think about how to reduce the number of parameters
DEFAULT_DETECTION_PARAMS = {
    "wave2": {
        "enabled": True,
        "wave": "wave2",

        # compare to resident detected cases only
        "observed_groups": ["B"],

        # routine resident PCR testing, roughly monthly
        "routine_testing": {
            "enabled": True,
            "period": 28.0,
            "first_test_time": None,  # None means random phase in [0, period)
            "test_type": "PCR",
            "coverage": 0.90,
            "sensitivity": {
                "E": 0.25,
                "Ip": 0.85,
                "Ia": 0.85,
                "Is": 0.95,
            },
        },

        # symptomatic residents are tested promtly
        "symptom_testing": {
            "enabled": True,
            "test_type": "PCR",
            "coverage": 0.95,
            "sensitivity": {
                "Is": 0.95,
            },
        },

        "outbreak_testing": {
            "enabled": True,

            # triggered based on recently detected resident positives
            "trigger_count": 2,
            "trigger_window": 14.0,

            # stop outbreak episode
            "cooldown": 14.0,

            # resident-only PCR testing
            "test_groups": ["B"],
            "test_type": "PCR",
            "coverage": 0.95,

            # round 1, repeat after 4 - 7 days.
            "rounds": [0.0, ("uniform", 4.0, 7.0)],

            "sensitivity": {
                "E": 0.25,
                "Ip": 0.85,
                "Ia": 0.85,
                "Is": 0.95,
            },
        }
    },
    "omicron": {
        "enabled": True,
        "wave": "omicron",

        # compare to resident detected cases only
        "observed_groups": ["B"],

        # LFD screening, roughly 1 - 2 times per week
        "routine_testing": {
            "enabled": True,
            "period": 3.5,
            "first_test_time": None,
            "test_type": "LFD",
            "coverage": 0.60,
            "sensitivity": {
                "E": 0.05,
                "Ip": 0.45,
                "Ia": 0.45,
                "Is": 0.65,
            },
        },

        # symptomatic residents often tested by LFD
        "symptom_testing": {
            "enabled": True,
            "test_type": "LFD",
            "coverage": 0.80,
            "sensitivity": {
                "Is": 0.70,
            },
        },
        
        # keep it simple, Omicron had less frequent whole-whome testing
        "outbreak_testing": {
            "enabled": False,
        }
    }
}


##############################
## Detection layer helpers
##############################

def deep_update(base, override):
    """
    Recursively update a nested dictionary
    """
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base

def build_detection_params(wave, detection_params=None):
    if wave is None:
        return {
            "enabled": False,
            "wave": None,
            "observed_groups": ["B"],
        }

    params = deepcopy(DEFAULT_DETECTION_PARAMS[wave])

    if detection_params is not None:
        params = deep_update(params, detection_params)

    return params

def detection_enabled(detection_params):
    return detection_params is not None and detection_params.get("enabled", False)

def detection_block(state):
    """
    Collapse disease state to the block used by the observation model

    Returns:
      "E", "Ip", "Ia", "Is", or None
    """
    if state == "S" or state == "R":  # S and R cannot test positive (assumption)
        return None
    if state.startswith("E:"):
        return "E"
    if state.startswith("Ip:"):
        return "Ip"
    if state.startswith("Ia:"):
        return "Ia"
    if state.startswith("Is:"):
        return "Is"

    raise ValueError(f"Unknown state in detection_block: {state}")

def initialize_detected_output(graph, detection_params):
    """
    Initialize auxiliary observed-case process
    """
    counts = compute_detected_counts(graph, detection_params)

    return {
        "wave": detection_params.get("wave"),
        "times": [0.0],
        "counts": [counts],
        "events": [],
        "outbreaks": [],
        "params": detection_params
    }

def compute_detected_counts(graph, detection_params):
    """
    Count observed/detected cases among observed groups

    Observed_groups is ["B"] because we compare to resident case counts
    """
    observed_groups = set(detection_params.get("observed_groups", ["B"]))

    D = 0
    D_by_group = {}

    for v in graph.nodes:
        g = graph.nodes[v].get("group")
        if g not in observed_groups:
            continue

        if graph.nodes[v].get("detected", False):
            D += 1
            D_by_group[g] = D_by_group.get(g, 0) + 1

    out = {"D": D}
    for g, d in D_by_group.items():
        out[f"D_{g}"] = d

    return out

def next_initial_routine_test_time(detection_params, rng):
    """
    Determine first routine testing time

    If first_test_time is None, choose a random phase in [0, period)

    NOTE: This avoids arbitrary aligning the epidemic start with the testing calendar
    """
    if not detection_enabled(detection_params):
        return np.inf

    routine = detection_params.get("routine_testing", {})
    if not routine.get("enabled", False):
        return np.inf

    period = float(routine["period"])
    first = routine.get("first_test_time", None)

    if first is None:
        return float(rng.uniform(0.0, period))

    return float(first)


def mark_detected(graph, v, t, detected_by, test_type):
    graph.nodes[v]["detected"] = True
    graph.nodes[v]["detected_time"] = float(t)
    graph.nodes[v]["detected_by"] = detected_by
    graph.nodes[v]["detected_test_type"] = test_type


def append_detected_event(detected, graph, detection_params, event):
    """
    Append one observation / testing event and the resulting cumulative D(t)
    """
    detected["events"].append(event)
    detected["times"].append(float(event["t"]))
    detected["counts"].append(compute_detected_counts(graph, detection_params))


def process_routine_testing_event(graph, t, detection_params, rng):
    """
    Routine testing of observed groups

    Testing is bookkeeping only, it does not alter disease state or rates
    """
    routine = detection_params["routine_testing"]

    observed_groups = set(detection_params.get("observed_groups", ["B"]))
    test_type = routine.get("test_type", None)
    coverage = float(routine.get("coverage", 1.0))
    sensitivity = routine.get("sensitivity", {})

    eligible_nodes = []
    tested_nodes = []
    positive_nodes = []

    for v in graph.nodes:
        if graph.nodes[v].get("group") not in observed_groups:
            continue

        # already counted: do not double count
        if graph.nodes[v].get("detected", False):
            continue

        eligible_nodes.append(v)

        # individual is included in this testing round with probability coverage
        if rng.random() >= coverage:
            continue

        tested_nodes.append(v)

        block = detection_block(graph.nodes[v]["state"])
        if block is None:
            continue

        p_pos = float(sensitivity.get(block, 0.0))

        if rng.random() < p_pos:
            mark_detected(graph, v, t, detected_by="routine", test_type=test_type)
            positive_nodes.append(v)

    counted_nodes = counted_positive_nodes(graph, positive_nodes, detection_params)
    return {
        "t": float(t),
        "kind": "routine",
        "test_type": test_type,
        "eligible_nodes": list(eligible_nodes),
        "tested_nodes": list(tested_nodes),
        "positive_nodes": list(positive_nodes),
        "counted_positive_nodes": list(counted_nodes),
        "node_states": {v: graph.nodes[v]["state"] for v in eligible_nodes},
        "n_eligible": len(eligible_nodes),
        "n_tested": len(tested_nodes),
        "n_positive": len(positive_nodes),
        "n_counted_positive": len(counted_nodes),
    }


def process_symptom_testing_event(graph, t, v, detection_params, rng):
    """
    Symptom-triggered testing for a resident who has just entered Is:1

    Called immediately after disease event Ip:K_pre -> Is:1
    """
    symptom = detection_params.get("symptom_testing", {})
    if not symptom.get("enabled", False):
        return None

    observed_groups = set(detection_params.get("observed_groups", ["B"]))

    if graph.nodes[v].get("group") not in observed_groups:
        return None

    if graph.nodes[v].get("detected", False):
        return None

    test_type = symptom.get("test_type", None)
    coverage = float(symptom.get("coverage", 1.0))
    sensitivity = symptom.get("sensitivity", {})

    eligible_nodes = [v]
    tested_nodes = []
    positive_nodes = []

    # resident develops symptoms but may or may not be tested
    if rng.random() < coverage:
        tested_nodes.append(v)

        block = detection_block(graph.nodes[v]["state"])
        p_pos = float(sensitivity.get(block, 0.0))

        if rng.random() < p_pos:
            mark_detected(
                graph=graph,
                v=v,
                t=t,
                detected_by="symptom",
                test_type=test_type,
            )
            positive_nodes.append(v)

    counted_nodes = counted_positive_nodes(graph, positive_nodes, detection_params)

    return {
        "t": float(t),
        "kind": "symptom",
        "test_type": test_type,
        "eligible_nodes": list(eligible_nodes),
        "tested_nodes": list(tested_nodes),
        "positive_nodes": list(positive_nodes),
        "counted_positive_nodes": list(counted_nodes),
        "node_states": {v: graph.nodes[v]["state"] for v in eligible_nodes},
        "n_eligible": len(eligible_nodes),
        "n_tested": len(tested_nodes),
        "n_positive": len(positive_nodes),
        "n_counted_positive": len(counted_nodes),
    }


def counted_positive_nodes(graph, positive_nodes, detection_params):
    """
    Returns positive nodes that contribute to the observed case count
    This means residents only, group B
    """
    observed_groups = set(detection_params.get("observed_groups", ["B"]))
    return [
        v for v in positive_nodes
        if graph.nodes[v].get("group") in observed_groups
    ]

# outbreak testing calendar
def next_outbreak_test_time(outbreak_test_times):
    if not outbreak_test_times:
        return np.inf
    return min(outbreak_test_times)

def pop_next_outbreak_test_time(outbreak_test_times):
    t_next = min(outbreak_test_times)
    outbreak_test_times.remove(t_next)
    return t_next

def resolve_outbreak_round_offset(round_spec, rng):
    """
    Convert a round specification to a time offset
    """
    if isinstance(round_spec, (int, float)):
        return float(round_spec)

    if isinstance(round_spec, (tuple, list)):
        if len(round_spec) == 3 and round_spec[0] == "uniform":
            return float(rng.uniform(float(round_spec[1]), float(round_spec[2])))

    raise ValueError(f"Unknown outbreak round specification: {round_spec!r}")

def maybe_trigger_outbreak_testing(
    t,
    detected,
    detection_params,
    rng,
    outbreak_state,
    outbreak_test_times):
    """
    Trigger outbreak testing if enough resident detections have occurred in the recent trigger window
    using observed / counted positives, not latent infections
    """
    outbreak = detection_params.get("outbreak_testing", {})
    if not outbreak.get("enabled", False):
        return False

    # Do not repeatedly trigger during the same outbreak episode.
    if t < outbreak_state.get("cooldown_until", -np.inf):
        return False

    trigger_count = int(outbreak.get("trigger_count", 2))
    trigger_window = float(outbreak.get("trigger_window", 14.0))

    recent_count = recent_counted_positive_count(
        detected=detected,
        t=t,
        window=trigger_window,
    )

    if recent_count < trigger_count:
        return False

    rounds = outbreak.get("rounds", [0.0, ("uniform", 4.0, 7.0)])

    scheduled_times = []
    for round_spec in rounds:
        offset = resolve_outbreak_round_offset(round_spec, rng)
        t_round = float(t + offset)
        outbreak_test_times.append(t_round)
        scheduled_times.append(t_round)

    cooldown = float(outbreak.get("cooldown", trigger_window))
    outbreak_state["cooldown_until"] = float(t + cooldown)
    outbreak_state["n_triggered"] = outbreak_state.get("n_triggered", 0) + 1

    detected["outbreaks"].append({
        "t_trigger": float(t),
        "recent_counted_positive_count": recent_count,
        "trigger_window": trigger_window,
        "scheduled_test_times": scheduled_times,
        "cooldown_until": outbreak_state["cooldown_until"],
    })

    return True

def recent_counted_positive_count(detected, t, window):
    """
    Count observed resident positives in the recent window
    """
    total = 0

    for ev in detected["events"]:
        ev_t = float(ev["t"])

        if t - ev_t <= window:
            total += int(ev.get("n_counted_positive", 0))

    return total

def process_outbreak_testing_event(graph, t, detection_params, rng):
    """
    Resident outbreak PCR testing
    """
    outbreak = detection_params.get("outbreak_testing", {})
    if not outbreak.get("enabled", False):
        return None

    test_groups = set(outbreak.get("test_groups", ["B"]))
    test_type = outbreak.get("test_type", "PCR")
    coverage = float(outbreak.get("coverage", 1.0))
    sensitivity = outbreak.get("sensitivity", {})

    eligible_nodes = []
    tested_nodes = []
    positive_nodes = []

    for v in graph.nodes:
        if graph.nodes[v].get("group") not in test_groups:
            continue

        # already counted: no double counting
        if graph.nodes[v].get("detected", False):
            continue

        eligible_nodes.append(v)

        if rng.random() >= coverage:
            continue

        tested_nodes.append(v)

        block = detection_block(graph.nodes[v]["state"])
        if block is None:
            continue

        p_pos = float(sensitivity.get(block, 0.0))

        if rng.random() < p_pos:
            mark_detected(
                graph=graph,
                v=v,
                t=t,
                detected_by="outbreak",
                test_type=test_type,
            )
            positive_nodes.append(v)

    return {
        "t": float(t),
        "kind": "outbreak",
        "test_type": test_type,
        "eligible_nodes": list(eligible_nodes),
        "tested_nodes": list(tested_nodes),
        "positive_nodes": list(positive_nodes),
        "counted_positive_nodes": list(positive_nodes),
        "node_states": {v: graph.nodes[v]["state"] for v in eligible_nodes},
        "n_eligible": len(eligible_nodes),
        "n_tested": len(tested_nodes),
        "n_positive": len(positive_nodes),
        "n_counted_positive": len(positive_nodes),
    }

def do_detection_checks(graph, detected, detection_params):
    """
    Basic checks for observed detection process
    """
    if not detection_enabled(detection_params):
        return

    observed_groups = set(detection_params.get("observed_groups", ["B"]))
    N_obs = sum(1 for v in graph.nodes if graph.nodes[v].get("group") in observed_groups)

    counts = detected["counts"]
    Ds = [c["D"] for c in counts]

    for D in Ds:
        if not (0 <= D <= N_obs):
            raise RuntimeError(f"Detected count out of bounds: D={D}, N_obs={N_obs}")

    for prev, curr in zip(Ds, Ds[1:]):
        if curr < prev:
            raise RuntimeError(f"Detected count decreased: prev={prev}, curr={curr}")
    
    # TODO: additional checks?


