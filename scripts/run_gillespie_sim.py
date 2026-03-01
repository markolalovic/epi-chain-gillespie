# scripts/run_gillespie_sim.py

# Test script
# for simple check that outputs are the same after modifications
# $ diff before after

from pathlib import Path
from pprint import pprint
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.contact_networks import example_contact_network
from src.gillespie_sim import gillespie_sim, get_node_order


if __name__ == "__main__":
    graph_g = example_contact_network()

    print("get_node_order(graph_g) =")
    print(get_node_order(graph_g))
    print()

    model_params = {
        "beta_AB": 2.0,  # transmission \beta_{A, B} (= \beta_{B, A})          
        "beta_AA": 2.0,  # transmission \beta_{A, A}
        "beta_BB": 2.0,  # transmission \beta_{B, B}
        "sigma": 1.0,    # E progression
        "mu": 1.0,       # I progression
        "K1": 2,         # number of E stages
        "K2": 2,         # number of I_asym stages
        "K3": 2,         # number of I_sym stages
        "alpha": 0.25,    # branching at E:K1 probability to go to I_asym:1
        "p_detect": 0.8,   # detection probability upon entry to Is:1
    }

    time_max = 50.0
    rng = np.random.default_rng(12345)
    n_init_infected = 2

    out = gillespie_sim(
        graph_g=graph_g,
        model_params=model_params,
        time_max=time_max,
        rng=rng,
        n_init_infected=n_init_infected,
    )

    print("out['init'] =")
    pprint(out['init'])
    print()

    print("out['node_order'] =")
    pprint(out['node_order'])
    print()

    print("out['events'][:5] =")
    pprint(out['events'][:5])
    print()

    tc = list(zip(out["times"], out["counts"]))

    print("first 10 time points with counts:")
    for t, c in tc[:10]:
        print(float(t), c)
    print()

    print("last 10 time points with counts:")
    for t, c in tc[-10:]:
        print(float(t), c)
    print()

    print("last 10 time points with D(t) only:")
    for t, c in tc[-10:]:
        print(float(t), c["D"])
    print()

