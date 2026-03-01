# scripts/run_gillespie_sim.py

# Test script
# for simple check that outputs are the same after modifications
# $ diff before after
# 

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
        "beta": 2.0,
        "mu": 1.0,
    }

    time_max = 5.0
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
    pprint(out["init"])
    print()

    print("out['node_order'] =")
    pprint(out["node_order"])
    print()

    print("out['events'][:5] =")
    pprint(out["events"][:5])

