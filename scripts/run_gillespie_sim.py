# scripts/run_gillespie_sim.py

from pathlib import Path
from pprint import pprint
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.contact_networks import example_contact_network
from src.gillespie_sim import gillespie_sim, get_node_order
from src.utils import *

if __name__ == "__main__":
    
    # 1. Contact network
    graph_g = example_contact_network()

    print("Example contact network")
    print("=======================")
    print("node_order =")
    print(get_node_order(graph_g))
    print()

    print("node groups =")
    pprint({v: graph_g.nodes[v]["group"] for v in get_node_order(graph_g)})
    print()

    # 2. Model parameters
    model_params = {
        "beta_AB": 2.0,      # transmission beta_{A,B} (= beta_{B,A})
        "beta_AA": 2.0,      # transmission beta_{A,A}
        "beta_BB": 2.0,      # transmission beta_{B,B}

        "sigma": 1.0,        # exposed progression rate
        "lambda_pre": 1.0,   # presymptomatic progression rate
        "mu": 1.0,           # asymptomatic/symptomatic recovery rate

        "K1": 2,             # number of exposed stages
        "K_pre": 2,          # number of presymptomatic infectious stages
        "K2": 2,             # number of asymptomatic infectious stages
        "K3": 2,             # number of symptomatic infectious stages

        "alpha": 0.25,       # P(Ip:K_pre -> Ia:1)
    }

    time_max = 50.0
    n_init_infected = 2
    seed = 12345
    wave = "wave2"

    print("Simulation inputs")
    print("=================")
    print(f"time_max = {time_max}")
    print(f"n_init_infected = {n_init_infected}")
    print(f"rng seed = {seed}")
    print(f"wave = {wave!r}")
    print()
    print("model_params =")
    pprint(model_params)
    print()

    # 3. Run simulation
    out = gillespie_sim(
        graph_g=graph_g,
        model_params=model_params,
        time_max=time_max,
        rng=np.random.default_rng(seed),
        n_init_infected=n_init_infected,
        wave=wave
    )

    # 4. Latent disease-process output
    print("Simulation output: latent disease process")
    print("=========================================")
    print("initial infected =")
    pprint(out["init"])
    print()

    print_time_counts(
        times=out["times"],
        counts=out["counts"],
        title="First time points with latent disease counts",
        n=10,
    )

    print_last_time_counts(
        times=out["times"],
        counts=out["counts"],
        title="Last time points with latent disease counts",
        n=10,
    )

    print_disease_events(
        events=out["events"],
        title="First disease events",
        n=10,
    )

    # 5. Detection process output
    print("Simulation output: observed resident cases")
    print("==========================================")
    print_detected_process(out["detected"])
    print_detection_events(out["detected"])
    print_outbreaks(out["detected"])

    print("Final summary")
    print("=============")
    print(f"final latent disease counts = {out['counts'][-1]}")
    print(f"final detected counts = {out['detected']['counts'][-1]}")
    print(f"detected resident nodes = {detected_resident_nodes(out)}")
