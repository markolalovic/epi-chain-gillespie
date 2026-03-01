# epi-chain-gillespie
This repository contains an implementation of Gillespie's algorithm for simulating a continuous-time Markov chain epidemic model on a contact network, with:

- weighted edges (interpreted as repeated contacts)

- two-type vertex dynamics (groups of vertices `A` and `B`)

- infection-chain tracing (recording parent infector on infection events)

- event-level output (time, event type, selected vertex, parent, counts, full state snapshots)

The simulation code is in `src/gillespie_sim.py`.

## Implemented model
Implemented model is a staged SEIR-type of model.

### State space
Node states are represented with strings:
- `"S"` (susceptible)
- `"E:1", ..., "E:K1"`  (stages of exposed block `E`)
- `"Ia:1", ..., "Ia:K2"` (stages of infectious asymptomatic block `Ia`)
- `"Is:1", ..., "Is:K3"` (stages of infectious symptomatic block `Is`)
- `"R"` (recovered)

(The infectious set is `Ia, Is`.)

### Exposed progression and branching
Exposed block has `K1` stages, with constant rate `sigma K1` at each stage:
- `E:k -> E:k+1` for `k < K1`
- `E:K1 -> Ia:1` with probability `alpha`
- `E:K1 -> Is:1` with probability `1 - alpha`

### Infectious progression and recovery
Asymptomatic block has `K2` stages, with rate `mu K2` at each stage:
- `Ia:k -> Ia:k+1` for `k < K2`
- `Ia:K2 -> R`

Symptomatic block has `K3` stages, with rate `mu K3` at each stage:
- `Is:k -> Is:k+1` for `k < K3`
- `Is:K3 -> R`

### Detection 
Detection is bookkeeping only (**not** a compartment):
- upon entry into `Is:1`, a Bernoulli(`p_detect`) trial is performed
- if successful, node attribute `detected=True` is set permanently
- and count `D(t)` = number of detected nodes at time `t` is incremented


## Example

### Example contact network
Load the example contact network
```python
from src.contact_networks import example_contact_network

graph_g = example_contact_network()

print(get_node_order(graph_g))
# ['a1', 'a2', 'a3', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6']
```

The example graph shown below has a partition into two groups (vertex types) `A` and `B`, with weighted edges between the groups. Edges with weight 1 are drawn without labels to reduce clutter:

![Example contact network](figures/contact_networks/g_weighted.svg)

### Example simulation
To run the simulation of staged SEIR with detection bookkeeping, using a fixed RNG seed, and store the output:
```python
from src.gillespie_sim import gillespie_sim, get_node_order

# set model parameters
model_params = {
    "beta_AB": 2.0,  # transmission \beta_{A, B}        
    "beta_AA": 2.0,  # transmission \beta_{A, B}
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
```

### Example output
The output is a dictionary with:

- `times` (event times starting at 0, ending at `time_max`)
- `counts` (aggregated counts of `S`, `E`, `I_asym`, `I_sym`, `R`, `D`)
- `events` (event list: `t`, `event_type`, `node`, `parent`, `detected_inc`, `states`)
- `node_order` (fixed node ordering used)
- `init` (initially infected sampled from group A)
- `graph` (final graph state includes `state`, `rate`, `detected` node attributes)

Example output:
```python
out['init']
# {'initial_infected_nodes': ['a2', 'a1']}

out["node_order"]
# ['a1', 'a2', 'a3', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6']

out['events'][:5]
# [{'detected_inc': 0,
#   'event_type': 'infection',
#   'node': 'b5',
#   'parent': 'a1',
#   'states': ['Ia:1', 'Ia:1', 'S', 'S', 'S', 'S', 'S', 'E:1', 'S'],
#   't': np.float64(0.0380083620196087)},
#  {'detected_inc': 0,
#   'event_type': 'infection',
#   'node': 'b4',
#   'parent': 'a2',
#   'states': ['Ia:1', 'Ia:1', 'S', 'S', 'S', 'S', 'E:1', 'E:1', 'S'],
#   't': np.float64(0.04812551960098333)},
#  {'detected_inc': 0,
#   'event_type': 'infection',
#   'node': 'b6',
#   'parent': 'a2',
#   'states': ['Ia:1', 'Ia:1', 'S', 'S', 'S', 'S', 'E:1', 'E:1', 'E:1'],
#   't': np.float64(0.08536049771657381)},
#  {'detected_inc': 0,
#   'event_type': 'infection',
#   'node': 'b3',
#   'parent': 'a1',
#   'states': ['Ia:1', 'Ia:1', 'S', 'S', 'S', 'E:1', 'E:1', 'E:1', 'E:1'],
#   't': np.float64(0.2205242250400909)},
#  {'detected_inc': 0,
#   'event_type': 'E_progress',
#   'node': 'b5',
#   'parent': None,
#   'states': ['Ia:1', 'Ia:1', 'S', 'S', 'S', 'E:1', 'E:1', 'E:2', 'E:1'],
#   't': np.float64(0.24967967610295402)}]

tc = list(zip(out["times"], out["counts"]))

print("first 10 time points with counts:")
for t, c in tc[:10]:
    print(float(t), c)
# 0.0 {'S': 7, 'E': 0, 'I_asym': 2, 'I_sym': 0, 'R': 0, 'D': 0}
# 0.0380083620196087 {'S': 6, 'E': 1, 'I_asym': 2, 'I_sym': 0, 'R': 0, 'D': 0}
# 0.04812551960098333 {'S': 5, 'E': 2, 'I_asym': 2, 'I_sym': 0, 'R': 0, 'D': 0}
# 0.08536049771657381 {'S': 4, 'E': 3, 'I_asym': 2, 'I_sym': 0, 'R': 0, 'D': 0}
# 0.2205242250400909 {'S': 3, 'E': 4, 'I_asym': 2, 'I_sym': 0, 'R': 0, 'D': 0}
# 0.24967967610295402 {'S': 3, 'E': 4, 'I_asym': 2, 'I_sym': 0, 'R': 0, 'D': 0}
# 0.30945569079928315 {'S': 2, 'E': 5, 'I_asym': 2, 'I_sym': 0, 'R': 0, 'D': 0}
# 0.3232687126991879 {'S': 2, 'E': 5, 'I_asym': 2, 'I_sym': 0, 'R': 0, 'D': 0}
# 0.3329481077232424 {'S': 1, 'E': 6, 'I_asym': 2, 'I_sym': 0, 'R': 0, 'D': 0}
# 0.3501603312215322 {'S': 1, 'E': 5, 'I_asym': 3, 'I_sym': 0, 'R': 0, 'D': 0}


print("last 10 time points with counts:")
for t, c in tc[-10:]:
    print(float(t), c)
# 1.711667568833028 {'S': 0, 'E': 1, 'I_asym': 1, 'I_sym': 2, 'R': 5, 'D': 4}
# 1.7356192593517656 {'S': 0, 'E': 1, 'I_asym': 0, 'I_sym': 2, 'R': 6, 'D': 4}
# 1.9280560667630224 {'S': 0, 'E': 0, 'I_asym': 0, 'I_sym': 3, 'R': 6, 'D': 5}
# 2.0886172548131707 {'S': 0, 'E': 0, 'I_asym': 0, 'I_sym': 3, 'R': 6, 'D': 5}
# 2.1214817426972568 {'S': 0, 'E': 0, 'I_asym': 0, 'I_sym': 2, 'R': 7, 'D': 5}
# 2.1804259208928687 {'S': 0, 'E': 0, 'I_asym': 0, 'I_sym': 2, 'R': 7, 'D': 5}
# 3.081834566581169 {'S': 0, 'E': 0, 'I_asym': 0, 'I_sym': 2, 'R': 7, 'D': 5}
# 3.200839788678001 {'S': 0, 'E': 0, 'I_asym': 0, 'I_sym': 1, 'R': 8, 'D': 5}
# 4.893070532827123 {'S': 0, 'E': 0, 'I_asym': 0, 'I_sym': 0, 'R': 9, 'D': 5}
# 50.0 {'S': 0, 'E': 0, 'I_asym': 0, 'I_sym': 0, 'R': 9, 'D': 5}

print("last 10 time points with D(t) only:")
for t, c in tc[-10:]:
    print(float(t), c["D"])
# 1.711667568833028 4
# 1.7356192593517656 4
# 1.9280560667630224 5
# 2.0886172548131707 5
# 2.1214817426972568 5
# 2.1804259208928687 5
# 3.081834566581169 5
# 3.200839788678001 5
# 4.893070532827123 5
# 50.0 5
```

See also `scripts/run_gillespie_sim.py` and `scripts/output.txt`.


### Infection-chain tracing 
The `parent` field records the sampled infector for an infection event, which allows infection-chain tracing.

For example, in the first event above at time `t = 0.038`:

- `node: 'b5'` means vertex `b5` became infected 
- `parent: 'a1'` means the simulation recorded `a1` as the infector of `b5`

So this event contributes the to infection-chain an edge:

- `a1 -> b5`

Thus, each infection event adds one directed edge `parent -> node` to the infection chain, which is a tree (forest in general) rooted at the initially infected vertices (`'a1', 'a2'` in this example).


## Setup

Minimal setup (tested with Python 3.14):
```zsh
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install numpy networkx
```

Optional: [SageMath](https://doc.sagemath.org/html/en/installation/index.html) for drawing and exporting graphs (tested with SageMath 10.8).



## References / Links

- [Wikipedia : Gillespie algorithm](https://en.wikipedia.org/wiki/Gillespie_algorithm) stochastic simulation algorithm (SSA)

- Tutorial [arXiv 2112.05293](https://arxiv.org/abs/2112.05293): "Gillespie algorithms for stochastic multiagent dynamics in populations and networks", by Naoki Masuda, Christian L. Vestergaard, (2021): includes SIR on networks, the first reaction method, Gillespie's direct method, and discussion of computational complexity.


### Note on terminology

According to the tutorial, this repository currently implements *Gillespie's direct method* for a continuous-time Markov-chain SIR model on a contact network, with **vertex-centric events**.

The simulation step is of the form:

1. Compute total rate $\Lambda = \sum_i \lambda_i$

2. Draw waiting time $\tau \sim \mathrm{Exp}(\Lambda)$

3. Select the next event with probability $\lambda_i / \Lambda$

4. Update the state and affected rates

### Additional notes

- The treatment of two-type vertex labels and infection-chain tracing is additional bookkeeping.

- Edge weights modify infection intensities and parent-sampling probabilities.
