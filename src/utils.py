# src/utils.py

def print_event_table(out):
    node_order = out["node_order"]

    # initial row at t = 0
    init_infected = out["init"]["initial_infected_nodes"]
    init_infected = [str(x) for x in init_infected]
    init_states = ["S"] * len(node_order)
    for v in init_infected:
        idx = node_order.index(v)
        init_states[idx] = "I"

    # if one initial infected, print that node, else print list
    init_node = init_infected[0] if len(init_infected) == 1 else init_infected

    print(f"{'time':<8} | {'event':<5} | {'i_selected':<10} | states")
    print(f"{0.0000:8.4f} | {'':<5} | {str(init_node):<10} | [{','.join(init_states)}]")

    for ev in out["events"]:
        event_code = "I" if ev["event_type"] == "infection" else "R"
        states_str = ",".join(ev["states"])
        print(f"{ev['t']:8.4f} | {event_code:<5} | {str(ev['node']):<10} | [{states_str}]")
