# src/utils.py

from src.gillespie_sim import get_node_order

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


def print_time_counts(times, counts, title, n=10):
    print(title)
    print("-" * len(title))

    for t, c in list(zip(times, counts))[:n]:
        print(f"t = {float(t):8.4f}   {c}")

    if len(times) > n:
        print("...")

    print()


def print_last_time_counts(times, counts, title, n=10):
    print(title)
    print("-" * len(title))

    for t, c in list(zip(times, counts))[-n:]:
        print(f"t = {float(t):8.4f}   {c}")

    print()


def print_disease_events(events, title, n=8):
    print(title)
    print("-" * len(title))

    for ev in events[:n]:
        print(
            f"t = {float(ev['t']):8.4f}   "
            f"type = {ev['event_type']:12s}   "
            f"node = {ev['node']}   "
            f"parent = {ev['parent']}"
        )

    if len(events) > n:
        print("...")

    print()


def print_detected_process(detected):
    print("Detected resident-case process")
    print("------------------------------")
    print(f"wave = {detected['wave']}")
    print()

    for t, c in zip(detected["times"], detected["counts"]):
        d_b = c.get("D_B", c["D"])
        print(f"t = {float(t):8.4f}   D = {c['D']:2d}   D_B = {d_b:2d}")

    print()


def print_detection_events(detected):
    print("Detection events")
    print("----------------")

    if not detected["events"]:
        print("<none>")
        print()
        return

    for ev in detected["events"]:
        print(
            f"t = {float(ev['t']):8.4f}   "
            f"kind = {ev['kind']:8s}   "
            f"test = {str(ev['test_type']):3s}   "
            f"n_tested = {ev['n_tested']:2d}   "
            f"n_pos = {ev['n_positive']:2d}   "
            f"n_counted = {ev['n_counted_positive']:2d}   "
            f"positive_nodes = {ev['positive_nodes']}"
        )

    print()


def print_outbreaks(detected):
    print("Outbreak triggers")
    print("-----------------")

    if not detected["outbreaks"]:
        print("<none>")
        print()
        return

    for ob in detected["outbreaks"]:
        scheduled = [round(float(x), 4) for x in ob["scheduled_test_times"]]
        print(
            f"t_trigger = {float(ob['t_trigger']):8.4f}   "
            f"recent_count = {ob['recent_counted_positive_count']}   "
            f"scheduled_test_times = {scheduled}"
        )

    print()


def detected_resident_nodes(out):
    graph = out["graph"]

    return [
        v for v in get_node_order(graph)
        if graph.nodes[v].get("group") == "B"
        and graph.nodes[v].get("detected", False)
    ]


