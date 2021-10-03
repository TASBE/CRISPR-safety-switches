import itertools

import workarounds

def transitive_closure(d: dict) -> dict:
    """Interpreting a dictionary as an acyclic directed graph, create a transitive closure of all k->v edges
    For example {1:[2,5], 2:[3,5], 3:[], 4:[5], 5:[6], 6:[]}
    returns {1:[2,3,5,6], 2:[3,5,6], 3:[], 4:[5,6], 5:[6], 6:[]}

    :param d: dictionary to close
    :return: closure dictionary
    """
    closure = d
    pending = closure
    while pending:
        # find out which are leaves
        resolvable = {k for k, v in pending.items() if not v}
        if not resolvable:
            raise ValueError(f'Cannot compute closure on cycle graph {d}')
        # union all the targets one step away
        for r in resolvable:
            closure[r] = sorted(set(closure[r]) | set(itertools.chain(*(closure[x] for x in closure[r]))))
        # remove everything that's been resolved
        pending = {k: [x for x in v if x not in resolvable] for k, v in pending.items() if k not in resolvable }
    return closure

