from typing import Tuple

import numpy as np


def neighbor_vertice(G: np.ndarray, v: int):
    """
    retunr neighbors of vertive v
    """
    result = set()
    result.add(v)
    for i in range(G.shape[1]):
        if G[v, i] != 0:
            result.add(i)
    return result


def neighbor_edge(G: np.ndarray, e: Tuple[int]):
    """
    return neighbors of edge e
    """
    return neighbor_vertice(G, e[0]).intersection(neighbor_vertice(G, e[1]))


def is_dominated_vertice(G: np.ndarray, e: Tuple[int], v: int):
    """
    return True is e is dominated by v
    """
    return neighbor_edge(G, e).issubset(neighbor_vertice(G, v)) and not (v in e)


def is_dominated(G: np.ndarray, e: Tuple[int]):
    """
    return True if e is dominated
    """
    for v in range(G.shape[1]):
        if is_dominated_vertice(G, e, v):
            return True
    return False


def exist_domination(G: np.ndarray):
    """
    return True is there is an edge dominated in G
    """
    for i in range(G.shape[0]):
        for j in range(i + 1, G.shape[1]):
            if G[i, j] != 0 and is_dominated(filtered(G, G[i, j]), (i, j)):
                return True
    return False


def filtered(G: np.ndarray, t: int):
    """
    return graph where edges that have larger filtration value than t is removed
    """
    return np.where(G <= t, G, 0)


def Optimize(G: np.ndarray):
    """
    :param G: original graph
    :return: optimized graph
    """
    while exist_domination(G):
        T = np.unique(np.sort(G.flatten()))[1:]
        for i in range(len(T)):
            for e in np.unique(np.argwhere(G == T[i]), axis=0):
                if is_dominated(filtered(G, T[i]), e):
                    if i == len(T) - 1:
                        G[e[0], e[1]] = 0
                        G[e[1], e[0]] = 0
                    else:
                        G[e[0], e[1]] = T[i + 1]
                        G[e[1], e[0]] = T[i + 1]
    return G


if __name__ == "__main__":
    # Complete Graph
    G = np.eye(10)
    G = -(G - 1)

    '''
    # Regular 2n-gon
    G = np.eye(4)
    G = -(G - 1)
    G[0,2] = np.sqrt(2)
    G[2, 0] = np.sqrt(2)
    G[3, 1] = np.sqrt(2)
    G[1, 3] = np.sqrt(2)
    '''
    '''
    # Random points
    points = [[1, 1], [7, 0], [4, 6], [9, 6], [0, 14], [2, 19], [9, 17]]
    rips_complex = gudhi.RipsComplex(points=points)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=1)
    G = np.zeros((len(points), len(points)))
    for filtered_value in simplex_tree.get_filtration():
        if len(filtered_value[0]) == 1:
            continue
        G[filtered_value[0][0], filtered_value[0][1]] = filtered_value[1]
        G[filtered_value[0][1], filtered_value[0][0]] = filtered_value[1]
    '''
    # store the original graph in the form compatible with Ripser
    with open("graph", "w") as f:
        for r in range(G.shape[0]):
            for c in range(r + 1, G.shape[1]):
                if G[r, c] != 0:
                    f.writelines(["{:} {:} {:}\n".format(r, c, G[r, c])])
    # able to read a Graph from a file
    with open("graph", "r") as f:
        lines = f.readlines()
    tmp = np.array([[float(i) for i in line.strip('\n').split(" ")] for line in lines])
    size = int(np.max(tmp[:, :2])) + 1
    G = np.zeros((size, size))
    for i in range(tmp.shape[0]):
        G[int(tmp[i, 0]), int(tmp[i, 1])] = tmp[i, 2]
        G[int(tmp[i, 1]), int(tmp[i, 0])] = tmp[i, 2]

    print(G)
    M = Optimize(G)
    print(M)
    # store the graph optimized in the form compatible with Ripser
    with open("result", "w") as f:
        for r in range(M.shape[0]):
            for c in range(r + 1, M.shape[1]):
                if M[r, c] != 0:
                    f.writelines(["{:} {:} {:}\n".format(r, c, M[r, c])])
