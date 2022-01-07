import numpy as np


def neighbor_vertice(G, v):
    result = set()
    result.add(v)
    for i in range(G.shape[1]):
        if G[v, i] != 0:
            result.add(i)
    return result


def neighbor_edge(G, e):
    return neighbor_vertice(G, e[0]).intersection(neighbor_vertice(G, e[1]))


def is_dominated_vertice(G, e, v):
    return neighbor_edge(G, e).issubset(neighbor_vertice(G, v)) and not (v in e)


def is_dominated(G, e):
    for v in range(G.shape[1]):
        if is_dominated_vertice(G, e, v):
            return True
    return False


def exist_domination(G):
    for i in range(G.shape[0]):
        for j in range(i + 1, G.shape[1]):
            if G[i, j] != 0 and is_dominated(filtered(G, G[i, j]), (i, j)):
                return True
    return False


def filtered(G, t):
    return np.where(G <= t, G, 0)


def Optimize(G):
    while exist_domination(G):
        T = np.unique(np.sort(G.flatten()))[1:]
        for i in range(len(T)):
            for e in np.unique(np.argwhere(G == T[i]), axis=0):
                if is_dominated(filtered(G, T[i]), e):
                    if i == len(T) - 1:
                        G[e[0], e[1]] = 0
                    else:
                        G[e[0], e[1]] = T[i + 1]
    return G


if __name__ == "__main__":
    G = np.eye(10)
    G = -(G - 1)
    with open("graph", "w") as f:
        for r in range(G.shape[0]):
            for c in range(r + 1, G.shape[1]):
                if G[r, c] != 0:
                    f.writelines(["{:} {:} {:}\n".format(r, c, G[r, c])])
    """
    with open("graph", "r") as f:
        lines = f.readlines()
    tmp = np.array([[float(i) for i in line.strip('\n').split(" ")] for line in lines])
    size = int(np.max(tmp[:, :2])) + 1
    G = np.zeros((size, size))
    for i in range(tmp.shape[0]):
        G[int(tmp[i, 0]), int(tmp[i, 1])] = tmp[i, 2]
        G[int(tmp[i, 1]), int(tmp[i, 0])] = tmp[i, 2]
    """
    print(G)
    M = Optimize(G)
    print(M)
    with open("result", "w") as f:
        for r in range(M.shape[0]):
            for c in range(r + 1, M.shape[1]):
                if M[r, c] != 0:
                    f.writelines(["{:} {:} {:}\n".format(r, c, M[r, c])])
