import time



def networkx_pagerank(file):
    import networkx as nx
    G = nx.DiGraph()
    with open(file) as f:
        for line in f:
            head, tail = [int(x) for x in line.split()]
            G.add_edge(head, tail)
    pr = nx.pagerank(G, alpha=0.85, max_iter=100, tol=1e-6)
    with open('result_networkx.txt', 'w') as f:
        for node, value in sorted(pr.items(), key=lambda x: x[1],
                                  reverse=True):
            f.write("{}\t{}\n".format(node, value))


def igraph_pagerank(file):
    from igraph import Graph as IGraph
    edges = []
    with open(file) as f:
        for line in f:
            head, tail = [int(x) for x in line.split()]
            edges.append((head, tail))
    g = IGraph.TupleList(edges, directed=True, vertex_name_attr='id')
    pg = g.pagerank(implementation='power')
    pgvs = []
    for p in zip(g.vs, pg):
        pgvs.append({'id': p[0]['id'], 'pg': p[1]})
    with open('result_igraph.txt', 'w') as f:
        for pgv in sorted(pgvs, key=lambda k: k['pg'], reverse=True):
            f.write("{}\t{}\n".format(pgv['id'], pgv['pg']))


def get_sparse_matrix(file, sparse_matrix):
    from collections import defaultdict
    m = defaultdict(lambda: [0, []])
    nodes = set()
    with open(file, 'r') as f:
        for line in f:
            from_, to = [int(x) for x in line.split()]
            m[from_][0] += 1
            m[from_][1].append(to)
    with open(sparse_matrix, 'w') as f:
        for from_, (degree, tos) in sorted(m.items(), key=lambda x: x[0]):
            f.write("{} {} {}\n".format(from_, degree,
                                        ' '.join(str(x) for x in sorted(tos))))
            nodes.add(from_)
            for to in tos:
                nodes.add(to)
    return sorted(nodes)


def base_pagerank(edges, sparse_matrix, r_old, beta, epsilon):
    '''
    ppt 中的2|r|+|M| 不知道咋做到的 
    计算err 和 写回r_new 好像没法同时完成
    所以这个应该是3|r|+|M|
    前三个参数 分别为数据集文件名、保存压缩矩阵的文件名和存储r的文件名
    '''
    nodes = get_sparse_matrix(edges, sparse_matrix)
    id2idx = {}
    for idx, id in enumerate(nodes):
        id2idx[id] = idx
    nodes_num = len(nodes)
    start = time.time()
    read_num = 0
    write_num = 0
    r = [1 / nodes_num] * nodes_num
    with open(r_old, 'w') as f:
        for r_ in r:
            f.write('{}\n'.format(r_))
            write_num += 1
    # while True:
    for i in range(100):
        r_file = open(r_old, 'r')
        r_file_idx = 0
        init = (1 - beta) / nodes_num
        r_new = [init] * nodes_num
        with open(sparse_matrix, 'r') as f:
            for line in f:
                tos = line.split(' ')
                read_num += len(tos)
                from_ = int(tos[0])
                degree = int(tos[1])
                tos = tos[2:]
                from_idx = id2idx[from_]
                while (r_file_idx != from_idx):
                    r_file_idx += 1
                    read_num += 1
                    _ = r_file.readline()
                r_ = float(r_file.readline())
                if r_ != r[from_idx]:
                    print('hi')
                r_file_idx += 1
                read_num += 1
                for to in tos:
                    idx = id2idx[int(to)]
                    r_new[idx] += beta * r_ / degree
        r_file.close()
        r_new_sum = sum(r_new)
        delta = (1 - r_new_sum) / nodes_num
        r_new = [r_ + delta for r_ in r_new]
        err = 0
        with open(r_old, 'r') as f:
            idx = 0
            for line in f:
                read_num += 1
                r_ = float(line)
                err += abs(r_ - r_new[idx])
                idx += 1
        with open(r_old, 'w') as f:
            for r_ in r_new:
                write_num += 1
                f.write('{}\n'.format(r_))
        r = r_new
        if err < epsilon:
            print('finish at iter {}'.format(i))
            break
    end = time.time()
    print('total time: {}s.'.format(end - start))
    print('total read number: {}'.format(read_num))
    print('total write number: {}'.format(write_num))
    dt = {}
    for idx, v in enumerate(r):
        dt[nodes[idx]] = v
    with open('result_base.txt', 'w') as f:
        num = 0
        for id, value in sorted(dt.items(), key=lambda x: x[1], reverse=True):
            f.write("{}\t{}\n".format(id, value))
            num += 1
            if (num == 100):
                break


if __name__ == "__main__":
    # networkx_pagerank('data.txt')
    # igraph_pagerank('data.txt')
    base_pagerank('data.txt', 'data_sparse.txt', 'r.txt', 0.85, 1e-6)