from base import *


def block_based_pagerank(edges, sparse_matrix, r1, r2, beta, epsilon, bsize):
    """
    bsize为块大小
    r1 r2 为保存r_old和r_new的临时文件
    """
    from collections import defaultdict

    nodes = get_sparse_matrix(edges, sparse_matrix)
    id2idx = {}
    for idx, id in enumerate(nodes):
        id2idx[id] = idx
    nodes_num = len(nodes)
    print(f'nodes num:{nodes_num}')
    start = time.time()
    read_num = 0
    write_num = 0
    r = [1 / nodes_num] * nodes_num
    with open(r1, 'w') as f:
        for r_ in r:
            f.write('{}\n'.format(r_))
            write_num += 1

    block_ranges = [
        range(start, min(start + bsize, nodes_num))
        for start in range(0, nodes_num, bsize)
    ]
    delta_new = 0
    # while True:
    for i in range(100):
        delta_old = delta_new
        delta_new = 0
        err = 0
        r_new_sum = 0

        #清空r2
        with open(r2, 'w') as f:
            f.truncate()

        #遍历每个子块
        for block_range in block_ranges:
            r_old = {}
            r_new = dict(zip(block_range, [0] * len(block_range)))
            # 遍历整个矩阵
            with open(sparse_matrix, 'r') as mfile, open(r1, 'r') as old_rfile:
                r_old_file_idx = 0
                for line in mfile:
                    tos = line.split(' ')
                    read_num += len(tos)
                    from_ = int(tos[0])  #源节点
                    degree = int(tos[1])  #源节点出度
                    tos = tos[2:]  #目的节点列表
                    from_idx = id2idx[from_]  #源节点编号

                    while r_old_file_idx <= from_idx:  #顺序读取上一次迭代的分数
                        # old_rfile保存的分数不包含delta_old，需加上
                        curr_old_r = float(old_rfile.readline()) + delta_old
                        read_num += 1
                        if r_old_file_idx in block_range:
                            r_old[r_old_file_idx] = curr_old_r
                        r_old_file_idx += 1

                    for to in tos:  #遍历所有目的节点
                        idx = id2idx[int(to)]  #目的节点编号
                        if idx in block_range:
                            r_new[idx] += beta * curr_old_r / degree  # 累加

            #计算err
            for idx in r_old.keys():
                err += abs((r_new[idx] - (r_old[idx] - delta_old)) / beta)

            #写回r_new到r2
            with open(r2, 'a') as f:
                for idx, score in r_new.items():
                    f.write(f'{score}\n')
                    write_num += 1

            r_new_sum += sum(r_new.values())

        #计算delta
        delta_new = (1 - r_new_sum) / nodes_num

        #交换r1 r2
        r1, r2 = r2, r1

        if err < epsilon:
            print('finish at iter {}'.format(i))
            break

    end = time.time()
    print('total time: {}s.'.format(end - start))
    print('total read number: {}'.format(read_num))
    print('total write number: {}'.format(write_num))

    r = []
    with open(r1, 'r') as f:
        for line in f:
            r += [float(line) + delta_new]
    # print(delta_new)
    dt = {}
    for idx, v in enumerate(r):
        dt[nodes[idx]] = v

    with open('result_blockbased.txt', 'w') as f:
        for id, value in sorted(dt.items(), key=lambda x: x[1],
                                reverse=True)[0:100]:
            f.write("{}\t{}\n".format(id, value))


if __name__ == "__main__":
    # networkx_pagerank('data.txt')
    # igraph_pagerank('data.txt')
    base_pagerank('data.txt', 'data_sparse.txt', 'r.txt', 0.85, 1e-6)
    block_based_pagerank('data.txt', 'data_sparse.txt', 'r1.txt', 'r2.txt',
                         0.85, 1e-6, 1000)
