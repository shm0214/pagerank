# 期中作业：Compute the PageRank scores on the given dataset

## 数据集说明

采用指定数据集——“data.txt”。文件每行表示一条有向边，格式为：$FromNodeID[空格]ToNodeID$。

## 关键代码细节

我们分别实现了基础算法（Base Algorithm, BA）、基础分块算法（Block-based Update Algorithm, BBUA）、条带分块算法（Block-Stripe Update Algorithm, BSUA），基础算法的实现考虑了dead ends 和spider trap节点、进行了优化稀疏矩阵，基础分块算法和条带分块算法进一步实现了两种分块运算。

### 基础算法

### 基础分块算法
基础的分块算法将新的分数列表分为若干子块，对每个子块的计算分别需要遍历一次矩阵和旧的分数。由于delta(即$\frac{1-\sum_{j} r'^{new}_{j}}{N} $)只能在所有分块遍历结束后获得，而每个分块在各自对应遍历结束后就需要写入文件。为了保证不额外进行读写，文件中写入的新的分数不包含delta项，下一次读取后再加上delta项。这样处理后，可以保证读写次数为$K|M|+(K+1)r$。

这段代码是为计算某个子块而遍历矩阵时，对矩阵某一行进行计算的代码，与基础算法稍有不同，仅对子块内的分数项进行更新。
```python
while r_old_file_idx <= from_idx: #顺序读取上一次迭代的分数
    # old_rfile保存的分数不包含delta项，需加上
    curr_old_r = float(old_rfile.readline())+delta_old
    read_num += 1
    if r_old_file_idx in block_range: #保存相关分数项，用于计算当前误差
        r_old[r_old_file_idx]=curr_old_r
    r_old_file_idx += 1

for to in tos:            #遍历所有目的节点
    idx = id2idx[int(to)] #目的节点编号
    if idx in block_range:
        r_new[idx] += beta * curr_old_r / degree # 累加
```

为了保证不额外进行读写，需要在遍历时就计算误差，但是此时新的分数还不包含delta项，因此做如下近似处理。
```python
#计算err       
for idx in r_old.keys():
    err+=abs((r_new[idx]-(r_old[idx]-delta_old))/beta)
```

### 条带分块算法


## 实验结果及分析

### 正确性验证

我们将结果与networkx和igraph中的pagerank实现进行了对比。我们的结果与igraph所得结果差异在误差范围内，与networkx所得结果差距较大，暂不清楚原因。

### 读写性能分析

我们统计了pagerank计算过程中的读写数量，以读写的浮点数数量为单位，用来比较BBUA与BSUA两种算法的读写性能差异。

| 算法 	| 设置分块大小 	| 运行时间 	| 读取浮点数数量 	| 写入浮点数数量 	|
|:----:	|-------------:	|---------:	|---------------:	|---------------:	|
|  BA  	|         6263 	|          	|        4723532 	|         281835 	|
| BBUA 	|         1000 	|          	|       31135720 	|         281835 	|
| BSUA 	|              	|          	|                	|                	|

…………