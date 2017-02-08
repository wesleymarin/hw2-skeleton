from hw2skeleton import cluster
from hw2skeleton import io
from hw2skeleton import utils
from hw2skeleton import cluster
import os
import matplotlib.pyplot as pp
import math


def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")

    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    utils.active_site_score([activesite_a, activesite_b])

    ## These active sites seem fairly different. I am not sure what this is meant for, since the scores
    #  are fairly arbitrary. I am asserting that these sites are different.
    assert cluster.compute_similarity(activesite_a, activesite_b) > 0.0


def test_partition_clustering():

    active_sites = io.read_active_sites("../Data")

    utils.active_site_score(active_sites)

    ## This section shows that less than 20 iterations are needed until the cluster score converges
    #  it tests up to 100 iterations.
    score_list = [0]
    x_list = []
    for i in range(0, 100):
        score = utils.cluster_score(cluster.cluster_by_partitioning(active_sites, 3, i+1))

        ## If it is the same for two iterations (only reason for two was to see the convergence line)
        if score == score_list[i-1]:
            score_list = score_list[1:]
            break
        score_list.append(score)
        x_list.append(i+1)

    pp.figure(figsize=(20,10))
    pp.plot(x_list, score_list, 'x-')
    pp.grid(True)
    pp.title("Iterations until cluster score convergence")
    pp.xlabel("number of iterations")
    pp.ylabel("cluster score")
    pp.show()

    ## This section shows the clustering quality based on the cluster_score function. The higher the y value
    #  the worse the cluster (many points far from the mean).
    y_part_vals = []
    x_part_vals = []
    step_size = math.floor(len(active_sites)/2/10)

    for i in range(0, 10):
        y_part_vals.append(utils.cluster_score(cluster.cluster_by_partitioning(active_sites, i*step_size+1, 20)))
        x_part_vals.append(i*step_size+1)

    pp.figure(figsize=(20,10))
    pp.semilogy(x_part_vals, y_part_vals, 'x-')
    pp.grid(True)
    pp.title("Partition clustering quality by number of clusters")
    pp.xlabel("number of clusters (k)")
    pp.ylabel("cluster score")
    pp.show()

    ## Another really large initial number to check against
    y_prev = 1000000000

    ## Assert that the overall cluster scores get smaller as the number of clusters increases
    for y in y_part_vals:
        assert y <= y_prev
        y_prev = y


def test_hierarchical_clustering():
    ## I am going for all active sites


    active_sites = io.read_active_sites("../Data")

    utils.active_site_score(active_sites)

    y_hier_vals = []
    x_hier_vals = []
    step_size = math.floor(len(active_sites)/2/10)
    for i in range(0,10):
        clusters = cluster.cluster_hierarchically(active_sites, i*step_size+1)
        y_hier_vals.append(utils.cluster_score(clusters))
        x_hier_vals.append(i*step_size+1)

    pp.figure(figsize=(20,10))
    pp.semilogy(x_hier_vals, y_hier_vals, 'x-')
    pp.grid(True)
    pp.title("Hierarchical clustering quality by number of clusters")
    pp.xlabel("number of clusters (k)")
    pp.ylabel("cluster score")
    pp.show()

    ## Assert that the overall cluster scores get smaller as the number of clusters increases
    y_prev = 1000000000
    for y in y_hier_vals:
        assert y <= y_prev
        y_prev = y

def test_cluster_differences():

    active_sites = io.read_active_sites("../Data")

    utils.active_site_score(active_sites)

    y_comb_vals = []
    x_vals = []
    step_size = math.floor(len(active_sites)/2/4)

    for i in range(0,4):
        y_hier_val = utils.cluster_score(cluster.cluster_hierarchically(active_sites, i*step_size+1))
        y_part_val = utils.cluster_score(cluster.cluster_by_partitioning(active_sites, i*step_size+1, 20))

        y_comb_vals.append(y_hier_val - y_part_val)

        x_vals.append(i*step_size+1)

    pp.figure(figsize=(20,10))
    pp.plot(x_vals, y_comb_vals, 'x-')
    pp.grid(True)
    pp.title("Hierarchical clustering score - Partition clustering score by k")
    pp.xlabel("number of clusters (k)")
    pp.ylabel("cluster score difference")
    pp.show()
