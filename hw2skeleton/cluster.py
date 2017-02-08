from .utils import Atom, Residue, ActiveSite
import matplotlib.pyplot as pp
import numpy as np

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    All of the heavy lifting is done by the active_site_scoring function, which computes a score
    for every active site. To check the similarity just find the absolute difference in the scores. A smaller
    difference means more similar.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """

    similarity = 0.0

    # Fill in your code here!

    similarity = abs(site_a.score - site_b.score)


    return similarity


def cluster_by_partitioning(active_sites, k, iterations):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.
    I am clustering by k_means.

    Input: a list of ActiveSite instances, k clusters, number of k_means iterations
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Fill in your code here!

    if len(active_sites) <= 1:
        return active_sites

    ## Initializing k_means dictionary with arbitrary k means
    k_means = {active_sites[i].score:[] for i in range(0,k)}

    ## Loop clustering for the desired number of iterations
    for i in range(0, iterations):

        ## Find the closest mean for every site in active sites, measured by score
        for site in active_sites:
            closest_mean = 1000000000 ## Initialize with massive value that will hopefully not cause bugs
                                      ## site.scores are centered around 0, so this number should never be met

            ## Check how far to every mean, save if closer mean is found
            for k_mean in k_means.keys():
                if abs(site.score - k_mean) < abs(site.score - closest_mean):
                    closest_mean = k_mean

            ## Use the closest found mean to append the active site to the k_means dictionary
            k_means[closest_mean].append(site)

        ## Placeholder for new k_means dictionary keys
        new_means = [sum([list_item.score for list_item in mean_list])/len(mean_list) for mean_list in k_means.values()]

        ## Save old k_means dictionary in case this is the last iteration
        return_means = k_means

        ## Create new k_means dictionary using newly generated means
        k_means = {new_mean:[] for new_mean in new_means}

    ## Once all iterations are complete, return the lists of the clustered active sites as a list
    return [value for value in return_means.values()]


def cluster_hierarchically(active_sites, k):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    if len(active_sites) <= 1:
        return active_sites

    ## Initialize starting clusters (every active site is its own cluster)
    clusters = {active_site:[active_site] for active_site in active_sites}

    ## Flag for if the clustering gets stuck, turns it from a reciprocal closest search, to a closest search
    super_add_mode = False

    ## While our total clusters are greater than k clusters
    while len(clusters) > k:

        ## Determine how many clusters are present (this becomes important at the end of the loop)
        cluster_length = len(clusters)

        ## Initialize new cluster dictionary
        new_clusters = {}

        ## Initialize best found cluster
        best_cluster = []

        ## Initialize best found cluster distance
        best_value = 1000000000

        ## Initialize best found cluster dictionary key
        best_key = []


        ## Loop over every cluster in the cluster dictionary
        for centroid_main in clusters.keys():

            ## Save the value of the current cluster (this is a list)
            centroid_cluster = clusters[centroid_main]

            ## Find the closest centroid to the current cluster centroid, excluding the current cluster
            ## Returns the closest cluster, as well as the score difference
            closest_centroid, closest_value = find_closest(centroid_main, clusters, [centroid_main])

            ## Find the closest centroid to the closest found centroid, excluding the closest foudn centroid
            reciprocal_centroid, reciprocal_value = find_closest(closest_centroid[0], clusters, closest_centroid)

            ## Check to see if both centroids are closest to each other, or see if clustering is stuck
            if centroid_main in reciprocal_centroid or super_add_mode:

                ## Is the found cluster score difference better than any other found cluster
                if abs(closest_value) < abs(best_value):

                    ## Set best found value to the cluster score difference
                    best_value = closest_value

                    ## Set the best cluster to this cluster
                    best_cluster = closest_centroid

                    ## Set the best cluster seed key to this cluster seed key
                    best_key = centroid_main

            ## Put cluster into new cluster dictionary (both cluster dicts are the same at this point)
            new_clusters[centroid_main] = centroid_cluster

        ## If a best cluster was found
        if len(best_cluster) > 0:

            ## Add the best cluster seed key to the list of clustered keys
            best_cluster = [best_key] + best_cluster

            ## Reset the cluster list
            centroid_cluster = []

            ## Remove all clustered elements from the new_clusters dictionary
            for centroid in best_cluster:
                centroid_cluster = centroid_cluster + clusters[centroid]
                new_clusters.pop(centroid, None)

            ## Add the clustered elements back into the new_clusters dictionary under the cluster seed key
            new_clusters[best_key] = centroid_cluster


        ## Replace the clusters dictionary with the new clusters dictionary, now ready for another round
        clusters = new_clusters

        ## Reset the stuck clustering flag
        super_add_mode = False

        ## If the original cluster length is the same as the new cluster length, turn on the stuck clustering flag
        if cluster_length == len(clusters):
            print("No new clusters found")
            super_add_mode = True

    ## Return a list of clusters
    return [cluster for cluster in clusters.values()]

def find_closest(check_key, check_dict, banned):
    """
    Find the closest cluster based on centroid, excluding any banned clusters                                                                #

    Input:  Dictionary key for cluster to check against
            Dictionary of clusters
            List of banned clusters (a list of keys)
    Output: The closest cluster, formatted as a list of keys
            The difference between the score of the checked cluster and the closest cluster
    """

    ## Find the centroid (average value) of the cluster to check against
    avg_value = find_centroid(check_key, check_dict)

    ## Create a list of all clusters in check_dict (a list of lists of active sites)
    cluster_list = [cluster for cluster in check_dict.values()]

    ## Checkpoint to make sure everything is working alright
    #  The test_for_validity value is a list of cluster averages
    #  The clust_len value is a total number of elements in cluster_list
    #  A flag can be set to true to display a graph of the clusters based on centroids
    test_for_validity, clust_len = visualize_h_cluster(cluster_list, False)

    ## Check to make sure the found centroid is also in the list of cluster averages
    if avg_value not in test_for_validity:
        print(avg_value)
        print(test_for_validity)

    ## Initializing a closest score
    closest_score = 1000000000 ## Going all in on massive number initialization

    ## Initializing the list of closest keys
    closest = []

    ## Initializing a closest difference
    closest_diff = 1000000000

    ## Loop over every key in the cluster dictionary to check the distance of every cluster
    for key in check_dict.keys():

        ## Make sure the key is not in the excluded list before continuing
        if key not in banned:

            ## Find the centroid for key in check_dict
            new_centroid = find_centroid(key, check_dict)

            ## If the new centroid is closer than the closest_score, keep the key, centroid, and difference
            if abs(avg_value - new_centroid) < abs(avg_value - closest_score):
                closest_score = new_centroid
                closest_diff = avg_value - new_centroid
                closest = [key]

            ## Also check to see if the new centroid is equal to the closest score, append the key
            elif abs(avg_value - new_centroid) == abs(avg_value - closest_score):
                closest.append(key)

    ## Return a list of the closest keys, and a float of the closest distance
    return (closest, closest_diff)

def find_centroid(check_key, check_dict):
    """
    Find the centroid of an active site (its average score)                                                              #

    Input:  Dictionary key for cluster to average
            Dictionary of clusters
    Output: The average score of the cluster
    """
    ## For every active site in the cluster, sum the score and divide by the number of sites
    avg_value = sum([value.score for value in check_dict[check_key]])/len(check_dict[check_key])
    return(avg_value)


def visualize_h_cluster(cluster_list, show_plot):

    """
    Find the cluster centroids by averaging, return averages and number of elements.
    Make a 1-d plot of the cluster centroids. The major differnece with this averaging method is that
    the input is a list instead of a dictionary. The output helps to check that the total cluster list is
    not growing or shrinking.

    This is not really needed anymore, but it really helped while building the hier clustering function
    to visualize that the right things were actually clustering. It also helped to determine if things
    were getting added or deleted to the overall active_site list.

    Input:  List of lists of active site clusters
            Boolean show plot
    Output: A list of the average values of the clusters
            An int of the number of total clustered elements
    """
    average_values = []
    clust_len = sum([len(cluster) for cluster in cluster_list])
    for cluster in cluster_list:
        average_values.append(sum([value.score for value in cluster])/len(cluster))
    if show_plot:
        pp.figure(figsize=(20,5))
        y_values = np.zeros_like(average_values)
        pp.plot(average_values, y_values, 'x')
        pp.show()
    return(average_values, clust_len)
