# Some utility classes to represent a PDB structure

class Atom:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type):
        self.type = type
        self.coords = (0.0, 0.0, 0.0)

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.type

class Residue:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type, number):
        self.type = type
        self.number = number
        self.atoms = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return "{0} {1}".format(self.type, self.number)

class ActiveSite:
    """
    A simple class for an active site
    """

    def __init__(self, name):
        self.name = name
        self.residues = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name

def active_site_score(active_sites):
    """
    Determine active_site.score. Normalize based on mean and standard deviation.

    Initially active sites are scored on multiple metrics, lysine positions, arginine positions,
    and overall hydrophobicity.

    These initial scores are then normalized around 0 and added together to get active_site.score

    This function directly writes to active_site.score.

    Input:  A list of active sites
    Output: None
    """

    if len(active_sites) <= 0:
        return

    ## The hydrophobicity dictionary
    score_metric_hydrophobicity = {'ARG': -4.5, 'HIS': -3.2, 'LYS': -3.9,
    'ASP': -3.5, 'GLU': -3.5, 'CYS': 2.5, 'GLY': -0.4, 'PRO': -1.6, 'ALA': 1.8,
    'VAL': 4.2, 'ILE': 4.5, 'LEU': 3.8, 'MET': 1.9, 'PHE': 2.8, 'TYR': -1.3,
    'TRP': -0.9, 'SER': -0.8, 'THR': -0.7, 'ASN': -3.5, 'GLN': -3.5}

    ## Determining lysine positions
    for active_site in active_sites:

        ## Initial lysine vector of 0
        lysine_coord_sum = [0.0, 0.0, 0.0]
        for residue in active_site.residues:
            if residue.type == "LYS":
                for atom in residue.atoms:
                    if atom.type == 'NZ':

                        ## Every time a lysine is found, add the vector
                        lysine_coord_sum[0] += atom.coords[0]
                        lysine_coord_sum[1] += atom.coords[1]
                        lysine_coord_sum[2] += atom.coords[2]

        ## Determine the magnitude of the summed lysine vector, and save to active site
        lysine_sum_magnitude = pow(sum([pow(pos, 2) for pos in lysine_coord_sum]), 0.5)
        active_site.lys_score = lysine_sum_magnitude

    ## Determine the mean and standard deviation of the lysine magnitudes
    lys_m = sum([active_site.lys_score for active_site in active_sites])/len(active_sites)
    lys_s = sum([abs(active_site.lys_score - lys_m) for active_site in active_sites])/len(active_sites)


    ## Repeat the same process for arginines
    for active_site in active_sites:
        arg_coord_sum = [0.0, 0.0, 0.0]
        for residue in active_site.residues:
            if residue.type == "ARG":
                for atom in residue.atoms:
                    if atom.type == 'CZ':
                        arg_coord_sum[0] += atom.coords[0]
                        arg_coord_sum[1] += atom.coords[1]
                        arg_coord_sum[2] += atom.coords[2]
        arg_sum_magnitude = pow(sum([pow(pos, 2) for pos in arg_coord_sum]), 0.5)
        active_site.arg_score = arg_sum_magnitude

    arg_m = sum([active_site.arg_score for active_site in active_sites])/len(active_sites)
    arg_s = sum([abs(active_site.arg_score - arg_m) for active_site in active_sites])/len(active_sites)


    ## Hydrophobic score determination usint the dictionary
    for active_site in active_sites:
        active_site.hyd_score = sum([score_metric_hydrophobicity[residue.type] for residue in active_site.residues])

    hyd_m = sum([active_site.hyd_score for active_site in active_sites])/len(active_sites)
    hyd_s = sum([abs(active_site.hyd_score - hyd_m) for active_site in active_sites])/len(active_sites)


    ## For every active site, normalize the scores, then add them together into the master score active_site.score
    for active_site in active_sites:
        active_site.lys_score = (active_site.lys_score - lys_m)/lys_s
        active_site.arg_score = (active_site.arg_score - arg_m)/arg_s
        active_site.hyd_score = (active_site.hyd_score - hyd_m)/hyd_s

        active_site.score = active_site.lys_score + active_site.arg_score + active_site.hyd_score


def cluster_score(cluster_list):
    """
    Determine a score for the cluster based on the active site scores inside of the cluster.
    The cluster score is determined by finding the distance of every active site to the cluster center,
    and summing those distances. Tighter clusters will have smaller sums.

    It is a very uncomplicated way of getting a cluster score, but it seems to work.

    Input: A list of clusters (the output from the clustering functions)
    Output: sum of distances from cluster mean
    """
    ## For every cluster, find the distance of each active_site.score from the mean
    dist_list = []
    for cluster in cluster_list:
        m = sum([active_site.score for active_site in cluster])/len(cluster)

        ## Append the sum of the distances to dist_list
        dist_list.append(sum([abs(active_site.score - m) for active_site in cluster]))

    ## Sum up the dist list
    return sum(dist_list)
