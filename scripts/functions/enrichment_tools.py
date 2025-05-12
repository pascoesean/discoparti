import numpy as np
import pandas as pd
rng = np.random.default_rng()
import seaborn as sns
import matplotlib.pyplot as plt

def calc_distance_to_archetype(archetype, cell):
    """
    Calculate the distance between a cell and an archetype.
    :param archetype: np.array of shape (n_pcs,)
    :param cell: np.array of shape (n_pcs,)
    :return: float
    """
    return np.linalg.norm(archetype - cell)

def calc_distance_to_archetypes(archetypes, cell):
    """
    Calculate the distance between a cell and each archetype.
    :param archetypes: np.array of shape (n_archetypes, n_pcs)
    :param cell: np.array of shape (n_pcs,)
    :return: np.array of shape (n_archetypes,)
    """
    distances = np.zeros(archetypes.shape[0])
    for i in range(archetypes.shape[0]):
        distances[i] = calc_distance_to_archetype(archetypes[i], cell)
    return distances

def calc_distance_to_archetypes_all(archetypes, cells):
    """
    Calculate the distance between each cell and each archetype.
    :param archetypes: np.array of shape (n_archetypes, n_pcs)
    :param cells: np.array of shape (n_cells, n_pcs)
    :return: np.array of shape (n_cells, n_archetypes)
    """
    distances = np.zeros((cells.shape[0], archetypes.shape[0]))
    for i in range(cells.shape[0]):
        distances[i] = calc_distance_to_archetypes(archetypes, cells[i])
    return distances


def calc_max_score(running_scores):
    # calculates the max (pos or negative) from a list of scores
    biggest_score = np.max(running_scores)
    smallest_score = np.min(running_scores)
    max_score = biggest_score
    if abs(smallest_score) > biggest_score:
        max_score = smallest_score

    return max_score


def do_pea(truth_series: np.ndarray, scores: np.ndarray, p: int = 1):
    # run phenotype enrichment to test enrichment of truth_series with corresponding scores 

    # returns: running_scores + max_enrichment_score
    
    abs_array = np.abs(scores)
    total_set_score = np.sum(np.power(abs_array[truth_series], p),axis=0)
    p_miss = 1 / (len(abs_array) - len(abs_array[truth_series]))
    
    # allocate array but add an extra one, start at zero, so we can backindex in loop
    running_scores = np.empty(len(truth_series) + 1) 
    running_scores[0] = 0

    for index, gene in enumerate(truth_series):
        # index will be one less than what we want
        scores_index = index + 1
        
        # if phenotype is true; add p_hit
        if truth_series[index]:
            p_hit = (abs_array[index]**p) / total_set_score
            # update next entry== last entry plus p_hit
            running_scores[scores_index] = p_hit + running_scores[scores_index - 1]

        else:
            running_scores[scores_index] = - p_miss + running_scores[scores_index - 1]


    # calc max (abs) score
    max_score = calc_max_score(running_scores)

    # we're done
    return running_scores[1:], max_score
        

# plot pea results
def plot_pea(running_scores, ax, phenotype: str = None, archetype: str = None, save_path=None):
    max_score = calc_max_score(running_scores)
    sns.lineplot(y = running_scores, x = range(len(running_scores)),
                 drawstyle='steps-post', color = "#39FF14", ax=ax)
    ax.hlines(y = 0, xmin=0, xmax = len(running_scores) - 0.9, linestyle = "--", color = "red")
    ax.hlines(y = max_score, xmin=0, xmax = len(running_scores) - 0.9, 
              linestyle = "dotted", color = "pink", label = f"max_score = {round(max_score,3)}")
    ax.set_ylim(top = abs(max_score) + 0.05, bottom = -abs(max_score) - 0.05)
    ax.set_xlim(left = 0, right = len(running_scores) - 0.9)
    ax.set_xlabel("Rank in distance to archetype " + str(archetype + 1))
    ax.set_ylabel("Enrichment Score")
    ax.set_title("Enrichment for Phenotype " + str(phenotype))
    ax.legend()
    if save_path is not None:
        plt.savefig(save_path + str(phenotype) + "_" + 'archetype_' + str(archetype + 1) + "_pea.png", dpi=300)
    

def pea_from_labels(labels, dists, label_of_interest, plot_save_path=None):
    # do phenotype enrichment analysis on the given label of interest
    truth = labels == label_of_interest
    dist_scores = dists - np.mean(dists, axis=0)  # NEED CENTERING!!

    # place to store the max scores
    pea_max_scores = np.zeros((dist_scores.shape[1],))

    # need to do separately for each archetype. idk how to do better than this </3
    for i in range(dist_scores.shape[1]):
        arch_i = dist_scores[:,i]
        indices = np.argsort(arch_i, axis=0) 
        sorted_truths = truth[indices]
        sorted_dists = arch_i[indices]
        running_score, max_score = do_pea(sorted_truths, sorted_dists, p=1)
        pea_max_scores[i] = max_score

        # 1-800 are we plotting?
        if plot_save_path is not None:
            fig, ax = plt.subplots(1, 1, figsize=(10, 5))
            plot_pea(running_score, ax, phenotype=label_of_interest, archetype=i, save_path=plot_save_path)
            plt.close(fig)

    return pea_max_scores

def calc_differential_enrichment(cell_labels, archetypes, scores, num_fake:int = 10000, plot_save_path=None):
    """
    Calculate the differential enrichment of each archetype in each cell type.
    :param cell_labels: np.array of shape (n_cells,) that has the labels of interest for each cell
    :return: np.array of shape (n_cell_types, n_archetypes)
    """
    # get unique labels
    unique_labs = np.unique(cell_labels)

    # get the distances
    dists = calc_distance_to_archetypes_all(archetypes, scores)

    # create a matrix to hold the differential enrichment scores: should be classes by archetypes
    true_enrich = np.zeros((len(unique_labs), archetypes.shape[0]))
    enrich_zscore = np.zeros((len(unique_labs), archetypes.shape[0]))

    # loop through each cell type and calculate the differential enrichment
    for i, label in enumerate(unique_labs):
        # which are the (true) relevant label right now?
        true_enrichment = pea_from_labels(cell_labels, dists, label, plot_save_path=plot_save_path)
        abs_true = np.abs(true_enrichment)

        # generate num_fake fake enrichments
        enrichments = np.zeros((num_fake, archetypes.shape[0]))
        for ii in range(num_fake):
            # generate a random label
            fake_labels = rng.permutation(cell_labels)
            # calculate the enrichment for the fake label
            enrichments[ii] = pea_from_labels(fake_labels, dists, label)
        
        # calculate the z score based on fake enrichments: how far away from the null mean are we?
        mus = np.mean(np.abs(enrichments), axis=0)
        sigmas = np.std(np.abs(enrichments), axis=0)
        # calculate the z score

        # save the S scores:
        enrich_zscore[i] = (abs_true - mus) / sigmas
        true_enrich[i] = true_enrichment
    
    return true_enrich, enrich_zscore, unique_labs, enrichments

def calc_enrichment_from_files(archetypes_filepath, scores_filepath, labels_filepath, num_fake=10000, plot_save_path=None, pandas=False):
    archs = pd.read_csv(archetypes_filepath, header=None)
    num_archetypes = len(archs)
    archs.columns = [f'PC_{i+1}' for i in range(num_archetypes - 1)]
    archs = archs.to_numpy()

    # once we know how many archetypes we want, we know how many PCs we're considering. lets load in those scores:
    scores = pd.read_csv(scores_filepath, usecols= list(range(num_archetypes - 1)),
                                names = [f'PC_{i+1}' for i in range(num_archetypes - 1)]).to_numpy()
    # lets smack some labels on there
    cell_labels = pd.read_csv(labels_filepath).to_numpy()[:,1]

    if pandas:
        # if we want to return a pandas dataframe
        nes_scores, zscores, labels, _ = calc_differential_enrichment(cell_labels, archs, scores, num_fake=num_fake, plot_save_path=plot_save_path)
        nes_pd = pd.DataFrame(nes_scores, index = labels, 
                      columns= [f'archetype_{i+1}' for i in range(nes_scores.shape[1])]
                      ).reset_index().melt(id_vars = 'index', 
                                           var_name = 'archetype', 
                                           value_name = 'nes_scores').rename(columns = {'index': 'labels'})
        zscores_pd = pd.DataFrame(zscores, index = labels, 
                      columns= [f'archetype_{i+1}' for i in range(zscores.shape[1])]
                      ).reset_index().melt(id_vars = 'index', 
                                           var_name = 'archetype', 
                                           value_name = 'zscores').rename(columns = {'index': 'labels'})
        
        return pd.merge(zscores_pd, nes_pd, on = ['labels', 'archetype'])


    return calc_differential_enrichment(cell_labels, archs, scores, num_fake=num_fake, plot_save_path=plot_save_path)


# helper function for getting distances for gsea:

def calc_distances_from_files(archetypes_filepath, scores_filepath, labels_filepath, num_fake=10000):
    archs = pd.read_csv(archetypes_filepath, header=None)
    num_archetypes = len(archs)
    archs.columns = [f'PC_{i+1}' for i in range(num_archetypes - 1)]
    archs = archs.to_numpy()

    # once we know how many archetypes we want, we know how many PCs we're considering. lets load in those scores:
    scores = pd.read_csv(scores_filepath, usecols= list(range(num_archetypes - 1)),
                                names = [f'PC_{i+1}' for i in range(num_archetypes - 1)]).to_numpy()
    # lets smack some labels on there
    cell_labels = pd.read_csv(labels_filepath).to_numpy()[:,1]

    unique_labs = np.unique(cell_labels)

    # get the distances
    dists = calc_distance_to_archetypes_all(archs, scores)

    return dists, cell_labels, unique_labs
