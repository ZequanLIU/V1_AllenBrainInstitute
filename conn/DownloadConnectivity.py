import pandas as pd
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
from utils_connectivity import *
import re

fit_with_distance=False # Whether to fit a Gaussian function to describe the probability of connection.
# Doesn't fit well for all cases. See the figures in plots/

# Gaussian model parameters:
# self.pmax * np.exp(-x**2 / (2 * self.size**2))
# pmax: float Maximum connection probability(at 0 intersomatic distance)
# size: float Gaussian sigma
db = SynphysDatabase.load_current('medium')  # small has connectivity properties

mouse_pairs = db.pair_query(project_name=["mouse V1 coarse matrix", "mouse V1 pre-production"],
                            synapse_probed=True  # filter: only cell pairs that were checked for connectivity
                            ).all()

print("loaded %d cell pairs" % len(mouse_pairs))

cell_class_criteria = {
    'VisL1_exc': {'cell_class': 'ex', 'cortical_layer': '1'},
    'VisL1_inh': {'cell_class': 'in', 'cortical_layer': '1'},
    'VisL2_3_exc': {'cell_class': 'ex', 'cortical_layer': '2/3'},
    'VisL2_3_inh': {'cell_class': 'in', 'cortical_layer': '2/3'},
    'VisL4_exc': {'cell_class': 'ex', 'cortical_layer': '4'},
    'VisL4_inh': {'cell_class': 'in', 'cortical_layer': '4'},
    'VisL5_exc': {'cell_class': 'ex', 'cortical_layer': '5'},
    'VisL5_inh': {'cell_class': 'in', 'cortical_layer': '5'},
    'VisL6_exc': {'cell_class': 'ex', 'cortical_layer': ('6a','6b')},
    'VisL6_inh': {'cell_class': 'in', 'cortical_layer': ('6a','6b')}
}

cell_classes = {name: CellClass(name=name, **criteria) for name, criteria in cell_class_criteria.items()}

# Group all cells by selected classes
cell_groups = classify_cells(cell_classes.values(), pairs=mouse_pairs)

# Group pairs into (pre_class, post_class) groups
pair_groups = classify_pairs(mouse_pairs, cell_groups)

# analyze matrix elements
results = pd.DataFrame(measure_connectivity(pair_groups, dist_measure='lateral_distance', fit_model=GaussianModel))

gaussianFits = {}
Pops = ['VisL1_exc','VisL2_3_exc','VisL4_exc','VisL5_exc','VisL6_exc',
        'VisL1_inh','VisL2_3_inh','VisL4_inh','VisL5_inh','VisL6_inh']

for i in results.columns:
    preSyn = str(i[0])
    postSyn = str(i[1])
    gaussianFits[preSyn + '->' + postSyn] = {}
    gaussianFits[preSyn+'->'+postSyn]['synapses'] = (results[i]['connectivity_fit'].pmax, results[i]['connectivity_fit'].size*1e6)
    gaussianFits[preSyn+'->'+postSyn]['gapJunctions'] = (results[i]['gap_fit'].pmax, results[i]['gap_fit'].size*1e6)
fitFile = pd.DataFrame(gaussianFits)
fitFile.to_csv('GaussianFits.txt')

if fit_with_distance:
    for i in results.columns:
        # Synapses
        x_probed = results[i]['probed_distances']
        # Created a mask to indicate which pairs are connected and which are not
        conn = [False] * len(x_probed)
        for ii in range(len(results[i]['probed_pairs'])):
            for jj in range(len(results[i]['connected_pairs'])):
                if results[i]['probed_pairs'][ii] == results[i]['connected_pairs'][jj]:
                    conn[ii] = True

        conn = np.asarray(conn)
        x_probed = np.asarray(x_probed)
        try:
            fig, ax = plt.subplots(1, 1, figsize=(15, 5))
            fit = results[i]['connectivity_fit']
            show_connectivity_profile(x_probed, conn, ax, fit)
            plt.title(str(i) + f" Gaussian fit pmax={fit.pmax:0.2f}, σ={fit.size*1e6:0.2f}μm")
            plt.savefig('plots/Synapses'+str(i) + '.png')
        except:
            1

        # Gap junctions
        x_probed = results[i]['gap_probed_distances']
        # Created a mask to indicate which pairs are connected and which are not
        conn = [False] * len(x_probed)
        for ii in range(len(results[i]['gaps_probed'])):
            for jj in range(len(results[i]['gap_pairs'])):
                if results[i]['gaps_probed'][ii] == results[i]['gap_pairs'][jj]:
                    conn[ii] = True

        conn = np.asarray(conn)
        x_probed = np.asarray(x_probed)
        try:
            fig, ax = plt.subplots(1, 1, figsize=(15, 5))
            fit = results[i]['gap_fit']
            show_connectivity_profile(x_probed, conn, ax, fit)
            plt.title(str(i) + f" Gaussian fit pmax={fit.pmax:0.2f}, σ={fit.size*1e6:0.2f}μm")
            plt.savefig('plots/GapJunctions'+str(i) + '.png')
        except:
            1

class_labels = {
    'VisL1_exc': 'L1 exc',
    'VisL1_inh': 'L1 inh',
    'VisL2_3_exc': 'L2/3 exc',
    'VisL2_3_inh': 'L2/3 inh',
    'VisL4_exc': 'L4 exc',
    'VisL4_inh': 'L4 inh',
    'VisL5_exc': 'L5 exc',
    'VisL5_inh': 'L5 inh',
    'VisL6_exc': 'L6 exc',
    'VisL6_inh': 'L6 inh',
    }

# analyze matrix elements
results = measure_connectivity(pair_groups, sigma=130e-6, dist_measure='lateral_distance')

ConnectionMatrix = {}
for i in results.keys():
    preSyn = str(i[0])
    postSyn = str(i[1])
    ConnectionMatrix[preSyn + '->' + postSyn] = {}
    ConnectionMatrix[preSyn+'->'+postSyn]['synapses'] = results[i]['adjusted_connectivity']
    ConnectionMatrix[preSyn+'->'+postSyn]['gapJunctions'] = results[i]['adjusted_gap_junction']
fitFile = pd.DataFrame(ConnectionMatrix)
fitFile.to_csv('ConnectionMatrix.txt')


# define a colormap and log normalization used to color the heatmap
norm = matplotlib.colors.Normalize(vmin=0.01, vmax=0.5)
cmap = matplotlib.cm.get_cmap('plasma')
fig, ax = plt.subplots(figsize=(15, 15))
# finally, draw the colormap using the provided function:
im, cbar, labels = show_connectivity_matrix(
    ax=ax,
    results=results,
    pre_cell_classes=cell_classes.values(),
    post_cell_classes=cell_classes.values(),
    class_labels=class_labels,
    cmap=cmap,
    norm=norm,
    distance_adjusted=True,
    ctype='chemical',
    alpha=True
)
plt.savefig('plots/'+"Connectivity_SynProbed.png", bbox_inches='tight')
# analyze matrix elements
results = measure_connectivity(pair_groups, sigma=75e-6, dist_measure='lateral_distance')
fig, ax = plt.subplots(figsize=(15, 15))
# finally, draw the colormap using the provided function:
im, cbar, labels = show_connectivity_matrix(
    ax=ax,
    results=results,
    pre_cell_classes=cell_classes.values(),
    post_cell_classes=cell_classes.values(),
    class_labels=class_labels,
    cmap=cmap,
    norm=norm,
    distance_adjusted=True,
    ctype='electrical',
    alpha=True
)
plt.savefig('plots/'+"Gaps_SynProbed.png", bbox_inches='tight')
