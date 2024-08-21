import numpy as np
from aisynphys.database import SynphysDatabase
from aisynphys.cell_class import CellClass, classify_cells, classify_pairs
from aisynphys.connectivity import measure_connectivity, pair_was_probed
from aisynphys.connectivity import GaussianModel
from aisynphys.ui.notebook import show_connectivity_profile

import matplotlib.colors, matplotlib.cm
import matplotlib.pyplot as plt
from aisynphys.ui.notebook import show_connectivity_matrix

def load_database():
    # Download and cache the sqlite file for the requested database
    #   (for available versions, see SynphysDatabase.list_versions())
    # Small (~170 MB) Medium (~7 GB) Full (~160 GB)
    db = SynphysDatabase.load_current('medium')  # medium in order to have gap junction info
    return db

def get_connectivity_values(db, cell_class_criteria, distance='lateral_distance', sigma=100e-6):
    # Load all cell pairs associated with mouse V1 projects
    mouse_pairs = db.pair_query(project_name=["mouse V1 coarse matrix", "mouse V1 pre-production"]).all()

    cell_classes = {name: CellClass(name=name, **criteria) for name, criteria in cell_class_criteria.items()}

    # Group all cells by selected classes
    cell_groups = classify_cells(cell_classes.values(), pairs=mouse_pairs)

    # Group pairs into (pre_class, post_class) groups
    pair_groups = classify_pairs(mouse_pairs, cell_groups)

    # analyze matrix elements
    results = measure_connectivity(pair_groups, sigma=sigma, dist_measure=distance)

    return results

def connectivity_model(db, distance='lateral_distance', cellMask='cell_class', preCell='ex', postCell='in',
                       lateralDistance=500e-6, connection_type='synapse'):
    # connection_type = 'synapse' or 'electrical'
    # The analysis above gives us an estimate of the relatives connectivities between cell types, but leaves out some
    # important details. In particular, we know that the probability of finding a connection between any two cells is
    # strongly related to the spatial relationship between the cells and their axo-dendritic morphology.
    # As an approximation, we think of cell morphology as being cylindrically symmetrical around the axis perpendicular to
    # the cortical surface. This means that the likelihood of two cells being connected by a synapse is strongly related to
    # the lateral distance between their cell bodies (the distance parallel to the cortical surface).

    # Gaussian model parameters:
    # self.pmax * np.exp(-x**2 / (2 * self.size**2))
    # pmax: float Maximum connection probability(at 0 intersomatic distance)
    # size: float Gaussian sigma

    mouse_pairs = db.pair_query(
        experiment_type='standard multipatch',  # filter: just multipatch experiments
        species='mouse',  # filter: only mouse data
        synapse_probed=True,  # filter: only cell pairs that were checked for connectivity
        preload=['cell']  # include tables that describe cell properties
    ).dataframe()

    mask = (
            (mouse_pairs['pre_cell.'+cellMask] == preCell) &
            (mouse_pairs['post_cell.'+cellMask] == postCell) &
            (mouse_pairs['pair.lateral_distance'] < lateralDistance) # To ensure we are looking at a cortical column
    )

    # distance = 'pair.distance' or 'pair.lateral_distance' or 'pair.vertical_distance'
    # [print(i) for i in mouse_pairs.keys()]
    # print('cortical_layer', mouse_pairs['post_cortical_cell_location.cortical_layer'][0])
    # print(mouse_pairs['post_cell.meta'][0])
    # print(mouse_pairs['post_cell.target_layer'][0])
    # print(mouse_pairs['pre_cell.target_layer'][0])
    # quit()
    x_probed = mouse_pairs[mask]['pair.'+distance].to_numpy(dtype=float)
    conn = mouse_pairs[mask]['pair.has_'+connection_type].to_numpy(dtype=bool)
    fit = GaussianModel.fit(x_probed, conn)

    return fit

