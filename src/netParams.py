from netpyne.batchtools import specs
from cfg import cfg
from neuron_morphology.features.layer.reference_layer_depths import DEFAULT_HUMAN_MTG_REFERENCE_LAYER_DEPTHS, DEFAULT_MOUSE_REFERENCE_LAYER_DEPTHS
import glob
import pandas as pd

cfg.update_cfg()
netParams = specs.NetParams()  # object of class NetParams to store the network parameters

###############################################################################
## NETWORK PARAMETERS
###############################################################################
# General network geometry
def network_geometry(netParams):
    #------------------------------------------------------------------------------
    # General network parameters
    #------------------------------------------------------------------------------
    netParams.scale = cfg.scale # Scale factor for number of cells
    netParams.sizeX = cfg.sizeX # x-dimension (horizontal length) size in um
    netParams.sizeY = cfg.sizeY # y-dimension (vertical height or cortical depth) size in um
    netParams.sizeZ = cfg.sizeZ # z-dimension (horizontal depth) size in um
    netParams.shape = 'cylinder' # cylindrical (column-like) volume
# General connectivity params
def connectivity_params(netParams):
    #------------------------------------------------------------------------------
    # General connectivity parameters
    #------------------------------------------------------------------------------
    netParams.defaultThreshold = 0.0
    netParams.defineCellShapes = True       # sets 3d geometry aligned along the y-axis
    netParams.defaultDelay = 2.0 # default conn delay (ms)
    netParams.propVelocity = 500.0 # propagation velocity (um/ms)
    netParams.probLambda = 100.0  # length constant (lambda) for connection probability decay (um)
# Definition of layers boundaries (normalized)
def layers_geometry():
    # Normalized layer distances
    # TODO: We may need to use the not normalized values when using a connectivity profile dependent on distance,
    #  or we will need to normalize the connectivity profile!
    layerNormalized = {}
    MaxPiaDistance = DEFAULT_MOUSE_REFERENCE_LAYER_DEPTHS['6b'][1]
    for layerLabel, object in DEFAULT_MOUSE_REFERENCE_LAYER_DEPTHS.items():
        layerNormalized[layerLabel] = [object.pia_side/MaxPiaDistance, object.wm_side/MaxPiaDistance]
    del layerNormalized['wm'] # White matter depth
    layerNormalized['LGN'] = [2.0, 3.0]
    # Normalize the key name in order to make it standard in the model
    layerNormalized['VisL1'] = layerNormalized.pop('1')
    layerNormalized['VisL2_3'] = layerNormalized.pop('2/3')
    layerNormalized['VisL4'] = layerNormalized.pop('4')
    layerNormalized['VisL5'] = layerNormalized.pop('5')
    layerNormalized['VisL6'] = layerNormalized.pop('6a')

    return layerNormalized

###############################################################################
## CELL TYPES AND POPULATION DEFINITIONS
###############################################################################

# Layer 1 definition. Components: i1Htr3a (2 models -> heterogeneous population)
def L1_pop(netParams, inh_cellNumbers={'Htr3a': 999}, model='Peri', layer='VisL1'):
    layerBoundaries=layers_geometry()
    for cellType, cellNumber in inh_cellNumbers.items():
        cellFiles = glob.glob('../cells/NetPyNE_cells/%s/%s%s*_%s_cellParams.pkl' % (layer,layer,cellType, model))
        k=0
        for neuronModel in cellFiles:
            label = layer + cellType
            netParams.loadCellParamsRule(label=label+'_'+str(k), fileName=neuronModel)
            netParams.renameCellParamsSec(label=label+'_'+str(k), oldSec='soma_0', newSec='soma')
            k+=1
        netParams.popParams[label]={'cellModel': model, 'cellType': label, 'numCells': cellNumber,
                                                   'diversity': True, 'ynormRange': layerBoundaries[layer]}
# Layer 2 definition. Components: i1Htr31 (8 models), Pvalb (3 models), i2/3Sst (4 models), i2/3Htr3a, E2/3 (3 models)
def L23_pop(netParams, inh_cellNumbers={'Pvalb': 640, 'Sst': 464,'Htr3a': 1107}, exc_cellNumbers={'E': 12689},
            model='Peri', layer='VisL2_3'):
    layerBoundaries = layers_geometry()
    for cellType, cellNumber in inh_cellNumbers.items():
        cellFiles = glob.glob('../cells/NetPyNE_cells/%s/%s%s*_%s_cellParams.pkl' % (layer,layer,cellType, model))
        k=0
        for neuronModel in cellFiles:
            label = layer + cellType
            netParams.loadCellParamsRule(label=label+'_'+str(k), fileName=neuronModel)
            netParams.renameCellParamsSec(label=label+'_'+str(k), oldSec='soma_0', newSec='soma')
            k+=1
        netParams.popParams[label]={'cellModel': model, 'cellType': label, 'numCells': cellNumber,
                                                   'diversity': True, 'ynormRange': layerBoundaries[layer]}

    for cellType, cellNumber in exc_cellNumbers.items():
        cellFiles = glob.glob('../cells/NetPyNE_cells/%s/%s%s*_%s_cellParams.pkl' % (layer,layer,cellType, model))
        k=0
        for neuronModel in cellFiles:
            label = layer + cellType
            netParams.loadCellParamsRule(label=label+'_'+str(k), fileName=neuronModel)
            netParams.renameCellParamsSec(label=label+'_'+str(k), oldSec='soma_0', newSec='soma')
            k+=1
        netParams.popParams[label]={'cellModel': model, 'cellType': label, 'numCells': cellNumber,
                                                   'diversity': True, 'ynormRange': layerBoundaries[layer]}

###############################################################################
## NETWORK CONNECTIVITY
###############################################################################
# Synaptic mechanisms definition
def initialize_synaptic_mechs(netParams):
    # Initialize synaptic mechanisms. We need to extract from bibliography which synaptic types are present between
    # each cell types. For now we use AMPA as excitatory and GABAA as inhibitory
    netParams.synMechParams['AMPA'] = {'mod': 'MyExp2SynBB', 'tau1': 0.05, 'tau2': 5.3, 'e': 0}
    netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB','tau1NMDA': 15, 'tau2NMDA': 150, 'r': 1, 'e': 0}
    netParams.synMechParams['GABAB'] = {'mod': 'MyExp2SynBB', 'tau1': 3.5, 'tau2': 260.9, 'e': -93}
    netParams.synMechParams['GABAA'] = {'mod': 'MyExp2SynBB', 'tau1': 0.07, 'tau2': 18.2, 'e': -80}
    netParams.synMechParams['GABAs'] = {'mod': 'MyExp2SynBB', 'tau1': 2, 'tau2': 100, 'e': -80}
    netParams.synMechParams['GABAss'] = {'mod': 'MyExp2SynBB', 'tau1': 200, 'tau2': 400, 'e': -80}
    netParams.synMechParams['gap'] = {'mod': 'ElectSyn', 'g': 1}

# LGN inputs definition
def LGN_pop(netParams, label='LGN', layer='LGN', numCells=2, spkTimes=[[0,100], [100,200]]):
    # Initialize stimuli. Is the LGN input to all layers
    layerBoundaries=layers_geometry()
    # Define the input populations
    # TODO: Run the LGN model from Allen Brain to extract spike times from a movie or image. Adjust conenctivity
    #  parameters
    netParams.popParams[label] = {'cellModel': 'VecStim', 'numCells': numCells, 'spkTimes': spkTimes,
                                  'ynormRange': layerBoundaries[layer]}
def LGN_stimuli(netParams, popInput='LGN', popTarget=['VisL1Htr3a']):
    for postPop in popTarget:
        netParams.connParams[popInput + '->' + postPop] = {
            'preConds': {'popLabel': 'LGN'},
            'postConds': {'popLabel': postPop},
            'weight': 0.05,
            'sec': 'soma',
            'probability': 1,
            'delay': 1,
            'loc': 0.5,
            'synMech': 'AMPA'}
# LGN inputs
def LGN_inputs(netParams, popInput='LGN', popTarget=['VisL1Htr3a', 'VisL2_3Pvalb', 'VisL2_3Sst', 'VisL2_3Htr3a', 'VisL2_3E']):
    for popPost in popTarget:
        LGN_stimuli(netParams, popInput=popInput, popTarget=popPost)
# Background activity in form of NetStim
# Population definition
def bkg_pop(netParams, popInput='Bkg', interval=150, number=int((1000 / 150) * cfg.duration), start=0, noise=0.5):
    # Random spike trains to population
    netParams.stimSourceParams[popInput] = {'type': 'NetStim', 'interval': interval,
                                            'number': number, 'start': start, 'noise': noise}
# General input definition
def bkg_input_def(netParams, popInput='Bkg', popTargetLabel='VisL1Htr3a', conds={'pop': 'VisL1Htr3a'}, sec='soma',
                  loc=0.5, weight=2e-3, delay=0.2, synMech='AMPA'):
    netParams.stimTargetParams[popInput+'->'+popTargetLabel] = {
        'source': popInput,
        'conds': conds,
        'sec': sec,
        'loc': loc,
        'weight': weight,
        'delay': delay,
        'synMech': synMech}
# Background inputs
def bkg_input(netParams, popTarget=['VisL1Htr3a', 'VisL2_3Pvalb', 'VisL2_3Sst', 'VisL2_3Htr3a', 'VisL2_3E'],
              sec='soma', loc=0.5, weight=10e-3, delay=0.2):
    for pop in popTarget:
        conds = {'popLabel': pop}
        bkg_input_def(netParams, popInput='Bkg', popTargetLabel=pop, conds=conds, sec=sec, loc=loc, weight=weight,
                      delay=delay, synMech='AMPA')

# Load connectivity matrix from Allen Brain dataset
def connection_probability():
    GaussianFit = pd.read_csv('../conn/GaussianFits.txt', index_col=0)
    ConnectionProbability = pd.read_csv('../conn/ConnectionMatrix.txt', index_col=0)
    return GaussianFit, ConnectionProbability
# Define the list of connections that are "well fitted" (eye-balled criterion)
def FittedList():
    FittedListSyn=[
        # L2/3-> rest of layers
        'VisL2_3_inh->VisL1_inh', 'VisL2_3_inh->VisL2_3_inh', 'VisL2_3_inh->VisL4_inh', 'VisL2_3_inh->VisL5_inh', # inh->inh
        'VisL2_3_inh->VisL2_3_exc', 'VisL2_3_inh->VisL5_exc',  # inh->exc
        'VisL2_3_exc->VisL2_3_inh', 'VisL2_3_exc->VisL4_inh',  # exc->inh
        'VisL2_3_exc->VisL4_exc', 'VisL2_3_exc->VisL5_exc',  # exc->exc
        # L4-> rest of layers
        'VisL4_inh->VisL2_3_inh', 'VisL4_inh->VisL4_inh', 'VisL4_inh->VisL5_inh', 'VisL4_inh->VisL6_inh',  # inh->inh
        'VisL4_inh->VisL2_3_exc', 'VisL4_inh->VisL4_exc', 'VisL4_inh->VisL5_exc',  # inh->exc
        'VisL4_exc->VisL2_3_inh', 'VisL4_exc->VisL4_inh', 'VisL4_exc->VisL5_inh',  # exc->inh
        'VisL4_exc->VisL2_3_exc', 'VisL4_exc->VisL4_exc', 'VisL4_exc->VisL5_exc'  # exc->exc
        # L5-> rest of layers
        'VisL5_inh->VisL2_3_inh', 'VisL5_inh->VisL4_inh', 'VisL5_inh->VisL5_inh',  # inh->inh
        'VisL5_inh->VisL2_3_exc', 'VisL5_inh->VisL4_exc', 'VisL5_inh->VisL5_exc',  # inh->exc
        'VisL5_exc->VisL2_3_inh', 'VisL5_exc->VisL4_inh', 'VisL5_exc->VisL5_inh', 'VisL5_exc->VisL6_inh',  # exc->inh
        'VisL5_exc->VisL4_exc', 'VisL5_exc->VisL5_exc',  # exc->exc
        # L6-> rest of layers
        'VisL6_inh->VisL2_3_inh', 'VisL6_inh->VisL4_inh', 'VisL6_inh->VisL6_inh',  # inh->inh
        'VisL6_inh->VisL5_exc',  # inh->exc
        'VisL6_exc->VisL5_inh', 'VisL6_exc->VisL5_inh',  # exc->inh
    ]
    # Gap junctions are bidirectional, thus we can use the same value for the symmetric connection
    FittedListGaps = ['VisL2_3_inh->VisL2_3_inh', 'VisL2_3_inh->VisL4_inh', 'VisL2_3_inh->VisL5_inh',
                      'VisL4_inh->VisL4_inh', 'VisL4_inh->VisL5_inh',
                      'VisL5_inh->VisL5_inh', 'VisL5_inh->VisL6_inh',
                      'VisL6_inh->VisL6_inh']
    return FittedListSyn, FittedListGaps
# General connection definition between two cell populations
def popConnection(netParams, popInput='VisL1Htr3a', popTarget='VisL1Htr3a', sec='soma', loc=0.5, weight=5e-3, delay=1,
              synMech='AMPA', probability=0.1, gap=False):
    netParams.connParams[popInput + '->' + popTarget] = {
        'preConds': {'popLabel': popInput},
        'postConds': {'popLabel': popTarget},
        'weight': weight,
        'sec': sec,
        'probability': probability, # Could be a number (connection probability) or a function if entered as string
        'delay': delay,
        'loc': loc,
        'synMech': synMech,
        'gapJunction': gap, # NetPyNE auto makes these junctions bidirectional so we don't want to read the direction in twice
    }
# Define excitatory connections
def exc_connections(netParams, popInput=['VisL2_3E'],
                    popTarget=['VisL1Htr3a', 'VisL2_3Pvalb', 'VisL2_3Sst', 'VisL2_3Htr3a', 'VisL2_3E'],
                    sec='soma', loc=0.5, weight=5e-3, delay=1, synMech='AMPA'):
    # Load excitatory connectivity profiles, either fitted with distance or probability of connection
    # For now, the distinction is between layers, and between excitatory and inhibitory neurons.
    # This could be changed in the future
    # TODO: Change the probability of connection to the one fitted (if in FittedList)
    #  or to the experimental probability of connection
    GaussianFit, ConnectionProbability = connection_probability()
    FittedListSyn, FittedListGaps = FittedList()

    for popsIn in popInput:
        for popsOut in popTarget:
            popConnection(netParams, popInput=popsIn, popTarget=popsOut, sec=sec, loc=loc, weight=weight,
                          delay=delay, synMech=synMech, probability=0.5)
# Define inhibitory connections
def inh_connections(netParams, popInput=['VisL1Htr3a', 'VisL2_3Pvalb', 'VisL2_3Sst', 'VisL2_3Htr3a'],
                    popTarget=['VisL1Htr3a', 'VisL2_3Pvalb', 'VisL2_3Sst', 'VisL2_3Htr3a', 'VisL2_3E'],
                    sec='soma', loc=0.5, weight=5e-3, delay=1, synMech='GABAA'):
    # Load inhibitory connectivity profiles, either fitted with distance or probability of connection
    # For now, the distinction is between layers, and between excitatory and inhibitory neurons.
    # This could be changed in the future
    # TODO: Change the probability of connection to the one fitted (if in FittedList)
    #  or to the experimental probability of connection
    GaussianFit, ConnectionProbability = connection_probability()
    FittedListSyn, FittedListGaps = FittedList()

    for popsIn in popInput:
        for popsOut in popTarget:
            popConnection(netParams, popInput=popsIn, popTarget=popsOut, sec=sec, loc=loc, weight=weight,
                          delay=delay, synMech=synMech, probability=0.5)
# Define gap junctions connections
def gap_connections(netParams, popInput=['VisL1Htr3a', 'VisL2_3Pvalb', 'VisL2_3Sst', 'VisL2_3Htr3a'],
                    popTarget=['VisL1Htr3a', 'VisL2_3Pvalb', 'VisL2_3Sst', 'VisL2_3Htr3a'],
                    sec='soma', loc=0.5, weight=1e-3, delay=0, synMech='gap'):
    # Load inhibitory connectivity profiles, either fitted with distance or probability of connection
    # For now, the distinction is between layers, and between excitatory and inhibitory neurons.
    # This could be changed in the future
    # TODO: Change the probability of connection to the one fitted (if in FittedList)
    #  or to the experimental probability of connection
    GaussianFit, ConnectionProbability = connection_probability()
    FittedListSyn, FittedListGaps = FittedList()

    for popsIn in popInput:
        for popsOut in popTarget:
            popConnection(netParams, popInput=popsIn, popTarget=popsOut, sec=sec, loc=loc, weight=weight,
                          delay=delay, synMech=synMech, probability=1.0, gap=True)

###############################################################################
## Printing network definitions for debugging
###############################################################################
def debug(netParams):
    replines=20
    print("-"*replines)
    print("cellParams content")
    [print(i, netParams.cellParams[i]['secs'].keys()) for i in netParams.cellParams.keys()]
    print("-"*replines)
    print("popParams content")
    [print(i, netParams.popParams[i]) for i in netParams.popParams.keys()]
    print("-"*replines)
    print("connParams content")
    [print(i, netParams.connParams[i]) for i in netParams.connParams.keys()]
    print("-"*replines)
    print("stimSourceParams content")
    [print(i, netParams.stimSourceParams[i]) for i in netParams.stimSourceParams.keys()]
    print("-"*replines)
    print("stimTargetParams content")
    [print(i, netParams.stimTargetParams[i]) for i in netParams.stimTargetParams.keys()]

###############################################################################
## netParams DEFINITION
###############################################################################
def initialize_netParams(netParams, model='Peri'):
    # Define network geometry and general parameters
    network_geometry(netParams)
    connectivity_params(netParams)
    # Define network populations
    L1_pop(netParams, model=model)
    L23_pop(netParams, model=model)
    # Define synaptic mechanisms
    initialize_synaptic_mechs(netParams)
    # Define LGN inputs
    LGN_pop(netParams, label='LGN', layer='LGN', numCells=800, spkTimes=[[0+i,10+i,20+i,30+i]
                                                                              for i in range(800)])
    # Connect LGN to pops
    LGN_inputs(netParams)
    # Define background input populations
    bkg_pop(netParams, popInput='Bkg', interval=1, number=int(1e6), start=0, noise=0.5)
    # Connect background inputs to V1 populations
    bkg_input(netParams)
    # Intra- and inter-layer connections
    inh_connections(netParams)
    exc_connections(netParams)
    gap_connections(netParams)

initialize_netParams(netParams, model='Peri') # 'AllActive' models are failing. Probably some error when trying to load the cell geometry
#debug(netParams)



