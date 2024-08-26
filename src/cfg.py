from netpyne.batchtools import specs

### config ###

cfg = specs.SimConfig()

# Simulation duration and time step
cfg.duration = 1000
cfg.dt = 0.05
cfg.hparams = {'v_init': -65.0}
cfg.verbose = False

# Data recording settings
#cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}
cfg.recordTraces = {}
cfg.recordStim = False
cfg.recordStep = 0.1            # Step size in ms to save data (eg. V traces, LFP, etc)
cfg.filename = '00'         # Set file output name
cfg.savePickle = False        # Save params, network and sim output to pickle file
cfg.saveDat = False
cfg.saveJson = True
cfg.createNEURONObj = True
cfg.createPyStruct = True
cfg.printRunTime = 0.1
cfg.simLabel = 'V1'
cfg.saveFolder = '../data'
cfg.sizeY = 910
cfg.sizeX = 400.0
cfg.sizeZ = 400.0
cfg.scale = 1/400. # Keep it between 1/200. and 1/400. so it can run in a laptop

cfg.checkErrors = False

allpops = ['VisL1Htr3a','VisL2_3Htr3a','VisL2_3Pvalb','VisL2_3Sst','VisL2_3E']
cfg.analysis['plotRaster'] = {'include': allpops, 'saveFig': True, 'legend': allpops} # Raster ok
cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}
cfg.analysis['plotTraces'] = {'include': [(i,0) for i in allpops], 'saveFig': True, 'legend': None} # Record one cell of each population

cfg.cache_efficient = True # better with MPI?
""" remove all of the unnecessary data """
cfg.saveCellSecs = True
cfg.saveCellConns = True

# Here we put the parameters to optimize or to run in a grid search
# Not used for now. In the future we can run grid search for several combinations of parameters
cfg.ampa={#AMPA search space
    "VisL2_3E->VisL1Htr3a": 0.36e-3,
    "VisL2_3E->VisL2_3Pvalb": 0.36e-3,
    "VisL2_3E->VisL2_3Sst": 0.02e-3,
    "VisL2_3E->VisL2_3Htr3a": 0.02e-3,
    "VisL2_3E->VisL2_3E": 0.02e-3,

    "VisL4E->VisL1Htr3a": 0.3e-3,  #  L4 
    "VisL4E->VisL2_3Pvalb": 0.3e-3,
    "VisL4E->VisL2_3Sst": 0.3e-3,
    "VisL4E->VisL2_3Htr3a": 0.3e-3,
    "VisL4E->VisL2_3E": 0.3e-3,
    "VisL4E->VisL4E": 0.3e-3,
    "VisL5E->VisL1Htr3a": 0.25e-3,  # L5 
    "VisL5E->VisL2_3Pvalb": 0.25e-3,
    "VisL5E->VisL2_3Sst": 0.25e-3,
    "VisL5E->VisL2_3Htr3a": 0.25e-3,
    "VisL5E->VisL2_3E": 0.25e-3,
    "VisL5E->VisL4E": 0.25e-3,
    "VisL5E->VisL5E": 0.25e-3,
    "VisL6E->VisL1Htr3a": 0.2e-3,   # L6 
    "VisL6E->VisL2_3Pvalb": 0.2e-3,
    "VisL6E->VisL2_3Sst": 0.2e-3,
    "VisL6E->VisL2_3Htr3a": 0.2e-3,
    "VisL6E->VisL2_3E": 0.2e-3,
    "VisL6E->VisL4E": 0.2e-3,
    "VisL6E->VisL5E": 0.2e-3,
    "VisL6E->VisL6E": 0.2e-3
}
cfg.gabaa = {#GABAA search space
    "VisL1Htr3a->VisL1Htr3a": 4.5e-3,
    "VisL1Htr3a->VisL2_3Pvalb": 0.72e-3,
    "VisL1Htr3a->VisL2_3Sst": 72e-3,
    "VisL1Htr3a->VisL2_3Htr3a": 72e-3,
    "VisL1Htr3a->VisL2_3E": 72e-3,
    "VisL2_3Pvalb->VisL1Htr3a": 72e-3,
    "VisL2_3Pvalb->VisL2_3Pvalb": 72e-3,
    "VisL2_3Pvalb->VisL2_3Sst": 72e-3,
    "VisL2_3Pvalb->VisL2_3Htr3a": 72e-3,
    "VisL2_3Pvalb->VisL2_3E": 72e-3,
    "VisL2_3Sst->VisL1Htr3a": 72e-3,
    "VisL2_3Sst->VisL2_3Pvalb": 72e-3,
    "VisL2_3Sst->VisL2_3Sst": 72e-3,
    "VisL2_3Sst->VisL2_3Htr3a": 72e-3,
    "VisL2_3Sst->VisL2_3E": 72e-3,
    "VisL2_3Htr3a->VisL1Htr3a": 72e-3,
    "VisL2_3Htr3a->VisL2_3Pvalb": 72e-3,
    "VisL2_3Htr3a->VisL2_3Sst": 72e-3,
    "VisL2_3Htr3a->VisL2_3Htr3a": 72e-3,
    "VisL2_3Htr3a->VisL2_3E": 72e-3,

    # L4 
    "VisL4Pvalb->VisL1Htr3a": 72e-3,
    "VisL4Pvalb->VisL2_3Pvalb": 72e-3,
    "VisL4Pvalb->VisL2_3Sst": 72e-3,
    "VisL4Pvalb->VisL2_3Htr3a": 72e-3,
    "VisL4Pvalb->VisL2_3E": 72e-3,
    "VisL4Pvalb->VisL4Pvalb": 72e-3,

    "VisL4Sst->VisL1Htr3a": 72e-3,
    "VisL4Sst->VisL2_3Pvalb": 72e-3,
    "VisL4Sst->VisL2_3Sst": 72e-3,
    "VisL4Sst->VisL2_3Htr3a": 72e-3,
    "VisL4Sst->VisL2_3E": 72e-3,
    "VisL4Sst->VisL4Sst": 72e-3,

    "VisL4Htr3a->VisL1Htr3a": 72e-3,
    "VisL4Htr3a->VisL2_3Pvalb": 72e-3,
    "VisL4Htr3a->VisL2_3Sst": 72e-3,
    "VisL4Htr3a->VisL2_3Htr3a": 72e-3,
    "VisL4Htr3a->VisL2_3E": 72e-3,
    "VisL4Htr3a->VisL4Htr3a": 72e-3,

    # L5
    "VisL5Pvalb->VisL1Htr3a": 72e-3,
    "VisL5Pvalb->VisL2_3Pvalb": 72e-3,
    "VisL5Pvalb->VisL2_3Sst": 72e-3,
    "VisL5Pvalb->VisL2_3Htr3a": 72e-3,
    "VisL5Pvalb->VisL2_3E": 72e-3,
    "VisL5Pvalb->VisL5Pvalb": 72e-3,

    "VisL5Sst->VisL1Htr3a": 72e-3,
    "VisL5Sst->VisL2_3Pvalb": 72e-3,
    "VisL5Sst->VisL2_3Sst": 72e-3,
    "VisL5Sst->VisL2_3Htr3a": 72e-3,
    "VisL5Sst->VisL2_3E": 72e-3,
    "VisL5Sst->VisL5Sst": 72e-3,

    "VisL5Htr3a->VisL1Htr3a": 72e-3,
    "VisL5Htr3a->VisL2_3Pvalb": 72e-3,
    "VisL5Htr3a->VisL2_3Sst": 72e-3,
    "VisL5Htr3a->VisL2_3Htr3a": 72e-3,
    "VisL5Htr3a->VisL2_3E": 72e-3,
    "VisL5Htr3a->VisL5Htr3a": 72e-3,

    # L6 
    "VisL6Pvalb->VisL1Htr3a": 72e-3,
    "VisL6Pvalb->VisL2_3Pvalb": 72e-3,
    "VisL6Pvalb->VisL2_3Sst": 72e-3,
    "VisL6Pvalb->VisL2_3Htr3a": 72e-3,
    "VisL6Pvalb->VisL2_3E": 72e-3,
    "VisL6Pvalb->VisL6Pvalb": 72e-3,

    "VisL6Sst->VisL1Htr3a": 72e-3,
    "VisL6Sst->VisL2_3Pvalb": 72e-3,
    "VisL6Sst->VisL2_3Sst": 72e-3,
    "VisL6Sst->VisL2_3Htr3a": 72e-3,
    "VisL6Sst->VisL2_3E": 72e-3,
    "VisL6Sst->VisL6Sst": 72e-3,

    "VisL6Htr3a->VisL1Htr3a": 72e-3,
    "VisL6Htr3a->VisL2_3Pvalb": 72e-3,
    "VisL6Htr3a->VisL2_3Sst": 72e-3,
    "VisL6Htr3a->VisL2_3Htr3a": 72e-3,
    "VisL6Htr3a->VisL2_3E": 72e-3,
    "VisL6Htr3a->VisL6Htr3a": 72e-3,
}

cfg.gap={#Gap junction search space
    "VisL1Htr3a->VisL1Htr3a": 4.5e-3,
    "VisL1Htr3a->VisL2_3Pvalb": 0.72e-3,
    "VisL1Htr3a->VisL2_3Sst": 72e-3,
    "VisL1Htr3a->VisL2_3Htr3a": 72e-3,
    "VisL2_3Pvalb->VisL1Htr3a": 72e-3,
    "VisL2_3Pvalb->VisL2_3Pvalb": 72e-3,
    "VisL2_3Pvalb->VisL2_3Sst": 72e-3,
    "VisL2_3Pvalb->VisL2_3Htr3a": 72e-3,
    "VisL2_3Sst->VisL1Htr3a": 72e-3,
    "VisL2_3Sst->VisL2_3Pvalb": 72e-3,
    "VisL2_3Sst->VisL2_3Sst": 72e-3,
    "VisL2_3Sst->VisL2_3Htr3a": 72e-3,
    "VisL2_3Htr3a->VisL1Htr3a": 72e-3,
    "VisL2_3Htr3a->VisL2_3Pvalb": 72e-3,
    "VisL2_3Htr3a->VisL2_3Sst": 72e-3,
    "VisL2_3Htr3a->VisL2_3Htr3a": 72e-3,

    # L4 
    "VisL4Pvalb->VisL4Pvalb": 72e-3,
    "VisL4Sst->VisL4Sst": 72e-3,
    "VisL4Htr3a->VisL4Htr3a": 72e-3,

    # L5 
    "VisL5Pvalb->VisL5Pvalb": 72e-3,
    "VisL5Sst->VisL5Sst": 72e-3,
    "VisL5Htr3a->VisL5Htr3a": 72e-3,

    # L6
    "VisL6Pvalb->VisL6Pvalb": 72e-3,
    "VisL6Sst->VisL6Sst": 72e-3,
    "VisL6Htr3a->VisL6Htr3a": 72e-3,
}
cfg.lgn={#LGN search space
    "LGN->VisL1Htr3a": 0.36e-3,
    "LGN->VisL2_3Htr3a": 0.36e-3,
    "LGN->VisL2_3Sst": 0.02e-3,
    "LGN->VisL2_3Pvalb": 0.02e-3,
    "LGN->VisL2_3E": 0.02e-3,

    # LGN to L4 
    "LGN->VisL4Htr3a": 0.36e-3,
    "LGN->VisL4Sst": 0.02e-3,
    "LGN->VisL4Pvalb": 0.02e-3,
    "LGN->VisL4E": 0.02e-3,

    # LGN to L5 
    "LGN->VisL5Htr3a": 0.36e-3,
    "LGN->VisL5Sst": 0.02e-3,
    "LGN->VisL5Pvalb": 0.02e-3,
    "LGN->VisL5E": 0.02e-3,

    # LGN to L6 
    "LGN->VisL6Htr3a": 0.36e-3,
    "LGN->VisL6Sst": 0.02e-3,
    "LGN->VisL6Pvalb": 0.02e-3,
    "LGN->VisL6E": 0.02e-3,
}
cfg.bkg = {#Background activity search space
    "Bkg->VisL1Htr3a": 0.36e-3,
    "Bkg->VisL2_3Htr3a": 0.36e-3,
    "Bkg->VisL2_3Sst": 0.02e-3,
    "Bkg->VisL2_3Pvalb": 0.02e-3,
    "Bkg->VisL2_3E": 0.02e-3,

    "Bkg->VisL4Htr3a": 0.36e-3,
    "Bkg->VisL4Sst": 0.02e-3,
    "Bkg->VisL4Pvalb": 0.02e-3,
    "Bkg->VisL4E": 0.02e-3,

    "Bkg->VisL5Htr3a": 0.36e-3,
    "Bkg->VisL5Sst": 0.02e-3,
    "Bkg->VisL5Pvalb": 0.02e-3,
    "Bkg->VisL5E": 0.02e-3,

    "Bkg->VisL6Htr3a": 0.36e-3,
    "Bkg->VisL6Sst": 0.02e-3,
    "Bkg->VisL6Pvalb": 0.02e-3,
    "Bkg->VisL6E": 0.02e-3,
}

# Path to LGN spike file
cfg.lgn_spike_file = '/spikes.h5'