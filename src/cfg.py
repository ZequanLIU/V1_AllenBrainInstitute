from netpyne.batchtools import specs
### config ###

cfg = specs.SimConfig()

cfg.duration = 1000
cfg.dt = 0.05
cfg.hparams = {'v_init': -65.0}
cfg.verbose = False
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
}
cfg.lgn={#LGN search space
    "LGN->VisL1Htr3a": 0.36e-3,
    "LGN->VisL2_3Htr3a": 0.36e-3,
    "LGN->VisL2_3Sst": 0.02e-3,
    "LGN->VisL2_3Pvalb": 0.02e-3,
    "LGN->VisL2_3E": 0.02e-3,
}
cfg.bkg = {#Background activity search space
    "Bkg->VisL1Htr3a": 0.36e-3,
    "Bkg->VisL2_3Htr3a": 0.36e-3,
    "Bkg->VisL2_3Sst": 0.02e-3,
    "Bkg->VisL2_3Pvalb": 0.02e-3,
    "Bkg->VisL2_3E": 0.02e-3,
}