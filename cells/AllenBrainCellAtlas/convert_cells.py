from utils import *
import os
import glob
import pandas as pd

cwd = os.getcwd()
cellModel = 'Peri' # Peri, AllActive
configfiles = glob.glob(cwd+'/CellModels%s/**/manifest.json' % cellModel, recursive=True)

# Load the models sued by Allen Brain
# Read the CSV file
specimen_model_id = pd.read_csv("bio_models_prop.csv", sep='\s+')
excitatory_cells = specimen_model_id.loc[specimen_model_id['ei'] == 'excitatory'][['model_id', 'cre_line']]

for i in range(len(configfiles)):
    path = os.path.normpath(configfiles[i])
    splitPath = path.split(os.sep)
    region = splitPath[-3]
    cellType, cell_id = splitPath[-2].split('_')
    workingDir = 'CellModels%s/%s/%s_%s/' % (cellModel,region,cellType,cell_id)
    if cellType in excitatory_cells['cre_line'].values and cell_id in str(excitatory_cells['model_id'].values):
        cellType = 'E'
    cell, cellRule = convert2NetPyNE(workingDir=workingDir, cellModel=cellModel, cellType=region+cellType, region=region,
                                     id=cell_id, SaveCell=True)
