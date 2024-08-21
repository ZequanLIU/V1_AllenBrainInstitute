import re
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
import airavata_cerebrum.atlas.data.abc_mouse as abc_mouse
import airavata_cerebrum.atlas.operations.netops as netops
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.api.queries.cell_types_api import CellTypesApi
import pandas as pd
import numpy as np
from allensdk.api.queries.biophysical_api import BiophysicalApi
import urllib
import os
from allensdk.model.biophys_sim.config import Config
from utils_AllenSDK import create_utils
from netpyne import specs

def compare_strings(string1, string2):
    pattern = re.compile(string2)
    match = re.findall(pattern, string1)
    return match


def loadMERFISHDatabase():
    # Downloading the necessary data
    download_base = "../cache/abc_mouse"
    abc_cache = AbcProjectCache.from_s3_cache(download_base)
    abc_cache.load_latest_manifest()
    abc_cache.get_directory_metadata("MERFISH-C57BL6J-638850")
    abc_cache.get_directory_metadata("MERFISH-C57BL6J-638850-CCF")
    abc_cache.get_directory_metadata("WMB-taxonomy")
    # Load Cell metadata and Gene metadata
    ABC_ATLAS_BASE = download_base
    manifest, file_meta = abc_mouse.merfish_files_meta()
    cell_meta = abc_mouse.cell_metadata(file_meta, ABC_ATLAS_BASE)
    gene_meta = abc_mouse.gene_metadata(file_meta, ABC_ATLAS_BASE)
    # Extend the cell metadata to include taxonomy data, class, subclass and types of all the cells
    cluster_details, cluster_colors = abc_mouse.taxonomy_cluster(manifest, ABC_ATLAS_BASE)
    cell_meta_ext = cell_meta.join(cluster_details, on="cluster_alias")
    cell_meta_ext = cell_meta_ext.join(cluster_colors, on="cluster_alias")

    return ABC_ATLAS_BASE, cell_meta, gene_meta, cell_meta_ext


def obtainCellCounts(ABC_ATLAS_BASE, numCells=300000, regionName="VISp", subregionName=["VISp1", "VISp2/3", "VISp4",
                     "VISp5","VISp6a","VISp6b"], modelName="v1"):
    ## Obtain Region Specific Ratios
    # Obtain region specific ratios based on the parcellation information of the Allen CCF.
    region_ei_df = abc_mouse.region_cell_type_ratios(regionName, ABC_ATLAS_BASE)

    ## Region Fractions to Counts
    modelNet_data = netops.atlasdata2regionfractions(region_ei_df, modelName)
    modelNet_data = netops.fractions2ncells(modelNet_data, numCells)
    modelNet_data = netops.subset_network(modelNet_data, subregionName)

    with open(modelName + "_stats.json", "w") as ofx:
        jsx = modelNet_data.model_dump_json(indent=4)
        ofx.write(jsx)

    return modelNet_data


def cellsModelSpecification(modelNet_data, regionList=["VISp1", "VISp2/3", "VISp4", "VISp5", "VISp6a", "VISp6b"],
                            structure_area_abbrev='VISp', structure_layer_name=["1", "2/3", "4", "5", "6a", "6b"]):
    '''
    Instantiate the CellTypesCache instance. The manifest_file argument tells it where to store the manifest,
    which is a JSON file that tracks file paths. If you supply a relative path (like this), it will go into your
    current working directory.
    TODO: There are too many models, how do we choose between them?
    '''
    ctc = CellTypesCache(manifest_file='cell_types/manifest.json')

    # mouse cells
    cells = ctc.get_cells(species=[CellTypesApi.MOUSE], require_reconstruction=True)
    cells = pd.DataFrame(cells)

    cells = cells.reset_index()  # make sure indexes pair with number of rows

    AreaList = []
    Model_specification = {}
    for layer in structure_layer_name:
        regionName = structure_area_abbrev + layer
        Model_specification[regionName] = {}
        for index, row in cells.iterrows():
            for neuronType in modelNet_data.locations[regionName].neurons:
                if (row['structure_layer_name'] == layer and row['structure_area_abbrev'] == structure_area_abbrev
                        and neuronType in compare_strings(row['name'], neuronType)):
                    Model_specification[regionName][str(neuronType)] = dict(
                        modelNet_data.locations[regionName].neurons[str(neuronType)])
                    AreaList.append([regionName, neuronType, row['id'], row['normalized_depth']])

    AreaList = np.stack(AreaList)
    for region in regionList:
        for neuronType in Model_specification[region].keys():
            maskRegion = AreaList[:, 0] == region
            maskNeuron = [neuronType in compare_strings(i, neuronType) for i in AreaList[:, 1]]
            TotalMask = maskRegion * maskNeuron
            Model_specification[region][neuronType]['id'] = np.array(AreaList[TotalMask, 2].flatten(),
                                                                     dtype=int).tolist()
            Model_specification[region][neuronType]['NormDepth'] = np.array(AreaList[TotalMask, 3].flatten(),
                                                                            dtype=float).tolist()

    return Model_specification


def downloadCellModels(Model_specification):
    # http://celltypes.brain-map.org
    bp = BiophysicalApi()
    bp.cache_stimulus = False  # change to False to not download the large stimulus NWB file

    CellModelsFolderPeri = 'CellModelsPeri/'
    CellModelsFolderAllActive = 'CellModelsAllActive/'

    NotLoadedPeri = []
    NotLoadedAllActive = []
    for i in Model_specification.keys():
        for j in Model_specification[i].keys():
            specimen_model_id = Model_specification[i][j]['id']  # get this from the website as above
            for id in specimen_model_id:
                workingDirPeri = CellModelsFolderPeri + '{:s}/{:s}_{:d}'.format(i.replace('/', '_'), j, id)
                workingDirAllActive = CellModelsFolderAllActive + '{:s}/{:s}_{:d}'.format(i.replace('/', '_'), j, id)

                pagePerisomatic = urllib.request.urlopen(
                    "http://api.brain-map.org/api/v2/data/query.xml?criteria=model::NeuronalModel,rma::critera,[specimen_id$eq{:d}],neuronal_model_template[name$il'*perisomatic*']".format(
                        id)).read()
                temp = re.findall(r"<id>\d+</id>", str(pagePerisomatic))
                temp = re.findall(r'\d+', str(temp))
                if len(temp) > 0:
                    resPerisomatic = int(temp[0])
                else:
                    resPerisomatic = []

                pageActive = urllib.request.urlopen(
                    "http://api.brain-map.org/api/v2/data/query.xml?criteria=model::NeuronalModel,rma::critera,[specimen_id$eq{:d}],neuronal_model_template[name$il'*active*']".format(
                        id)).read()
                temp = re.findall(r"<id>\d+</id>", str(pageActive))
                temp = re.findall(r"\d+", str(temp))
                if len(temp) > 0:
                    resActive = int(temp[0])
                else:
                    resActive = []
                # If model does not exist
                if not os.path.isdir(workingDirPeri):
                    try:
                        bp.cache_data(resPerisomatic, working_directory=workingDirPeri)
                    except:
                        NotLoadedPeri.append([workingDirPeri])
                        with open('Errors_Perisomatic.txt', 'a') as f:
                            f.write("No data available for {} \n".format(workingDirPeri))
                if not os.path.isdir(workingDirAllActive):
                    try:
                        bp.cache_data(resActive, working_directory=workingDirAllActive)
                    except:
                        NotLoadedAllActive.append([workingDirAllActive])
                        with open('Errors_AllActive.txt', 'a') as f:
                            f.write("No data available for {} \n".format(workingDirAllActive))


# ------------------------------------------------------------------------------------------------------------
# Set cell dynamic params into a cell rule (netParams.cellParams) from Json
# ------------------------------------------------------------------------------------------------------------
def setCellRuleDynamicParamsFromJson(cell_dynamic_params, cellRule, cellModel='Peri'):

    passive = cell_dynamic_params['passive'][0]
    conditions = cell_dynamic_params['conditions'][0]
    genome = cell_dynamic_params['genome']

    # Set passive properties
    if cellModel=='Peri':
        cm_dict = dict([(c['section'], c['cm']) for c in passive['cm']])
        for secName, sec in cellRule['secs'].items():
            sec['geom']['Ra'] = passive['ra']
            sec['geom']['cm'] = cm_dict[secName.split('_')[0]]
            sec['mechs'] = {'pas': {'e': passive["e_pas"]}}

        # Insert channels and set parameters
        for p in genome:
            sections = [s for s in cellRule['secs'] if s.split('_')[0] == p["section"]]
            for sec in sections:
                if p["mechanism"] != "":
                    cellRule['secs'][sec]['mechs'][p['mechanism']] = {p['name'].split('_')[0]: p['value']}
    elif cellModel=='AllActive':
        for secName, sec in cellRule['secs'].items():
            sec['geom']['Ra'] = passive['ra']

        # Insert channels and set parameters
        for p in genome:
            sections = [s for s in cellRule['secs'] if s.split('_')[0] == p["section"]]
            for sec in sections:
                if p["mechanism"] != "":
                    cellRule['secs'][sec]['mechs'][p['mechanism']] = {p['name'].split('_')[0]: p['value']}
                elif p["mechanism"]=="" and p['name']=='g_pas':
                    cellRule['secs'][sec]['mechs']['pas'] = {p['name'].split('_')[0]: p['value']}
                elif p["mechanism"]=="" and p['name']=='e_pas':
                    cellRule['secs'][sec]['mechs']['pas'] = {p['name'].split('_')[0]: p['value']}
                elif p["mechanism"]=="" and p['name']=='cm':
                    cellRule['secs'][sec]['geom']['cm'] = {p['name'].split('_')[0]: p['value']}


    # Set reversal potentials
    for erev in conditions['erev']:
        sections = [s for s in cellRule['secs'] if s.split('_')[0] == erev["section"]]
        for sec in sections:
            for eion in erev:
                if eion.startswith('e'):
                    if 'ions' not in cellRule['secs'][sec]:
                        #print(sec, eion)
                        cellRule['secs'][sec]['ions'] = {}
                    cellRule['secs'][sec]['ions'][eion[1:]] = {'e': erev[eion]}

    if 'v_init' in conditions:
        for sec in cellRule['secs'].values():
            sec['vinit'] = conditions['v_init']

    return cellRule


def convert2NetPyNE(workingDir, cellType, region, id, cellModel='Peri', SaveCell=True):
    '''Single-process main function for simulating sweeps in a biophysical experiment.

    Parameters
    ----------
    description : Config
        All information needed to run the experiment.
    sweeps : list
        list of experiment sweep numbers to simulate.  If None, simulate all sweeps.
    '''
    cwd = os.getcwd()
    #os.system('nrnivmodl ' + workingDir + 'modfiles/')
    os.chdir(workingDir)
    description = Config().load('manifest.json')
    # configure NEURON
    utils = create_utils(description)
    #h = utils.h

    # configure model
    manifest = description.manifest
    morphology_path = description.manifest.get_path('MORPHOLOGY')
    #utils.generate_morphology(morphology_path.encode('ascii', 'ignore'))
    #utils.load_cell_parameters()

    netParams = specs.NetParams()  # object of class NetParams to store the network parameters
    full_morpho_path = utils.manifest.path_info['BASEDIR']['spec'] + morphology_path[1:]
    if cellModel=='Peri':
        cellRuleAux = netParams.importCellParams(label=cellType+'_'+region+'_'+id,
                                          conds={'cellType': cellType, 'cellModel': cellModel},
                                          fileName=cwd+'/templates/BioAxonStub.hoc',
                                          cellName='BioAxonStub',
                                          cellArgs=[0, full_morpho_path])  # 0 in cellArgs is just placeholder
    elif cellModel=='AllActive':
        cellRuleAux = netParams.importCellParams(label=cellType+'_'+region+'_'+id,
                                          conds={'cellType': cellType, 'cellModel': cellModel},
                                          fileName=cwd+'/templates/Biophys1.hoc',
                                          cellName='Biophys1',
                                          cellArgs=[0, full_morpho_path])  # 0 in cellArgs is just placeholder

    netParams.renameCellParamsSec(cellType+'_'+region+'_'+id, 'soma_0', 'soma')  # rename imported section 'soma_0' to 'soma'
    cellRule = setCellRuleDynamicParamsFromJson(utils.description.data, cellRuleAux, cellModel=cellModel)
    os.chdir(cwd)

    cellRule['secLists']['alldend'] = [sec for sec in cellRule.secs if ('dend' in sec or 'apic' in sec)]  # basal+apical
    cellRule['secLists']['basal'] = [sec for sec in cellRule.secs if ('dend' in sec and 'apic' not in sec)]  # basal
    cellRule['secLists']['apic'] = [sec for sec in cellRule.secs if ('apic' in sec)]  # apical
    cellRule['secLists']['axon'] = [sec for sec in cellRule.secs if ('axon' in sec)]  # axons
    cellRule['secLists']['somatic'] = [sec for sec in cellRule.secs if ('soma' in sec)]  # soma
    if SaveCell:
        netParams.saveCellParamsRule(label=cellType+'_'+region+'_'+id, fileName=cwd+'/../NetPyNE_cells/%s/%s_%s_cellParams.pkl' % (region, cellType+'_'+id, cellModel))

    return utils, cellRule
