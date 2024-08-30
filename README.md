# V1_AllenBrainInstitute
 V1 model based on the Allen Brain V1 model, and using data extracted from the synaptic physiology and connectivity databases \\

## Environment Setup

**For conda environment:**

1. `conda create -n <ENV> python=3.10`
2. `conda activate <ENV>`
3. Install dependencies: `conda install -c conda-forge mpi4py mpich`, `pip install neuron`
4. Install raytune package: `pip install -U "ray[default]"`
5. Install NetPyNE from developer branch: 
```
git clone https://github.com/Neurosim-lab/netpyne.git
cd netpyne
git checkout batch
pip install -e .
```
6. Install batchtools
```
git clone https://github.com/jchen6727/batchtk.git
cd batchtk
pip install -e .
``` 

## Preparing simulation

Move to the dir `src/` and type `nrnivmodl ../mod/` to compile the NEURON mechanisms. A folder with compiled C code 
should be created, depending on the architecture.

## Running the simulation

From the same folder, type `mpiexec -n <Node>  nrniv -mpi -python init.py`


## TO DO in the future

1. Add realistic LGN inputs
2. Tune connectivity following experimental data from the synaptic physiology database. Change parameters in cfg.py
3. Add synaptic dynamics (LTP, STD) adn synaptic strength from same database
4. Add input normalization to all neurons, to reproduce *synaptic democracy*
5. Tune unknown parameters for reproducing *in-vivo* activity.

## Bugs

1. The models with `'AllActive'` segments are not currently working.
2. Others?
