# SDinSlice
## Overview
A tissue-scale model of spreading depolarization (SD) in brain slices.
We used the NEURON simulator's reaction-diffusion framework to implement embed thousands of neurons 
(based on the the model from Wei et al. 2014)
in the extracellular space of a brain slice, which is itself embedded in an bath solution.
We initiate SD in the slice by elevating extracellular K+ in a spherical region at the center of the slice.
Effects of hypoxia and propionate on the slice were modeled by appropriate changes to the volume fraction 
and tortuosity of the extracellular space and oxygen/chloride concentrations.
Users need to install [NEURON](https://neuron.yale.edu/neuron/), and we recommend using 
[MPI](https://www.open-mpi.org/) to parallelize simulations.

## Code
**SpatialModel.py** -- Simulation of SD with user specification of slice and cell properties via a json configuration file.

**genCfgs.py** -- Generates json configuration files that specifies slice dimensions, cell density, neuronal volume fraction,
neuronal surface area to volume ratio, slice oxygenation, etc.

**analyzeNeuromorpho.py** -- Computes average neuronal surface to volume ratios for various neuronal cell types from 
different brain regions in rats, mice, and humans using data from [NeuroMorpho](http://neuromorpho.org/).

**SpatialModelDynAlpha.py** -- Simulation of SD with dynamic changes in volume fraction of the extracellular 
space in perfused slice.

**analysis.py** -- Functions for analyzing output from SD simulations.

**figures.py** -- Functions for plotting output from SD simulations.

## Basic Usage 
### SD in small, perfused slice for 2s
The following runs a simulation of SD in a small (500 um x 500 um x 200 um), perfused 
slice for 2 seconds using a pre-made configuration file.  MPI is used for parallelizing 
and is highly recommended, especially for larger simulations.
```
mpiexec -n 6 nrniv -python -mpi SpatialModel.py cfgs/small_sim.json
```
This simulation can take over 30 minutes to run.  If MPI is not installed, it may 
be run with:
```
nrniv -python SpatialModel.py cfgs/small_sim.json
```
but this will take much longer.

To plot the output from that simulation:
```
python3 basicPlots.py Data/small_sim/
```

### SD in small, hypoxic slice for 2s
The following runs a simulation of SD in a small, hypoxic slice for 2 seconds using 
another pre-made configuration file.
```
mpiexec -n 6 nrniv -python -mpi SpatialModel.py cfgs/small_hypoxic_sim.json
```
Similarly, to plot the out for that simulation:
```
python3 basicPlots.py Data/small_hypoxic_sim/
```

### Comparing results from those two simulations
We have included a script for comparing the output of the previous two 
simulations shown above.  After both simulations, the following will create
a figure (**small_sim_comparison.png**) comparing the radial trajectories of the K+ waves and reduced raster plots
(plotting only the first spike rather than all for visualization purposes).
```
python3 compareSims.py
```

### SD in larger, hypoxic slice for 10 s (recommend running on HPC)
The following uses **genCfgs.py** to create a configuration file for 
simulating SD in a larger (1 mm x 1 mm x 400 um), hypoxic slice for 10 s, 
then runs the simulation with MPI.  Because of the size and duration of the 
simulation, we recommend only running this on an HPC.  The number of threads (nthreads) 
and number of processes passed to MPI (-n) can be changed depending on the available 
resources.  
```
python3 genCfgs.py --tstop=10000 --ox=anoxic --k0=70 --r0=100 --pas=-70.0 --uniformRec=True \
--nthreads=40 --nrec=40 --dir=Data/hypox_1mmmx1mmx400um_10s/ --sa2v=3.0 --O2consume=True \
cfgs/hypox_1mmx1mmx400um_10s.json
mpiexec -n 40 nrniv -python -mpi SpatialModel.py cfgs/hypox_1mmx1mmx400um_10s.json
```

## References
Wei, Yina, Ghanim Ullah, and Steven J. Schiff. "Unification of neuronal spikes, seizures, and spreading depression." Journal of Neuroscience 34, no. 35 (2014): 11733-11743.
https://doi.org/10.1523/JNEUROSCI.0516-14.2014



