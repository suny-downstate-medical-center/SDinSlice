from netpyne import specs
import numpy as np
#------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
#------------------------------------------------------------------------------

# Run parameters
cfg = specs.SimConfig()       # object of class cfg to store simulation configuration
cfg.duration = 2e3        # Duration of the simulation, in ms
cfg.hParams['v_init'] = -70.0   # set v_init to -65 mV
cfg.hParams['celsius'] = 37.0
cfg.dt = 0.1 #0.025              # Internal integration timestep to use
cfg.verbose = False            # Show detailed messages 
cfg.recordStep = 1             # Step size in ms to save data (eg. V traces, LFP, etc)
cfg.filename = 'normox_network_noelk_O2consume/'   # Set file output name
cfg.Kceil = 15.0

 # Network dimensions
cfg.sizeX = 250.0 #250.0 #1000
cfg.sizeY = 900.0 #250.0 #1000
cfg.sizeZ = 250.0 #200.0
cfg.density = 90000.0
cfg.Vtissue = cfg.sizeX * cfg.sizeY * cfg.sizeZ

## densities and E/I proporions based loosely on M1 paper 
cfg.L23density = 110e3
cfg.L4density = 111e3 
cfg.L5density = 100e3 
cfg.L23_eprop = 1766 / (1766 + 315)
cfg.L23_iprop = 1 - cfg.L23_eprop
cfg.L4_eprop = 1641 / (1641 + 89 + 182)
cfg.L4_iprop = 1 - cfg.L4_eprop 
cfg.L5_eprop = 1/2 #2310 / (2310 + 613)
cfg.L5_iprop = 1 - cfg.L5_eprop

cfg.N_L23_E = int(cfg.L23density * (cfg.sizeX * (1/3) * cfg.sizeY * cfg.sizeZ * 1e-9) * cfg.L23_eprop) 
cfg.N_L23_I = int(cfg.L23density * (cfg.sizeX * (1/3) * cfg.sizeY * cfg.sizeZ * 1e-9) * cfg.L23_iprop) 
cfg.N_L4_E  = int(cfg.L4density  * (cfg.sizeX * (1/3) * cfg.sizeY * cfg.sizeZ * 1e-9) * cfg.L4_eprop) 
cfg.N_L4_I  = int(cfg.L4density  * (cfg.sizeX * (1/3) * cfg.sizeY * cfg.sizeZ * 1e-9) * cfg.L4_iprop) 
cfg.N_L5_E  = int(cfg.L5density  * (cfg.sizeX * (1/3) * cfg.sizeY * cfg.sizeZ * 1e-9) * cfg.L5_eprop) 
cfg.N_L5_I  = int(cfg.L5density  * (cfg.sizeX * (1/3) * cfg.sizeY * cfg.sizeZ * 1e-9) * cfg.L5_iprop) 

# slice conditions 
cfg.ox = 'normoxic' #'perfused'
if cfg.ox == 'perfused':
    cfg.o2_bath = 0.1
    cfg.alpha_ecs = 0.2 
    cfg.tort_ecs = 1.6
if cfg.ox == 'normoxic':
    cfg.o2_bath = 0.04
    cfg.alpha_ecs = 0.2 
    cfg.tort_ecs = 1.6
elif cfg.ox == 'hypoxic':
    cfg.o2_bath = 0.01
    cfg.alpha_ecs = 0.07 
    cfg.tort_ecs = 1.8

cfg.O2consume = True 
cfg.gliaFactor = 1.0
# cfg.prep = 'invitro' #'invivo'

# neuron params 
cfg.betaNrn = 0.24
cfg.Ncell = int(cfg.density*(cfg.sizeX*cfg.sizeY*cfg.sizeZ*1e-9)) # default 90k / mm^3
if cfg.density == 90000.0:
    cfg.rs = ((cfg.betaNrn * cfg.Vtissue) / (2 * np.pi * cfg.Ncell)) ** (1/3)
else:
    cfg.rs = 7.52

cfg.epas = -70 # False
cfg.gpas = 0.0001
cfg.sa2v = 3.0 # False
if cfg.sa2v:
    cfg.somaR = (cfg.sa2v * cfg.rs**3 / 2.0) ** (1/2)
else:
    cfg.somaR = cfg.rs
cfg.cyt_fraction = cfg.rs**3 / cfg.somaR**3

# sd init params 
cfg.k0 = 3.5
cfg.r0 = 100.0

cfg.nRec = 240

# Recording/plotting parameters
# cfg.recordCells = [('E', [0,50,100,150,200,250,300,350,400,450])]
# cfg.recordTraces = {'V_soma':{'sec': 'soma','loc': 0.5,'var': 'v'}}#,
                        #   'ik_soma': {'sec': 'soma', 'loc': 0.5, 'var': 'ik'},
                        #   'cai_soma': {'sec': 'soma', 'loc':0.5, 'var': 'cai'},
                        #   'cao_soma': {'sec': 'soma', 'loc':0.5, 'var': 'cao'}}
# cfg.analysis['plotTraces'] = {'include': [('E', i*25) for i in range(20)], 'oneFigPer': 'cell', 'figSize': (10,4), 'saveFig': True, 'showFig': False}
# cfg.analysis['plotRaster'] = {'orderBy': 'y', 'orderInverse': True, 'saveFig': True, 'popNumCells' : [50 for i in range(6)]}         # Plot a raster

# cfg.recordLFP = [[-15, y, 1.0*cfg.sizeZ] for y in range(int(cfg.sizeY/3), int(cfg.sizeY), int(cfg.sizeY/3))]

# cfg.anacfg.sizeYsis['plotTraces']={'include': random.sample(range(cfg.Ncell),100), 'saveFig' : True}
# cfg.anacfg.sizeYsis['plotTraces']={'include': random.sample(range(cfg.Ncell),20), 'saveFig' : True, 'showFig' : False}
# cfg.anacfg.sizeYsis['plotRaster'] = {'orderBy': 'y', 'orderInverse': True, 'saveFig': True, 'figSize': (9,3), 'showFig' : False}      # Plot a raster
# cfg.anacfg.sizeYsis['plot2Dnet'] = {'saveFig': True, 'showFig' : False}   
#cfg.anacfg.sizeYsis['plotLFP'] = {'includeAxon': False, 'figSize': (6,10), 'NFFT': 256, 'noverlap': 48, 'nperseg': 64, 'saveFig': True} 
# cfg.anacfg.sizeYsis['plotRxDConcentration'] = {'speciesLabel': 'k', 'regionLabel': 'ecs', 'showFig' : False}
 
