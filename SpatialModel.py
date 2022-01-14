"""Simulation of spreading depression"""
# from mpi4py import MPI
from numpy import random
from neuron import h, rxd
#h.nrnmpi_init()
from neuron.rxd import v
from neuron.rxd.rxdmath import exp, log, tanh
from neuron.units import sec, mM
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot, colors, colorbar
import math
from math import pi
from matplotlib import pyplot
import numpy
import os
import sys
import pickle
import json 
import itertools
from scipy.signal import find_peaks 

# when using multiple processes get the relevant id and number of hosts 
pc = h.ParallelContext()
pcid = pc.id()
nhost = pc.nhost()
# pcid = 0
# nhost = 20
pc.timeout(0)

# load sim configuration
with open(sys.argv[-1],'rb') as fileObj:
    args = json.load(fileObj)

# set the save directory 
outdir = os.path.abspath(args['dir'])
if pcid == 0 and not os.path.exists(outdir):
    try:
        os.makedirs(outdir)
        os.makedirs(outdir+'/gifFigs/')
    except:
        print("Unable to create the directory %r for the data and figures"
              % outdir)
        os._exit(1)

# save sim state functions 
def saveRxd():
    for sp in rxd.species._all_species:
        s = sp()
        numpy.save(os.path.join(outdir, s.name + '_concentrations_' + str(pcid) + '.npy'), s.nodes.concentration)

def runSS():
    svst = h.SaveState()
    svst.save()
    f = h.File(os.path.join(outdir,'save_test_' + str(pcid) + '.dat'))
    svst.fwrite(f)

rxd.nthread(args['nthreads'])  # set the number of rxd threads - original 4
rxd.options.enable.extracellular = True # enable extracellular rxd

h.load_file('stdrun.hoc')
h.celsius = 37
# h.dt = 0.025
h.dt = 0.1 # original above, going for longer runs

e_charge = 1.60217662e-19
scale = 1e-14/e_charge
gnabar = (30/1000) * scale     # molecules/um2 ms mV 
gnabar_l = (0.0247/1000) * scale
gkbar = (25/1000) * scale
gkbar_l = (0.05/1000) * scale
gclbar_l = (0.1/1000) * scale
ukcc2 = 0.3 * mM/sec 
unkcc1 = 0.1 * mM/sec 
alpha = 5.3
epsilon_k_max = 0.25/sec
epsilon_o2 = 0.17/sec
vtau = 1/250.0
g_gliamax = 5 * mM/sec
beta0 = 7.0
avo = 6.0221409*(10**23)

# modified parameter
p_max = args['p_max'] #8 #0.8 * mM/sec

nao_initial = 144.0
nai_initial = 18.0 #22.5 #
gnai_initial = 18.0
gki_initial = 80.0

#MODIFIED -- original 142.5, 7.8 ---
ko_initial = 3.5
ki_initial = 140.0 #133.0 #
clo_initial = 130
cli_initial = 6.0 #3.8 #

if args['ox'] == 'primed':
    clo_initial = clo_initial / 2
    cli_initial = cli_initial / 2

if args['varCl']:
    factor = args['varCl'] / clo_initial
    clo_initial = clo_initial * factor 
    cli_initial = cli_initial * factor

if args['ox'] == 'anoxic':
    oa_bath = 0.01
elif args['ox'] == 'orig' or args['ox'] == 'primed' or args['ox'] == 'mannitol':
    oa_bath = 0.1
else:
    oa_bath = 0.04 #args.bathO2/alpha

if args['varO2']:
    oa_bath = args['varO2']

# oa_bath = 0.1 # original value from adam's code
v_initial = -70 #-74.7 #-70

#sodium activation 'm'
alpha_m = (0.32 * (v + 54.0))/(1.0 - exp(-(v + 54)/4.0))
beta_m = (0.28 * (v + 27.0))/(exp((v + 27.0)/5.0) - 1.0)

alpha_m0 =(0.32 * (v_initial + 54.0))/(1.0 - math.exp(-(v_initial + 54)/4.0))
beta_m0 = (0.28 * (v_initial + 27.0))/(math.exp((v_initial + 27.0)/5.0) - 1.0)
m_initial = alpha_m0/(beta_m0 + 1.0)

#sodium inactivation 'h'
alpha_h = 0.128 * exp(-(v + 50.0)/18.0)
beta_h = 4.0/(1.0 + exp(-(v + 27.0)/5.0))
alpha_h0 = 0.128 * math.exp(-(v_initial + 50.0)/18.0)
beta_h0 = 4.0/(1.0 + math.exp(-(v_initial + 27.0)/5.0))
h_initial = alpha_h0/(beta_h0 + 1.0)

#potassium activation 'n'
alpha_n = (0.032 * (v + 52.0))/(1.0 - exp(-(v + 52.0)/5.0))
beta_n = 0.5 * exp(-(v + 57.0)/40.0)
alpha_n0 = (0.032 * (v_initial + 52.0))/(1.0 - math.exp(-(v_initial + 52.0)/5.0))
beta_n0 = 0.5 * math.exp(-(v_initial + 57.0)/40.0)
n_initial = alpha_n0/(beta_n0 + 1.0)

numpy.random.seed(args['randSeed']+pcid)    # use a difference seed for each process

# simulation parameters
if args['size'] == 'small':
    Lx, Ly = 500.0, 500.0
elif args['size'] == 'mini':
    Lx, Ly = 200.0, 200.0
elif args['size'] == 'big':
    Lx, Ly = 2000.0, 2000.0
elif args['size'] == 'bigger':
    Lx, Ly = 3000.0, 3000.0
else:
    Lx, Ly = 1000.0, 1000.0#170 # value fr fovea. 1000 #250 # 750, 750, 375 # size of the extracellular space mu m^3 - original 500, 500, 250
    
Lz = args['Lz']

if args['size'] == 'cube':
    Lx, Ly, Lz = 500.0, 500.0, 500.0

Kceil = 15.0                     # threshold used to determine wave speed
Vtissue = Lx*Ly*Lz 
# cell numbers 
Ncell = int(args['density']*(Lx*Ly*Lz*1e-9)) # default 90k / mm^3
Nrec = args['nrec']
if args['density'] == 90000:
    rs = ((args['alphaNrn']*Vtissue)/(2*numpy.pi*Ncell)) ** (1/3) # defaults to 7.52: appropriate radius for neuronal volume fraction of 0.24 given cylinders whose height is the diameter
else:
    rs = 7.52

# compute appropriate radius for given surface area to volume ratio
if args['sa2v']:
    somaR = (args['sa2v'] * rs**3 / 2.0) ** (1/2)
else:
    somaR = rs #10      # larger than in the paper - original 15. 

# ECS params
alpha0, alpha1, alpha2, alpha3, alpha4 = 0.07, 0.2, 0.12, 0.3, 0.32  # anoxic and normoxic volume fractions 
tort0, tort1, tort2, tort3 = 1.8, 1.6, 2.0, 1.4    # anoxic and normoxic tortuosities 

r0 = args['r0']  # radius for initial elevated K+

# allows for cmd line switching normox vs anox
if args['ox'] == 'normoxic' or args['ox'] ==  'orig':
    alpha0 = alpha1 # calls are made w/ initial anoxic vars
    tort0 = tort1
elif args['ox'] == 'brainstem':
    alpha0 = alpha3 
    tort0 = tort1
elif args['ox'] == 'primed':
    alpah0 = alpha2 
    tort0 = tort1
elif args['ox'] == 'mannitol':
    alpha0 = alpha4
    tort0 = tort3

if args['alphaECS']:
    alpha0 = args['alphaECS']

if args['lambdaECS']:
    tort0 = args['lambdaECS']

soma_list = h.SectionList()
dend_list = h.SectionList()

class Neuron:
    """ A neuron with soma and dendrite with; fast and persistent sodium
    currents, potassium currents, passive leak and potassium leak and an
    accumulation mechanism. """
    def __init__(self, x, y, z, rec=False):
        self.x = x
        self.y = y
        self.z = z

        self.soma = h.Section(name='soma', cell=self)
        # add 3D points to locate the neuron in the ECS  
        self.soma.pt3dadd(x, y, z + somaR, 2.0*somaR)
        self.soma.pt3dadd(x, y, z - somaR, 2.0*somaR)
        if args['pas']:
            self.soma.insert('pas')
            self.soma(0.5).pas.e = args['pas']
            self.soma(0.5).pas.g = args['gpas']

        soma_list.append(self.soma)
    
        if rec: # record membrane potential (shown in figure 1C)
            self.somaV = h.Vector()
            self.somaV.record(self.soma(0.5)._ref_v)

class NeuronD:
    """ A neuron with soma and dendrite with; fast and persistent sodium
    currents, potassium currents, passive leak and potassium leak and an
    accumulation mechanism. """
    def __init__(self, x, y, z, height, rec=False):
        self.x = x
        self.y = y
        self.z = z

        self.soma = h.Section(name='soma', cell=self)
        # add 3D points to locate the neuron in the ECS  
        self.soma.pt3dadd(x, y, z + somaR, 2.0*somaR)
        self.soma.pt3dadd(x, y, z - somaR, 2.0*somaR)
        self.dend = h.Section(name='dend', cell=self)
        self.dend.pt3dadd(x, y, z + somaR, 4)
        self.dend.pt3dadd(x, y, z + somaR + height, 4)
        self.dend.diam = 4
        self.dend.L = height
        self.dend.connect(self.soma(1))

        soma_list.append(self.soma)
        dend_list.append(self.dend)
    
        if rec: # record membrane potential (shown in figure 1C)
            self.somaV = h.Vector()
            self.somaV.record(self.soma(0.5)._ref_v)

# Randomly distribute 1000 neurons which we record the membrane potential
# every 50ms
if args['dendL'] > 0:
    if pcid == 0:
        print("Neuron with dend")
    rec_neurons = [NeuronD(
       (numpy.random.random()*2.0 - 1.0) * (Lx/2.0 - somaR), 
       (numpy.random.random()*2.0 - 1.0) * (Ly/2.0 - somaR), 
       (numpy.random.random()*2.0 - 1.0) * (Lz/2.0 - somaR - args['dendL']), args['dendL'], 50)
       for i in range(0, int(Nrec/nhost))]
    all_neurons = [NeuronD(
        (numpy.random.random()*2.0 - 1.0) * (Lx/2.0 - somaR),
        (numpy.random.random()*2.0 - 1.0) * (Ly/2.0 - somaR),
        (numpy.random.random()*2.0 - 1.0) * (Lz/2.0 - somaR - args['dendL']), args['dendL'])
        for i in range(int(Nrec/nhost), int(Ncell/nhost))]    
else:
    if pcid == 0:
        print('point neuron')
    if args['uniformRec']:
        center_rec_neurons = [Neuron(
            (numpy.random.random()*2.0 - 1.0) * (Lx/2.0 - somaR), 
            (numpy.random.random()*2.0 - 1.0) * (Ly/2.0 - somaR), 
            (numpy.random.random()*2.0 - 1.0) * (Lz/2.0 - somaR), 50)
            for i in range(0, int(Nrec/nhost))]
        periph_rec_neurons = []
    else:
        center_rec_neurons = [Neuron(
            (numpy.random.random()*2.0 - 1.0) * (Lx/2.0 - somaR), 
            (numpy.random.random()*2.0 - 1.0) * (Ly/2.0 - somaR), 
            0, 50)
            for i in range(0, int(Nrec/2/nhost))]
        #(numpy.random.random()*2.0 - 1.0) * 0.5
        periph_rec_neurons = [Neuron(
            (numpy.random.random()*2.0 - 1.0) * (Lx/2.0 - somaR), 
            (numpy.random.random()*2.0 - 1.0) * (Ly/2.0 - somaR), 
            numpy.random.choice((-1,1)) * (Lz/2.0 - somaR), 50)
            for i in range(0, int(Nrec/2/nhost))]
    # Randomly distribute the remaining neurons
    all_neurons = [Neuron(
        (numpy.random.random()*2.0 - 1.0) * (Lx/2.0 - somaR),
        (numpy.random.random()*2.0 - 1.0) * (Ly/2.0 - somaR),
        (numpy.random.random()*2.0 - 1.0) * (Lz/2.0 - somaR))
        for i in range(int(Nrec/nhost), int(Ncell/nhost))]


# Where? -- define the extracellular space
if args['edemaCore'] or args['ischemEdemaCore']: # need args['ox'] == 'anoxic'
    def alphaecs(x, y, z) :
        if x**2 + y**2 + z**2 < r0**2:
            return alpha0 
        else:
            return min(alpha1, alpha0 + (alpha1-alpha0) *((x**2+y**2+z**2)**0.5-r0)/(Lx/2))

    def tortecs(x, y, z) :
        if x**2 + y**2 + z**2 < r0**2:
            return tort0 
        else:
            return max(tort1, tort0 - (tort0-tort1) *((x**2+y**2+z**2)**0.5-r0)/(Lx/2))

    ecs = rxd.Extracellular(-Lx/2.0, -Ly/2.0,
                        -Lz/2.0, Lx/2.0, Ly/2.0, Lz/2.0, dx=25,
                        volume_fraction=alphaecs, tortuosity=tortecs) # switched to ischemic
else:
    ecs = rxd.Extracellular(-Lx/2.0, -Ly/2.0,
                        -Lz/2.0, Lx/2.0, Ly/2.0, Lz/2.0, dx=25,
                        volume_fraction=alpha0, tortuosity=tort0) # switched to ischemic

## separate ecs for o2 
ecs_o2 = rxd.Extracellular(-Lx/2.0, -Ly/2.0,
                        -Lz/2.0, Lx/2.0, Ly/2.0, Lz/2.0, dx=25,
                        volume_fraction=1.0, tortuosity=1.0)

if args['sa2v']:
    cyt_frac = rs**3 / somaR**3
    cyt = rxd.Region(h.allsec(), name='cyt', nrn_region='i', geometry=rxd.FractionalVolume(cyt_frac, surface_fraction=1))
else:
    cyt = rxd.Region(h.allsec(), name='cyt', nrn_region='i')

mem = rxd.Region(h.allsec(), name='mem', geometry=rxd.membrane())

# What? -- define the species
def concentration(i, o):
    return lambda nd: i if isinstance(nd, rxd.node.Node1D) else o 

# if args['BC'] == 'invivo':
k_bc = ko_initial
na_bc = nao_initial
cl_bc = clo_initial
o2_bc = oa_bath
# else:
#     k_bc = None 
#     na_bc = None 
#     cl_bc = None 
#     o2_bc = None 

k = rxd.Species([cyt, mem, ecs], name='k', d=2.62, charge=1,
                initial=lambda nd: ki_initial if
                isinstance(nd, rxd.node.Node1D) else ( args['k0']
                if nd.x3d**2 + nd.y3d**2 + nd.z3d**2 <= r0**2 else ko_initial),
                ecs_boundary_conditions=k_bc)

na = rxd.Species([cyt, mem, ecs], name='na', d=1.78, charge=1,
                 initial=concentration(nai_initial, nao_initial),
                 ecs_boundary_conditions=na_bc)

cl = rxd.Species([cyt, mem, ecs], name='cl', d=2.1, charge=-1,
                 initial=concentration(cli_initial, clo_initial),
                 ecs_boundary_conditions=cl_bc)

# rescale mM/ms to molecules/um**2/ms 
volume = cyt.geometry.volumes1d(center_rec_neurons[0].soma)[0]
area = cyt.geometry.surface_areas1d(center_rec_neurons[0].soma)[0]

volume_scale = 1e-18 * avo * volume / area

#extracellular oxygen concentration
if args['ischemCore'] or args['ischemEdemaCore']:
    o2_extracellular = rxd.Species([ecs_o2], name='o2', d=3.3, initial = lambda nd: 0.01
                if nd.x3d**2 + nd.y3d**2 + nd.z3d**2 <= r0**2 else 0.04, ecs_boundary_conditions=0.04) # changed for separate ecs for o2 
else:
    o2_extracellular = rxd.Species([ecs_o2], name='o2', d=3.3, initial=oa_bath, ecs_boundary_conditions=o2_bc) # changed for separate ecs for o2 

o2ecs = o2_extracellular[ecs_o2]
o2switch = (1.0 + tanh(1e4*(o2ecs-5e-4)))/2.0

#volume ratio
vol_ratio = rxd.State([cyt, ecs], name='volume', initial=1.0)
vir = vol_ratio[cyt] # intracellular ratio of volume at time t to initial volume
vor = vol_ratio[ecs] # extracellular ratio of volume at time t to initial volume
                     # vor(t) == beta0*( 1.0 - vir(t) ) + 1.0)
                     # but two states are needed to support the regions.

# boundary conditions - as is like in vitro
def bc(nd):
    if (abs(nd.x3d - ecs._xlo) < ecs._dx[0] or
        abs(nd.x3d - ecs._xhi) < ecs._dx[0] or
        abs(nd.y3d - ecs._ylo) < ecs._dx[1] or
        abs(nd.y3d - ecs._yhi) < ecs._dx[1] or
        abs(nd.z3d - ecs._zlo) < ecs._dx[2] or
        abs(nd.z3d - ecs._zhi) < ecs._dx[2]):
        return 1.0
    return 0.0

# in vivo - only restrict in one z-direction, can diffuse out
def bcvivo(nd):
    if abs(nd.z3d - ecs._zlo) >= ecs._dx[2]:
        return 0.0
    return 1.0

# core conditions
def core(nd):
    if nd.x3d**2 + nd.y3d**2 + nd.z3d**2 <= r0**2:
        return 1.0
    return 0.0

def anticore(nd):
    if nd.x3d**2 + nd.y3d**2 + nd.z3d**2 <= r0**2:
        return 0.0
    return 1.0

iscore = rxd.Parameter([ecs_o2, mem], name='iscore', value = lambda nd: core(nd))

notcore = rxd.Parameter([ecs, ecs_o2, mem], name='notcore', value = lambda nd: anticore(nd))

dump = rxd.Parameter([cyt, ecs, ecs_o2], name='dump')

ecsbc = rxd.Parameter([ecs, ecs_o2], name='ecsbc', value = lambda nd: bc(nd))

ki, ko, nai, nao, cli, clo = k[cyt], k[ecs], na[cyt], na[ecs], cl[cyt], cl[ecs]
#STATES-------------------------------------------------------------------------

#gating variables (m, h, n)
mgate = rxd.State([cyt, mem], name='mgate', initial=m_initial) 
hgate = rxd.State([cyt, mem], name='hgate', initial=h_initial) 
ngate = rxd.State([cyt, mem], name='ngate', initial=n_initial) 

#ALL EQUATIONS------------------------------------------------------------------
gna = gnabar*mgate**3*hgate
gk = gkbar*ngate**4
fko = 1.0 / (1.0 + exp(16.0 - ko / vor))
nkcc1 = unkcc1*fko*(log((ki * cli / vir**2) / (ko * clo / vor**2)) + log((nai * cli / vir**2) / (nao * clo / vor**2)))
kcc2 = ukcc2 * log((ki * cli / vir**2) / (ko * clo / vor**2))

#Nerst equation - reversal potentials
ena = 26.64 * log(nao*vir/(nai*vor))
ek = 26.64 * log(ko*vir/(ki*vor))
ecl = 26.64 * log(cli*vor/(clo*vir))

p = o2switch * p_max / (1.0 + exp((args['pparam'] - (o2ecs/vor) * alpha)/3.0))
pump = args['nrnPumpFactor'] * (p / (1.0 + exp((25.0 - nai / vir)/3.0))) * (1.0 / (1.0 + exp(3.5 - ko / vor)))
gliapump = args['glialPumpFactor'] * (1.0/3.0) * (p / (1.0 + exp((25.0 - gnai_initial) / 3.0))) * (1.0 / (1.0 + exp(3.5 - ko/vor)))
g_glia = g_gliamax / (1.0 + exp(-(o2ecs*alpha/vor - 2.5)/0.2))
glia12 = g_glia / (1.0 + exp((18.0 - ko / vor)/2.5))

epsilon_k = (epsilon_k_max/(1.0 + exp(-((o2ecs/vor) * alpha - 2.5)/0.2))) * (1.0/(1.0 + exp((-20 + ((1.0+1.0/beta0 -vor)/vor) /2.0))))
           
#RATES--------------------------------------------------------------------------

#dm/dt
m_gate = rxd.Rate(mgate, (alpha_m * (1.0 - mgate)) - (beta_m * mgate))

#dh/dt
h_gate = rxd.Rate(hgate, (alpha_h * (1.0 - hgate)) - (beta_h * hgate))

#dn/dt
n_gate = rxd.Rate(ngate, (alpha_n * (1.0 - ngate)) - (beta_n * ngate))

#Diffusion
o2diff = rxd.Rate(o2ecs, ecsbc*(epsilon_o2 * (oa_bath - o2ecs/vor))) 
kdiff = rxd.Rate(ko, ecsbc*(epsilon_k * (ko_initial - ko/vor))) 
nadiff = rxd.Rate(nao, ecsbc*(epsilon_k * (nao_initial - nao/vor))) 
cldiff = rxd.Rate(clo, ecsbc*(epsilon_k * (clo_initial - clo/vor)))

# K+ infusion 
if args['infuse']:
    kinfuse = rxd.Rate(ko, iscore * (epsilon_k_max * (args['k0'] - ko)))

#change in volume
osm = (1.1029 - 0.1029*exp( ( (nao + ko + clo + 18.0)/vor - 
                             (nai + ki + cli + 132.0)/vir)/20.0))
scalei = (avo*1e-18)
scaleo = (avo*1e-18)

vol_dyn = rxd.MultiCompartmentReaction(vir, dump[ecs],
                                       -scalei*vtau*(osm-vir),
                                       mass_action=False,
                                       membrane=mem,
                                       scale_by_area=False,
                                       membrane_flux=False)

vol_dyn_ecs = rxd.MultiCompartmentReaction(dump[cyt], vor,
                                       -scaleo*vtau*(osm-vir),
                                       mass_action=False,
                                       membrane=mem,
                                       scale_by_area=False,
                                       membrane_flux=False)


#CURRENTS/LEAKS ----------------------------------------------------------------
#sodium (Na) current
if args['dendL'] > 0:
    if pcid == 0:
        print('L = ' + str(args['dendL']))
    def difDend(nd):
        if nd.sec.name().split('.')[-1] == 'soma':
            return 1.0
        return 0.0
    dendP = rxd.Parameter([mem], name='dendP', value = lambda nd: difDend(nd))
    na_current = rxd.MultiCompartmentReaction(nai, nao, dendP * gna * (v - ena),
                                            mass_action=False, membrane=mem,
                                            membrane_flux=True)
else:
    if pcid == 0:
        print('sphere')
    na_current = rxd.MultiCompartmentReaction(nai, nao, gna * (v - ena),
                                            mass_action=False, membrane=mem,
                                            membrane_flux=True)
#potassium (K) current
k_current = rxd.MultiCompartmentReaction(ki, ko, gk * (v - ek),
                                         mass_action=False, membrane=mem,
                                         membrane_flux=True)
#nkcc1 (Na+/K+/2Cl- cotransporter)
nkcc1_current1 = rxd.MultiCompartmentReaction(cli, clo, 2.0 * nkcc1 * volume_scale,
                                              mass_action=False, membrane=mem,  
                                              membrane_flux=True)
nkcc1_current2 = rxd.MultiCompartmentReaction(ki, ko, nkcc1*volume_scale,
                                              mass_action=False, membrane=mem,
                                              membrane_flux=True)
nkcc1_current3 = rxd.MultiCompartmentReaction(nai, nao, nkcc1*volume_scale,
                                              mass_action=False, membrane=mem,
                                              membrane_flux=True)
#kcc2 (K+/Cl- cotransporter)
kcc2_current1 = rxd.MultiCompartmentReaction(cli, clo, kcc2*volume_scale,
                                             membrane=mem, mass_action=False,
                                             membrane_flux=True)
kcc2_current2 = rxd.MultiCompartmentReaction(ki, ko, kcc2*volume_scale,
                                             membrane=mem, mass_action=False,
                                             membrane_flux=True)
#sodium leak
na_leak = rxd.MultiCompartmentReaction(nai, nao, gnabar_l*(v - ena),
                                       mass_action=False, membrane=mem,
                                       membrane_flux=True)
#potassium leak
k_leak = rxd.MultiCompartmentReaction(ki, ko, gkbar_l*(v - ek),
                                      mass_action=False, membrane=mem,
                                      membrane_flux=True)
#chlorine (Cl) leak                              
cl_current = rxd.MultiCompartmentReaction(cli, clo, gclbar_l * (ecl - v),
                                          mass_action=False, membrane=mem,
                                          membrane_flux=True)


if args['ouabain']:
    #Na+/K+ pump current in neuron (2K+ in, 3Na+ out)
    pump_current = rxd.MultiCompartmentReaction(ki, ko, -2.0*pump*volume_scale*notcore,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True) 
    pump_current_na = rxd.MultiCompartmentReaction(nai, nao, 3.0*pump*volume_scale*notcore,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True)
    #Na+/K+ pump current in glia (2K+ in, 3Na+ out)
    gpump_current_na = rxd.Rate(nao, 3.0*gliapump*notcore)
    #Glia K+ current 
    glia_k_current = rxd.Rate(ko, -glia12 - 2*gliapump*notcore)
    # O2 dynamics
    o2_pump = rxd.Rate(o2ecs, -gliapump * notcore)
    oxygen = rxd.MultiCompartmentReaction(o2ecs, dump[cyt], pump * volume_scale * notcore,
                                        mass_action=False, membrane=mem)
elif args['ischemCore']:
    #Na+/K+ pump current in neuron (2K+ in, 3Na+ out)
    pump_current = rxd.MultiCompartmentReaction(ki, ko, -2.0*pump*volume_scale,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True) 
    pump_current_na = rxd.MultiCompartmentReaction(nai, nao, 3.0*pump*volume_scale,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True)
    #Na+/K+ pump current in glia (2K+ in, 3Na+ out)
    gpump_current_na = rxd.Rate(nao, 3.0*gliapump)
    #Glia K+ current 
    glia_k_current = rxd.Rate(ko, -glia12 - 2*gliapump)
    # O2 dynamics
    o2_pump = rxd.Rate(o2ecs, -gliapump * iscore)
    oxygen = rxd.MultiCompartmentReaction(o2ecs, dump[cyt], pump * volume_scale * iscore,
                                        mass_action=False, membrane=mem)
elif args['edemaCore']:
    #Na+/K+ pump current in neuron (2K+ in, 3Na+ out)
    pump_current = rxd.MultiCompartmentReaction(ki, ko, -2.0*pump*volume_scale,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True) 
    pump_current_na = rxd.MultiCompartmentReaction(nai, nao, 3.0*pump*volume_scale,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True)
    #Na+/K+ pump current in glia (2K+ in, 3Na+ out)
    gpump_current_na = rxd.Rate(nao, 3.0*gliapump)
    #Glia K+ current 
    glia_k_current = rxd.Rate(ko, -glia12 - 2*gliapump)
elif args['ischemEdemaCore']:
    #Na+/K+ pump current in neuron (2K+ in, 3Na+ out)
    pump_current = rxd.MultiCompartmentReaction(ki, ko, -2.0*pump*volume_scale,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True) 
    pump_current_na = rxd.MultiCompartmentReaction(nai, nao, 3.0*pump*volume_scale,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True)
    #Na+/K+ pump current in glia (2K+ in, 3Na+ out)
    gpump_current_na = rxd.Rate(nao, 3.0*gliapump)
    #Glia K+ current 
    glia_k_current = rxd.Rate(ko, -glia12 - 2*gliapump)
    # #O2 dynamics
    o2_pump = rxd.Rate(o2ecs, -gliapump * iscore)
    oxygen = rxd.MultiCompartmentReaction(o2ecs, dump[cyt], pump * volume_scale * iscore,
                                        mass_action=False, membrane=mem)
elif args['O2consume']:
    #Na+/K+ pump current in neuron (2K+ in, 3Na+ out)
    pump_current = rxd.MultiCompartmentReaction(ki, ko, -2.0*pump*volume_scale,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True) 
    pump_current_na = rxd.MultiCompartmentReaction(nai, nao, 3.0*pump*volume_scale,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True)
    #Na+/K+ pump current in glia (2K+ in, 3Na+ out)
    gpump_current_na = rxd.Rate(nao, 3.0*gliapump)
    #Glia K+ current 
    glia_k_current = rxd.Rate(ko, -glia12 - 2*gliapump)
    o2_pump = rxd.Rate(o2ecs, -gliapump)
    oxygen = rxd.MultiCompartmentReaction(o2ecs, dump[cyt], pump * volume_scale,
                                          mass_action=False, membrane=mem)
else:
    #Na+/K+ pump current in neuron (2K+ in, 3Na+ out)
    pump_current = rxd.MultiCompartmentReaction(ki, ko, -2.0*pump*volume_scale,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True) 
    pump_current_na = rxd.MultiCompartmentReaction(nai, nao, 3.0*pump*volume_scale,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True)
    #Na+/K+ pump current in glia (2K+ in, 3Na+ out)
    gpump_current_na = rxd.Rate(nao, 3.0*gliapump)
    #Glia K+ current 
    glia_k_current = rxd.Rate(ko, -glia12 - 2*gliapump)

pc.set_maxstep(100) # required when using multiple processes

t = h.Vector().record(h._ref_t)
soma_v = []
soma_ki = []
soma_nai = []
soma_cli = []
soma_nao = []
soma_ko = []
soma_clo = []
soma_o2 = []
soma_vir = []
soma_vor = []
rpos = []

cell_positions = [(sec.x3d(0)**2 + sec.y3d(0)**2 + sec.z3d(0)**2)**(0.5) for sec in soma_list]

def saveconc():
    numpy.save(os.path.join(outdir,"k_%i.npy" % int(h.t)), k[ecs].states3d)
    numpy.save(os.path.join(outdir,"na_%i.npy" % int(h.t)), na[ecs].states3d)
    numpy.save(os.path.join(outdir,"cl_%i.npy" % int(h.t)), cl[ecs].states3d)
    numpy.save(os.path.join(outdir,'o2_%i.npy' % int(h.t)), o2ecs.states3d)

for i in range(int(Lx//10)):
    # for r, soma in zip(cell_positions, h.allsec()):
    for r, soma in zip(cell_positions, soma_list):
        if (10.0*i-2.5) < r < (10.0*i+2.5):
            print(i,r)
            rpos.append((soma.x3d(0), soma.y3d(0), soma.z3d(0)))
            soma_v.append(h.Vector().record(soma(0.5)._ref_v))
            soma_nai.append(h.Vector().record(soma(0.5)._ref_nai))
            soma_ki.append(h.Vector().record(soma(0.5)._ref_ki))
            soma_cli.append(h.Vector().record(soma(0.5)._ref_cli))

            soma_nao.append(h.Vector().record(soma(0.5)._ref_nao))
            soma_ko.append(h.Vector().record(soma(0.5)._ref_ko))
            soma_clo.append(h.Vector().record(soma(0.5)._ref_clo))


            soma_o2.append(h.Vector().record(o2ecs.node_by_location(soma.x3d(0),soma.y3d(0),soma.z3d(0))._ref_concentration))
            soma_vir.append(h.Vector().record(soma(0.5)._ref_volumei))
            soma_vor.append(h.Vector().record(vor.node_by_location(soma.x3d(0),soma.y3d(0),soma.z3d(0))._ref_value))
            break

recs = {'v':soma_v, 'ki':soma_ki, 'nai':soma_nai, 'cli':soma_cli,
        't':t,      'ko':soma_ko, 'nao':soma_nao, 'clo':soma_clo,
        'pos':rpos, 'o2':soma_o2, 'vir':soma_vir, 'vor':soma_vor,
        'rad':cell_positions}

# initialize and set the intracellular concentrations

def progress_bar(tstop, size=40):
    """ report progress of the simulation """
    prog = h.t/float(tstop)
    fill = int(size*prog)
    empt = size - fill
    progress = '#' * fill + '-' * empt
    sys.stdout.write('[%s] %2.1f%% %6.1fms of %6.1fms\r' % (progress, 100*prog, h.t, tstop))
    sys.stdout.flush()

def plot_rec_neurons():
    """ Produces plots of record neurons membrane potential (shown in figure 1C) """
    # load the recorded neuron data
    somaV, pos = [], []
    for i in range(nhost):
        fin = open(os.path.join(outdir,'membrane_potential_%i.pkl' % i),'rb')
        [sV, p] = pickle.load(fin)
        fin.close()
        somaV.extend(sV)
        pos.extend(p)

        for idx in range(somaV[0].size()): 
            # create a plot for each record (100ms)

            fig = pyplot.figure()
            ax = fig.add_subplot(111,projection='3d')
            ax.set_position([0.0,0.05,0.9,0.9])
            ax.set_xlim([-Lx/2.0, Lx/2.0])
            ax.set_ylim([-Ly/2.0, Ly/2.0])
            ax.set_zlim([-Lz/2.0, Lz/2.0])
            ax.set_xticks([int(Lx*i/4.0) for i in range(-2,3)])
            ax.set_yticks([int(Ly*i/4.0) for i in range(-2,3)])
            ax.set_zticks([int(Lz*i/4.0) for i in range(-2,3)])

            cmap = pyplot.get_cmap('jet')
            for ii in range(Nrec):
                x = pos[ii]
                soma_z = [x[2]-somaR,x[2]+somaR]
                cell_x = [x[0],x[0]]
                cell_y = [x[1],x[1]]
                scolor = cmap((somaV[ii].get(idx)+70.0)/70.0)
                # plot the soma
                ax.plot(cell_x, cell_y, soma_z, linewidth=2, color=scolor, 
                        alpha=0.5)
    
            norm = colors.Normalize(vmin=-70,vmax=0)
            pyplot.title('Neuron membrane potentials; t = %gms' % (idx * 100))

            # add a colorbar 
            ax1 = fig.add_axes([0.88,0.05,0.04,0.9])
            cb1 = colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,
                                                orientation='vertical')
            cb1.set_label('mV')
            
            # save the plot
            filename = 'neurons_{:05d}.png'.format(idx)
            pyplot.savefig(os.path.join(outdir,filename))
            pyplot.close()

def plot_image_data(data, min_val, max_val, filename, title):
    """Plot a 2d image of the data"""
    # sb = scalebar.ScaleBar(1e-6)
    # sb.location='lower left'
    pyplot.imshow(data, extent=k[ecs].extent('xy'), vmin=min_val,
                  vmax=max_val, interpolation='nearest', origin='lower')
    pyplot.colorbar()
    # sb = scalebar.ScaleBar(1e-6)
    # sb.location='lower left'
    ax = pyplot.gca()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    # ax.add_artist(sb)
    pyplot.title(title)
    pyplot.xlim(k[ecs].extent('x'))
    pyplot.ylim(k[ecs].extent('y'))
    pyplot.savefig(os.path.join(outdir,filename))
    pyplot.close()


def boxoff(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def plotVm(name):
    fig =  pyplot.figure(dpi=200)
    ax = pyplot.subplot(111)
    # pyplot.plot(cell_positions, [sec.v for sec in h.allsec()], '.')
    pyplot.plot(cell_positions, [sec.v for sec in soma_list], '.')
    pyplot.xlabel('distance (μm)')
    pyplot.ylabel('membrane potential (mV)')
    pyplot.ylim([-80,40])
    pyplot.title('t = %6.0fms' % h.t)
    boxoff(ax)
    fig.savefig(name)
    pyplot.close()

def run(tstop):
    """ Run the simulations saving figures every 100ms and recording the wave progression every time step"""
    if pcid == 0:
        # record the wave progress (shown in figure 2)
        name = ''
        fout = open(os.path.join(outdir,'wave_progress%s.txt' % name),'a')
    last_plot = 0
    last_print = 0
    plotnum = 0
    bwinsz = 10
    time = []
    saveint = 100
    dumpint = 1000
    last_dump = 0

    while h.t < tstop:
        time.append(h.t)
        if int(h.t) % saveint == 0:
            # plot extracellular concentrations averaged over depth every 100ms 
            if pcid == 0:
                plot_image_data(k[ecs].states3d.mean(2), 3.5, 40,
                                'k_mean_%05d' % int(h.t/100),
                                'Potassium concentration; t = %6.0fms'
                                % h.t)
                                
                plot_image_data(o2ecs.states3d.mean(2), oa_bath * 0.99, oa_bath, 
                                'o2_mean_%05d' % int(h.t/100),
                                'Oxygen concentration; t = %6.0fms'
                                % h.t)

                saveconc()

            if pcid == nhost - 1:
                plot_image_data(na[ecs].states3d.mean(2), 120, 150,
                                'na_mean_%05d' % int(h.t/100),
                                'Sodium concentration; t = %6.0fms'
                                % h.t)
                
                plot_image_data(cl[ecs].states3d.mean(2), 100, 150,
                                'cl_mean_%05d' % int(h.t/100),
                                'Chloride concentration; t = %6.0fms'
                                % h.t)

        if pcid == 0: progress_bar(tstop)
        pc.psolve(pc.t(0)+h.dt)  # run the simulation for 1 time step

        # h.fadvance()
        # determine the furthest distance from the origin where
        # extracellular potassium exceeds Kceil (dist)
        # And the shortest distance from the origin where the extracellular
        # extracellular potassium is below Kceil (dist1)
        if pcid == 0 and h.t - last_print > 1.0:
            last_print = h.t
            dist = 0
            dist1 = 1e9
            for nd in k.nodes:
                if str(nd.region).split('(')[0] == 'Extracellular':
                    r = (nd.x3d**2+nd.y3d**2+nd.z3d**2)**0.5
                    if nd.concentration>Kceil and r > dist:
                        dist = r
                    if nd.concentration<=Kceil and r < dist1:
                        dist1 = r
            fout.write("%g\t%g\t%g\n" %(h.t, dist, dist1))
            fout.flush()

    if pcid == 0:
        progress_bar(tstop)
        fout.close()
        with open(os.path.join(outdir,"recs.pkl"),'wb') as fout:
            pickle.dump(recs,fout)
        print("\nSimulation complete. Plotting membrane potentials")

    # # save membrane potentials
    soma, pos = [], []
    for n in center_rec_neurons:
        soma.append(n.somaV)
        pos.append([n.x,n.y,n.z])
    with open(os.path.join(outdir,"centermembrane_potential_%i.pkl" % pcid),'wb') as pout:
        pickle.dump([soma,pos,time],pout)
    if periph_rec_neurons:
        soma, pos = [], []
        for n in periph_rec_neurons:
            soma.append(n.somaV)
            pos.append([n.x,n.y,n.z])
        with open(os.path.join(outdir,"periphmembrane_potential_%i.pkl" % pcid),'wb') as pout:
            pickle.dump([soma,pos,time],pout)    

    pc.barrier()    # wait for all processes to save

## restore from previous sim 
if args['restoredir']:
    restoredir = args['restoredir']

    # restore sim state functions 
    def restoreSS():
        svst = h.SaveState()
        f = h.File(os.path.join(restoredir, 'save_test_'+str(pcid) + '.dat'))
        svst.fread(f)
        svst.restore()

    def restoreSim():
        restoreSS()
        for sp in rxd.species._all_species:
            s = sp()
            print(s.name)
            temp = numpy.load(os.path.join(restoredir, s.name + '_concentrations_' + str(pcid) + '.npy'))
            s.nodes.concentration = list(temp)
            print('PCID ' + str(pcid) + ': resotred ' + s.name)

    fih = h.FInitializeHandler(1, restoreSim)
    h.finitialize()
else:
    h.finitialize(v_initial)

# run the simulation
run(args['tstop'])

# save final sim state 
runSS()
saveRxd()

pc.barrier()  
h.quit()

# v0.0 - realisitc cell dendity for cortex, increased SA to V ratio, cmdline specify anoxic vs normoxic
# v0.1 - fixed dend diameter, dend length now cmdline specified
# v0.2 - expanded dimensions and number recorded cells, fixed vel smoothing
# v0.3 - back to original dimensions, plotting K concentration and membrane potential together to generate a gif, doubled number of threads
# v0.4 - r0 = 75, Lz = 500, initial elevate K+ 70mM
# v0.5 - fixed [K+] plotting, changed method for calculating SD wave position
# v0.6 - changing wave position method again, including firing back into elevated K+ radius, allow user to change between invivo and invitro boundary conds
# v0.7 - trying parallel context, still working on wave front postition
# v0.8 - fixed so cell position and voltages are only somatic, no more dendritic Vs
# v0.9 - reverted back to original calculations of wavefront position and velocity. changed r0 to 100um
# v0.10 - cut down dt, adding code for raster plots, potential fix for K wave progress issues
# v0.11 - fixing typos, turn off timeout
# v0.12 - 1mm^3, save raster data
# v0.13 - back down to former volume, 24 threads, and ecs [k+] in core up to 6
# v0.14 - try out recording lfps with LFPsimpy
# v0.15 - drop ecs [k+] in core back down to 40, save raster plot
# v0.16 - more user specified cmd line args, double volume
# v0.17 - reduced oxygen by 95% in core, try keeping [k+] initially uniform
# v0.18 - user specified o2 factor
# v0.19 - trying fovea-like dimensions
# v0.20 - changes to boundary conds for in vivo, diff in Lz for invivo vs invitro, save concs every 100ms
# v0.21 - resolved overwriting oxygen rate
# v0.22 - added user specification of O2 bath value
# v0.23 - remove Na flux from dendrites
# v0.24 - remove mg/mL*s conversion factor alpha
# v0.25 - switched alpha back for now
# v0.26 - translated o2switch by 3e-4 mM O2 hoping to resolve negative [O2]
# v0.27 - applied o2switch to o2 consumption by Na/K pumps rather than pump activity to account for anaerobic glycolysis
# v0.28 - switch back on pump, upped translation to 5e-4
# v0.29 - trying constant infusion of K+, turn of saving figs for gifs
# v0.30 - added option for 2mm x 2mm x 170 um
# v0.31 - removed anox oa_bath
# v0.32 - removed dependance on o2switch altogether
# v0.33 - neuronal pump x100
# v0.34 - back to original oa_bath, no increase in pump activity
# v0.35 - user speecifies whether to infuse and factors for neural and glial Na/K pumps
# v0.36 - user specification of cell density and option for brainstem-like volume fraction
# v0.37 - attempt at 10x pump use of O2 effeciency (as if 10x more o2ecs) 
# v0.38 - adding stimulation
# v0.39 - attempting state saving
# v0.40 - only save rxd using collect... use savestate for ephys
# v0.41 - cleaned up saving, not useing SaveState
# v0.42 - reinstated oa_bath differences, ditched saving, no stims, mini size
# v0.43 - reduce saving interval for mini sized / long run sims 
# v0.44 - user specifies surface area to vol ratio, calculates cell radius and fractional cell volume to keep neuron volume fraction 0.24
# v0.45 - save state at end of the sim, option to restore state 
# v0.46 - ditch mpi4py, run with nrniv rather than python3
# v0.46.1 - ditched argparse, loading arguments from json file 
# v0.46.2 - tried moving reinstation to just before calling run()
# v0.47 - user specification of p_max
# v0.47.1 - looking for no depol for normox by increasing pump activity - changed parameter in equation for p
# v0.47.2 - added interval saving every 5 seconds
# v0.47.3 - user specifies param in equation for p, toggles O2 consumption
# v0.48 - allow for ischemic core where O2 is reduced and not replenished and ECS propertied are altered
# v0.48.2 - allow for core with reduced O2 but replenished in addition to previous configuration
# v0.48.3 - interval saving of membrane potential files as well
# v0.49 - separate ecs for o2 reflecting its a gas 
# v0.49.1 - still working on spontaneous depol issue, kcc2 was being overwritten, trying first version 
# v0.49.2 - back to original kcc2 and inserted pas
# v0.49.3 - make sure inufuse is just true, not necessarily yes
# v0.49.4 - user specifies gpas 
# v0.49.5 - user specifies Lz 
# v0.50 - abandoned interval saving, recs_#.pkl files were corrupted
# v0.51 - placing rec_neurons at the same depth (z=0)
# v0.51.1 - two sets of rec neurons, one at center another near the margin
# v0.51.2 - toggle between uniform recordings and at specific depths
# v0.51.3 - centermembrane needs to be reflected for raster plotting
# v0.51.4 - evenly distributed the center and peripheral recorded neurons
# v0.51.5 - includes alphas for primed with propionate and brainstem, changed periph vs center layout
# v0.52 - original concentrations, option for original o2 bath, user specified random seeds
# v0.53 - setup for cores of ouabain, ischemia, and edema 
# v0.54 - user specifies neuronal volume fraction, sim computes appropriate radius for cell volume 
# v0.55 - nonuniform recording setup now has center as random pos in middle 50%, periph as random pos in top/bottom 25%
# v0.56 - fixed lambda for mannitol, added user specification of tortuosity and O2 with otherwise perfused parameters
# v1.00 - cleaned up version for fresh repo associated w/ the paper 