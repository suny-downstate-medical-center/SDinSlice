"""Simulation of spreading depression"""
# from mpi4py import MPI
import sys
sys.path.insert(0, '../nrn/')
from numpy import random
from neuron import h, rxd
#h.nrnmpi_init()
from neuron.rxd import v
from neuron.rxd.rxdmath import exp, log, tanh
from neuron.units import sec, mM, ms
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot, colors, colorbar
# from matplotlib_scalebar import scalebar
from mpl_toolkits.mplot3d import Axes3D
import math
from math import pi
from matplotlib import pyplot
import numpy
import os
import pickle
from pylab import convolve
import json 
import itertools
from scipy.signal import find_peaks 
import datetime 

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
# numpy.random.seed(2549637+pcid) # different random seed for some low s:v sims

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
        print(datetime.datetime.now())
    if args['uniformRec']:
        center_rec_neurons = [Neuron(
            (numpy.random.random()*2.0 - 1.0) * (Lx/2.0 - somaR), 
            (numpy.random.random()*2.0 - 1.0) * (Ly/2.0 - somaR), 
            (numpy.random.random()*2.0 - 1.0) * (Lz/2.0 - somaR), 50)
            for i in range(0, int(Nrec/nhost))]
        periph_rec_neurons = []
    else:
        # xy = numpy.linspace(0, Lx/2 - somaR - 1, nhost)
        # center_rec_neurons = [Neuron(
        #     xy[pcid], 
        #     xy[pcid], 
        #     0, 50)
        #     for i in range(0, int(Nrec/nhost))]
        # periph_rec_neurons = [Neuron(
        #     xy[pcid], 
        #     xy[pcid], 
        #     Lz/2 - somaR - 1, 50)
        #     for i in range(0, int(Nrec/nhost))]
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
    #(numpy.random.random()*0.5 + numpy.random.choice((-1, 0.5))) * (Lz/2.0 - somaR), 50)
    # Randomly distribute the remaining neurons
    all_neurons = [Neuron(
        (numpy.random.random()*2.0 - 1.0) * (Lx/2.0 - somaR),
        (numpy.random.random()*2.0 - 1.0) * (Ly/2.0 - somaR),
        (numpy.random.random()*2.0 - 1.0) * (Lz/2.0 - somaR))
        for i in range(int(Nrec/nhost), int(Ncell/nhost))]

## stimulation of recorded cells 
# stims = [] 
# amps = numpy.linspace(0.005, 0.03, 15, endpoint=True)
# amps = numpy.logspace(numpy.log10(0.005), numpy.log10(0.1), 15, endpoint=True)
# amps = numpy.logspace(numpy.log10(0.01), numpy.log10(0.5), 15, endpoint=True)

# for nrn in rec_neurons:
#     # delay = 1000
#     delay = 100000
#     for amp in amps:
#         stims.append(h.IClamp(nrn.soma(0.5)))
#         stims[-1].dur = 250
#         stims[-1].delay = delay 
#         stims[-1].amp = amp
#         delay = delay + 500


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

# o2_extracellular = rxd.Species([ecs], name='o2', d=1,
#     initial=lambda nd: oa_bath * args.o2factor
#     if nd.x3d**2 + nd.y3d**2 + nd.z3d**2 <= r0**2 else oa_bath,
#     ecs_boundary_conditions=o2_bc)
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
dumpi = dump[cyt]
dumpo = dump[ecs]

# if args['BC'] == 'invitro':
ecsbc = rxd.Parameter([ecs, ecs_o2], name='ecsbc', value = lambda nd: bc(nd))
# else:
#     ecsbc = rxd.Parameter([ecs], name='ecsbc', value = lambda nd: bcvivo(nd))

ki, ko, nai, nao, cli, clo = k[cyt], k[ecs], na[cyt], na[ecs], cl[cyt], cl[ecs]

## dynamic ecs
initial_alpha = alpha1
alpha_ecs = rxd.State([ecs], name='alpha_ecs', initial=initial_alpha)
dalpha_ecs = rxd.Rate(alpha_ecs, 5e-7 * (ko-3.5) * (ko-70.0) * alpha_ecs)
ecs.alpha = alpha_ecs

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
# nkcc1 = unkcc1*fko*(log((ki * cli * vor**2) / (ko * clo * vir**2)) + log((nai * cli * vor**2) / (nao * clo * vir**2)))
# kcc2 = ukcc2 * log((ki * cli * vor**2) / (ko * clo * vir**2))

#Nerst equation - reversal potentials
ena = 26.64 * log(nao*vir/(nai*vor))
ek = 26.64 * log(ko*vir/(ki*vor))
ecl = 26.64 * log(cli*vor/(clo*vir))

# p = p_max / (1.0 + exp((20.0 - (o2ecs/vor) * alpha)/3.0))
# p = o2switch * p_max / (1.0 + exp((20.0 - (o2ecs/vor) * alpha)/3.0))
p = o2switch * p_max / (1.0 + exp((args['pparam'] - (o2ecs/vor) * alpha)/3.0))
# p = o2switch * p_max / (1.0 + exp((20.0 - (o2ecs/vor))/3.0))
pump = args['nrnPumpFactor'] * (p / (1.0 + exp((25.0 - nai / vir)/3.0))) * (1.0 / (1.0 + exp(3.5 - ko / vor)))
gliapump = args['glialPumpFactor'] * (1.0/3.0) * (p / (1.0 + exp((25.0 - gnai_initial) / 3.0))) * (1.0 / (1.0 + exp(3.5 - ko/vor)))
g_glia = g_gliamax / (1.0 + exp(-(o2ecs*alpha/vor - 2.5)/0.2))
# g_glia = g_gliamax / (1.0 + exp(-(o2ecs/vor - 2.5)/0.2))
glia12 = g_glia / (1.0 + exp((18.0 - ko / vor)/2.5))

epsilon_k = (epsilon_k_max/(1.0 + exp(-((o2ecs/vor) * alpha - 2.5)/0.2))) * (1.0/(1.0 + exp((-20 + ((1.0+1.0/beta0 -vor)/vor) /2.0))))
# epsilon_k = (epsilon_k_max/(1.0 + exp(-((o2ecs/vor) - 2.5)/0.2))) * (1.0/(1.0 + exp((-20 + ((1.0+1.0/beta0 -vor)/vor) /2.0))))
           
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
                                       -scalei*vtau*(osm-vir) * initial_alpha / alpha_ecs[ecs],
                                       mass_action=False,
                                       membrane=mem,
                                       scale_by_area=False,
                                       membrane_flux=False)

vol_dyn_ecs = rxd.MultiCompartmentReaction(dump[cyt], vor,
                                       -scaleo*vtau*(osm-vir) * initial_alpha / alpha_ecs[ecs],
                                       mass_action=False,
                                       membrane=mem,
                                       scale_by_area=False,
                                       membrane_flux=False)

if pcid == 0:
    print(datetime.datetime.now())
    print('Vol Dyn Complete\n')


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
        print(datetime.datetime.now())
        print('sphere')
    # na_current = rxd.MultiCompartmentReaction(nai, nao, gna * (v - ena),
    #                                         mass_action=False, membrane=mem,
    #                                         membrane_flux=True)
    na_current_1 = rxd.MultiCompartmentReaction(dumpi, nao, gna * (v - ena) * initial_alpha / alpha_ecs[ecs],
                                            mass_action=False, membrane=mem,
                                            membrane_flux=True)
    na_current_2 = rxd.MultiCompartmentReaction(nai, dumpo, gna * (v - ena),
                                            mass_action=False, membrane=mem,
                                            membrane_flux=True)

if pcid == 0:
    print(datetime.datetime.now())
    print('Na Complete\n')
#potassium (K) current
# k_current = rxd.MultiCompartmentReaction(ki, ko, gk * (v - ek),
#                                          mass_action=False, membrane=mem,
#                                          membrane_flux=True)
k_current_1 = rxd.MultiCompartmentReaction(dumpi, ko, gk * (v - ek) * initial_alpha / alpha_ecs[ecs],
                                         mass_action=False, membrane=mem,
                                         membrane_flux=True)
k_current_2 = rxd.MultiCompartmentReaction(ki, dumpo, gk * (v - ek),
                                         mass_action=False, membrane=mem,
                                         membrane_flux=True)
if pcid == 0:
    print(datetime.datetime.now())
    print('K Complete\n')
#nkcc1 (Na+/K+/2Cl- cotransporter)
# nkcc1_current1 = rxd.MultiCompartmentReaction(cli, clo, 2.0 * nkcc1 * volume_scale,
#                                               mass_action=False, membrane=mem,  
#                                               membrane_flux=True)
nkcc1_current1_1 = rxd.MultiCompartmentReaction(dumpi, clo, 2.0 * nkcc1 * volume_scale * initial_alpha / alpha_ecs[ecs],
                                              mass_action=False, membrane=mem,  
                                              membrane_flux=True)
nkcc1_current1_2 = rxd.MultiCompartmentReaction(cli, dumpo, 2.0 * nkcc1 * volume_scale,
                                              mass_action=False, membrane=mem,  
                                              membrane_flux=True)
if pcid == 0:
    print(datetime.datetime.now())
    print('NKCC1 Current 1 Complete\n')

# nkcc1_current2 = rxd.MultiCompartmentReaction(ki, ko, nkcc1*volume_scale,
#                                               mass_action=False, membrane=mem,
#                                               membrane_flux=True)
nkcc1_current2_1 = rxd.MultiCompartmentReaction(dumpi, ko, nkcc1*volume_scale*initial_alpha/alpha_ecs[ecs],
                                              mass_action=False, membrane=mem,
                                              membrane_flux=True)
nkcc1_current2_2 = rxd.MultiCompartmentReaction(ki, dumpo, nkcc1*volume_scale,
                                              mass_action=False, membrane=mem,
                                              membrane_flux=True)
if pcid == 0:
    print(datetime.datetime.now())
    print('NKCC1 Current 2 Complete\n')
    
# nkcc1_current3 = rxd.MultiCompartmentReaction(nai, nao, nkcc1*volume_scale,
#                                               mass_action=False, membrane=mem,
#                                               membrane_flux=True)
nkcc1_current3_1 = rxd.MultiCompartmentReaction(dumpi, nao, nkcc1*volume_scale*initial_alpha/alpha_ecs[ecs],
                                              mass_action=False, membrane=mem,
                                              membrane_flux=True)
nkcc1_current3_2 = rxd.MultiCompartmentReaction(nai, dumpo, nkcc1*volume_scale,
                                              mass_action=False, membrane=mem,
                                              membrane_flux=True)
if pcid == 0:
    print(datetime.datetime.now())
    print('NKCC1 Current 3 Complete\n')
#kcc2 (K+/Cl- cotransporter)
# kcc2_current1 = rxd.MultiCompartmentReaction(cli, clo, kcc2*volume_scale,
#                                              membrane=mem, mass_action=False,
#                                              membrane_flux=True)
kcc2_current1_1 = rxd.MultiCompartmentReaction(dumpi, clo, kcc2*volume_scale*initial_alpha/alpha_ecs[ecs],
                                             membrane=mem, mass_action=False,
                                             membrane_flux=True)
kcc2_current1_2 = rxd.MultiCompartmentReaction(cli, dumpo, kcc2*volume_scale,
                                             membrane=mem, mass_action=False,
                                             membrane_flux=True)
if pcid == 0:
    print(datetime.datetime.now())
    print('KCC1 Current 1 Complete\n')
# kcc2_current2 = rxd.MultiCompartmentReaction(ki, ko, kcc2*volume_scale,
#                                              membrane=mem, mass_action=False,
#                                              membrane_flux=True)
kcc2_current2_1 = rxd.MultiCompartmentReaction(dumpi, ko, kcc2*volume_scale*initial_alpha/alpha_ecs[ecs],
                                             membrane=mem, mass_action=False,
                                             membrane_flux=True)
kcc2_current2_2 = rxd.MultiCompartmentReaction(ki, dumpo, kcc2*volume_scale,
                                             membrane=mem, mass_action=False,
                                             membrane_flux=True)
if pcid == 0:
    print(datetime.datetime.now())
    print('KCC1 Current 2 Complete\n')
#sodium leak
# na_leak = rxd.MultiCompartmentReaction(nai, nao, gnabar_l*(v - ena),
#                                        mass_action=False, membrane=mem,
#                                        membrane_flux=True)
na_leak_1 = rxd.MultiCompartmentReaction(dumpi, nao, gnabar_l*(v - ena)*initial_alpha/alpha_ecs[ecs],
                                       mass_action=False, membrane=mem,
                                       membrane_flux=True)
na_leak_2 = rxd.MultiCompartmentReaction(nai, nao, gnabar_l*(v - ena)*initial_alpha/alpha_ecs,
                                       mass_action=False, membrane=mem,
                                       membrane_flux=True)
if pcid == 0:
    print(datetime.datetime.now())
    print('Na Leak Complete\n')
#potassium leak
# k_leak = rxd.MultiCompartmentReaction(ki, ko, gkbar_l*(v - ek),
#                                       mass_action=False, membrane=mem,
#                                       membrane_flux=True)
k_leak_1 = rxd.MultiCompartmentReaction(dumpi, ko, gkbar_l*(v - ek) * initial_alpha/alpha_ecs[ecs],
                                      mass_action=False, membrane=mem,
                                      membrane_flux=True)
k_leak_2 = rxd.MultiCompartmentReaction(ki, dumpo, gkbar_l*(v - ek),
                                      mass_action=False, membrane=mem,
                                      membrane_flux=True)
if pcid == 0:
    print(datetime.datetime.now())
    print('K Leak Complete\n')
#chlorine (Cl) leak                              
# cl_current = rxd.MultiCompartmentReaction(cli, clo, gclbar_l * (ecl - v),
#                                           mass_action=False, membrane=mem,
#                                           membrane_flux=True)
cl_current_1 = rxd.MultiCompartmentReaction(dumpi, clo, gclbar_l * (ecl - v)*initial_alpha/alpha_ecs[ecs],
                                          mass_action=False, membrane=mem,
                                          membrane_flux=True)
cl_current_2 = rxd.MultiCompartmentReaction(cli, dumpo, gclbar_l * (ecl - v),
                                          mass_action=False, membrane=mem,
                                          membrane_flux=True)
if pcid == 0:
    print(datetime.datetime.now())
    print('Cl Leak Complete\n')

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
    # #O2 dynamics
    # o2_pump = rxd.Rate(o2ecs, -gliapump * iscore)
    # oxygen = rxd.MultiCompartmentReaction(o2ecs, dump[cyt], pump * volume_scale * iscore,
    #                                     mass_action=False, membrane=mem)
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
    # pump_current = rxd.MultiCompartmentReaction(ki, ko, -2.0*pump*volume_scale,
    #                                             mass_action=False, membrane=mem,
    #                                             membrane_flux=True)
    pump_current_1 = rxd.MultiCompartmentReaction(dumpi, ko, -2.0*pump*volume_scale*initial_alpha/alpha_ecs[ecs],
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True)
    pump_current_2 = rxd.MultiCompartmentReaction(ki, dumpo, -2.0*pump*volume_scale,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True) 

    if pcid == 0:
        print(datetime.datetime.now())
        print('K pump current Complete\n')
    # pump_current_na = rxd.MultiCompartmentReaction(nai, nao, 3.0*pump*volume_scale,
    #                                             mass_action=False, membrane=mem,
    #                                             membrane_flux=True)
    pump_current_na_1 = rxd.MultiCompartmentReaction(dumpi, nao, 3.0*pump*volume_scale*initial_alpha/alpha_ecs[ecs],
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True)
    pump_current_na_2 = rxd.MultiCompartmentReaction(nai, dumpo, 3.0*pump*volume_scale,
                                                mass_action=False, membrane=mem,
                                                membrane_flux=True)
    if pcid == 0:
        print(datetime.datetime.now())
        print('K pump current Complete\n')
    #Na+/K+ pump current in glia (2K+ in, 3Na+ out)
    gpump_current_na = rxd.Rate(nao, 3.0*gliapump)
    #Glia K+ current 
    glia_k_current = rxd.Rate(ko, -glia12 - 2*gliapump)
    o2_pump = rxd.Rate(o2ecs, -gliapump)
    oxygen = rxd.MultiCompartmentReaction(o2ecs, dump[cyt], pump * volume_scale,
                                          mass_action=False, membrane=mem)
    if pcid == 0:
        print(datetime.datetime.now())
        print('O2 Complete\n')
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

# o2_pump = rxd.Rate(o2ecs, -gliapump * o2switch)
# oxygen = rxd.MultiCompartmentReaction(o2ecs, dump[cyt], pump * volume_scale * o2switch,
#                                       mass_action=False, membrane=mem)

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

# cell_positions = [(sec.x3d(0)**2 + sec.y3d(0)**2 + sec.z3d(0)**2)**(0.5) for sec in h.allsec()]
cell_positions = [(sec.x3d(0)**2 + sec.y3d(0)**2 + sec.z3d(0)**2)**(0.5) for sec in soma_list]

def saveconc():
    numpy.save(os.path.join(outdir,"k_%i.npy" % int(h.t)), k[ecs].states3d)
    numpy.save(os.path.join(outdir,"na_%i.npy" % int(h.t)), na[ecs].states3d)
    numpy.save(os.path.join(outdir,"cl_%i.npy" % int(h.t)), cl[ecs].states3d)
    numpy.save(os.path.join(outdir,'o2_%i.npy' % int(h.t)), o2ecs.states3d)
    numpy.save(os.path.join(outdir,'alpha_%i.npy' % int(h.t)), alpha_ecs[ecs].states3d)


# def setup_events(): # ditched, doing fixed time step
#     cv = h.CVode()
#     for i in range(60000//50):
#         cv.event(50*i, saveconc)

# fih = h.FInitializeHandler(2, setup_events)

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
print(datetime.datetime.now())
recs = {'v':soma_v, 'ki':soma_ki, 'nai':soma_nai, 'cli':soma_cli,
        't':t,      'ko':soma_ko, 'nao':soma_nao, 'clo':soma_clo,
        'pos':rpos, 'o2':soma_o2, 'vir':soma_vir, 'vor':soma_vor,
        'rad':cell_positions}

# initialize and set the intracellular concentrations
# h.finitialize(-70)

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

        # if (int(h.t) % dumpint == 0) and (h.t - last_dump > 500):
        #     ## interval saving for restoring 
        #     last_dump = h.t
        #     ### save state
        #     runSS()
        #     saveRxd()
            ### save recs 
            # with open(os.path.join(outdir,"recs" + str(numpy.round(h.t)) + ".pkl"),'wb') as fileObj:
            #     pickle.dump(recs,fileObj)
            # ### save membrane potentials
            # soma, pos = [], []
            # for n in rec_neurons:
            #     soma.append(n.somaV)
            #     pos.append([n.x,n.y,n.z])
            # with open(os.path.join(outdir,"membrane_potential_%i_%i.pkl" % (pcid, numpy.round(h.t))),'wb') as pout:
            #     pickle.dump([soma,pos,time],pout)

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

        # Plot concentration and membrane potential vs radius every 10 ms
        # if int(h.t) % 10 == 0 and pcid == 0:
        #     k_rs = []
        #     concs = []
        #     for nd in k.nodes:
        #         if str(nd.region).split('(')[0] == 'Extracellular':
        #             k_rs.append((nd.x3d**2+nd.y3d**2+nd.z3d**2)**0.5)
        #             concs.append(nd.concentration)
        #     fig = pyplot.figure(dpi=200)
        #     pyplot.subplot(211)
        #     pyplot.scatter(k_rs, concs)
        #     pyplot.title(str(int(h.t))+' ms')
        #     pyplot.ylabel('[K+]')
        #     pyplot.ylim(0,70)
        #     pyplot.subplot(212)
        #     # vs = [sec.v for sec in h.allsec()]
        #     vs = [sec.v for sec in soma_list]
        #     pyplot.scatter(cell_positions, vs)
        #     pyplot.ylabel('Membrane Potential (mV)')
        #     pyplot.ylim(-80,40)
        #     pyplot.xlabel('Radius (microns)')
        #     fig.savefig(outdir+'/gifFigs/k_conc_'+str(int(h.t))+'.png')
        #     pyplot.close()
        #     del(fig)

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
    print('Start Init')
    print(datetime.datetime.now())
    h.finitialize(v_initial)
    print('Finished Init')
    print(datetime.datetime.now())

# run the simulation
run(args['tstop'])

# save final sim state 
runSS()
saveRxd()

# plot raster
if pcid == 0:
    raster = {}
    files = os.listdir(outdir)
    mem_files = [file for file in files if (file[:8] == 'centermembrane')]
    for file in mem_files:
        fileObj = open(os.path.join(outdir,file), 'rb')
        data = pickle.load(fileObj)
        fileObj.close()
        for v, pos in zip(data[0],data[1]):
            pks, _ = find_peaks(v.as_numpy(), 0)
            if len(pks):
                r = (pos[0]**2 + pos[1]**2 + pos[2]**2)**(0.5)
                raster[r] = [data[2][ind] for ind in pks]
    raster_fig = pyplot.figure(figsize=(16,8))
    for key in raster.keys():
        pyplot.plot(numpy.divide(raster[key],1000), [key for i in range(len(raster[key]))], 'b.')
    pyplot.xlabel('Time (s)', fontsize=16)
    pyplot.ylabel('Distance (microns)', fontsize=16)
    pyplot.xticks(fontsize=14)
    pyplot.yticks(fontsize=14)
    raster_fig.savefig(os.path.join(outdir,'raster_plot.png'))

pc.barrier()  
h.quit()

# v1.0 - Port of SpatialModel.py with phenomenological implementation of drop in ECS volume fraction during SD 