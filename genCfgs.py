# uses command line arguments to generate sim config files for 
# SpatialModelRealistic.py

import json 
import argparse
import os 

# try:
parser = argparse.ArgumentParser(description = '''Run the spreading
                                depression simulation''')
parser.add_argument('--tstop', nargs='?', type=float, default=60,
                    help='''duration of the simulation in ms (defaults
                    to 60ms)''')
parser.add_argument('--dir', type=str, default='Data/default_sim_cfg/',
                    help='a directory to save the figures and data')
parser.add_argument('--ox', nargs='?', type=str, default='normoxic',
                    help='Either normoxic or anoxic, defaults to normoxic')
parser.add_argument('--dendL', nargs='?', type=float, default=0,
                    help='length of dendritic compartment, defaults to no dendrite')
parser.add_argument('--BC', nargs='?', type=str, default='invitro',
                    help='boundary conditions, options are invitro (default) or invivo')
parser.add_argument('--k0', nargs='?', type=float, default=40, 
                    help='initial K+ conc in ecs in core, defaults to 40')
parser.add_argument('--r0', nargs='?', type=float, default=50,
                    help='radius of the core, defaults to 50')
parser.add_argument('--o2factor', nargs='?', type=float, default=1.0,
                    help='factor to reduce o2 in core, defaults to 1')
parser.add_argument('--nthreads', nargs='?', type=int, default=1, 
                    help='number of rxd threads, defaults to 1')
parser.add_argument('--nrec', nargs='?', type=int, default=100,
                    help='number of recorded cells, defaults to 100')
parser.add_argument('--size', nargs='?', type=str, default='mm', 
                    help='size of Lx and Ly, defualts to large')
parser.add_argument('--infuse', nargs='?', type=str, default=False,
                    help='infuse K: yes')
parser.add_argument('--nrnPumpFactor', nargs='?', type=float, default=1)
parser.add_argument('--glialPumpFactor', nargs='?', type=float, default=1)
parser.add_argument('--density', nargs='?', type=int, default=90000)
parser.add_argument('--sa2v', nargs='?', type=float, default=None)
parser.add_argument('--restoredir', nargs='?', type=str, default=None)
parser.add_argument('--p_max', nargs='?', type=float, default=0.8)
parser.add_argument('--pparam', nargs='?', type=float, default=20.0)
parser.add_argument('--O2consume', nargs='?', type=str, default=False)
parser.add_argument('--core', nargs='?', type=str, default=False)
parser.add_argument('--pas', nargs='?', type=float, default=None)
parser.add_argument('--gpas', nargs='?', type=float, default=0.0001)
parser.add_argument('--Lz', nargs='?', type=float, default=400.0)
parser.add_argument('--uniformRec', nargs='?', type=str, default=False)
parser.add_argument('--randSeed', nargs='?', type=int, default=6324555)
parser.add_argument('--ischemCore',nargs='?', type=str, default=False)
parser.add_argument('--ischemEdemaCore', nargs='?', type=str, default=False)
parser.add_argument('--edemaCore', nargs='?', type=str, default=False)
parser.add_argument('--ouabain', nargs='?', type=str, default=False)
parser.add_argument('--alphaNrn', nargs='?', type=float, default=0.24)
parser.add_argument('--alphaECS', nargs='?', type=float, default=None)
parser.add_argument('--lambdaECS', nargs='?', type=float, default=None)
parser.add_argument('--varO2', nargs='?', type=float, default=None)
parser.add_argument('--varCl', nargs='?', type=float, default=None)
parser.add_argument('filename', metavar='filename', type=str)
args = parser.parse_args()
# except:
#     os._exit(1)

out = {'tstop' : args.tstop, 
         'dir' : args.dir, 
         'ox' : args.ox, 
         'dendL' : args.dendL, 
         'BC' : args.BC, 
         'k0' : args.k0, 
         'r0' : args.r0,  
         'o2factor' : args.o2factor, 
         'nthreads' : args.nthreads, 
         'nrec' : args.nrec, 
         'size' : args.size, 
         'infuse' : False, 
         'nrnPumpFactor' : args.nrnPumpFactor, 
         'glialPumpFactor' : args.glialPumpFactor, 
         'density' : args.density, 
         'sa2v' : args.sa2v, 
         'p_max' : args.p_max,
         'pparam' : args.pparam,
         'O2consume' : args.O2consume,
         'core' : args.core,
         'pas' : args.pas,
         'gpas' : args.gpas,
         'Lz' : args.Lz,
         'uniformRec' : args.uniformRec,
         'restoredir' : args.restoredir,
         'randSeed' : args.randSeed,
         'ouabain' : args.ouabain,
         'ischemCore' : args.ischemCore,
         'ischemEdemaCore' : args.ischemEdemaCore,
         'edemaCore' : args.edemaCore,
         'alphaNrn' : args.alphaNrn,
         'alphaECS' : args.alphaECS,
         'lambdaECS' : args.lambdaECS,
         'varO2' : args.varO2,
         'varCl' : args.varCl}

with open(args.filename, 'w') as fileObj:
    json.dump(out, fileObj)

# v1.0 - script for generating sim configs for SpatialModel.py and SpatialModelDynAlpha.py