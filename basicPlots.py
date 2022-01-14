from analysis import traceExamples, allSpeciesMov
import sys 

datadir = sys.argv[-1]

# plot example traces 
cell_inds = [0, 6, 8, 9, 10 ]
traceExamples(datadir, datadir + 'example_traces.png', iss=cell_inds)

# generate movie for all ions and o2 over time 
## movie params 
vmins = [3.5, 100.0, 30.0, 0.01] # [k, cl, na, o2]
vmaxes = [40.0, 130.0, 140.0, 0.10] 
dur = 2 # sim duration in seconds
extent= (-250, 250, -250, 250) # (xmin, xmax, ymin, ymax)
allSpeciesMov(datadir, datadir + '/movFigs/', vmins, vmaxes, 
                datadir + 'species_movie.mp4', dur=dur, extent=extent)

# v1.0 - generates basic plots for sim specified by its output data dir as cmd line arg