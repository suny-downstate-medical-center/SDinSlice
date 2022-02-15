import numpy as np
from matplotlib import pyplot as plt 
import os
from scipy.signal import find_peaks 
import pickle
import imageio

def rasterPlot(datadir, center = [125, -450, 125], uniform=True, figname='raster.png', orderBy='y', position='center'):
    files = os.listdir(datadir)
    mem_files = [file for file in files if (file.startswith(position + 'membrane'))]
    raster = {}
    for file in mem_files:
        with open(os.path.join(datadir,file), 'rb') as fileObj:
            data = pickle.load(fileObj)
        for v, pos, pop in zip(data[0], data[1], data[2]):
            pks, _ = find_peaks(v.as_numpy(), 0)
            if len(pks):
                if uniform:
                    r = ((pos[0]-center[0])**2 + (pos[1]-center[1])**2 + (pos[2]-center[2])**2)**(0.5)
                else:
                    r = (pos[0]**2 + pos[1]**2)**(0.5)
                if orderBy == 'y':
                    raster[pos[1]] = {'t': [data[3][ind] for ind in pks],
                        'pop' : pop}
                elif orderBy == 'x':
                    raster[pos[0]] = {'t': [data[3][ind] for ind in pks],
                        'pop' : pop}
                elif orderBy == 'z':
                    raster[pos[3]] = {'t': [data[3][ind] for ind in pks],
                        'pop' : pop}
                elif orderBy == 'r':
                    raster[r] = {'t': [data[3][ind] for ind in pks],
                        'pop' : pop}
                else:
                    print("Invalid oderBy. Using y")
                    raster[pos[1]] = {'t': [data[3][ind] for ind in pks],
                        'pop' : pop}
    pops = np.array(['E2', 'I2', 'E4', 'I4', 'E5', 'I5'])
    cols = ['blue', 'red', 'yellow', 'purple', 'green', 'black']
    for key in raster.keys():
        c = cols[np.argwhere(pops==raster[key]['pop'])[0][0]]
        plt.plot(np.divide(raster[key]['t'],1000), [key for i in range(len(raster[key]['t']))], '.', color=c)
    plt.savefig(figname)

def allSpeciesMov(datadir, outpath, vmins, vmaxes, figname, condition='Perfused', dur=10, extent=None, includeSpks=False):
    """"Generates an mp4 video with heatmaps for K+, Cl-, Na+, and O2 overlaid with spiking data"""
    try:
        os.mkdir(outpath)
    except:
        pass
    plt.ioff()
    specs = ['k', 'cl', 'na', 'o2']
    k_files = [specs[0]+'_'+str(i)+'.npy' for i in range(int(dur*1000)) if (i%100)==0]    
    cl_files = [specs[1]+'_'+str(i)+'.npy' for i in range(int(dur*1000)) if (i%100)==0]    
    na_files = [specs[2]+'_'+str(i)+'.npy' for i in range(int(dur*1000)) if (i%100)==0]    
    o2_files = [specs[3]+'_'+str(i)+'.npy' for i in range(int(dur*1000)) if (i%100)==0]    
    for k_file, cl_file, na_file, o2_file in zip(k_files, cl_files, na_files, o2_files):
        t = int(k_file.split('.')[0].split('_')[1])
        ttl = 't = ' + str(float(k_file.split('.')[0].split('_')[1]) / 1000) + ' s'
        fig = plt.figure(figsize=(12,7.6))
        fig.text(.45, 0.825, ttl, fontsize=20)
        fig.text(0.45, 0.9, condition, fontsize=20)
        if includeSpks:
            posBySpkTime = xyOfSpikeTime(datadir)
            spkTimes = [key for key in posBySpkTime if (t-50 < key <= t+50)]
        ## K+ plot
        ax1 = fig.add_subplot(141)
        data = np.load(datadir+k_file)
        im = plt.imshow(np.transpose(data.mean(2)), vmin=vmins[0], vmax=vmaxes[0], interpolation='nearest', origin='lower', extent=extent)
        plt.colorbar(im, fraction=0.046, pad=0.04)
        if includeSpks:
            if len(spkTimes):
                for spkTime in spkTimes:
                    plt.plot(posBySpkTime[spkTime]['x'], posBySpkTime[spkTime]['y'], 'w*')
        plt.title(r'[K$^+$]$_{ECS}$ ', fontsize=20)
        ## Cl plot
        ax2 = fig.add_subplot(142)
        data = np.load(datadir+cl_file)
        im = plt.imshow(np.transpose(data.mean(2)), vmin=vmins[1], vmax=vmaxes[1], interpolation='nearest', origin='lower', extent=extent)
        plt.colorbar(im, fraction=0.046, pad=0.04)
        if includeSpks:
            if len(spkTimes):
                for spkTime in spkTimes:
                    plt.plot(posBySpkTime[spkTime]['x'], posBySpkTime[spkTime]['y'], 'w*')
        plt.title(r'[Cl$^-$]$_{ECS}$ ', fontsize=20)
        ## Na plot
        ax3 = fig.add_subplot(143)
        data = np.load(datadir+na_file)
        im = plt.imshow(np.transpose(data.mean(2)), vmin=vmins[2], vmax=vmaxes[2], interpolation='nearest', origin='lower', extent=extent)
        plt.colorbar(im, fraction=0.046, pad=0.04)
        if includeSpks:
            if len(spkTimes):
                for spkTime in spkTimes:
                    plt.plot(posBySpkTime[spkTime]['x'], posBySpkTime[spkTime]['y'], 'w*')
        plt.title(r'[Na$^+$]$_{ECS}$ ', fontsize=20)
        ## O2 plot 
        ax4 = fig.add_subplot(144)
        data = np.load(datadir+o2_file)
        im = plt.imshow(np.transpose(data.mean(2)), vmin=vmins[3], vmax=vmaxes[3], interpolation='nearest', origin='lower', extent=extent)
        plt.colorbar(im, fraction=0.046, pad=0.04)
        if includeSpks:
            if len(spkTimes):
                for spkTime in spkTimes:
                    plt.plot(posBySpkTime[spkTime]['x'], posBySpkTime[spkTime]['y'], 'w*')
        plt.title(r'[O$_2$]$_{ECS}$ ', fontsize=20)
        plt.tight_layout()
        
        # plt.tight_layout()
        fig.savefig(outpath + k_file[:-4] + '.png')
        plt.close()
    times = []
    filenames = os.listdir(path=outpath)
    for file in filenames: 
        times.append(float(file.split('_')[-1].split('.')[0])) 
    inds = np.argsort(times)
    filenames_sort = [filenames[i] for i in inds]
    imagesc = []
    for filename in filenames_sort: 
        imagesc.append(imageio.imread(outpath+filename)) 
    imageio.mimsave(figname, imagesc)

def xyOfSpikeTime(datadir):
    pass