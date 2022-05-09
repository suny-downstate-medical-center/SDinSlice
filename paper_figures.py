from matplotlib import pyplot as plt
from analysis import *
from sklearn import linear_model 
from sklearn.metrics import r2_score
from matplotlib_scalebar.scalebar import ScaleBar
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits import axisartist
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes

def centerVsPeriphThick():
    # Figure 7A
    fig = plt.figure(figsize=(12, 5))
    ax = plt.subplot(111)
    datadirs = ['Data/perfuse_lz100/',
                'Data/perfuse_lz200/',
                'Data/perfuse_lz300/',
                'Data/perfuse_standard/',
                'Data/perfuse_lz500/',
                'Data/perfuse_lz600/',
                'Data/SD_Data/perfuse_lz700/',
                'Data/perfuse_lz800/']
    rmaxs = [300, 300, 400, 500, 600, 600, 600, 600]
    speeds_core_perfuse = []
    speeds_periph_perfuse = []
    for datadir, rmax in zip(datadirs, rmaxs):
        core, periph = centerVsPeriphKspeed(datadir, 10, rmax=rmax)
        speeds_core_perfuse.append(core)
        speeds_periph_perfuse.append(periph)
    speeds_periph_perfuse[1] = 0
    datadirs = ['Data/SD_Data/hypox_lz100/',
                'Data/SD_Data/hypox_lz200/',
                'Data/SD_Data/hypox_lz300/',
                'Data/hypox_standard/',
                'Data/SD_Data/hypox_lz500/',
                'Data/SD_Data/hypox_lz600/',
                'Data/SD_Data/hypox_lz700/',
                'Data/hypox_lz800/']
    speeds_core_hypox = []
    speeds_periph_hypox = []
    for datadir in datadirs:
        core, periph = centerVsPeriphKspeed(datadir, 10)
        speeds_core_hypox.append(core)
        speeds_periph_hypox.append(periph)
    thick = [100, 200, 300, 400, 500, 600, 700, 800]
    ax.plot(thick, speeds_core_perfuse, 'b*-', label='Perfused - Core', linewidth=4, markersize=8)
    ax.plot(thick, speeds_periph_perfuse, 'b*--', label='Perfused - Periphery', linewidth=4, markersize=8)
    ax.plot(thick, speeds_core_hypox, 'r*-', label='Hypoxic - Core', linewidth=4, markersize=8)
    ax.plot(thick, speeds_periph_hypox, 'r*--', label='Hypoxic - Periphery', linewidth=4, markersize=8)
    ax.set_xlabel(r'Slice Thickness ($\mu$m)', fontsize=16)
    ax.set_ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=16) 
    ax.set_ylim(-0.1, 11)
    plt.setp(ax.get_xticklabels(), fontsize=14)
    plt.setp(ax.get_yticklabels(), fontsize=14)
    leg = plt.legend(fontsize=14, title='Slice Condition - Depth')
    plt.setp(leg.get_title(), fontsize=16)
    ax.set_title('Influence of Slice Thickness & Depth', fontsize=20, fontweight='bold')
    ax.text(-0.075, 1.1, 'A)', transform=ax.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')

def dkdt():
    # datadir = 'Data/SD_Data/perfuse_standard_highNrec/'
    datadir = 'Data/SD_Data/hypox_standard_highNrec/'
    filename = 'recs.pkl'
    with open(datadir+filename, 'rb') as fileObj:
        data = pickle.load(fileObj)
    posish = []
    for i in range(len(data['v'])):
        pos = (data['pos'][i][0] ** 2 + data['pos'][i][1] ** 2 + data['pos'][i][2] ** 2)**(0.5)
        posish.append(pos)
    posish = np.divide(posish, np.max(posish))
    cols = plt.cm.jet(posish)
    fig, axs = plt.subplots(3, 1, sharex=True)
    fig.set_figheight(9)
    fig.set_figwidth(18)
    for i, c in enumerate(cols):
        dkidt = np.multiply(np.diff(data['ki'][i]),1000)
        dnaidt = np.multiply(np.diff(data['nai'][i]),1000)
        dclidt = np.multiply(np.diff(data['cli'][i]),1000)
        # v = data['v'][i]
        time = np.divide(list(data['t'])[1:],1000)
        ind = 0
        has_spiked = False
        dkidt_cut, dnaidt_cut, dclidt_cut = [], [], []
        time_cut = []
        while not has_spiked and ind < len(dkidt):
            if dkidt[ind] > -8:
                dkidt_cut.append(dkidt[ind])
                dnaidt_cut.append(dnaidt[ind])
                dclidt_cut.append(dclidt[ind])
                time_cut.append(time[ind])
                ind = ind + 1
            else:
                has_spiked = True
        axs[0].plot(time_cut, dkidt_cut, color=c)
        axs[1].plot(time_cut, dnaidt_cut, color=c)
        axs[2].plot(time_cut, dclidt_cut, color=c)
    axs[0].set_ylabel(r'$\frac{d[K^+]_i}{dt}$ (mM/s)', fontsize=16)
    axs[1].set_ylabel(r'$\frac{d[Na^+]_i}{dt}$ (mM/s)', fontsize=16)
    axs[2].set_ylabel(r'$\frac{d[Cl^-]_i}{dt}$ (mM/s)', fontsize=16)
    axs[2].set_xlabel('Time (s)', fontsize=16)
    plt.subplot(311)
    plt.ylim(-0.05, 0.2)
    plt.plot([0,10],[0,0],'k--')
    plt.xlim(-0.05, 10)
    sm = plt.cm.ScalarMappable(cmap=plt.cm.jet)
    cbar = plt.colorbar(sm)


def exampleTraces():
    iss = [0, 7, 10, 17, 25, 30, 35]
    datadir = 'Data/SD_Data/perfuse_standard_highNrec/'
    filename = 'recs.pkl'
    
    if isinstance(datadir, list):
        data = combineRecs(datadir)
    else:
        with open(datadir+filename, 'rb') as fileObj:
            data = pickle.load(fileObj)
    # fig = plt.figure(figsize=(18,9))
    fig, axs = plt.subplots(4, 1, sharex=True)
    fig.set_figheight(9)
    fig.set_figwidth(18)
    # axs3 = plt.subplot(4,3,10)
    # axs4 = plt.subplot(4,3,11)
    # axs5 = plt.subplot(4,3,12)
    for i in iss:
        l = r'%s $\mu$m' % str(np.round((data['pos'][i][0] ** 2 + data['pos'][i][1] ** 2 + data['pos'][i][2] ** 2)**(0.5),1))
        axs[0].plot(np.divide(data['t'],1000), data['v'][i], label=l)
        axs[1].plot(np.divide(data['t'],1000), data['ko'][i])
        axs[2].plot(np.divide(list(data['t'])[1:],1000), np.multiply(np.diff(data['ki'][i]),1000))
        axs[3].plot(np.divide(data['t'],1000), data['o2'][i])
        # axs[1].plot(np.divide(data['t'],1000), data['ki'][i])
        # axs[2].plot(np.divide(data['t'],1000), data['ko'][i])
        # axs4.plot(np.divide(data['t'],1000), data['nai'][i])
        # axs5.plot(np.divide(data['t'],1000), data['cli'][i])
    plt.subplot(413)
    plt.ylim(-0.1, 0.2)
    plt.plot([0,10],[0,0],'k--')
    plt.xlim(-0.05, 3.0)

    leg = axs[0].legend(title='Radial Position', fontsize=11, bbox_to_anchor=(-0.15, 1.05))
    plt.setp(leg.get_title(), fontsize=15)
    axs[0].set_ylabel(r'V$_{memb}$ (mV)', fontsize=16)
    plt.setp(axs[0].get_xticklabels(), fontsize=14)
    plt.setp(axs[0].get_yticklabels(), fontsize=14)
    axs[0].text(-0.1, 1.1, 'A)', transform=axs[0].transAxes,
        fontsize=16, fontweight='bold', va='top', ha='right')

    # axs[1][0].set_ylabel(r'Extracellular [O$_{2}$] (mM)', fontsize=16)
    # plt.setp(axs[1][0].get_xticklabels(), fontsize=14)
    # plt.setp(axs[1][0].get_yticklabels(), fontsize=14)
    # axs[1][0].text(-0.1, 1.1, 'E)', transform=axs[1][0].transAxes,
    #     fontsize=16, fontweight='bold', va='top', ha='right')
    
    axs[1].set_ylabel(r'[K$^{+}$]$_{ECS}$ (mM)', fontsize=16)
    plt.setp(axs[1].get_xticklabels(), fontsize=14)
    plt.setp(axs[1].get_yticklabels(), fontsize=14)
    axs[1].text(-0.1, 1.1, 'B)', transform=axs[1].transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')

    axs[2].set_ylabel(r'$\frac{d[K^+]_i}{dt}$ (mM/s)', fontsize=16)
    plt.setp(axs[2].get_xticklabels(), fontsize=14)
    plt.setp(axs[2].get_yticklabels(), fontsize=14)
    axs[2].text(-0.1, 1.1, 'C)', transform=axs[2].transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')

    axs[3].set_ylabel(r'[O$_2$]$_{ECS}$ (mM)', fontsize=16)
    plt.setp(axs[3].get_xticklabels(), fontsize=14)
    plt.setp(axs[3].get_yticklabels(), fontsize=14)
    axs[3].text(-0.1, 1.1, 'D)', transform=axs[3].transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')
    
    # axs[2].set_ylabel(r'Extracellular [K$^{+}$] (mM)', fontsize=16)
    # plt.setp(axs[2].get_xticklabels(), fontsize=14)
    # plt.setp(axs[2].get_yticklabels(), fontsize=14)
    # axs[2]].text(-0.1, 1.1, 'F)', transform=axs[1][1].transAxes,
    #   fontsize=18, fontweight='bold', va='top', ha='right')
    
    # axs[0][2].set_ylabel(r'Intracellular [Na$^{+}$] (mM)', fontsize=16)
    # plt.setp(axs[0][2].get_xticklabels(), fontsize=14)
    # plt.setp(axs[0][2].get_yticklabels(), fontsize=14)
    # axs[0][2].text(-0.1, 1.1, 'C)', transform=axs[0][2].transAxes,
    #   fontsize=18, fontweight='bold', va='top', ha='right')

    # axs[1][2].set_ylabel(r'Extracellular [Na$^{+}$] (mM)', fontsize=16)
    # plt.setp(axs[1][2].get_xticklabels(), fontsize=14)
    # plt.setp(axs[1][2].get_yticklabels(), fontsize=14)
    # axs[1][2].text(-0.1, 1.1, 'G)', transform=axs[1][2].transAxes,
    #   fontsize=18, fontweight='bold', va='top', ha='right')
    
    # axs[0][3].set_ylabel(r'Intracellular [Cl$^{-}$] (mM)', fontsize=16)
    # plt.setp(axs[0][3].get_xticklabels(), fontsize=14)
    # plt.setp(axs[0][3].get_yticklabels(), fontsize=14)
    # axs[0][3].text(-0.1, 1.1, 'D)', transform=axs[0][3].transAxes,
    #   fontsize=18, fontweight='bold', va='top', ha='right')
    
    # axs[1][3].set_ylabel(r'Extracellular [Cl$^{-}$] (mM)', fontsize=16)
    # plt.setp(axs[1][3].get_xticklabels(), fontsize=14)
    # plt.setp(axs[1][3].get_yticklabels(), fontsize=14)
    # axs[1][3].text(-0.1, 1.1, 'H)', transform=axs[1][3].transAxes,
    #   fontsize=18, fontweight='bold', va='top', ha='right')

    fig.text(0.55, 0.01, 'Time (s)', fontsize=16)
    plt.tight_layout()
    # plt.savefig(figname)

def SD():
    datadir = 'Data/SD_Data/perfuse_standard_highNrec/'
    k_files = ['k_0.npy', 'k_2500.npy', 'k_5000.npy', 'k_7500.npy']
    
    fig = plt.figure(figsize=(14,8))
    ax0 = plt.subplot(241)
    with open(datadir + k_files[0], 'rb') as fileObj:
        data = np.load(fileObj)
    im = ax0.imshow(data.mean(2), vmin=3.5, vmax=40, interpolation='nearest', origin='lower', cmap='inferno')#, extent=(-500,500,-500,500))
    pos = ax0.get_position()
    pos.x0 = pos.x0 - 0.1
    ax0.set_position(pos)
    ax0.set_xticks([])
    ax0.set_yticks([])
    ax0.set_title('t = 0.0 s', fontsize=18)
    scalebar = ScaleBar(25, 'um')
    scalebar.location = 'lower left'
    ax0.add_artist(scalebar)
    ax1 = plt.subplot(242)
    with open(datadir + k_files[1], 'rb') as fileObj:
        data = np.load(fileObj)
    im = ax1.imshow(data.mean(2), vmin=3.5, vmax=40, interpolation='nearest', origin='lower', cmap='inferno')#, extent=(-500,500,-500,500))
    pos = ax1.get_position()
    pos.x0 = pos.x0 - 0.1
    ax1.set_position(pos)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_title('t = 2.5 s', fontsize=18)
    ax2 = plt.subplot(243)
    with open(datadir + k_files[2], 'rb') as fileObj:
        data = np.load(fileObj)
    im = ax2.imshow(data.mean(2), vmin=3.5, vmax=40, interpolation='nearest', origin='lower', cmap='inferno')#, extent=(-500,500,-500,500))    
    pos = ax2.get_position()
    pos.x0 = pos.x0 - 0.1
    ax2.set_position(pos)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title('t = 5.0 s', fontsize=18)
    ax3 = plt.subplot(244)
    with open(datadir + k_files[3], 'rb') as fileObj:
        data = np.load(fileObj)
    im = ax3.imshow(data.mean(2), vmin=3.5, vmax=40, interpolation='nearest', origin='lower', cmap='inferno')#, extent=(-500,500,-500,500))
    pos = ax3.get_position()
    pos.x0 = pos.x0 - 0.1
    ax3.set_position(pos)
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.set_title('t = 7.5 s', fontsize=18)
    cbar_ax = fig.add_axes([0.87, 0.56, 0.035, 0.29])
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=[3.5, 40])
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(14)
    cbar.set_label(r'[K$^+$] (mM)', fontsize=16, rotation=270)
    ax4 = plt.subplot(212)
    datadir = 'Data/SD_Data/perfuse_nonuniform_highNrec_v4/'
    raster = getRaster(datadir)
    for key in raster.keys():
        ax4.plot(np.divide(raster[key],1000), [key for i in range(len(raster[key]))], 'k.', markersize=1)
    ax4.set_xlabel('Time (s)', fontsize=16)
    ax4.set_ylabel(r'Cell Radial Distance ($\mu$m)', fontsize=16)
    plt.setp(ax4.get_xticklabels(), fontsize=14)
    plt.setp(ax4.get_yticklabels(), fontsize=14)
    ax4.set_title('Spreading Depolarization Raster', fontsize=18)
    fig.text(0.415, 0.9, r'Spreading K$^+$ Wave', fontsize=18)
    ax0.text(-0.075, 1.1, 'A)', transform=ax0.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')
    ax4.text(-0.075, 1.1, 'B)', transform=ax4.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')

def sliceDepthKwave():
    # Figure 7B
    datadir = 'Data/SD_Data/perfuse_standard_highNrec/'
    k_files = ['k_0.npy', 'k_2500.npy', 'k_5000.npy', 'k_7500.npy']
    times = ['0', '2.5', '5.0', '7.5']
    fig, axs = plt.subplots(2,len(k_files))
    fig.set_figwidth(14)
    fig.set_figheight(8)
    for ind, k_file in enumerate(k_files):
        with open(datadir + k_file, 'rb') as fileObj:
            data = np.load(fileObj)
        center = int(data.shape[2]/2)
        data_core = data[:,:,center-1:center+2]
        data_periph = np.concatenate((data[:,:,:2], data[:,:,-3:]), axis=2)
        im0 = axs[0][ind].imshow(data_periph.mean(2), vmin=3.5, vmax=40, interpolation='nearest', origin='lower', cmap='inferno')
        pos = axs[0][ind].get_position()
        pos.x0 = pos.x0 - 0.1
        axs[0][ind].set_position(pos)
        axs[0][ind].set_xticks([])
        axs[0][ind].set_yticks([])
        axs[0][ind].set_title('t = ' + times[ind] + ' s', fontsize=18)
        im1 = axs[1][ind].imshow(data_core.mean(2), vmin=3.5, vmax=40, interpolation='nearest', origin='lower', cmap='inferno')
        pos = axs[1][ind].get_position()
        pos.x0 = pos.x0 - 0.1
        axs[1][ind].set_position(pos)
        axs[1][ind].set_xticks([])
        axs[1][ind].set_yticks([])
        axs[1][ind].set_title('t = ' + times[ind] + ' s', fontsize=18)
    cbar_ax = fig.add_axes([0.87, 0.13, 0.035, 0.72])
    cbar = fig.colorbar(im1, cax=cbar_ax, ticks=[3.5, 40])
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(14)
    cbar.set_label(r'[K$^+$] (mM)', fontsize=16, rotation=270)
    fig.text(0.415, 0.91, 'Periphery', fontsize=20, fontweight='bold')
    fig.text(0.435, 0.48, 'Core', fontsize=20, fontweight='bold')
    axs[0][0].text(-0.075, 1.1, 'B)', transform=axs[0][0].transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')
    axs[1][0].text(-0.075, 1.1, 'C)', transform=axs[1][0].transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')
    scalebar = ScaleBar(25, 'um')
    scalebar.location = 'lower left'
    axs[0][0].add_artist(scalebar)
    scalebar2 = ScaleBar(25, 'um')
    scalebar2.location = 'lower left'
    axs[1][0].add_artist(scalebar2)

def dynamicAlpha():
    fig, axs = plt.subplots(1,5, sharex=True, sharey=True)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    for ind, t in enumerate(['0', '1000', '2000', '3000', '4000']):
        data = np.load('Data/dyn_alpha_5s/alpha_'+t+'.npy')
        im = axs[ind].imshow(data.mean(2), vmin=0, vmax=0.2, interpolation='nearest', origin='lower', cmap='inferno')
        axs[ind].set_title(t[0]+' s', fontsize=18)
        axs[ind].set_xticks([])
        axs[ind].set_yticks([])
        pos = axs[ind].get_position()
        pos.x0 = pos.x0 - 0.1
        axs[ind].set_position(pos)
    cbar_ax = fig.add_axes([pos.x1, pos.y0, 0.03, pos.y1-pos.y0])
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=[0, 0.2])
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(14)
    cbar.set_label(r'$\alpha_{ECS}$', fontsize=16, rotation=270)
    scalebar = ScaleBar(25, 'um')
    scalebar.location = 'lower left'
    axs[0].add_artist(scalebar)
    axs[0].text(-0.1, 1.1, 'B)', fontsize=16, fontweight='bold', transform=axs[0].transAxes,
                va='top', ha='right')


def sliceThickness():
    datadirs = ['Data/perfuse_lz100/',
                'Data/perfuse_lz200/',
                'Data/perfuse_lz300/',
                'Data/perfuse_standard/',
                'Data/perfuse_lz500/',
                'Data/perfuse_lz600/',
                'Data/SD_Data/perfuse_lz700/',
                'Data/perfuse_lz800/']
    speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=620)
    speeds_normox[0] = 0
    labels_normox = [r'100 $\mu$m', r'200 $\mu$m', r'300 $\mu$m', r'400 $\mu$m', r'500 $\mu$m', r'600 $\mu$m', r'700 $\mu$m', r'800 $\mu$m']
    thick_normox = [100, 200, 300, 400, 500, 600, 700, 800]
    legendTitle = 'Slice Thickness'
    compareKwaves(datadirs, labels_normox, legendTitle, sbplt=131)
    
    datadirs = ['Data/SD_Data/hypox_lz100/',
                'Data/SD_Data/hypox_lz200/',
                'Data/SD_Data/hypox_lz300/',
                'Data/hypox_standard/',
                'Data/SD_Data/hypox_lz500/',
                'Data/SD_Data/hypox_lz600/',
                'Data/SD_Data/hypox_lz700/',
                'Data/hypox_lz800/']
    speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=620)
    labels_anox = [r'100 $\mu$m', r'200 $\mu$m', r'300 $\mu$m', r'400 $\mu$m', r'500 $\mu$m', r'500 $\mu$m', r'600 $\mu$m', r'700 $\mu$m', r'800 $\mu$m']
    thick_anox = [100, 200, 300, 400, 500, 600, 700, 800]
    compareKwaves(datadirs, labels_anox, legendTitle, sbplt=132)
    plt.subplot(133)
    plt.plot(thick_normox, speeds_normox, '*-', label='Perfused')
    plt.plot(thick_anox, speeds_anox, '*-', label='Hypoxic')
    legend = plt.legend(title='Tissue Oxygenation', fontsize=12)
    plt.setp(legend.get_title(), fontsize=14)
    plt.xlabel(r'Slice Thickness ($\mu$m$^{-2}$)', fontsize=14)
    plt.ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.subplot(131)
    plt.xlim(0,10)
    plt.ylim(-10,710)
    plt.title('Perfused', fontsize=18)
    plt.subplot(132)
    plt.xlim(0,10)
    plt.ylim(-10,710)
    plt.title('Hypoxic', fontsize=18)
    plt.subplot(133)
    plt.title('SD Wave Speed', fontsize=18)
    fig = plt.gca()
    fig.figure.set_figheight(10)
    fig.figure.set_figwidth(14)
    plt.tight_layout()

def perfusePrimedHypoxic():
    datadirs = ['Data/perfuse_standard/',
                'Data/SD_Data/primed_standard/',
                'Data/hypox_standard/']
    speeds = getKwaveSpeed(datadirs, r0=100, tcut=10,rmax=600)
    labels = [cond + r': %.2f mm/min' % (s) for cond, s in zip(['Perfused','Primed','Hypoxic'],speeds)]
    legendTitle = 'Slice Condition'
    fig = plt.figure()
    compareKwaves(datadirs, labels, legendTitle, sbplt=211)
    plt.title(r'[K$^+$]$_{bolus}$ = 70 mM', fontsize=18)
    datadirs = ['Data/SD_Data/perfuse_ouabain/',
                'Data/SD_Data/primed_ouabain/',
                'Data/SD_Data/hypox_ouabain/']
    speeds = getKwaveSpeed(datadirs, r0=100, tcut=10,rmax=600)
    labels = [cond + r': %.2f mm/min' % (s) for cond, s in zip(['Perfused','Primed','Hypoxic'],speeds)]
    compareKwaves(datadirs, labels, legendTitle, sbplt=212)
    fig = plt.gca()
    fig.figure.set_figheight(9)
    fig.figure.set_figwidth(7)
    plt.title(r'Ouabain', fontsize=18)
    plt.tight_layout()

def surfaceToVolume():
    datadirs = ['Data/SD_Data/hypox_sv02/', 
                'Data/SD_Data/hypox_sv1/',
                'Data/SD_Data/hypox_sv2/',  
                'Data/hypox_standard/',
                'Data/SD_Data/hypox_sv4/', 
                'Data/SD_Data/hypox_sv5/',
                'Data/SD_Data/hypox_sv6/']
    speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=600)
    labels_anox = [r'0.02 $\mu$m$^{-1}$', r'1.0 $\mu$m$^{-1}$', r'2.0 $\mu$m$^{-1}$', r'3.0 $\mu$m$^{-1}$',
            r'4.0 $\mu$m$^{-1}$', r'5.0 $\mu$m$^{-1}$',  r'6.0 $\mu$m$^{-1}$']
    ratios_anox = [0.02, 1, 2, 3, 4, 5, 6]
    # ratios_anox = [0.02, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0]
    legendTitle = 'S:V'
    compareKwaves(datadirs, labels_anox, legendTitle, sbplt=132)
    datadirs = ['Data/SD_Data/perufse_sv02/', 
                'Data/SD_Data/perufse_sv1/',
                'Data/SD_Data/perufse_sv2/',  
                'Data/perfuse_standard/',
                'Data/SD_Data/perufse_sv4/', 
                'Data/SD_Data/perufse_sv5/',
                'Data/SD_Data/perufse_sv6/']
    ratios_normox = ratios_anox
    labels_normox = [r'0.02 $\mu$m$^{-1}$', r'1.0 $\mu$m$^{-1}$', r'2.0 $\mu$m$^{-1}$', 
                    r'3.0 $\mu$m$^{-1}$', r'4.0 $\mu$m$^{-1}$', r'5.0 $\mu$m$^{-1}$',
                    r'6.0 $\mu$m$^{-1}$']
    speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=400)
    speeds_normox[0] = 0
    speeds_normox[1] = 0
    compareKwaves(datadirs, labels_normox, legendTitle, sbplt=131)
    plt.subplot(133)
    plt.plot(ratios_normox, speeds_normox, '*-', label='Perfused')
    plt.plot(ratios_anox, speeds_anox, '*-', label='Hypoxic')
    legend = plt.legend(title='Tissue Oxygenation', fontsize=12)
    plt.setp(legend.get_title(), fontsize=14)
    plt.xlabel(r'S:V ($\mu$m$^{-1}$)', fontsize=14)
    plt.ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.subplot(131)
    plt.xlim(0,10)
    plt.ylim(-10,650)
    plt.title('Perfused', fontsize=18)
    plt.subplot(132)
    plt.xlim(0,10)
    plt.ylim(-10,650)
    plt.title('Hypoxic', fontsize=18)
    plt.subplot(133)
    plt.title('SD Wave Speed', fontsize=18)
    fig = plt.gca()
    fig.figure.set_figheight(7.6)
    fig.figure.set_figwidth(14)
    plt.tight_layout()

def cellDensity():    
    datadirs = ['Data/SD_Data/hypox_d45/',
                'Data/SD_Data/hypox_d675/',  
                'Data/hypox_standard/',
                'Data/SD_Data/hypox_d1125/', 
                'Data/SD_Data/hypox_d120/']
    speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=600)
    labels_anox = ['45k', '67.5k', '90k', '112.5k', '120k']
    densities_anox = [45, 67.5, 90, 112.5, 120]
    legendTitle = r'Neuron Density (cells/mm$^3$)'
    compareKwaves(datadirs, labels_anox, legendTitle, sbplt=132)
    datadirs = ['Data/SD_Data/perfuse_d45/',
                'Data/SD_Data/perfuse_d675/',  
                'Data/perfuse_standard/',
                'Data/SD_Data/perfuse_d1125/', 
                'Data/SD_Data/perfuse_d120/']
    speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=400)
    densities_normox = [45, 67.5, 90, 112.5, 120]
    compareKwaves(datadirs, labels_anox, legendTitle, sbplt=131)
    plt.subplot(133)
    plt.plot(densities_normox, speeds_normox, '*-', label='Perfused')
    plt.plot(densities_anox, speeds_anox, '*-', label='Hypoxic')
    legend = plt.legend(title='Tissue Oxygenation', fontsize=12)
    plt.setp(legend.get_title(), fontsize=14)
    plt.xlabel(r'Neuron Density (x 1000 cells/mm$^3$)', fontsize=14)
    plt.ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.subplot(131)
    plt.xlim(0,10)
    plt.ylim(-10,650)
    plt.title('Perfused', fontsize=18)
    plt.subplot(132)
    plt.xlim(0,10)
    plt.ylim(-10,650)
    plt.title('Hypoxic', fontsize=18)
    plt.subplot(133)
    plt.title('SD Wave Speed', fontsize=18)
    fig = plt.gca()
    fig.figure.set_figheight(7.6)
    fig.figure.set_figwidth(14)
    plt.tight_layout()

def cellDensityCompareThickness():    
    datadirs = ['Data/SD_Data/hypox_d45/',
                'Data/SD_Data/hypox_d675/',  
                'Data/hypox_standard/',
                'Data/SD_Data/hypox_d1125/', 
                'Data/SD_Data/hypox_d120/']
    speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=600)
    datadirs = ['Data/SD_Data/perfuse_d45/',
                'Data/SD_Data/perfuse_d675/',  
                'Data/perfuse_standard/',
                'Data/SD_Data/perfuse_d1125/', 
                'Data/SD_Data/perfuse_d120/']
    speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=400)
    datadirs = ['Data/SD_Data/hypox_d45_lz200/',
                'Data/SD_Data/hypox_d675_lz200/',  
                'Data/SD_Data/hypox_lz200/',
                'Data/SD_Data/hypox_d1125_lz200/', 
                'Data/SD_Data/hypox_d120_lz200/']
    speeds_anox_lz200 = getKwaveSpeed(datadirs, r0=100, rmax=600)
    datadirs = ['Data/SD_Data/perfuse_d45_lz200/',
                'Data/SD_Data/perfuse_d675_lz200/',  
                'Data/perfuse_lz200/',
                'Data/SD_Data/perfuse_d1125_lz200/', 
                'Data/SD_Data/perfuse_d120_lz200/']
    speeds_normox_lz200 = getKwaveSpeed(datadirs, r0=100, rmax=400)
    densities = [45, 67.5, 90, 112.5, 120]
    plt.subplot(121)
    plt.plot(densities, speeds_normox_lz200, '*-', label=r'200 $\mu$m')
    plt.plot(densities, speeds_normox, '*-', label=r'400 $\mu$m')
    legend = plt.legend(title='Slict Thickness', fontsize=12)
    plt.setp(legend.get_title(), fontsize=14)
    plt.xlabel(r'Neuron Density (x 1000 cells/mm$^3$)', fontsize=14)
    plt.ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title('Perfused')
    plt.subplot(122)
    plt.plot(densities, speeds_anox_lz200, '*-', label=r'200 $\mu$m')
    plt.plot(densities, speeds_anox, '*-', label=r'400 $\mu$m')
    legend = plt.legend(title='Slict Thickness', fontsize=12)
    plt.setp(legend.get_title(), fontsize=14)
    plt.xlabel(r'Neuron Density (x 1000 cells/mm$^3$)', fontsize=14)
    plt.ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title('Hypoxic')
    fig = plt.gca()
    fig.figure.set_figheight(7.6)
    fig.figure.set_figwidth(14)
    plt.tight_layout()

def bolusProps():
    datadirs = [['Data/SD_Data/perfuse_r050/',
                'Data/SD_Data/r050_v0/',
                'Data/SD_Data/r050_v1/',
                'Data/SD_Data/r050_v2/',
                'Data/SD_Data/r050_v3/'],
                ['Data/SD_Data/perfuse_r075/',
                'Data/SD_Data/r075_v0/',
                'Data/SD_Data/r075_v1/',
                'Data/SD_Data/r075_v2/',
                'Data/SD_Data/r075_v3/'],
                ['Data/SD_Data/perfuse_standard_highNrec/',
                'Data/SD_Data/r0100_v0/',
                'Data/SD_Data/r0100_v1/',
                'Data/SD_Data/r0100_v2/',
                'Data/SD_Data/r0100_v3/',
                'Data/SD_Data/k070_v0/',
                'Data/SD_Data/k070_v1/',
                'Data/SD_Data/k070_v2/',
                'Data/SD_Data/k070_v3/'],
                ['Data/SD_Data/perfuse_r0125/',
                'Data/SD_Data/r0125_v0/',
                'Data/SD_Data/r0125_v1/',
                'Data/SD_Data/r0125_v2/',
                'Data/SD_Data/r0125_v3/'],
                ['Data/SD_Data/perfuse_r0150/',
                'Data/SD_Data/r0150_v0/',
                'Data/SD_Data/r0150_v1/',
                'Data/SD_Data/r0150_v2/',
                'Data/SD_Data/r0150_v3/']]
    speeds_avg = []
    speeds_err = []
    for d in datadirs:
        s = getKwaveSpeed(d, rmax=475)
        speeds_avg.append(np.mean(s))
        speeds_err.append(np.std(s))
    rad = [100, 150, 200, 250, 300]
    fig = plt.figure()
    ax1 = fig.add_axes([0.15, 0.1, 0.65, 0.8], axes_class=HostAxes)
    ax2 = ParasiteAxes(ax1, sharey=ax1)
    ax1.parasites.append(ax2)
    ax1.errorbar(rad, speeds_avg, yerr=speeds_err, color='blue', label='Bolus Diameter')
    ax1.set_ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=16)
    ax1.axis['left'].label.set_fontsize(16)
    ax1.set_xlabel(r'Bolus Diameter ($\mu$m)', fontsize=16)
    ax1.axis['bottom'].label.set_fontsize(16)
    fig.text(0.3, 0.97, r'Influence of K$^+$ Bolus Properties', fontsize=18, fontweight='bold')
    ax1.set_ylim(0,3)
    # plt.setp(ax1.get_xticklabels(), fontsize=14)
    # plt.setp(ax1.get_yticklabels(), fontsize=14)
    ax1.axis['bottom'].major_ticklabels.set_fontsize(14)
    ax1.axis['left'].major_ticklabels.set_fontsize(14)
    # ax1.text(-0.15, 1.1, 'A)', transform=ax1.transAxes,
    #         fontsize=18, fontweight='bold', va='top', ha='right')
    datadirs = [['Data/SD_Data/perfuse_k020/',
                'Data/SD_Data/perfuse_k020_v0/',
                'Data/SD_Data/perfuse_k020_v1/',
                'Data/SD_Data/perfuse_k020_v2/',
                'Data/SD_Data/perfuse_k020_v3/'],
                ['Data/SD_Data/perufse_k030/',
                'Data/SD_Data/k030_v0/',
                'Data/SD_Data/k030_v1/',
                'Data/SD_Data/k030_v2/',
                'Data/SD_Data/k030_v3/'],
                ['Data/SD_Data/perufse_k040/',
                'Data/SD_Data/k040_v0/',
                'Data/SD_Data/k040_v1/',
                'Data/SD_Data/k040_v2/',
                'Data/SD_Data/k040_v3/'],
                ['Data/SD_Data/perufse_k050/',
                'Data/SD_Data/k050_v0/',
                'Data/SD_Data/k050_v1/',
                'Data/SD_Data/k050_v2/',
                'Data/SD_Data/k050_v3/'],
                ['Data/SD_Data/perufse_k060/',
                'Data/SD_Data/k060_v0/',
                'Data/SD_Data/k060_v1/',
                'Data/SD_Data/k060_v2/',
                'Data/SD_Data/k060_v3/'],
                ['Data/SD_Data/perfuse_standard_highNrec/',
                'Data/SD_Data/k070_v0/',
                'Data/SD_Data/k070_v1/',
                'Data/SD_Data/k070_v2/',
                'Data/SD_Data/k070_v3/',
                'Data/SD_Data/r0100_v0/',
                'Data/SD_Data/r0100_v1/',
                'Data/SD_Data/r0100_v2/',
                'Data/SD_Data/r0100_v3/',]]
    speeds_avg = []
    speeds_err = []
    for d in datadirs:
        s = getKwaveSpeed(d, rmax=475)
        speeds_avg.append(np.mean(s))
        speeds_err.append(np.std(s))
    ks = [20, 30, 40, 50, 60, 70]
    ax2.axis['top'].set_visible(True)
    ax2.axis['top'].major_ticklabels.set_visible(True)
    ax2.axis['top'].label.set_visible(True)
    ax2.errorbar(ks, speeds_avg, yerr=speeds_err, linestyle='--', color='blue', label=r'Bolus [K$^+$]$_{ECS}$')
    # ax2.set_ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=16)
    ax2.set_xlabel(r'Bolus [K$^+$]$_{ECS}$ (mM)', fontsize=16)
    ax2.axis['top'].label.set_fontsize(16)
    # plt.setp(ax2.get_xticklabels(), fontsize=14)
    # plt.setp(ax2.get_yticklabels(), fontsize=14)
    ax2.axis['top'].major_ticklabels.set_fontsize(14)
    ax2.set_ylim(0,6)
    # ax2.text(-0.15, 1.1, 'B)', transform=ax2.transAxes,
    #         fontsize=18, fontweight='bold', va='top', ha='right')
    plt.legend(fontsize=16)

def bolusRadius():
    datadirs = ['Data/SD_Data/hypox_r050/',
                'Data/SD_Data/hypox_r075/',
                'Data/hypox_standard/',
                'Data/SD_Data/hypox_r0125/',
                'Data/SD_Data/hypox_r0150/']
    speeds_anox = getKwaveSpeed(datadirs, r0=100, tcut=10)
    labels = ['50 $\mu$m', '75 $\mu$m', '100 $\mu$m', '125 $\mu$m', '150 $\mu$m']
    radii = [50, 75, 100, 125, 150]
    legendTitle = 'Bolus Radius'
    compareKwaves(datadirs, labels, legendTitle, sbplt=132)
    datadirs = ['Data/SD_Data/perfuse_r050/',
                'Data/SD_Data/perfuse_r075/',
                'Data/perfuse_standard/',
                'Data/SD_Data/perfuse_r0125/',
                'Data/SD_Data/perfuse_r0150/']
    speeds_normox = getKwaveSpeed(datadirs, r0=100, tcut=10)
    compareKwaves(datadirs, labels, legendTitle, sbplt=131)
    plt.subplot(133)
    plt.plot(radii, speeds_normox, '*-', label='Perfused')
    plt.plot(radii, speeds_anox, '*-', label='Hypoxic')
    plt.xlabel(r'Initial Bolus Radius ($\mu$m)', fontsize=14)
    plt.ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    legend = plt.legend(title='Tissue Oxygenation', fontsize=12)
    plt.setp(legend.get_title(), fontsize=14)
    plt.subplot(131)
    plt.xlim(0,10)
    plt.ylim(-10,710)
    plt.title('Perfused', fontsize=18)
    plt.subplot(132)
    plt.xlim(0,10)
    plt.ylim(-10,710)
    plt.title('Hypoxic', fontsize=18)
    plt.subplot(133)
    plt.title('SD Wave Speed', fontsize=18)
    fig = plt.gca()
    fig.figure.set_figheight(7.6)
    fig.figure.set_figwidth(14)
    plt.tight_layout()

def bolusK():
    datadirs = ['Data/SD_Data/perufse_k030/',
                'Data/SD_Data/perufse_k040/',
                'Data/SD_Data/perufse_k050/',
                'Data/SD_Data/perufse_k060/',
                'Data/perfuse_standard/']
    speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=475)
    labels_normox = ['30 mM', '40 mM', '50 mM', '60 mM', '70 mM']
    densities = [30, 40, 50, 60, 70]
    legendTitle = r'Initial Bolus [K$^+$]'
    compareKwaves(datadirs, labels_normox, legendTitle, sbplt=131)
    datadirs = ['Data/SD_Data/hypox_k030/',
                'Data/SD_Data/hypox_k040/',
                'Data/SD_Data/hypox_k050/',
                'Data/SD_Data/hypox_k060/',
                'Data/hypox_standard/']
    speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=600)
    labels_anox = ['30 mM', '40 mM', '50 mM', '60 mM', '70 mM']
    compareKwaves(datadirs, labels_anox, legendTitle, sbplt=132)
    plt.subplot(133)
    plt.plot(densities, speeds_normox, '*-', label='Perfused')
    plt.plot(densities, speeds_anox, '*-', label='Hypoxic')
    plt.xlabel(r'Initial Bolus [K$^+$] (mM)', fontsize=14)
    plt.ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    legend = plt.legend(title='Tissue Oxygenation', fontsize=12)
    plt.setp(legend.get_title(), fontsize=14)
    plt.subplot(131)
    plt.xlim(0,10)
    plt.ylim(-10,710)
    plt.title('Perfused', fontsize=18)
    plt.subplot(132)
    plt.xlim(0,10)
    plt.ylim(-10,710)
    plt.title('Hypoxic', fontsize=18)
    plt.subplot(133)
    plt.title('SD Wave Speed', fontsize=18)
    fig = plt.gca()
    fig.figure.set_figheight(7.6)
    fig.figure.set_figwidth(14)
    plt.tight_layout()

def alphaNrn():
    datadirs = ['SD_Data/perfuse_nv165/',
                'Data/perfuse_standard/',
                'SD_Data/perfuse_nv315/',
                'SD_Data/perfuse_nv39/',
                'SD_Data/perfuse_nv465/',
                'SD_Data/perfuse_nv54/']
    speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=475)
    labels = ['0.165', '0.240', '0.315', '0.390', '0.465', '0.540']
    legendTitle = r'$\alpha_{nrn}$'
    compareKwaves(datadirs, labels, legendTitle, sbplt=131)
    datadirs = ['SD_Data/hypox_nv165/',
                'Data/hypox_standard/',
                'SD_Data/hypox_nv315/',
                'SD_Data/hypox_nv39/',
                'SD_Data/hypox_nv465/',
                'SD_Data/hypox_nv54/']
    speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=600)
    compareKwaves(datadirs, labels, legendTitle, sbplt=132)
    alphas = [0.165, 0.24, 0.315, 0.39, 0.465, 0.54]
    plt.subplot(133)
    plt.plot(alphas, speeds_normox, '*-', label='Perfused')
    plt.plot(alphas, speeds_anox, '*-', label='Hypoxic')
    plt.xlabel(r'$\alpha_{nrn}$', fontsize=14)
    plt.ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    legend = plt.legend(title='Tissue Oxygenation', fontsize=12)
    plt.setp(legend.get_title(), fontsize=14)
    plt.subplot(131)
    plt.xlim(0,10)
    plt.ylim(-10,710)
    plt.title('Perfused', fontsize=18)
    plt.subplot(132)
    plt.xlim(0,10)
    plt.ylim(-10,710)
    plt.title('Hypoxic', fontsize=18)
    plt.subplot(133)
    plt.title('SD Wave Speed', fontsize=18)
    fig = plt.gca()
    fig.figure.set_figheight(7.6)
    fig.figure.set_figwidth(14)
    plt.tight_layout()

def centerAndThick():
    ## Figure 8
    # datadirs = ['Data/perfuse_lz100/',
    #             'Data/perfuse_lz200/',
    #             'Data/perfuse_lz300/',
    #             'Data/perfuse_standard/',
    #             'Data/perfuse_lz500/',
    #             'Data/perfuse_lz600/',
    #             'Data/SD_Data/perfuse_lz700/',
    #             'Data/perfuse_lz800/']
    # speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=620)
    # speeds_normox[0] = 0
    # thick = [100, 200, 300, 400, 500, 600, 700, 800]
    # legendTitle = 'Slice Thickness'
    # datadirs = ['Data/SD_Data/hypox_lz100/',
    #             'Data/SD_Data/hypox_lz200/',
    #             'Data/SD_Data/hypox_lz300/',
    #             'Data/hypox_standard/',
    #             'Data/SD_Data/hypox_lz500/',
    #             'Data/SD_Data/hypox_lz600/',
    #             'Data/SD_Data/hypox_lz700/',
    #             'Data/hypox_lz800/']
    # speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=620)
    # fig = plt.figure(figsize=(12, 5))
    # ax0 = plt.subplot(111)
    # ax0.plot(thick, speeds_normox, '*-', linewidth=4, markersize=8, color='blue', label='Perfused')
    # ax0.plot(thick, speeds_anox, '*-', linewidth=4, markersize=8, color='red', label='Hypoxic')
    # ax0.set_xlabel(r'Slice Thickness ($\mu$m)', fontsize=16)
    # ax0.set_ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=16)
    # ax0.set_title('Influence of Slice Thickness', fontsize=18)
    # ax0.text(-0.075, 1.1, 'A)', transform=ax0.transAxes,
    #         fontsize=18, fontweight='bold', va='top', ha='right')
    # plt.setp(ax0.get_xticklabels(), fontsize=14)
    # plt.setp(ax0.get_yticklabels(), fontsize=14)
    datadir = 'Data/SD_Data/perfuse_nonuniform_highNrec_v4/'
    spatialBin = 50 
    noverlap = 10 
    poolsize=16 
    rmax = 680 
    duration = 10
    fig = plt.figure(figsize=(12,10))
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)
    hmap1, cbar1 = comboVmembPlot(datadir, duration, fig, ax1, -80, -15, spatialbin=spatialBin, noverlap=noverlap, poolsize=poolsize, rmax=rmax, top=True, left=True)
    hmap2, cbar2 = comboVmembPlot(datadir, duration, fig, ax2, -80, -15, spatialbin=spatialBin, noverlap=noverlap, poolsize=poolsize, rmax=rmax, position='periph', top=True)
    datadir = 'Data/SD_Data/hypox_nonuniform_highNrec_v4/'
    hmap3, cbar3 = comboVmembPlot(datadir, duration, fig, ax3, -80, -15, spatialbin=spatialBin, noverlap=noverlap, poolsize=poolsize, rmax=rmax, left=True)
    hmap4, cbar4 = comboVmembPlot(datadir, duration, fig, ax4, -80, -15, spatialbin=spatialBin, noverlap=noverlap, poolsize=poolsize, rmax=rmax, position='periph')
    ax1.set_ylim(0,600)
    ax1.set_title('Perfused - Core', fontsize=18)
    ax1.text(-0.25, 1.1, 'A)', transform=ax1.transAxes,
            fontsize=18, fontweight='bold', va='top', ha='right')
    ax2.set_ylim(0,600)
    ax2.set_title('Perfused - Periphery', fontsize=18)
    ax2.text(-0.25, 1.1, 'B)', transform=ax2.transAxes,
            fontsize=18, fontweight='bold', va='top', ha='right')
    ax3.set_ylim(0,600)
    ax3.set_title('Hypoxic - Core', fontsize=18)
    ax3.text(-0.25, 1.1, 'C)', transform=ax3.transAxes,
            fontsize=18, fontweight='bold', va='top', ha='right')
    ax4.set_ylim(0,600)
    ax4.set_title('Hypoxic - Periphery', fontsize=18)
    ax4.text(-0.25, 1.1, 'D)', transform=ax4.transAxes,
            fontsize=18, fontweight='bold', va='top', ha='right')
    fig.text(0.3, 0.95, 'Depth Dependent SD Propagation', fontsize=20, fontweight='bold')


    # fig.text(0.35, 0.95, 'Depth-Dependent SD Propagation', fontsize=18)

    # plt.tight_layout()
    

def centerVsPeriphSpkProp():
    datadir = 'Data/SD_Data/perfuse_nonuniform/'
    normox_center = getSpkMetrics(datadir)
    normox_periph = getSpkMetrics(datadir, position='periph')
    datadir = 'Data/SD_Data/hypox_nonuniform/'
    anox_center = getSpkMetrics(datadir)
    anox_periph = getSpkMetrics(datadir, position='periph')
    anox_center_r = [k for k in anox_center.keys()]
    anox_periph_r = [k for k in anox_periph.keys()]
    normox_center_r = [k for k in normox_center.keys()]
    normox_periph_r = [k for k in normox_periph.keys()]
    anox_center_spkFreq = [anox_center[k]['spkFreq'] for k in anox_center.keys()]
    anox_periph_spkFreq = [anox_periph[k]['spkFreq'] for k in anox_periph.keys()]
    normox_center_spkFreq = [normox_center[k]['spkFreq'] for k in normox_center.keys()]
    normox_periph_spkFreq = [normox_periph[k]['spkFreq'] for k in normox_periph.keys()]
    anox_center_spkDur = [anox_center[k]['spkDur'] for k in anox_center.keys()]
    anox_periph_spkDur = [anox_periph[k]['spkDur'] for k in anox_periph.keys()]
    normox_center_spkDur = [normox_center[k]['spkDur'] for k in normox_center.keys()]
    normox_periph_spkDur = [normox_periph[k]['spkDur'] for k in normox_periph.keys()]
    anox_center_nSpks = [anox_center[k]['nSpks'] for k in anox_center.keys()]
    anox_periph_nSpks = [anox_periph[k]['nSpks'] for k in anox_periph.keys()]
    normox_center_nSpks = [normox_center[k]['nSpks'] for k in normox_center.keys()]
    normox_periph_nSpks = [normox_periph[k]['nSpks'] for k in normox_periph.keys()]
    anox_center_t2firstSpk = [anox_center[k]['t2firstSpk'] for k in anox_center.keys()]
    anox_periph_t2firstSpk = [anox_periph[k]['t2firstSpk'] for k in anox_periph.keys()]
    normox_center_t2firstSpk = [normox_center[k]['t2firstSpk'] for k in normox_center.keys()]
    normox_periph_t2firstSpk = [normox_periph[k]['t2firstSpk'] for k in normox_periph.keys()]
    fig, big_axes = plt.subplots( figsize=(14, 7.6) , nrows=1, ncols=1, sharey=True)
    # for row, big_ax in enumerate(zip(big_axes,['Perfused']), start=1):
    big_axes.set_title('Perfused', fontsize=22)
    # Turn off axis lines and ticks of the big subplot 
    # obs alpha is 0 in RGBA string!
    big_axes.tick_params(labelcolor=(1.,1.,1., 0.0), top='off', bottom='off', left='off', right='off')
    # removes the white frame
    big_axes._frameon = False
    ax1 = fig.add_subplot(1,4,1)
    plt.scatter(normox_center_r, normox_center_t2firstSpk, label='Center')
    plt.scatter(normox_periph_r, normox_periph_t2firstSpk, label='Periphery')
    plt.xlabel(r'Radius ($\mu$m)', fontsize=14)
    plt.ylabel('Time to First Spike (s)', fontsize=14)
    plt.legend(fontsize=14)
    ax2 = fig.add_subplot(1,4,2)
    plt.scatter(normox_center_r, normox_center_spkFreq, label='Center')
    plt.scatter(normox_periph_r, normox_periph_spkFreq, label='Periphery')
    plt.xlabel(r'Radius ($\mu$m)', fontsize=14)
    plt.ylabel('Spike Frequency (Hz)', fontsize=14)
    ax3 = fig.add_subplot(1,4,3)
    plt.scatter(normox_center_r, normox_center_spkDur, label='Center')
    plt.scatter(normox_periph_r, normox_periph_spkDur, label='Periphery')
    plt.xlabel(r'Radius ($\mu$m)', fontsize=14)
    plt.ylabel('Spiking Duration (s)', fontsize=14)
    # plt.ylim(-.05, .9)
    ax4 = fig.add_subplot(1,4,4)
    plt.scatter(normox_center_r, normox_center_nSpks, label='Center')
    plt.scatter(normox_periph_r, normox_periph_nSpks, label='Periphery')
    plt.xlabel(r'Radius ($\mu$m)', fontsize=14)
    plt.ylabel('Number of Spikes', fontsize=14)
    plt.tight_layout()
    fig2, big_axes2 = plt.subplots( figsize=(14, 7.6) , nrows=1, ncols=1, sharey=True)
    # for row, big_ax in enumerate(zip(big_axes,['Perfused']), start=1):
    big_axes2.set_title('Hypoxic', fontsize=22)
    # Turn off axis lines and ticks of the big subplot 
    # obs alpha is 0 in RGBA string!
    big_axes2.tick_params(labelcolor=(1.,1.,1., 0.0), top='off', bottom='off', left='off', right='off')
    # removes the white frame
    big_axes2._frameon = False
    ax5 = fig2.add_subplot(1,4,1)
    plt.scatter(anox_center_r, anox_center_t2firstSpk, label='Center')
    plt.scatter(anox_periph_r, anox_periph_t2firstSpk, label='Periphery')
    plt.xlabel(r'Radius ($\mu$m)', fontsize=14)
    plt.ylabel('Time to First Spike (s)', fontsize=14)
    plt.legend()
    ax6 = fig2.add_subplot(1,4,2)
    plt.scatter(anox_center_r, anox_center_spkFreq, label='Center')
    plt.scatter(anox_periph_r, anox_periph_spkFreq, label='Periphery')
    plt.xlabel(r'Radius ($\mu$m)', fontsize=14)
    plt.ylabel('Spike Frequency (Hz)', fontsize=14)
    ax7 = fig2.add_subplot(1,4,3)
    plt.scatter(anox_center_r, anox_center_spkDur, label='Center')
    plt.scatter(anox_periph_r, anox_periph_spkDur, label='Periphery')
    plt.xlabel(r'Radius ($\mu$m)', fontsize=14)
    plt.ylabel('Duration of Spiking (s)', fontsize=14)
    ax8 = fig2.add_subplot(1,4,4)
    plt.scatter(anox_center_r, anox_center_nSpks, label='Center')
    plt.scatter(anox_periph_r, anox_periph_nSpks, label='Periphery')
    plt.xlabel(r'Radius ($\mu$m)', fontsize=14)
    plt.ylabel('Number of Spikes', fontsize=14)
    # fig.set_figheight(7.6)
    # fig.set_figwidth(14)
    plt.tight_layout()

def hypoxicSDlikeDepol():
    ## Figure 5A,B
    # datadirs = ['Data/SD_Data/perfuse_standard_highNrec/',
    #             'Data/SD_Data/mannitol_standard_higNrec_v2/',
    #             'Data/SD_Data/primed_standard_highNrec_v2/',
    #             'Data/SD_Data/hypox_standard_highNrec/']
    datadirs = ['Data/SD_Data/perfuse_standard_highNrec/',
                '/u/craig/SDinSlice/Data/pad_k015_o2bc_0.1_o2bath_0.03_10s/']
    # colors = ['blue', 'green', 'red']
    colors = ['blue', 'red']
    pos = ['center' for d in datadirs]
    norm_speed = getKwaveSpeed(datadirs[0], r0=100, tcut=8)
    anox_speed = getKwaveSpeed(datadirs[1], r0=0, tcut=8)
    speeds = [norm_speed, anox_speed]
    # labels = [cond + r': %.2f mm/min' % (s) for cond, s in zip(['Perfused', 'Mannitol Treated', 'Propionate Treated','Hypoxic'],speeds)]
    # labels = [cond for cond, s in zip(['Perfused', 'Mannitol', 'Propionate','Hypoxic'],speeds)]
    # labels = [cond for cond, s in zip(['Perfused', 'Propionate','Hypoxic'],speeds)]
    labels = [cond for cond, s in zip(['SD in Perfused Slice', 'Hypoxic SD-like Depolarization'],speeds)]
    legendTitle = 'Depolarization Type'
    fig = plt.figure()
    fig.set_figheight(7)
    fig.set_figwidth(10)
    ax0 = plt.subplot(211)
    compareKwaves(datadirs, labels, legendTitle, colors=colors, sbplt=211)
    ax0.set_xlim(0.0, 10)
    ax0.set_ylim(0,700)
    metrics = [getSpkMetrics(d, uniform=True, position=p) for d, p in zip(datadirs,pos)]
    ax0.set_title(r'K$^+$ Wave', fontsize=18)
    ax0.text(-0.125, 1.25, 'A)', transform=ax0.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')
    spkSpeeds = getSpkWaveSpeed(datadirs, pos, r0=100)
    # labels = [cond + r': %.2f mm/min' % (s) for cond, s in zip(['Perfused', 'Mannitol Treated', 'Propionate Treated','Hypoxic'],spkSpeeds)]
    # labels = [cond + r': %.2f mm/min' % (s) for cond, s in zip(['Perfused', 'Propionate Treated','Hypoxic'],spkSpeeds)]
    labels = [cond + r': %.2f mm/min' % (s) for cond, s in zip(['Perfused', r'Dynamic $\alpha_{ECS}$', 'Propionate Treated','Hypoxic'],spkSpeeds)]
    ax2 = plt.subplot(212)
    for m, l, c in zip(metrics, labels, colors):
        r = [k for k in m.keys()]
        spkDur = [m[k]['spkDur'] for k in m.keys()]
        spkFreq = [m[k]['spkFreq'] for k in m.keys()]
        nSpks = [m[k]['nSpks'] for k in m.keys()]
        t2firstSpk = [m[k]['t2firstSpk'] for k in m.keys()]
        ax2.scatter(t2firstSpk, r, label=l, color=c)
    ax2.set_ylabel(r'Radial Cell Position ($\mu$m)', fontsize=16)
    ax2.set_xlabel('Time to First Spike (s)', fontsize=16)
    ax2.set_title('Depolarization Wave', fontsize=18)
    plt.setp(ax2.get_xticklabels(), fontsize=14)
    plt.setp(ax2.get_yticklabels(), fontsize=14)
    ax2.set_xlim(0,10)
    ax2.set_ylim(0,700)
    ax2.text(-0.125, 1.25, 'B)', transform=ax2.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')
    # legend = plt.legend(title=legendTitle, fontsize=12, bbox_to_anchor=(-0.2, 1.05))
    # plt.setp(legend.get_title(), fontsize=14)

    # host = host_subplot(212, axes_class=axisartist.Axes)
    # par1 = host.twiny()
    # par2 = host.twiny()
    # par1.axis['bottom'] = par1.new_fixed_axis(loc='bottom', offset=(0,-50))
    # par2.axis['bottom'] = par2.new_fixed_axis(loc='bottom', offset=(0,-100))

    # ax3 = plt.subplot(313)
    # datadirs = ['Data/SD_Data/perfuse_alpha07/',
    #             'Data/SD_Data/perfuse_alpha13/',
    #             'Data/SD_Data/perfuse_alpha20/',
    #             'Data/SD_Data/perfuse_alpha26/',
    #             'Data/SD_Data/perfuse_alpha32/',
    #             'Data/alpha42/']
    # speeds_alpha = getKwaveSpeed(datadirs, r0=100, tcut=8)
    # alphas = [0.07, 0.13, 0.20, 0.26, 0.32, 0.42]
    # # p0, = host.plot(alphas, speeds, '*-', linewidth=4, markersize=8, color='darkorange', label=r'$\alpha_{ECS}$')
    # datadirs = ['Data/SD_Data/perfuse_lambda14/',
    #             'Data/SD_Data/perfuse_lambda155/',
    #             'Data/SD_Data/perfuse_lambda17/',
    #             'Data/SD_Data/perfuse_lambda185/',
    #             'Data/SD_Data/perfuse_lambda2/']
    # lambdas = [1.4, 1.55, 1.7, 1.85, 2.0]
    # speeds_lambda = getKwaveSpeed(datadirs, r0=100, tcut=8)
    # # p1, = par1.plot(lambdas, speeds, '*-', linewidth=4, markersize=8, color='aqua', label=r'$\lambda_{ECS}$')
    # datadirs = ['Data/SD_Data/varO201/',
    #             'Data/SD_Data/varO20325/',
    #             'Data/SD_Data/varO2055/',
    #             'Data/SD_Data/varO20775/',
    #             'Data/SD_Data/varO21/']
    # o2 = [0.01, 0.0325, 0.055, 0.0775, 0.1]
    # speeds_o2 = getKwaveSpeed(datadirs, r0=100, tcut=8)
    # datadirs = ['Data/SD_Data/cl65/',
    #             'Data/SD_Data/cl81/',
    #             'Data/SD_Data/cl98/',
    #             'Data/SD_Data/cl114/',
    #             'Data/SD_Data/perfuse_standard_highNrec/']
    # speeds_cl = getKwaveSpeed(datadirs, r0=100, tcut=8)
    # ax3.plot([i for i in range(len(speeds_alpha))], speeds_alpha, '*-', linewidth=4, markersize=8, color='orange', label=r'$\alpha_{ECCS}$')
    # ax3.plot([i for i in range(len(speeds_lambda))], speeds_lambda, '*-', linewidth=4, markersize=8, color='aqua', label=r'$\lambda_{ECS}$')
    # ax3.plot([i for i in range(len(speeds_o2))], speeds_o2, '*-', linewidth=4, markersize=8, color='black', label=r'[O$_2$]')
    # ax3.plot([i for i in range(len(speeds_cl))], speeds_cl, '*-', linewidth=4, markersize=8, color='lime', label=r'[Cl$^-$]')
    # ax3.set_ylabel(r'K$^+$ Wave Speed', fontsize=16)
    # ax3.set_xticks([])
    # ax3.set_title(r'Influence of $\alpha_{ECS}$, $\lambda_{ECS}$, [Cl$^-$], [O$_2$]', fontsize=18)
    # ax3.set_xlim(0,len(speeds_cl)-1)
    # plt.setp(ax3.get_yticklabels(), fontsize=14)
    # ax3.text(-0.125, 1.25, 'C)', transform=ax3.transAxes,
    #     fontsize=18, fontweight='bold', va='top', ha='right')
    # # legend = plt.legend(title='Variable', fontsize=12)#, bbox_to_anchor=(-0.2, 1.05))
    # legend = plt.legend(fontsize=12)#, bbox_to_anchor=(-0.2, 1.05))
    # plt.setp(legend.get_title(), fontsize=14)

    # host.set_ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=16)
    # host.set_xlabel(r'$\alpha_{ECS}$', fontsize=16)
    # par1.set_xlabel(r'$\lambda_{ECS}$', fontsize=16)
    # par2.set_xlabel(r'[O$_2$]$_{bath}$', fontsize=16)
    # host.axis['bottom'].label.set_color(p0.get_color())
    # par1.axis['bottom'].label.set_color(p1.get_color())
    # par2.axis['bottom'].label.set_color(p2.get_color())
    # host.axis['bottom'].label.set_fontsize(16)
    # par1.axis['bottom'].label.set_fontsize(16)
    # par2.axis['bottom'].label.set_fontsize(16)
    # host.axis['bottom'].label.set_fontweight('bold')
    # par1.axis['bottom'].label.set_fontweight('bold')
    # par2.axis['bottom'].label.set_fontweight('bold')
    # host.axis['left'].label.set_fontsize(16)
    # host.text(-0.05, 1.25, 'C)', transform=host.transAxes,
    #     fontsize=18, fontweight='bold', va='top', ha='right')
    # host.axis['bottom'].major_ticklabels.set_fontsize(14)
    # par1.axis['bottom'].major_ticklabels.set_fontsize(14)
    # par2.axis['bottom'].major_ticklabels.set_fontsize(14)
    # host.axis['left'].major_ticklabels.set_fontsize(14) 
    # ttl = r'Influence of $\alpha_{ECS}$, $\lambda_{ECS}$, and [O$_2$]$_{bath}$'
    # host.set_title(ttl)
    # host.title.set_fontsize(18)  

    # plt.setp(host.get_xticklabels(), fontsize=14)
    # plt.setp(host.get_yticklabels(), fontsize=14)
    # plt.setp(par1.get_xticklabels(), fontsize=14)
    # # plt.setp(par1.get_yticklabels(), fontsize=14)
    # plt.setp(par2.get_xticklabels(), fontsize=14)
    # plt.setp(par2.get_yticklabels(), fontsize=14)

    # ax3.plot(alphas, speeds, '*-', linewidth=4, markersize=8, color='blue')
    # ax3.set_xlabel(r'$\alpha_{ECS}$', fontsize=16)
    # ax3.set_ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=16)
    # ax3.set_title(r'Influence of $\alpha_{ECS}$', fontsize=18)
    # ax3.text(-0.05, 1.1, 'C)', transform=ax3.transAxes,
    #     fontsize=18, fontweight='bold', va='top', ha='right')
    # plt.setp(ax3.get_xticklabels(), fontsize=14)
    # plt.setp(ax3.get_yticklabels(), fontsize=14)
    # fig.set_figheight(9)
    # fig.set_figwidth(14)
    plt.draw()
    plt.tight_layout()

def sliceConds():
    ## Figure 5
    # datadirs = ['Data/SD_Data/perfuse_standard_highNrec/',
    #             'Data/SD_Data/mannitol_standard_higNrec_v2/',
    #             'Data/SD_Data/primed_standard_highNrec_v2/',
    #             'Data/SD_Data/hypox_standard_highNrec/']
    datadirs = ['Data/SD_Data/perfuse_standard_highNrec/',
                'Data/dyn_alpha_10s/',
                'Data/SD_Data/primed_standard_highNrec_v2/',
                '/u/craig/SDinSlice/Data/pad_k015_o2bc_0.1_o2bath_0.03_10s/']
                # 'Data/SD_Data/hypox_standard_highNrec/']
    # colors = ['blue', 'green', 'red']
    colors = ['blue', 'purple', 'green', 'red']
    pos = ['center' for d in datadirs]
    speeds = getKwaveSpeed(datadirs, r0=100, tcut=8)
    # labels = [cond + r': %.2f mm/min' % (s) for cond, s in zip(['Perfused', 'Mannitol Treated', 'Propionate Treated','Hypoxic'],speeds)]
    # labels = [cond for cond, s in zip(['Perfused', 'Mannitol', 'Propionate','Hypoxic'],speeds)]
    # labels = [cond for cond, s in zip(['Perfused', 'Propionate','Hypoxic'],speeds)]
    labels = [cond for cond, s in zip(['Perfused', r'Dynamic $\alpha_{ECS}$', 'Propionate','Hypoxic'],speeds)]
    legendTitle = 'Slice Condition'r'Dynamic $\alpha_{ECS}$'
    fig = plt.figure()
    fig.set_figheight(9)
    fig.set_figwidth(18)
    ax0 = plt.subplot(311)
    compareKwaves(datadirs, labels, legendTitle, colors=colors, sbplt=311)
    ax0.set_xlim(0.0, 10)
    ax0.set_ylim(0,700)
    metrics = [getSpkMetrics(d, uniform=True, position=p) for d, p in zip(datadirs,pos)]
    ax0.set_title(r'K$^+$ Wave', fontsize=18)
    ax0.text(-0.125, 1.25, 'A)', transform=ax0.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')
    spkSpeeds = getSpkWaveSpeed(datadirs, pos, r0=100)
    # labels = [cond + r': %.2f mm/min' % (s) for cond, s in zip(['Perfused', 'Mannitol Treated', 'Propionate Treated','Hypoxic'],spkSpeeds)]
    # labels = [cond + r': %.2f mm/min' % (s) for cond, s in zip(['Perfused', 'Propionate Treated','Hypoxic'],spkSpeeds)]
    labels = [cond + r': %.2f mm/min' % (s) for cond, s in zip(['Perfused', r'Dynamic $\alpha_{ECS}$', 'Propionate Treated','Hypoxic'],spkSpeeds)]
    ax2 = plt.subplot(312)
    for m, l, c in zip(metrics, labels, colors):
        r = [k for k in m.keys()]
        spkDur = [m[k]['spkDur'] for k in m.keys()]
        spkFreq = [m[k]['spkFreq'] for k in m.keys()]
        nSpks = [m[k]['nSpks'] for k in m.keys()]
        t2firstSpk = [m[k]['t2firstSpk'] for k in m.keys()]
        ax2.scatter(t2firstSpk, r, label=l, color=c)
    ax2.set_ylabel(r'Radial Cell Position ($\mu$m)', fontsize=16)
    ax2.set_xlabel('Time to First Spike (s)', fontsize=16)
    ax2.set_title('Depolarization Wave', fontsize=18)
    plt.setp(ax2.get_xticklabels(), fontsize=14)
    plt.setp(ax2.get_yticklabels(), fontsize=14)
    ax2.set_xlim(0,10)
    ax2.set_ylim(0,700)
    ax2.text(-0.125, 1.25, 'B)', transform=ax2.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')
    # legend = plt.legend(title=legendTitle, fontsize=12, bbox_to_anchor=(-0.2, 1.05))
    # plt.setp(legend.get_title(), fontsize=14)

    # host = host_subplot(212, axes_class=axisartist.Axes)
    # par1 = host.twiny()
    # par2 = host.twiny()
    # par1.axis['bottom'] = par1.new_fixed_axis(loc='bottom', offset=(0,-50))
    # par2.axis['bottom'] = par2.new_fixed_axis(loc='bottom', offset=(0,-100))

    ax3 = plt.subplot(313)
    datadirs = ['Data/SD_Data/perfuse_alpha07/',
                'Data/SD_Data/perfuse_alpha13/',
                'Data/SD_Data/perfuse_alpha20/',
                'Data/SD_Data/perfuse_alpha26/',
                'Data/SD_Data/perfuse_alpha32/',
                'Data/alpha42/']
    speeds_alpha = getKwaveSpeed(datadirs, r0=100, tcut=8)
    alphas = [0.07, 0.13, 0.20, 0.26, 0.32, 0.42]
    # p0, = host.plot(alphas, speeds, '*-', linewidth=4, markersize=8, color='darkorange', label=r'$\alpha_{ECS}$')
    datadirs = ['Data/SD_Data/perfuse_lambda14/',
                'Data/SD_Data/perfuse_lambda155/',
                'Data/SD_Data/perfuse_lambda17/',
                'Data/SD_Data/perfuse_lambda185/',
                'Data/SD_Data/perfuse_lambda2/']
    lambdas = [1.4, 1.55, 1.7, 1.85, 2.0]
    speeds_lambda = getKwaveSpeed(datadirs, r0=100, tcut=8)
    # p1, = par1.plot(lambdas, speeds, '*-', linewidth=4, markersize=8, color='aqua', label=r'$\lambda_{ECS}$')
    datadirs = ['Data/SD_Data/varO201/',
                'Data/SD_Data/varO20325/',
                'Data/SD_Data/varO2055/',
                'Data/SD_Data/varO20775/',
                'Data/SD_Data/varO21/']
    o2 = [0.01, 0.0325, 0.055, 0.0775, 0.1]
    speeds_o2 = getKwaveSpeed(datadirs, r0=100, tcut=8)
    datadirs = ['Data/SD_Data/cl65/',
                'Data/SD_Data/cl81/',
                'Data/SD_Data/cl98/',
                'Data/SD_Data/cl114/',
                'Data/SD_Data/perfuse_standard_highNrec/']
    speeds_cl = getKwaveSpeed(datadirs, r0=100, tcut=8)
    ax3.plot([i for i in range(len(speeds_alpha))], speeds_alpha, '*-', linewidth=4, markersize=8, color='orange', label=r'$\alpha_{ECCS}$')
    ax3.plot([i for i in range(len(speeds_lambda))], speeds_lambda, '*-', linewidth=4, markersize=8, color='aqua', label=r'$\lambda_{ECS}$')
    ax3.plot([i for i in range(len(speeds_o2))], speeds_o2, '*-', linewidth=4, markersize=8, color='black', label=r'[O$_2$]')
    ax3.plot([i for i in range(len(speeds_cl))], speeds_cl, '*-', linewidth=4, markersize=8, color='lime', label=r'[Cl$^-$]')
    ax3.set_ylabel(r'K$^+$ Wave Speed', fontsize=16)
    ax3.set_xticks([])
    ax3.set_title(r'Influence of $\alpha_{ECS}$, $\lambda_{ECS}$, [Cl$^-$], [O$_2$]', fontsize=18)
    ax3.set_xlim(0,len(speeds_cl)-1)
    plt.setp(ax3.get_yticklabels(), fontsize=14)
    ax3.text(-0.125, 1.25, 'C)', transform=ax3.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')
    # legend = plt.legend(title='Variable', fontsize=12)#, bbox_to_anchor=(-0.2, 1.05))
    legend = plt.legend(fontsize=12)#, bbox_to_anchor=(-0.2, 1.05))
    plt.setp(legend.get_title(), fontsize=14)

    # host.set_ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=16)
    # host.set_xlabel(r'$\alpha_{ECS}$', fontsize=16)
    # par1.set_xlabel(r'$\lambda_{ECS}$', fontsize=16)
    # par2.set_xlabel(r'[O$_2$]$_{bath}$', fontsize=16)
    # host.axis['bottom'].label.set_color(p0.get_color())
    # par1.axis['bottom'].label.set_color(p1.get_color())
    # par2.axis['bottom'].label.set_color(p2.get_color())
    # host.axis['bottom'].label.set_fontsize(16)
    # par1.axis['bottom'].label.set_fontsize(16)
    # par2.axis['bottom'].label.set_fontsize(16)
    # host.axis['bottom'].label.set_fontweight('bold')
    # par1.axis['bottom'].label.set_fontweight('bold')
    # par2.axis['bottom'].label.set_fontweight('bold')
    # host.axis['left'].label.set_fontsize(16)
    # host.text(-0.05, 1.25, 'C)', transform=host.transAxes,
    #     fontsize=18, fontweight='bold', va='top', ha='right')
    # host.axis['bottom'].major_ticklabels.set_fontsize(14)
    # par1.axis['bottom'].major_ticklabels.set_fontsize(14)
    # par2.axis['bottom'].major_ticklabels.set_fontsize(14)
    # host.axis['left'].major_ticklabels.set_fontsize(14) 
    # ttl = r'Influence of $\alpha_{ECS}$, $\lambda_{ECS}$, and [O$_2$]$_{bath}$'
    # host.set_title(ttl)
    # host.title.set_fontsize(18)  

    # plt.setp(host.get_xticklabels(), fontsize=14)
    # plt.setp(host.get_yticklabels(), fontsize=14)
    # plt.setp(par1.get_xticklabels(), fontsize=14)
    # # plt.setp(par1.get_yticklabels(), fontsize=14)
    # plt.setp(par2.get_xticklabels(), fontsize=14)
    # plt.setp(par2.get_yticklabels(), fontsize=14)

    # ax3.plot(alphas, speeds, '*-', linewidth=4, markersize=8, color='blue')
    # ax3.set_xlabel(r'$\alpha_{ECS}$', fontsize=16)
    # ax3.set_ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=16)
    # ax3.set_title(r'Influence of $\alpha_{ECS}$', fontsize=18)
    # ax3.text(-0.05, 1.1, 'C)', transform=ax3.transAxes,
    #     fontsize=18, fontweight='bold', va='top', ha='right')
    # plt.setp(ax3.get_xticklabels(), fontsize=14)
    # plt.setp(ax3.get_yticklabels(), fontsize=14)
    fig.set_figheight(9)
    fig.set_figwidth(14)
    plt.draw()
    plt.tight_layout()

def waveSpeedVsSurfaceArea(sbplt):
    # if sbplt:
    #     plt.subplot(sbplt)
    # else:
    #     fig = plt.figure(figsize=(10,10))

    ## s:v
    datadirs = ['/u/craig/spreadingdepression/Data/SD_Data/perufse_sv02/', 
                '/u/craig/spreadingdepression/Data/SD_Data/perufse_sv1/',
                '/u/craig/spreadingdepression/Data/SD_Data/perufse_sv2/',  
                '/u/craig/spreadingdepression/Data/perfuse_standard/',
                '/u/craig/spreadingdepression/Data/SD_Data/perufse_sv4/', 
                '/u/craig/spreadingdepression/Data/SD_Data/perufse_sv5/',
                '/u/craig/spreadingdepression/Data/SD_Data/perufse_sv6/']
    sa2vs = [0.02, 1.0, 2.0, 3.0, 4.0 , 5.0, 6.0]
    sa_normox = [totalSurfaceArea(1000, 1000, 400, 90000, 0.24, sa2v) for sa2v in sa2vs]
    speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=400)
    datadirs = ['Data/pad_data/s2v02/', 
                'Data/pad_data/s2v1/',
                'Data/pad_data/s2v2/',  
                'Data/pad_data/pad_standard/',
                'Data/pad_data/s2v4/', 
                'Data/pad_data/s2v5/',
                'Data/pad_data/s2v6/']
    speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=400, dur=6)

    ## cell density - constant morphology
    densities = [45, 67.5, 90, 112.5, 120]
    datadirs = ['/u/craig/spreadingdepression/Data/SD_Data/perfuse_d45/',
                '/u/craig/spreadingdepression/Data/SD_Data/perfuse_d675/',  
                '/u/craig/spreadingdepression/Data/perfuse_standard/',
                '/u/craig/spreadingdepression/Data/SD_Data/perfuse_d1125/', 
                '/u/craig/spreadingdepression/Data/SD_Data/perfuse_d120/']
    sa_normox.extend([totalSurfaceArea(1000, 1000, 400, d*1000, 0.24, 3.0, rs=7.52) for d in densities])
    speeds_normox.extend(getKwaveSpeed(datadirs, r0=100, rmax=400))
    datadirs = ['Data/pad_data/d45000_constBeta/',
                'Data/pad_data/d67500_constBeta/',  
                'Data/pad_data/pad_standard/',
                'Data/pad_data/d112500_constBeta/', 
                'Data/pad_data/d120000_constBeta/']
    speeds_anox.extend(getKwaveSpeed(datadirs, r0=100, rmax=400, dur=6))

    ## cell density - constant neuronal volume fraction 
    densities = [45, 67.5, 90, 112.5, 120]
    datadirs = ['/u/craig/spreadingdepression/Data/SD_Data/perfuse_d45_volfrac/',
                '/u/craig/spreadingdepression/Data/SD_Data/perfuse_d675_volfrac/',
                '/u/craig/spreadingdepression/Data/perfuse_standard/',
                '/u/craig/spreadingdepression/Data/SD_Data/perfuse_d1125_volfrac/', 
                '/u/craig/spreadingdepression/Data/SD_Data/perfuse_d120_volfrac/']
    speeds_normox.extend(getKwaveSpeed(datadirs, r0=100, rmax=475))
    sa_normox.extend([totalSurfaceArea(1000, 1000, 400, d*1000, 0.24, 3.0) for d in densities])
    datadirs = ['Data/pad_data/d45000/',
                'Data/pad_data/d67500/',  
                'Data/pad_data/pad_standard/',
                'Data/pad_data/d112500/', 
                'Data/pad_data/d120000/']
    speeds_anox.extend(getKwaveSpeed(datadirs, r0=100, rmax=400, dur=6))

    ## neuronal volume fraction 
    datadirs = ['/u/craig/spreadingdepression/SD_Data/perfuse_nv165/',
                '/u/craig/spreadingdepression/Data/perfuse_standard/',
                '/u/craig/spreadingdepression/SD_Data/perfuse_nv315/',
                '/u/craig/spreadingdepression/SD_Data/perfuse_nv39/',
                '/u/craig/spreadingdepression/SD_Data/perfuse_nv465/',
                '/u/craig/spreadingdepression/SD_Data/perfuse_nv54/']
    speeds_normox.extend(getKwaveSpeed(datadirs, r0=100, rmax=475))
    alphas = [0.165, 0.24, 0.315, 0.39, 0.465, 0.54]
    sa_normox.extend([totalSurfaceArea(1000, 1000, 400, 90000, a, 3.0) for a in alphas])
    datadirs = ['Data/pad_data/betaNrn_165/',
                'Data/pad_data/pad_standard/',
                'Data/pad_data/betaNrn_315/',
                'Data/pad_data/betaNrn_39/',
                'Data/pad_data/betaNrn_465/',
                'Data/pad_data/betaNrn_54/']
    speeds_anox.extend(getKwaveSpeed(datadirs, r0=100, rmax=400, dur=6))

    ## linear regression
    perfuse_train = np.array(speeds_normox[:-8]).reshape(-1,1)
    perfuse_test = np.array(speeds_normox[-8:]).reshape(-1,1)
    hypoxic_train = np.array(speeds_anox[:-8]).reshape(-1,1)
    hypoxic_test = np.array(speeds_anox[-8:]).reshape(-1,1)
    sa_train = np.array(sa_normox[:-8]).reshape(-1,1)
    sa_test = np.array(sa_normox[-8:]).reshape(-1,1)

    regr_perfuse = linear_model.LinearRegression()
    regr_hypoxic = linear_model.LinearRegression()
    regr_perfuse.fit(sa_train, perfuse_train)
    regr_hypoxic.fit(sa_train, hypoxic_train)
    perfuse_pred = regr_perfuse.predict(sa_test)
    hypoxic_pred = regr_hypoxic.predict(sa_test)
    perfuse_r2 = r2_score(perfuse_test, perfuse_pred)
    hypoxic_r2 = r2_score(hypoxic_test, hypoxic_pred)

    sa_pred = np.linspace(np.min(sa_normox), np.max(sa_normox), 10).reshape(-1,1)
    p_pred = regr_perfuse.predict(sa_pred)
    h_pred = regr_hypoxic.predict(sa_pred)

    perfuse_label = r'Perfused: r$^2$ = %0.2f' % (perfuse_r2)
    hypoxic_label = r'Hypoxic: r$^2$ = %0.2f' % (hypoxic_r2)
    sbplt.plot(sa_normox, speeds_normox, '*', color='blue', markersize=8)
    sbplt.plot(sa_normox, speeds_anox, '*', color='red', markersize=8)
    sbplt.plot(sa_pred, p_pred, '--', color='blue', label=perfuse_label, linewidth=4)
    sbplt.plot(sa_pred, h_pred, '--', color='red', label=hypoxic_label, linewidth=4)
    sbplt.set_xlabel(r'Total Neuronal Surface Area ($\mu$m$^2$)', fontsize=16)
    if not sbplt:
        plt.ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=16)
    sbplt.set_title('Total Neuronal Surface Area', fontsize=18)
    # sbplt.set_ylim(-0.1, 15.5)
    sbplt.set_ylim(-0.1, 7.3)
    leg = plt.legend(fontsize=14, title='Slice Condition')
    plt.setp(leg.get_title(), fontsize=16)
    plt.setp(sbplt.get_xticklabels(), fontsize=14)
    plt.setp(sbplt.get_yticklabels(), fontsize=14)

# def bolusProps():
#     fig = plt.figure(figsize=(16,8))
#     ## bolus radius 
#     datadirs = ['Data/SD_Data/hypox_r050/',
#                 'Data/SD_Data/hypox_r075/',
#                 'Data/hypox_standard/',
#                 'Data/SD_Data/hypox_r0125/',
#                 'Data/SD_Data/hypox_r0150/']
#     speeds_anox = getKwaveSpeed(datadirs, r0=100, tcut=8)
#     labels = ['50 $\mu$m', '75 $\mu$m', '100 $\mu$m', '125 $\mu$m', '150 $\mu$m']
#     radii = [50, 75, 100, 125, 150]
#     datadirs = ['Data/SD_Data/perfuse_r050/',
#                 'Data/SD_Data/perfuse_r075/',
#                 'Data/perfuse_standard/',
#                 'Data/SD_Data/perfuse_r0125/',
#                 'Data/SD_Data/perfuse_r0150/']
#     speeds_normox = getKwaveSpeed(datadirs, r0=100, tcut=8)
#     ax0 = plt.subplot(121)
#     ax0.plot(radii, speeds_normox, '*-', linewidth=4, markersize=8, color='blue', label='Perfused')
#     ax0.plot(radii, speeds_anox, '*-', linewidth=4, markersize=8, color='red', label='Hypoxic')
#     ax0.set_xlabel(r'Initial Bolus Radius ($\mu$m)', fontsize=16)
#     ax0.set_ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=16)
#     plt.setp(ax0.get_xticklabels(), fontsize=14)
#     plt.setp(ax0.get_yticklabels(), fontsize=14)
#     legend = plt.legend(fontsize=14, title='Slice Condition')
#     plt.setp(legend.get_title(), fontsize=16)
#     ax0.text(-0.05, 1.1, 'A)', transform=ax0.transAxes,
#         fontsize=18, fontweight='bold', va='top', ha='right')
#     ax0.set_title(r'K$^+$ Bolus Radius', fontsize=18)
#     ax0.set_ylim(0, 7.3)
    
#     datadirs = ['Data/SD_Data/perufse_k030/',
#                 'Data/SD_Data/perufse_k040/',
#                 'Data/SD_Data/perufse_k050/',
#                 'Data/SD_Data/perufse_k060/',
#                 'Data/perfuse_standard/']
#     speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=475)
#     densities = [30, 40, 50, 60, 70]
#     datadirs = ['Data/SD_Data/hypox_k030/',
#                 'Data/SD_Data/hypox_k040/',
#                 'Data/SD_Data/hypox_k050/',
#                 'Data/SD_Data/hypox_k060/',
#                 'Data/hypox_standard/']
#     speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=600)
#     ax1 = plt.subplot(122)
#     ax1.plot(densities, speeds_normox, '*-', linewidth=4, markersize=8, color='blue', label='Perfused')
#     ax1.plot(densities, speeds_anox, '*-', linewidth=4, markersize=8, color='red', label='Hypoxic')
#     ax1.set_xlabel(r'Initial Bolus [K$^+$] (mM)', fontsize=16)
#     ax1.set_ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=16)
#     ax1.set_ylim(0, 7.3)
#     plt.setp(ax1.get_xticklabels(), fontsize=14)
#     plt.setp(ax1.get_yticklabels(), fontsize=14)
#     ax1.text(-0.05, 1.1, 'B)', transform=ax1.transAxes,
#         fontsize=18, fontweight='bold', va='top', ha='right')
#     ax1.set_title(r'Initial Bolus [K$^+$]', fontsize=18)
#     plt.tight_layout()

def tissueProps():
    # Figure 6
    fig = plt.figure(figsize=(14,12))
    ## surface to volume ratio 
    datadirs = ['/u/craig/spreadingdepression/Data/SD_Data/perufse_sv02/', 
                '/u/craig/spreadingdepression/Data/SD_Data/perufse_sv1/',
                '/u/craig/spreadingdepression/Data/SD_Data/perufse_sv2/',  
                '/u/craig/spreadingdepression/Data/perfuse_standard/',
                '/u/craig/spreadingdepression/Data/SD_Data/perufse_sv4/', 
                '/u/craig/spreadingdepression/Data/SD_Data/perufse_sv5/',
                '/u/craig/spreadingdepression/Data/SD_Data/perufse_sv6/']
    speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=400)
    ratios_anox = [0.02, 1, 2, 3, 4, 5, 6]
    datadirs = ['Data/pad_data/s2v02/', 
                'Data/pad_data/s2v1/',
                'Data/pad_data/s2v2/',  
                'Data/pad_data/pad_standard/',
                'Data/pad_data/s2v4/', 
                'Data/pad_data/s2v5/',
                'Data/pad_data/s2v6/']
    speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=400, dur=6)
    ax0 = plt.subplot(321)
    ax0.plot(ratios_anox, speeds_normox, '*-', color='blue', linewidth=4, markersize=8, label='Perfused')
    ax0.plot(ratios_anox, speeds_anox, '*-', color='red', linewidth=4, markersize=8, label='Hypoxic')
    plt.setp(ax0.get_xticklabels(), fontsize=14)
    plt.setp(ax0.get_yticklabels(), fontsize=14)
    ax0.set_xlabel(r'S:V (mm$^{-1}$)', fontsize=16)
    ax0.set_title('S:V', fontsize=18)
    # ax0.set_ylim(-0.1, 15.5)
    ax0.set_ylim(-0.1, 7.3)
    # leg = plt.legend(fontsize=14, title='Slice Condition')
    # plt.setp(leg.get_title(), fontsize=16)
    ax0.text(-0.1, 1.1, 'A)', transform=ax0.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')

    ## cell volume fraction 
    datadirs = ['/u/craig/spreadingdepression/SD_Data/perfuse_nv165/',
                '/u/craig/spreadingdepression/Data/perfuse_standard/',
                '/u/craig/spreadingdepression/SD_Data/perfuse_nv315/',
                '/u/craig/spreadingdepression/SD_Data/perfuse_nv39/',
                '/u/craig/spreadingdepression/SD_Data/perfuse_nv465/',
                '/u/craig/spreadingdepression/SD_Data/perfuse_nv54/']
    speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=400)
    datadirs = ['Data/pad_data/betaNrn_165/',
                'Data/pad_data/pad_standard/',
                'Data/pad_data/betaNrn_315/',
                'Data/pad_data/betaNrn_39/',
                'Data/pad_data/betaNrn_465/',
                'Data/pad_data/betaNrn_54/']
    speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=400, dur=6)
    alphas = [0.165, 0.24, 0.315, 0.39, 0.465, 0.54]
    ax1 = plt.subplot(322)
    ax1.plot(alphas, speeds_normox, '*-', color='blue', linewidth=4, markersize=8, label='Perfused')
    ax1.plot(alphas, speeds_anox, '*-', color='red', linewidth=4, markersize=8, label='Hypoxic')
    plt.setp(ax1.get_xticklabels(), fontsize=14)
    plt.setp(ax1.get_yticklabels(), fontsize=14)
    ax1.set_xlabel(r'$\beta_{nrn}$', fontsize=16)
    ax1.set_title(r'Neuronal Volume Fraction ($\beta_{nrn}$)', fontsize=18)
    # ax1.set_ylim(-0.1, 15.5)
    ax1.set_ylim(-0.1, 7.3)
    ax1.text(-0.1, 1.1, 'B)', transform=ax1.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')

    ## cell density - constant morph 
    datadirs = ['Data/pad_data/d45000_constBeta/',
                'Data/pad_data/d67500_constBeta/',  
                'Data/pad_data/pad_standard/',
                'Data/pad_data/d112500_constBeta/', 
                'Data/pad_data/d120000_constBeta/']
    speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=400, dur=6)
    densities_anox = [45, 67.5, 90, 112.5, 120]
    datadirs = ['/u/craig/spreadingdepression/Data/SD_Data/perfuse_d45/',
                '/u/craig/spreadingdepression/Data/SD_Data/perfuse_d675/',  
                '/u/craig/spreadingdepression/Data/perfuse_standard/',
                '/u/craig/spreadingdepression/Data/SD_Data/perfuse_d1125/', 
                '/u/craig/spreadingdepression/Data/SD_Data/perfuse_d120/']
    speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=400)
    ax2 = plt.subplot(323)
    ax2.plot(densities_anox, speeds_normox, '*-', color='blue', linewidth=4, markersize=8, label='Perfused')
    ax2.plot(densities_anox, speeds_anox, '*-', color='red', linewidth=4, markersize=8, label='Hypoxic')
    plt.setp(ax2.get_xticklabels(), fontsize=14)
    plt.setp(ax2.get_yticklabels(), fontsize=14)
    ax2.set_xlabel(r'Cell Density (neurons/mm$^3$)', fontsize=16)
    ax2.set_title(r'Cell Density - Constant S$_{nrn}$, vol$_{nrn}$', fontsize=18)
    ax2.set_ylabel(r'K$^+$ Wave Speed (mm/min)', fontsize=16)
    # # ax2.set_ylim(-0.1, 15.5)
    ax2.set_ylim(-0.1, 7.3)
    ax2.text(-0.1, 1.4, 'C)', transform=ax2.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')

    ## cell density - constant neuronal volume fraction 
    ## cell density - constant neuronal volume fraction 
    densities = [45, 67.5, 90, 112.5, 120]
    datadirs = ['/u/craig/spreadingdepression/Data/SD_Data/perfuse_d45_volfrac/',
                '/u/craig/spreadingdepression/Data/SD_Data/perfuse_d675_volfrac/',  
                '/u/craig/spreadingdepression/Data/perfuse_standard/',
                '/u/craig/spreadingdepression/Data/SD_Data/perfuse_d1125_volfrac/', 
                '/u/craig/spreadingdepression/Data/SD_Data/perfuse_d120_volfrac/']
    speeds_normox = getKwaveSpeed(datadirs, r0=100, rmax=400)
    datadirs = ['Data/pad_data/d45000/',
                'Data/pad_data/d67500/',  
                'Data/pad_data/pad_standard/',
                'Data/pad_data/d112500/', 
                'Data/pad_data/d120000/']
    speeds_anox = getKwaveSpeed(datadirs, r0=100, rmax=400, dur=6)
    ax3 = plt.subplot(324)
    ax3.plot(densities_anox, speeds_normox, '*-', color='blue', linewidth=4, markersize=8, label='Perfused')
    ax3.plot(densities_anox, speeds_anox, '*-', color='red', linewidth=4, markersize=8, label='Hypoxic')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax3.set_xlabel(r'Cell Density (neurons/mm$^3$)', fontsize=16)
    ax3.set_title(r'Cell Density - Constant $\beta_{nrn}$, S:V', fontsize=18)
    # ax3.set_ylim(-0.1, 15.5)
    ax3.set_ylim(-0.1, 7.3)
    ax3.text(-0.1, 1.4, 'D)', transform=ax3.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')

    ax4 = plt.subplot(313)
    waveSpeedVsSurfaceArea(sbplt=ax4)
    ax4.text(-0.05, 1.1, 'E)', transform=ax4.transAxes,
        fontsize=18, fontweight='bold', va='top', ha='right')

    plt.tight_layout()