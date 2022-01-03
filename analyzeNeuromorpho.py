from numpy.core.numeric import True_
import requests 
import numpy as np
from matplotlib import pyplot as plt 
plt.ion()

brainstem_rxvs = ['Williams_Wilson', 'Mena-Segovia', 'Cameron', 'Irintchev']
neo_pyramidal_rxvs = ['BICCN-MOp-miniatlas-anatomy', 'Dendritica', 'Suter_Shepherd', 'Groh', 'Buchs', 'DeKock', 'Kawaguchi', 'Jacobs', 'Meyer']
neo_int_rxvs = ['Yuste', 'Sun', 'Bikson']
hipp_pyr_rxvs = ['Arnold_Johnston', 'Soltesz']
granule_rxvs = ['Beining', 'Bausch']
hipp_int_rxvs = ['Anstoetz', 'Gulyas', 'Hajos', 'Turner']
out = {}

base_url = 'http://neuromorpho.org/api/neuron/select?q=archive:'

all_rxvs = brainstem_rxvs + neo_pyramidal_rxvs + neo_int_rxvs + hipp_pyr_rxvs + granule_rxvs + hipp_int_rxvs

for rxv in all_rxvs:
    full_url = base_url + rxv
    r = requests.get(full_url)
    data = r.json()
    out[rxv] = {}
    for cell in data['_embedded']['neuronResources']:
        if (cell['attributes'] == 'Diameter, 3D, Angles') and (cell['domain'].startswith('Dendrites')):
            if 'neocortex' in cell['brain_region'] and (('principal cell' in cell['cell_type']) or ('pyramidal' in cell['cell_type'])):
                if 'neo_pyr' not in list(out[rxv].keys()):
                    out[rxv]['neo_pyr'] = {'sv' : [float(cell['surface']) / float(cell['volume'])],
                                            'min_age' : [cell['min_age']],
                                            'max_age' : [cell['max_age']],
                                            'reference' : cell['reference_doi'],
                                            'avg_sv' : None,
                                            'std_sv' : None,
                                            'brain_region' : cell['brain_region'],
                                            'species' : cell['species'],
                                            'cell_type' : cell['cell_type'],
                                            'n' : None}
                else:
                    out[rxv]['neo_pyr']['sv'].append(float(cell['surface']) / float(cell['volume']))
                    out[rxv]['neo_pyr']['min_age'].append(cell['min_age'])
                    out[rxv]['neo_pyr']['max_age'].append(cell['max_age'])
                    out[rxv]['neo_pyr']['avg_sv'] = np.mean(out[rxv]['neo_pyr']['sv'])
                    out[rxv]['neo_pyr']['std_sv'] = np.std(out[rxv]['neo_pyr']['sv'])
                    out[rxv]['neo_pyr']['n'] = len(out[rxv]['neo_pyr']['sv'])
            if ('neocortex' in cell['brain_region']) and ('interneuron' in cell['cell_type']):
                if 'neo_int' not in list(out[rxv].keys()):
                    out[rxv]['neo_int'] = {'sv' : [float(cell['surface']) / float(cell['volume'])],
                                            'min_age' : [cell['min_age']],
                                            'max_age' : [cell['max_age']],
                                            'reference' : cell['reference_doi'],
                                            'avg_sv' : None,
                                            'std_sv' : None,
                                            'brain_region' : cell['brain_region'],
                                            'species' : cell['species'],
                                            'cell_type' : cell['cell_type'],
                                            'n' : None}
                else:
                    out[rxv]['neo_int']['sv'].append(float(cell['surface']) / float(cell['volume']))
                    out[rxv]['neo_int']['min_age'].append(cell['min_age'])
                    out[rxv]['neo_int']['max_age'].append(cell['max_age'])
                    out[rxv]['neo_int']['avg_sv'] = np.mean(out[rxv]['neo_int']['sv'])
                    out[rxv]['neo_int']['std_sv'] = np.std(out[rxv]['neo_int']['sv'])
                    out[rxv]['neo_int']['n'] = len(out[rxv]['neo_int']['sv'])
            if ('brainstem' in cell['brain_region']):
                if 'brainstem' not in list(out[rxv].keys()):
                    out[rxv]['brainstem'] = {'sv' : [float(cell['surface']) / float(cell['volume'])],
                                            'min_age' : [cell['min_age']],
                                            'max_age' : [cell['max_age']],
                                            'reference' : cell['reference_doi'],
                                            'avg_sv' : None,
                                            'std_sv' : None,
                                            'brain_region' : cell['brain_region'],
                                            'species' : cell['species'],
                                            'cell_type' : [cell['cell_type']],
                                            'n' : None}
                else:
                    out[rxv]['brainstem']['sv'].append(float(cell['surface']) / float(cell['volume']))
                    out[rxv]['brainstem']['min_age'].append(cell['min_age'])
                    out[rxv]['brainstem']['max_age'].append(cell['max_age'])
                    out[rxv]['brainstem']['avg_sv'] = np.mean(out[rxv]['brainstem']['sv'])
                    out[rxv]['brainstem']['std_sv'] = np.std(out[rxv]['brainstem']['sv'])
                    out[rxv]['brainstem']['n'] = len(out[rxv]['brainstem']['sv'])
                    out[rxv]['brainstem']['cell_type'].append(cell['cell_type'])
            if ('hippocampus' in cell['brain_region']) and ('Control' in cell['experiment_condition']) and (('pyramidal' in cell['cell_type']) or ('pyramidal cell' in cell['cell_type'])):
                if 'hipp_pyr' not in list(out[rxv].keys()):
                    out[rxv]['hipp_pyr'] = {'sv' : [float(cell['surface']) / float(cell['volume'])],
                                            'min_age' : [cell['min_age']],
                                            'max_age' : [cell['max_age']],
                                            'reference' : cell['reference_doi'],
                                            'avg_sv' : None,
                                            'std_sv' : None,
                                            'brain_region' : cell['brain_region'],
                                            'species' : cell['species'],
                                            'cell_type' : cell['cell_type'],
                                            'n' : None}
                else:
                    out[rxv]['hipp_pyr']['sv'].append(float(cell['surface']) / float(cell['volume']))
                    out[rxv]['hipp_pyr']['min_age'].append(cell['min_age'])
                    out[rxv]['hipp_pyr']['max_age'].append(cell['max_age'])
                    out[rxv]['hipp_pyr']['avg_sv'] = np.mean(out[rxv]['hipp_pyr']['sv'])
                    out[rxv]['hipp_pyr']['std_sv'] = np.std(out[rxv]['hipp_pyr']['sv'])
                    out[rxv]['hipp_pyr']['n'] = len(out[rxv]['hipp_pyr']['sv'])
            if ('hippocampus' in cell['brain_region']) and ('Control' in cell['experiment_condition']) and ('granule' in cell['cell_type']):
                if 'hipp_gran' not in list(out[rxv].keys()):
                    out[rxv]['hipp_gran'] = {'sv' : [float(cell['surface']) / float(cell['volume'])],
                                            'min_age' : [cell['min_age']],
                                            'max_age' : [cell['max_age']],
                                            'reference' : cell['reference_doi'],
                                            'avg_sv' : None,
                                            'std_sv' : None,
                                            'brain_region' : cell['brain_region'],
                                            'species' : cell['species'],
                                            'cell_type' : cell['cell_type'],
                                            'n' : None}
                else:
                    out[rxv]['hipp_gran']['sv'].append(float(cell['surface']) / float(cell['volume']))
                    out[rxv]['hipp_gran']['min_age'].append(cell['min_age'])
                    out[rxv]['hipp_gran']['max_age'].append(cell['max_age'])
                    out[rxv]['hipp_gran']['avg_sv'] = np.mean(out[rxv]['hipp_gran']['sv'])
                    out[rxv]['hipp_gran']['std_sv'] = np.std(out[rxv]['hipp_gran']['sv'])
                    out[rxv]['hipp_gran']['n'] = len(out[rxv]['hipp_gran']['sv'])
            if ('hippocampus' in cell['brain_region']) and ('Control' in cell['experiment_condition']) and ('interneuron' in cell['cell_type']):
                if 'hipp_inter' not in list(out[rxv].keys()):
                    out[rxv]['hipp_inter'] = {'sv' : [float(cell['surface']) / float(cell['volume'])],
                                            'min_age' : [cell['min_age']],
                                            'max_age' : [cell['max_age']],
                                            'reference' : cell['reference_doi'],
                                            'avg_sv' : None,
                                            'std_sv' : None,
                                            'brain_region' : cell['brain_region'],
                                            'species' : cell['species'],
                                            'cell_type' : cell['cell_type'],
                                            'n' : None}
                else:
                    out[rxv]['hipp_inter']['sv'].append(float(cell['surface']) / float(cell['volume']))
                    out[rxv]['hipp_inter']['min_age'].append(cell['min_age'])
                    out[rxv]['hipp_inter']['max_age'].append(cell['max_age'])
                    out[rxv]['hipp_inter']['avg_sv'] = np.mean(out[rxv]['hipp_inter']['sv'])
                    out[rxv]['hipp_inter']['std_sv'] = np.std(out[rxv]['hipp_inter']['sv'])
                    out[rxv]['hipp_inter']['n'] = len(out[rxv]['hipp_inter']['sv'])

# plt.plot(sv, '*')
# plt.scatter(age, sv)

for archive in out.keys(): 
    print(archive) 
    for cell_type in out[archive].keys():
        print(cell_type) 
        print('s:v ' + str(out[archive][cell_type]['avg_sv'])) 
        print('std ' + str(out[archive][cell_type]['std_sv'])) 
        print('species ' + out[archive][cell_type]['species']) 
        print('cell type' + str(out[archive][cell_type]['cell_type']))
        print('n ' + str(out[archive][cell_type]['n']))
        print('\n')

bstem = []
adlt_bstem = []
bstem_refs = []
for rxv in brainstem_rxvs:
    adlt = []
    for sv, age in zip(out[rxv]['brainstem']['sv'], out[rxv]['brainstem']['min_age']):
        if age != 'Not Reported' and age != '4.5':
            if float(age) > 14:
                adlt.append(sv)
        else:
            adlt.append(sv)
    # adlt = [sv for sv, age in zip(out[rxv]['brainstem']['sv'], out[rxv]['brainstem']['min_age']) if (age is not 'Not Reported') and float(age) > 14]
    adlt_bstem = adlt_bstem + adlt
    bstem = bstem + out[rxv]['brainstem']['sv']
    bstem_refs.append(out[rxv]['brainstem']['reference'])

rat_neopyr = []
rat_neopyr_rxvs = []
rat_neopyr_refs = []
for rxv in all_rxvs:
    if 'neo_pyr' in (out[rxv].keys()):
        if out[rxv]['neo_pyr']['species'] == 'rat':
            rat_neopyr_rxvs.append(rxv)
            rat_neopyr = rat_neopyr + out[rxv]['neo_pyr']['sv']
            rat_neopyr_refs.append(out[rxv]['neo_pyr']['reference'])

mouse_neopyr_rxvs = []
mouse_neopyr_refs = []
for rxv in all_rxvs:
    if 'neo_pyr' in (out[rxv].keys()):
        if out[rxv]['neo_pyr']['species'] == 'mouse':
            mouse_neopyr_rxvs.append(rxv)
            mouse_neopyr = rat_neopyr + out[rxv]['neo_pyr']['sv']
            mouse_neopyr_refs.append(out[rxv]['neo_pyr']['reference'])

rat_neoint = []
rat_neoint_rxvs = []
rat_neoint_refs = []
for rxv in all_rxvs:
    if 'neo_int' in (out[rxv].keys()):
        if out[rxv]['neo_int']['species'] == 'rat':
            rat_neoint_rxvs.append(rxv)
            rat_neoint = rat_neopyr + out[rxv]['neo_int']['sv']
            rat_neoint_refs.append(out[rxv]['neo_int']['reference'])

mouse_neoint = []
mouse_neoint_rxvs = []
mouse_neoint_refs = []
for rxv in all_rxvs:
    if 'neo_int' in (out[rxv].keys()):
        if out[rxv]['neo_int']['species'] == 'mouse':
            mouse_neoint_rxvs.append(rxv)
            mouse_neoint = rat_neopyr + out[rxv]['neo_int']['sv']
            mouse_neoint_refs.append(out[rxv]['neo_int']['reference'])


print('Brainstem: ' + str(round(np.mean(bstem),2)) + ' +/- ' + 
        str(round(np.std(bstem),2)) + ' (n=' + str(len(bstem)) + ')\n')
print('Mature Brainstem: ' + str(round(np.mean(adlt_bstem),2)) + ' +/- ' + 
        str(round(np.std(adlt_bstem),2)) + ' (n=' + str(len(adlt_bstem)) + ')\n')
refstr = ''
for ref in bstem_refs:
    if ref:
        refstr = refstr + ref[0] + ', '
    else:
        refstr = refstr + 'None, '
refstr = refstr[:-2] + '\n'
print(refstr)

print('Rat Nctx PYR: ' + str(round(np.mean(rat_neopyr),2)) + ' +/- ' + 
        str(round(np.std(rat_neopyr),2)) + ' (n=' + str(len(rat_neopyr)) + ')\n')
refstr = ''
for ref, rxv in zip(rat_neopyr_refs, rat_neopyr_rxvs):
    if ref:
        refstr = refstr + ref[0] + ', '
    else:
        refstr = refstr + rxv + ', '
refstr = refstr[:-2] + '\n'
print(refstr)

print('Rat Nctx INT: ' + str(round(np.mean(rat_neoint),2)) + ' +/- ' + 
        str(round(np.std(rat_neoint),2)) + ' (n=' + str(len(rat_neoint)) + ')\n')
for ref, rxv in zip(rat_neoint_refs, rat_neoint_rxvs):
    if ref:
        refstr = refstr + ref[0] + ', '
    else:
        refstr = refstr + rxv + ', '
refstr = refstr[:-2] + '\n'
print(refstr)

print('Mouse Nctx PYR: ' + str(round(np.mean(mouse_neopyr),2)) + ' +/- ' + 
        str(round(np.std(mouse_neopyr),2)) + ' (n=' + str(len(mouse_neopyr)) + ')\n')
for ref, rxv in zip(mouse_neopyr_refs, mouse_neopyr_rxvs):
    if ref:
        refstr = refstr + ref[0] + ', '
    else:
        refstr = refstr + rxv + ', '
refstr = refstr[:-2] + '\n'
print(refstr)

print('Mouse Nctx INT: ' + str(round(np.mean(mouse_neoint),2)) + ' +/- ' + 
        str(round(np.std(mouse_neoint),2)) + ' (n=' + str(len(mouse_neoint)) + ')\n')
for ref, rxv in zip(mouse_neoint_refs, mouse_neoint_rxvs):
    if ref:
        refstr = refstr + ref[0] + ', '
    else:
        refstr = refstr + rxv + ', '
refstr = refstr[:-2] + '\n'
print(refstr)

# v1.0 - analysis of S:V for cells from NeuroMorpho with largely complete 3D dendritic reconstructions