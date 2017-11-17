from pyteomics import mgf
import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def func(x, l):
    return l * np.exp(-l * x)

# read_mgf will read in an mgf file and return a list of the
# spectra contained within that file
def read_mgf(filename):
    spectra = mgf.read(filename)
    return spectra

'''
def filter_noise(spectrum):
    intensity = np.array(spectrum['intensity array'])
    is_noise = [True] * len(intensity)
    if len(intensity) == 0:
        return []

    while True:
        old_avg = np.mean(intensity[is_noise])
        old_std = np.std(intensity[is_noise])
        for i in range(0, len(intensity)):
            if intensity[i] > old_avg + 3 * old_std:
                is_noise[i] = False
        avg = np.mean(intensity[is_noise])
        if old_avg == avg:
            break

    return is_noise
'''

def filter_noise(spectrum):
    intensities = spectrum['intensity array']
    n_steps = max(intensities)/4
    step_size = max(intensities)/n_steps
    print 'step size',step_size
    bins = np.arange(0,max(intensities)+step_size,step_size)
    inds = np.digitize(intensities, bins)
    #print bins
    unique, counts = np.unique(inds, return_counts=True)
    x = []
    y = []
    n = len(intensities)
    for i in range(1, len(bins)):
        x.append(int(i))
        y.append(float((inds == i).sum())/n)
    #print y[1:]
    #plt.plot(x[1:50], y[1:50], 'b-', label='data')
    sub_n = n/4
    l = (sub_n - 2)/np.sum(y[1:sub_n])
    for thing in x:
        calc = func(thing, l)
        if calc * n < 1:
            cutoff = bins[thing]
            print l,cutoff
            return cutoff
    return 0

filename = sys.argv[1]
print 'Processing',filename

outfile_name = './' + filename.split('.')[0] + '_filtered.' + filename.split('.')[1]
outfile = open(outfile_name, 'w')

processed_spectra = []
spectra = read_mgf(filename)
for spectrum in spectra:
    print len(processed_spectra)
    new_entry = {}
    mz = np.array(spectrum['m/z array'])
    intensity = np.array(spectrum['intensity array'])
    #is_noise = filter_noise(spectrum)
    #new_mz = mz[np.invert(is_noise)]
    #new_intensity = intensity[np.invert(is_noise)]
    cutoff = filter_noise(spectrum)
    new_mz = []
    new_intensity = []
    for i in range(0, len(mz)):
        if intensity[i] > cutoff:
            new_mz.append(mz[i])
            new_intensity.append(intensity[i])
    new_entry['m/z array'] = new_mz
    new_entry['intensity array'] = new_intensity
    new_entry['params'] = spectrum['params']
    new_entry['charge array'] = spectrum['charge array']
    processed_spectra.append(new_entry)

mgf.write(processed_spectra, output=outfile)
outfile.close()
