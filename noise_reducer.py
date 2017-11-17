from pyteomics import mgf
import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from multiprocessing import Pool as ThreadPool

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
    mz = spectrum['m/z array']
    intensities = spectrum['intensity array']
    sorted_intensities = sorted(intensities, reverse=True)
    running_variance = []
    for i in range(0, len(sorted_intensities)):
        running_variance.append(np.var(sorted_intensities[i:len(sorted_intensities)]))
    
    relative_variance_change = []
    for i in range(1, len(running_variance)):
        relative_variance_change.append(-(running_variance[i] - running_variance[i-1])/running_variance[i-1]*100)
    
    cutoff_index = -1
    for i in range(0,len(relative_variance_change)):
        change = relative_variance_change[i]
        if change < 1:
            cutoff_index = i+1
            break
    if cutoff_index > -1:
        cutoff = sorted_intensities[cutoff_index]
    else:
        cutoff = 2*np.median(intensities)

    spectrum['params']['noise_level'] = cutoff

    new_entry = {}
    new_mz = []
    new_intensity = []
    for i in range(0, len(mz)):
        if intensities[i] > cutoff:
            new_mz.append(mz[i])
            new_intensity.append(intensities[i])
    new_entry['m/z array'] = new_mz
    new_entry['intensity array'] = new_intensity
    new_entry['params'] = spectrum['params']
    new_entry['charge array'] = spectrum['charge array']

    return new_entry


if __name__ == '__main__':
    filename = sys.argv[1]
    print('Processing',filename)

    outfile_name = '.\\' + filename.split('.')[1] + '_filtered.' + filename.split('.')[2]

    outfile = open(outfile_name, 'w')

    spectra = read_mgf(filename)
    pool = ThreadPool()

    processed_spectra = pool.map(filter_noise, spectra)
    pool.close()
    pool.join()    

    print('Writing to file...')
    mgf.write(processed_spectra, output=outfile)
    outfile.close()
