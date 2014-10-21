"""
Test case for the Spike-Train Communities algorithm (Humphries, 2011).

Run this file ($ python test_cases.py) for a demonstration.

Thomas Sharp, 2013
thomas.sharp@riken.jp
"""
import itertools
import numpy
import pylab

import __init__ as assembler
import utils



def demonstrate():
    """
    Demonstrates the assembly algorithm.
    """
    from test_data import n_trains, spikes, train_period
    # Load and sanitise test data, and run grouping algorithm
    spikes[:,0] -= 1 # To establish zero-indexing
    spikes[:,1] *= 1e3 # Convert spike times into milliseconds
    spikes[:,[0,1]] = spikes[:,[1,0]] # Switch column order
    train_period *= 1e3 # Convert period into milliseconds
    Q, L, C, T = assembler.find_assemblies(n_trains, train_period, spikes)

    # Plot results
    fig, ax = pylab.subplots(1,1,1)
    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'])
    count = 0
    # Scatter each group of spikes separately, for colouring
    spk = utils.spikes_decompress(spikes)
    for i in range(L.max() + 1):
        idx = numpy.argwhere(L == i).flatten()
        sp = utils.spikes_recompress(spk[:,idx])
        ax.scatter(sp[:,0], sp[:,1] + count, color=next(colors), s=5.)
        count += idx.size
    # Set axis details
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Cell ID')
    ax.set_xlim(-50,1050)
    ax.set_ylim(-5,110)
    pylab.show()


if __name__ == '__main__':
    demonstrate()
