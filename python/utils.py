"""
Utilities for plotting the test case.

Thomas Sharp, 2013
thomas.sharp@riken.jp
"""
import numpy



def spikes_decompress(spikes, shape=None):
    """
    Converts spike array from sparse to dense encoding.

    :param numpy.ndarray spikes:
        spikes in two columns corresponding to firing times in milliseconds and
        cell IDs in unique integers.
    :param tuple shape:
        the dimensions of the output matrix, corresponding to the number of
        cells and the period of time that the matrix covers. If no parameter is
        provided, the maximum cell ID and spike time is read from the inputs.

    :returns:
        2d array with neurons on columns, time on rows and spikes represented
        by ones.
    """
    # Find the output array shape if necessary
    if not shape:
        shape = spikes.max(axis=0) + 1
    # Set up the output array
    decompressed = numpy.zeros(shape)
    # Ensure that we have array-addressable values
    spikes = spikes.astype(numpy.int)
    # Write spikes into decompressed format
    decompressed[spikes[:,0], spikes[:,1]] = 1

    return decompressed
    
    
def spikes_recompress(spikes):
    """
    Converts spikes from dense to sparse encoding.

    :param numpy.ndarray spikes:
        2d array with neurons on columns, time on rows and spikes represented
        by ones.

    :returns:
        spikes in two columns corresponding to firing times in milliseconds and
        cell IDs in unique integers.
    """
    # Return spikes to compressed format
    output = numpy.argwhere(spikes)

    return output
