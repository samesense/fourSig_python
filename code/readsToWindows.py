def calcWindows(data, windowSize):
    """Return a list of read counts in overlapping windows
       of size windowSize.
       data is a list of reads for this chrom."""
    ret = []
    for idx in range(data):
        ret = sum(data[idx:idx+windowSize])
    return ret


    
