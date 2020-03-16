import pandas as pd
import numpy as np

'''
plot variants on axis
'''
def plot(variants, axis, name):
    axis.plot(variants.location / 1000000, variants.n_events, color="black")
    axis.set_ylabel(name)

    return axis


'''
prepare variants data to be plotted at a chrom
'''
def prepare_at_chrom(variants, chrom):
    variants = variants[variants["chr"] == str(chrom)]
    return bin_frequencies(variants.pos, 200, variants.pos.min(),
                           variants.pos.max())


'''
bin variant data
'''
def bin_frequencies(locations, n_bins, start, extent):
    bins = np.linspace(start, extent, n_bins)
    digitized = np.digitize(locations, bins)

    binned_loc = [locations[digitized == i].mean() for i in range(1, len(bins))]
    n_events = [len(locations[digitized == i]) for i in range(1, len(bins))]

    return pd.DataFrame({"location": binned_loc,
                       "n_events": n_events})


'''
read in
'''
def read(f):
    data= pd.DataFrame(parse(open(f), "\t"), columns=["chr", "pos"])
    data = data.astype({"pos": np.int64})
    return data


def parse(f, sep):
    return [line.split(sep)[:2] for line in f
            if not line.startswith('#')]
