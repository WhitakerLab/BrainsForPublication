import numpy as np
import nibabel as nb
from scipy.ndimage import label as sci_label
from matplotlib import pylab
from matplotlib import pyplot as plt
from nipy.labs import viz
from os.path import join as opj
from os.path import basename as opb
from glob import glob as gg
from IPython.display import Image, display
import shutil


def get_labels(data, threshold=0, min_extent=0):
    """
    Get number of clusters in dataset as well as labeled volume
    Minimal extend of cluster and voxel-vise threshold can be specified
    """
    binarized_data = abs(data) > threshold

    # corrects directly the loaded data file
    data[np.invert(binarized_data)] = 0

    labels, nlabels = sci_label(binarized_data)
    for idx in range(1, nlabels + 1):
        if np.sum(labels == idx) < min_extent:
            labels[labels == idx] = 0
    nlabels = len(np.setdiff1d(np.unique(labels), [0]))
    return sci_label(labels)


def get_cluster_info(img, affine, data):
    """
    Returns peak coordinations of a dataset, if labeled file,
        data set and affine matrix is provided
    """
    coords = []
    labelID = np.setdiff1d(np.unique(img.ravel()), [0])
    cs = []
    maxcoords = []
    valuesInfo = []
    for lab in labelID:
        sumval = np.sum(img == lab)
        cs.append(sumval)
        maxval = np.max(data[img == lab])
        minval = np.min(data[img == lab])
        meanval = np.mean(data[img == lab])
        stdval = np.std(data[img == lab])
        maxidx = np.nonzero(np.multiply(data, img == lab) == maxval)
        maxcoords.append([m[0] for m in maxidx])
        valuesInfo.append([sumval, maxval, minval, meanval, stdval])

    maxcoords = np.asarray(maxcoords)
    valuesInfo = np.asarray(valuesInfo)
    maxcoords = maxcoords[np.argsort(cs)[::-1], :]
    valuesInfo = valuesInfo[np.argsort(cs)[::-1], :]
    for i, lab in enumerate(labelID[np.argsort(cs)[::-1]]):
        coords.append(np.dot(affine,
                             np.hstack((maxcoords[i], 1)))[:3].tolist())
    valuesInfo = np.hstack((np.array(coords), valuesInfo))
    return coords, valuesInfo


def show_slices(data, affine, coords=None, cmap=None, prefix=None,
                show_colorbar=None, formatter='%.2f', showCross=False):

    if cmap is None:
        cmap = pylab.cm.hot

    anatimg = nb.load('../templates/MNI152_T1_1mm_brain.nii.gz')
    anatdata, anataff = anatimg.get_data(), anatimg.get_affine()
    anatdata = anatdata.astype(np.float)
    anatdata[anatdata < 10.] = np.nan

    for idx, coord in enumerate(coords):
            outfile = 'cluster%02d' % idx
            if prefix:
                outfile = '_'.join((prefix, outfile))
            osl = viz.plot_map(
                np.asarray(data), affine, anat=anatdata, anat_affine=anataff,
                threshold=0.1, cmap=cmap, annotate=True, black_bg=False,
                cut_coords=coord, draw_cross=showCross, slicer='ortho')
            if show_colorbar:
                im = plt.gca().get_images()[1]
                cb = plt.colorbar(im, cax=plt.axes([0.4, 0.05, 0.2, 0.025]),
                                  orientation='horizontal', format=formatter)
                cb.set_ticks([cb._values.min(), cb._values.max()])

            osl.frame_axes.figure.savefig(
                opj('figures', outfile + '.png'), dpi=300,
                bbox_inches='tight', transparent=True)
            plt.close()
