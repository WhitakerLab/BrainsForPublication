import os
import sys
import numpy as np
import nibabel as nb
from scipy.ndimage import label as sci_label
from matplotlib import pylab
from matplotlib import pyplot as plt
from nipy.labs import viz
from os.path import join as opj
from os.path import basename as opb
from glob import glob as gg


def get_labels(data, threshold=0, min_extent=0):
    """
    Get number of clusters in dataset as well as a labeled volume
    Minimal extend of cluster and voxel-vise threshold can be specified
    """
    binarized_data = abs(data) > threshold

    labels, nlabels = sci_label(binarized_data)
    for idx in range(1, nlabels + 1):
        if np.sum(labels == idx) < min_extent:
            labels[labels == idx] = 0
    labels, nlabels = sci_label(labels)

    # overwrites the input data with the thresholded data
    binarized_data = labels.astype('bool')
    data[np.invert(binarized_data)] = 0

    return labels, nlabels


def get_cluster_info(img, affine, data):
    """
    Returns peak coordinations and cluster information of a given dataset,
    if labeled file and affine is provided
    """
    coords = []
    labelID = np.setdiff1d(np.unique(img.ravel()), [0])
    cs = []
    maxcoords = []
    clusterInfo = []
    for lab in labelID:
        # Computes Max, Min, Mean, Std and peak coordinate of each cluster
        sumval = np.sum(img == lab)
        cs.append(sumval)
        maxval = np.max(data[img == lab])
        minval = np.min(data[img == lab])
        meanval = np.mean(data[img == lab])
        stdval = np.std(data[img == lab])
        maxidx = np.nonzero(np.multiply(data, img == lab) == maxval)
        maxcoords.append([m[0] for m in maxidx])
        clusterInfo.append([sumval, maxval, minval, meanval, stdval])

    maxcoords = np.asarray(maxcoords)
    clusterInfo = np.asarray(clusterInfo)
    maxcoords = maxcoords[np.argsort(cs)[::-1], :]
    clusterInfo = clusterInfo[np.argsort(cs)[::-1], :]
    for i, lab in enumerate(labelID[np.argsort(cs)[::-1]]):
        coords.append(np.dot(affine,
                             np.hstack((maxcoords[i], 1)))[:3].tolist())
    clusterInfo = np.hstack((np.array(coords), clusterInfo))

    # Returns peak coordination and additional cluster infos
    return coords, clusterInfo


def show_slices(data, affine, coords=None, cmap=None, prefix='cluster',
                foldername=None, show_colorbar=None, showCross=False,
                threshold=0, annotate=True, orientation='horizontal',
                template='../scripts/templates/MNI152_T1_1mm_brain.nii.gz',
                dpiRes=300, suffix='png', showTitle=False):

    # Specify colormap
    if cmap is None:
        cmap = pylab.cm.hot

    # Prepare background image
    anatimg = nb.load(template)
    anatdata, anataff = anatimg.get_data(), anatimg.get_affine()
    anatdata = anatdata.astype(np.float)
    anatdata[anatdata < 10.] = np.nan

    # Create output figure for each coordinate
    for idx, coord in enumerate(coords):
        outfile = '%02d' % idx
        if prefix:
            outfile = '_'.join((prefix, outfile))
        if showTitle:
            title = outfile + ' %s' % coord
        else:
            title = ''
        osl = viz.plot_map(
            np.asarray(data), affine, anat=anatdata, anat_affine=anataff,
            threshold=threshold, cmap=cmap, annotate=annotate,
            black_bg=False, cut_coords=coord, draw_cross=showCross,
            slicer='ortho', title=title)
        if show_colorbar:
            if orientation == 'horizontal':
                cbarLocation = [0.4, 0.01, 0.2, 0.025]
            elif orientation == 'vertical':
                cbarLocation = [-0.1, 0.2, 0.015, 0.6]
            im = plt.gca().get_images()[1]
            cb = plt.colorbar(im, cax=plt.axes(cbarLocation),
                              orientation=orientation, format='%.2f')
            cb.set_ticks([cb._values.min(), cb._values.max()])

        osl.frame_axes.figure.savefig(
            opj('figures', foldername, outfile + '.%s' % suffix),
            dpi=dpiRes, bbox_inches='tight', transparent=True)
        plt.close()


def create_output(filePath, cluster_extend, threshold, template, create_CSV,
                  show_cross, annotate_figure, show_colorbar,
                  colorbar_orientation, show_title, dpi, imageType, prefix):

    img = nb.load(filePath)
    data = img.get_data()
    affine = img.get_affine()
    labels, nlabels = get_labels(data, threshold=threshold,
                                 min_extent=cluster_extend)

    # Catch if nlabels is 0, i.e. no clusters survived thresholding
    if nlabels == 0:
        print 'No clusters survive the thresholds in %s' % filePath
    else:
        print '%s clusters were found in %s' % (nlabels, filePath)
        coords, clusterInfo = get_cluster_info(labels, affine, data)

        # Get file prefix
        if filePath.endswith('.nii'):
            filename = opb(filePath)[:-4]
        elif filePath.endswith('.nii.gz'):
            filename = opb(filePath)[:-7]

        # Create output folder
        if not os.path.exists(opj('figures', filename)):
            os.makedirs(opj('figures', filename))

        # Create figures
        show_slices(
            data, affine, coords, cmap=pylab.cm.hot, dpiRes=dpi, prefix=prefix,
            show_colorbar=show_colorbar, showCross=show_cross,
            suffix=imageType, orientation=colorbar_orientation,
            threshold=threshold, annotate=annotate_figure,
            showTitle=show_title, foldername=filename,)

        # Create CSV output
        if create_CSV:
            header = 'X,Y,Z,Size,Max,Min,Mean,Std'
            np.savetxt(
                opj('figures', filename, 'cluster_info.csv'), clusterInfo,
                delimiter=',', fmt='%.8f', header=header, comments='')

            # Print cluster info in terminal
            row_format = "{:>8}{:>8}{:>8}{:>10}{:>16}{:>16}{:>16}{:>16}"
            print row_format.format(
                *['X', 'Y', 'Z', 'Size', 'Max', 'Min', 'Mean', 'Std'])
            for c in clusterInfo:
                print row_format.format(*c)
            print '\n'


if __name__ == "__main__":

    cluster_extend = int(sys.argv[1])
    threshold = float(sys.argv[2])
    template = str(sys.argv[3])
    create_CSV = bool(sys.argv[4])
    show_cross = bool(sys.argv[5])
    annotate_figure = bool(sys.argv[6])
    show_colorbar = bool(sys.argv[7])
    colorbar_orientation = str(sys.argv[8])
    show_title = bool(sys.argv[9])
    dpi = int(sys.argv[10])
    imageType = str(sys.argv[11])
    prefix = str(sys.argv[12])

    # Go through all the files in the data folder
    fileList = gg('data/*')
    for fpath in fileList:
        create_output(
            fpath, cluster_extend, threshold, template, create_CSV, show_cross,
            annotate_figure, show_colorbar, colorbar_orientation, show_title,
            dpi, imageType, prefix)
