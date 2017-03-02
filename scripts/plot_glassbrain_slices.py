import os
import sys
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from glob import glob as gg
from os.path import join as opj
from os.path import basename as opb
from nilearn.image import smooth_img
from nilearn.plotting import plot_roi
from nilearn.plotting.find_cuts import find_cut_slices
from nibabel.nifti1 import Nifti1Image


def getEqualSpacing(dirMin, dirMax, ortho, nCuts):
    """
    Computes cut coordinates with equal spacing in a given direction
    """

    if ortho == 'x':
        idx = 0
    elif ortho == 'y':
        idx = 1
    elif ortho == 'z':
        idx = 2

    sign = -1.0 * np.sign(dirMin[idx])
    stepsize = sign * int(np.abs(dirMin[idx] - dirMax[idx]) / (nCuts + 1))
    cut_order = np.arange(dirMin[idx], dirMax[idx], stepsize)
    if len(cut_order) == nCuts + 2:
        cut_coords = cut_order[1:-1]
    else:
        cut_coords = np.delete(cut_order, np.argmax(np.abs(cut_order)))

    return cut_coords


def plotGlassbrainSlices(niftipath, mnipath, ortho='z', nRows=2, nCuts=6,
                         threshpos=0, threshneg=0, figLayout='Both',
                         showLRannot=True, findOptimalCut=True,
                         imageType='svg'):
    """
    Creates nice glassbrain slice figures in the direction x, y and z
    """

    # Initiation of relevant parameters
    img = nb.load(niftipath)
    lineW = 2. / (nRows + int((figLayout == 'Brain' or figLayout == 'Both')))

    # Reduce 4D volume to 3D
    if len(img.shape) == 4:
        data4D = img.get_data()
        data4D = data4D.reshape(data4D.shape[:-1])
        img = Nifti1Image(data4D, img.get_affine())

    # Get voxel extend in all directions
    dirMin = np.dot(img.get_affine(), [0, 0, 0, 1])[:3]
    dirMax = np.dot(img.get_affine(),
                    np.array(img.shape).tolist() + [1])[:3]

    if findOptimalCut:
        # Find cuts automatically
        cut_coords = find_cut_slices(img, direction=ortho, n_cuts=nCuts)
    else:
        # Split orientation in x-equal parts
        cut_coords = getEqualSpacing(dirMin, dirMax, ortho, nCuts)

    # Split cuts according nRows
    cut_coords = [cut_coords[int(i * len(cut_coords) / np.float(nRows)):
                             int((i + 1) * len(cut_coords) / np.float(nRows))]
                  for i in range(nRows)]

    # Create Slices
    for i in range(nRows):

        # Create axes for plotting
        ax = plt.subplot(nRows + int((figLayout == 'Brain' or
                                      figLayout == 'Both')),
                         1, i + 1)

        # Plot the white background for all slices as a zeros value brain
        # (without it, the view focuses around the first area plotted)
        zerobrain = Nifti1Image(img.get_data() * 0, img.get_affine())
        brain = plot_roi(
            zerobrain, zerobrain, colorbar=False, cut_coords=cut_coords[i],
            display_mode=ortho, alpha=1, draw_cross=False, cmap=plt.cm.gray,
            black_bg=False, axes=ax, annotate=False)

        # Plot positive values
        posdata = np.copy(img.get_data())
        posdata[posdata <= threshpos] = 0.001  # = 0 crashes contour function
        posbrain = Nifti1Image(posdata, img.get_affine())
        brain.add_contours(
            posbrain, filled=False, cmap=plt.cm.hot, alpha=1, linewidths=lineW)

        # Plot negative values
        negdata = np.copy(img.get_data())
        negdata[negdata >= -threshneg] = 0.001  # = 0 crashes contour function
        negbrain = Nifti1Image(negdata, img.get_affine())
        brain.add_contours(
            negbrain, filled=False, cmap=plt.cm.winter, alpha=1,
            linewidths=lineW)

        # Plot outer MNI contours
        brain.add_contours(
            smooth_img(mnipath, 4), alpha=1, filled=False,
            levels=[100], linewidths=lineW, cmap=plt.cm.gray)

        # Plot inner MNI contours
        brain.add_contours(
            nb.load(mnipath), alpha=0.8, levels=[5000], linewidths=lineW,
            cmap=plt.cm.gray)

        # Add annotation if requested
        if figLayout == 'Both' or figLayout == 'Number':
            brain.annotate(left_right=showLRannot, size=int(12 * lineW))

    # Plot overview Brain at the bottom
    if figLayout == 'Brain' or figLayout == 'Both':

        # Create axes for overview brain
        ax = plt.subplot(nRows + 1, 1, nRows + 1)

        # Find overview view direction
        if ortho == 'z':
            direction = 'x'
        elif ortho == 'x':
            direction = 'z'
        elif ortho == 'y':
            direction = 'z'

        # Plot the white backgroundas a zeros value brain
        brain = plot_roi(
            zerobrain, zerobrain, colorbar=False, cut_coords=[0],
            display_mode=direction, alpha=1, draw_cross=False,
            cmap=plt.cm.gray, black_bg=False, axes=ax, annotate=False)

        # Plot positive values
        brain.add_contours(
            posbrain, filled=False, cmap=plt.cm.hot, alpha=1, linewidths=lineW)

        # Plot negative values
        brain.add_contours(
            negbrain, filled=False, cmap=plt.cm.winter, alpha=1,
            linewidths=lineW)

        # Plot outer MNI contours
        brain.add_contours(
            smooth_img(mnipath, 4), alpha=1, filled=False,
            levels=[100], linewidths=lineW, cmap=plt.cm.gray)

        # Plot inner MNI contours
        brain.add_contours(
            nb.load(mnipath), alpha=0.8, levels=[5000], linewidths=lineW,
            cmap=plt.cm.gray)

        # Plot the line indicating the cut
        for i in np.array(cut_coords).flatten():
            if ortho == 'z' or ortho == 'y':
                ax.plot([-100, 100], [i, i], 'k-', lw=lineW)
            elif ortho == 'x':
                ax.plot([i, i], [-100, 100], 'k-', lw=lineW)

        if ortho == 'z':
            ax.axis((-300.0, 300.0, dirMin[2], dirMax[2]))
        elif ortho == 'y':
            ax.axis((-300.0, 300.0, dirMin[1], dirMax[1]))
        elif ortho == 'x':
            stretcher = (nRows + 1) / 2.
            ax.axis((-300.0 * stretcher, 300.0 * stretcher, -100.0, 100.0))

        # Add annotation if requested
        if figLayout == 'Both' or figLayout == 'Number':
            brain.annotate(left_right=showLRannot, size=int(12 * lineW))

    # Get file prefix
    if niftipath.endswith('.nii'):
        filename = opb(niftipath)[:-4]
    elif niftipath.endswith('.nii.gz'):
        filename = opb(niftipath)[:-7]

    # Create output folder
    path2Figure = opj(os.path.split(os.path.realpath(niftipath))[0], 'figures')
    if not os.path.exists(opj(path2Figure)):
        os.makedirs(opj(path2Figure))

    # Save figure
    figname = '_'.join([filename, '%s-cut' % ortho])
    plt.savefig(opj(path2Figure, '%s.%s' % (figname, imageType)))
    plt.clf()


if __name__ == "__main__":

    niftipath = str(sys.argv[1])
    mnipath = str(sys.argv[2])
    ortho = str(sys.argv[3])
    nRows = int(sys.argv[4])
    nCuts = int(sys.argv[5])
    showLRannot = bool(int(sys.argv[6]))
    figLayout = str(sys.argv[7])
    threshpos = int(sys.argv[8])
    threshneg = int(sys.argv[9])
    findOptimalCut = bool(int(sys.argv[10]))
    imageType = str(sys.argv[11])

    # Go through all the files in the data folder if requested
    if niftipath == 'data':
        fileList = gg('data/*.nii*')
        for fpath in fileList:
            for o in list(ortho):
                plotGlassbrainSlices(
                    fpath, mnipath, o, nRows, nCuts, threshpos, threshneg,
                    figLayout, showLRannot, findOptimalCut, imageType)
    else:
        for o in list(ortho):
            plotGlassbrainSlices(
                niftipath, mnipath, o, nRows, nCuts, threshpos, threshneg,
                figLayout, showLRannot, findOptimalCut, imageType)
