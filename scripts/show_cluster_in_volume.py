#!/usr/bin/env python

#=============================================================================
# Created by Michael Notter
# at OHBM 2016 Brainhack in Lausanne, June 2016
# Edited with more comments by Kirstie Whitaker
# at Cambridge Brainhack-Global 2017, March 2017
# Contact: kw401@cam.ac.uk
#=============================================================================

#=============================================================================
# IMPORTS
#=============================================================================
import argparse
import future        # pip install future
from glob import glob as gg
import os
from os.path import join as opj
from os.path import basename as opb
import sys
import textwrap

import numpy as np
from matplotlib import pylab
from matplotlib import pyplot as plt
import nibabel as nb
import nilearn
from nipy.labs import viz
from scipy.ndimage import label as sci_label

#=============================================================================
# FUNCTIONS
#=============================================================================
def setup_argparser():
    '''
    Code to read in arguments from the command line
    Also allows you to change some settings
    '''
    # Build a basic parser.
    help_text = ('Show the locations of clusters in a statistical map in MNI space.')

    sign_off = 'Author: Kirstie Whitaker <kw401@cam.ac.uk>'

    parser = argparse.ArgumentParser(description=help_text,
                                     epilog=sign_off,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # Now add the arguments
    parser.add_argument(dest='stats_file',
                            type=str,
                            metavar='stats_file',
                            help=textwrap.dedent('3D nifti file in MNI space containing the statistical values\n ' +
                                                  'you want to visualise.\n' +
                                                  'Note that this file can be pre-thresholded or not.\n' +
                                                  'If your file is already thresholded then you will need to\n ' +
                                                  'pass an argument to the -t option otherwise it will default to 2.3.\n ' +
                                                  'A suggested value is 0.01.' ))

    parser.add_argument('-ce, --cluster_extent',
                            type=str,
                            metavar='cluster_extent',
                            help=textwrap.dedent("Minimum cluster extent for a region to be included in the visualisation\n (integer)\n Default: 20"),
                            default=20)

    parser.add_argument('-t, --cluster_thr',
                            type=str,
                            metavar='threshold',
                            help=textwrap.dedent("Minimum statistical value for a region to be included in the visualisation\n (float)\n Default: 2.3"),
                            default=2.3)

    parser.add_argument('--csv',
                            action='store_true',
                            help=textwrap.dedent('Create a csv file with cluster information.\n Default: False'),
                            default=False)

    parser.add_argument('--cluster_title',
                            action='store_true',
                            help=textwrap.dedent('Show cluster information in the title of the plot.\n Default: False'),
                            default=False)

    parser.add_argument('-c', '--cmap',
                            type=str,
                            metavar='cmap',
                            help=textwrap.dedent('Any matplotlib colormap listed at\n  http://matplotlib.org/examples/color/colormaps_reference.html\n Default: RdBu_r'),
                            default='hot')

    parser.add_argument('-cb', '--cbar',
                            action='store_true',
                            help=textwrap.dedent('Display a colorbar on the right of the plots\n Default: False'),
                            default=False)

    parser.add_argument('--black_bg',
                            action='store_true',
                            help=textwrap.dedent('Set the background to black.\n Default: White'),
                            default=False)
    """
    parser.add_argument('--thr_abs',
                            type=float,
                            metavar='thr_abs',
                            help=textwrap.dedent('Mask the input image such that all values\n  which have an absolute value less than this threshold\n  are not shown.\nIf None then no thresholding is undertaken.\nDefault: None'),
                            default=None)

    parser.add_argument('--thr_pos',
                            type=float,
                            metavar='thr_pos',
                            help=textwrap.dedent('Mask the input image such that all values\n  less than this threshold are not shown.\nIf None then no thresholding is undertaken.\nDefault: None'),
                            default=None)
    parser.add_argument('-l', '--lower',
                            type=float,
                            metavar='lowerthr',
                            help='Lower limit for colormap',
                            default=None)

    parser.add_argument('-u', '--upper',
                            type=float,
                            metavar='upperthr',
                            help='Upper limit for colormap',
                            default=None)
    """

    parser.add_argument('--dpi',
                            type=float,
                            metavar='dpi',
                            help='DPI of output png file\n Default: 300',
                            default=300)

    parser.add_argument('--format',
                            type=float,
                            metavar='format',
                            help=textwrap.dedent('Format of the output image file.\n Eg: png, pdf, tif, jpeg, svg. \n Default: png'),
                            default='png')

    arguments = parser.parse_args()

    return arguments, parser


def get_labels(data, cluster_thr=0, min_extent=0):
    """
    Get number of clusters in dataset as well as a labeled volume
    Minimal extent of each cluster and voxel-vise threshold can be specified
    """

    # Threshold the data by zeroing all voxels with values that have an
    # absolute value less than the cluster threshold
    thr_data = abs(data) > cluster_thr

    # Find all the clusters in the thresholded data
    labels, nlabels = sci_label(thr_data)

    # Now loop through all the clusters
    # and if a cluster is smaller than the minimum cluster extent
    # exclude it from the list (set the values to zero)
    for idx in range(1, nlabels + 1):
        if np.sum(labels == idx) < min_extent:
            labels[labels == idx] = 0

    # Re-run the clustering command to get only the clusters
    # that are larger than the minimum extent
    labels, nlabels = sci_label(labels)

    # overwrites the input data with the thresholded data
    binarized_data = labels.astype('bool')
    data[np.invert(binarized_data)] = 0

    return labels, nlabels, data, binarized_data



def get_cluster_info(img, affine, data):
    """
    Returns peak coordinations and cluster information of a given dataset,
    if labeled file and affine is provided
    """

    # Set up some variables we're going to need
    coords = []          #
    cs = []              # cluster sum values
    maxcoords = []       # peak coordinate locations
    clusterInfo = []     # list of lists containing max, min,
                         # mean and std of the cluster

    # Find all the label ids (that aren't 0!)
    labelID = np.setdiff1d(np.unique(img.ravel()), [0])

    # Loop through the labels
    for lab in labelID:
        # Calculate the voume of the cluster
        sumval = np.sum(img == lab)
        cs.append(sumval)

        # Calculate the max, min, mean and std of the cluster
        maxval = np.max(data[img == lab])
        minval = np.min(data[img == lab])
        meanval = np.mean(data[img == lab])
        stdval = np.std(data[img == lab])

        # Save these values in a list
        clusterInfo.append([sumval, maxval, minval, meanval, stdval])

        # Get the location of the peak coordinate
        maxidx = np.nonzero(np.multiply(data, img == lab) == maxval)
        maxcoords.append([m[0] for m in maxidx])

    # Transform the lists into numpy arrays
    maxcoords = np.asarray(maxcoords)
    clusterInfo = np.asarray(clusterInfo)

    # Sort the lists by the volume of the clusters
    maxcoords = maxcoords[np.argsort(cs)[::-1], :]
    clusterInfo = clusterInfo[np.argsort(cs)[::-1], :]

    # Loop through the clusters and put the peak coordinates
    # in MNI space
    for i, lab in enumerate(labelID[np.argsort(cs)[::-1]]):
        coords.append(np.dot(affine,
                             np.hstack((maxcoords[i], 1)))[:3].tolist())

    # Add the peak coordinate information to the clusterInfo array
    clusterInfo = np.hstack((np.array(coords), clusterInfo))

    # Returns peak coordination and additional cluster infos
    return coords, clusterInfo

def show_slices(data, affine,
                    coords=None,
                    cmap=None,
                    show_colorbar=None,
                    showCross=False,
                    cluster_thr=0,
                    annotate=True,                              ###### KW DOCUMENT
                    template='../scripts/templates/MNI152_T1_1mm_brain.nii.gz', ####### KW DOCUMENT
                    dpiRes=300,
                    suffix='png',
                    show_title=False):

    # Prepare background image
    anatimg = nb.load(template)
    anatdata, anataff = anatimg.get_data(), anatimg.affine()
    anatdata = anatdata.astype(np.float)
    anatdata[anatdata < 10.] = np.nan

    # Create output figure for each peak coordinate
    # (so a different figure for each cluster)
    for idx, coord in enumerate(coords):

        # Name the output file to include the cluster id,
        # the cluster threshold and the minimum cluster extent
        outfile = 'Cluster_{}_thr{:04.2f}_minext{:03:0f}'.format(idx, cluster_thr, cluster_extent)

        # If show_title argument has been set to true then print the file name
        # and the peak coordinates in the title of the figure
        if show_title:
            title = '{} {}'.format(outfile + coord)
        else:
            title = ''

        # Woooo plot three orthogonal views of the cluster sliced through the
        # peak coordinate
        osl = viz.plot_map(
                    np.asarray(data), affine, anat=anatdata, anat_affine=anataff,
                    threshold=cluster_thr, cmap=cmap, annotate=annotate,
                    black_bg=False, cut_coords=coord, draw_cross=showCross,
                    slicer='ortho', title=title)

        # If the show colorbar option is true then show the color bar on the
        # right hand side of the image
        if show_colorbar:
            cbarLocation = [-0.1, 0.2, 0.015, 0.6]
            im = plt.gca().get_images()[1]
            cb = plt.colorbar(im, cax=plt.axes(cbarLocation),
                              orientation='horizontal', format='%.2f')
            cb.set_ticks([cb._values.min(), cb._values.max()])

        # Save the figure!
        osl.frame_axes.figure.savefig(
            opj(output_folder, '{}.{}'.format(outfile, suffix)),
            dpi=dpiRes, bbox_inches='tight', transparent=True)

        # DONE! Close the plot
        plt.close()


#=============================================================================
# SET SOME VARIABLES
#=============================================================================
# Read in the arguments from argparse
arguments, parser = setup_argparser()

stats_file = arguments.stats_file
cluster_extent = arguments.cluster_extent
cluster_thr = arguments.cluster_thr
store_csv = arguments.csv
cluster_title = arguments.cluster_title
cmap = arguments.cmap
show_colorbar = arguments.cbar
#thr_abs = arguments.thr_abs
#thr_pos = arguments.thr_pos
#black_bg = arguments.black_bg
#lower_thresh = arguments.lower
#upper_thresh = arguments.upper
dpi = arguments.dpi
image_format = arguments.format

# Set a couple of hard coded options
#symmetric_cbar=False

#===============================================================================
# Get the colormap from nilearn
#===============================================================================
if hasattr(cm, cmap):
    cmap = getattr(cm, cmap)

#===============================================================================
# Create the output folder
#===============================================================================
output_folder = '{}_CLUSTERS'.format(stats_file.rsplit('.nii', 1)[0])

if not os.isdir(output_folder):
    os.path.makedirs(output_folder)


def create_output(stats_file, cluster_extent, threshold, template, create_CSV,
                  show_cross, annotate_figure, cmap, show_colorbar,
                  show_title, dpi, imageType):

    # Read in the stats file
    img = nb.load(stats_file)
    data = img.get_data()
    affine = img.affine()

    # Find the clusters
    labels, nlabels, data, binarized_data = get_labels(data,
                                                       cluster_thr=cluster_thr,
                                                       min_extent=cluster_extent)

    # Catch if nlabels is 0, i.e. no clusters survived thresholding
    if nlabels == 0:
        print('No clusters survive the thresholds in {}'.format(stats_file))
        return

    # If there *are* cluster though, then get the cluster information
    # for each of them
    print('{} clusters were found in {}'.format(nlabels, stats_file))
    coords, clusterInfo = get_cluster_info(labels, affine, data)

    """
    # Get file prefix
    if filePath.endswith('.nii'):
        filename = opb(filePath)[:-4]
    elif filePath.endswith('.nii.gz'):
        filename = opb(filePath)[:-7]
    """

    # Create output folder
    output_folder = '{}_CLUSTERS'.format(stats_file.rsplit('.nii', 1)[0])
    if not os.isdir(output_folder):
        os.path.makedirs(output_folder)

    # Create figures
    show_slices(data, affine,
                        coords=coords,
                        cmap=cmap,
                        show_colorbar=show_colorbar,
                        showCross=False,                            ####### KW THINK ABOUT
                        cluster_thr=cluster_thr,
                        annotate=True,                              ###### KW DOCUMENT
                        template='../scripts/templates/MNI152_T1_1mm_brain.nii.gz', ####### KW DOCUMENT
                        dpiRes=dpi,
                        suffix=image_format,
                        show_title=show_title)

    # Create CSV output
    if create_CSV:
        header = 'X,Y,Z,Size,Max,Min,Mean,Std'
        np.savetxt(
            opj('figures', filename, 'cluster_info.csv'), clusterInfo,
            delimiter=',', fmt='%.8f', header=header, comments='')

        # Print cluster info in terminal
        row_format = "{:>8}{:>8}{:>8}{:>10}{:>16}{:>16}{:>16}{:>16}"
        print(row_format.format(
            *['X', 'Y', 'Z', 'Size', 'Max', 'Min', 'Mean', 'Std']))
        for c in clusterInfo:
            print(row_format.format(*c))
        print('\n')

#===============================================================================
# Save the figure
#===============================================================================
output_folder = '{}_CLUSTERS'.format(stats_file.rsplit('.nii', 1)[0])

if not os.isdir(output_folder):
    os.path.makedirs(output_folder)


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

    #=========================================================================
    # SET SOME VARIABLES
    #=========================================================================
    # Read in the arguments from argparse
    arguments, parser = setup_argparser()

    stats_file = arguments.stats_file
    cluster_extent = arguments.cluster_extent
    cluster_thr = arguments.cluster_thr
    store_csv = arguments.csv
    cluster_title = arguments.cluster_title
    cmap = arguments.cmap
    show_colorbar = arguments.cbar
    #thr_abs = arguments.thr_abs
    #thr_pos = arguments.thr_pos
    #black_bg = arguments.black_bg
    #lower_thresh = arguments.lower
    #upper_thresh = arguments.upper
    dpi = arguments.dpi
    image_format = arguments.format


    #===============================================================================
    # Get the colormap from nilearn
    #===============================================================================
    if hasattr(cm, cmap):
        cmap = getattr(cm, cmap)


    #===============================================================================
    # Create the output folder
    #===============================================================================
    output_folder = '{}_CLUSTERS'.format(stats_file.rsplit('.nii', 1)[0])

    if not os.isdir(output_folder):
        os.path.makedirs(output_folder)

    #===============================================================================
    # Create the figures and CSV output
    #===============================================================================
    create_output(stats_file, cluster_extent, threshold, template, create_CSV,
                      show_cross, annotate_figure, cmap, show_colorbar,
                      show_title, dpi, imageType)
