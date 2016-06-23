#!/usr/bin/env python

#=============================================================================
# Created by Kirstie Whitaker
# September 2014
# Contact: kw401@cam.ac.uk
#=============================================================================

#=============================================================================
# IMPORTS
#=============================================================================
import os
import sys
import argparse
import textwrap

import numpy as np

import nibabel as nib
from nilearn import plotting

#=============================================================================
# FUNCTIONS
#=============================================================================
def setup_argparser():
    '''
    Code to read in arguments from the command line
    Also allows you to change some settings
    '''
    # Build a basic parser.
    help_text = ('Plot a statistical volume in MNI space within a glass brain')

    sign_off = 'Author: Kirstie Whitaker <kw401@cam.ac.uk>'

    parser = argparse.ArgumentParser(description=help_text,
                                     epilog=sign_off,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # Now add the arguments
    parser.add_argument(dest='stats_file',
                            type=str,
                            metavar='stats_file',
                            help='Nifti file in MNI space containing the \n  statistical values you want to visualise')

    parser.add_argument('--display_mode',
                            type=str,
                            metavar='display_mode',
                            help=textwrap.dedent("The layout for the glass brain projections.\nOptions are:\n  x: sagittal\n  y: coronal\n  z: axial\n  l: sagittal left hemisphere only\n  r: sagittal right hemisphere only\n  ortho: three cuts in orthogonal directions.\nYou can combine them too, for example:\n  xz: sagittal & axial views\n  lr: left and right sagittal views\n  lyrz: left hemi sagittal, coronal,\n        right hemi sagittal & axial views\nDefault: ortho"),
                            default='ortho')

    parser.add_argument('-c', '--cmap',
                            type=str,
                            metavar='cmap',
                            help=textwrap.dedent('Any matplotlib colormap listed at\n  http://matplotlib.org/examples/color/colormaps_reference.html\nDefault: RdBu_r'),
                            default='RdBu_r')

    parser.add_argument('-cb', '--cbar',
                            action='store_true',
                            help=textwrap.dedent('Display a colorbar on the right of the plots\nDefault: False'),
                            default=False)

    parser.add_argument('--black_bg',
                            action='store_true',
                            help=textwrap.dedent('Set the background to black.\nDefault: White'),
                            default=False)

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

    parser.add_argument('--dpi',
                            type=float,
                            metavar='dpi',
                            help='DPI of output png file\nDefault: 300',
                            default=300)

    arguments = parser.parse_args()

    return arguments, parser


#=============================================================================
# SET SOME VARIABLES
#=============================================================================
# Read in the arguments from argparse
arguments, parser = setup_argparser()

stats_file = arguments.stats_file
display_mode = arguments.display_mode
cmap = arguments.cmap
colorbar = arguments.cbar
thr_abs = arguments.thr_abs
thr_pos = arguments.thr_pos
black_bg = arguments.black_bg
lower_thresh = arguments.lower
upper_thresh = arguments.upper
dpi = arguments.dpi

# Set a couple of hard coded options
symmetric_cbar=False

#===============================================================================
# Sort out the thresholding
#===============================================================================
# Read in the data using nibabel
img = nib.load(stats_file)
data = img.get_data()

# If you want to show both positive and negative values
# on either side of a threshold value then you want to use
# thr_abs.
if thr_abs is not None:

    threshold = thr_abs

    # We want to set the upper and lower boundaries of the colormap
    # to be centered around 0.
    symmetric_cbar = True

    # If the absolute value of the upper threshold is larger than the
    # absolute value of the lower threshold then the lower threshold
    # should be changed to -1 * upper threshold. If the absolute value
    # of the lower threshold is larger than the absolute value of the
    # upper threshold then we should set the upper threshold to
    # -1 * lower_threshold.
    if np.abs(lower_thresh) < np.abs(upper_thresh):
        lower_thresh = np.abs(upper_thresh) * -1
    else:
        upper_thresh = np.abs(lower_thresh) * -1

# If you want to only show positive values (or values larger than a certain
# threshold value then you want to use thr_pos.
elif thr_pos is not None:

    threshold = thr_pos

    # We need to set all the data that we don't want to see to equal thr_pos
    data[data<thr_pos] = thr_pos

#===============================================================================
# Plot away!
#===============================================================================
slicer = plotting.plot_glass_brain(img,
                                    threshold=threshold,
                                    plot_abs=False,
                                    symmetric_cbar=symmetric_cbar,
                                    vmin=lower_thresh,
                                    vmax=upper_thresh,
                                    cmap=cmap,
                                    colorbar=colorbar,
                                    display_mode=display_mode,
                                    black_bg=black_bg)

#===============================================================================
# Save the figure
#===============================================================================
output_file = os.path.join(stats_file.rsplit('.nii', 1)[0]
                           + '_GlassBrain_{}.png'.format(display_mode))

slicer.savefig(output_file, dpi=dpi)
