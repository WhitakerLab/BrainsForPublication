#!/usr/bin/env python

#=============================================================================
# Created by Kirstie Whitaker
# at OHBM 2016 Brainhack in Lausanne, June 2016
# Contact: kw401@cam.ac.uk
#=============================================================================

#=============================================================================
# IMPORTS
#=============================================================================
import os
import sys
import argparse
import textwrap
from glob import glob

import numpy as np

import nibabel as nib
from nilearn import plotting
from nilearn.plotting import cm
from nilearn.image import reorder_img
from nilearn.image.resampling import coord_transform

import imageio

#=============================================================================
# FUNCTIONS
#=============================================================================
def setup_argparser():
    '''
    Code to read in arguments from the command line
    Also allows you to change some settings
    '''
    # Build a basic parser.
    help_text = ('Plot an anatomical volume in subject space,\noptionally overlay another image in the same space,\nand make into a gif')

    sign_off = 'Author: Kirstie Whitaker <kw401@cam.ac.uk>'

    parser = argparse.ArgumentParser(description=help_text,
                                     epilog=sign_off,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # Now add the arguments
    parser.add_argument(dest='anat_file',
                            type=str,
                            metavar='anat_file',
                            help='Nifti or mgz file in subject space that you want to visualise')

    parser.add_argument('-o', '--overlay_file',
                            type=str,
                            metavar='overlay_file',
                            help='Nifti or mgz file in subject space that you want to overlay',
                            default=None)

    parser.add_argument('-a,', '--axis',
                            type=str,
                            metavar='axis',
                            help=textwrap.dedent("The axis you'd like to project.\nOptions are:\n  x: sagittal\n  y: coronal\n  z: axial\n\nDefault: ortho"),
                            default='x')

    parser.add_argument('-c', '--cmap',
                            type=str,
                            metavar='cmap',
                            help=textwrap.dedent('Any matplotlib colormap listed at\n  http://matplotlib.org/examples/color/colormaps_reference.html\nDefault: gray'),
                            default='gray')

    parser.add_argument('-oc', '--overlay_cmap',
                            type=str,
                            metavar='overlay_cmap',
                            help=textwrap.dedent('Any matplotlib colormap listed at\n  http://matplotlib.org/examples/color/colormaps_reference.html\nDefault: prism'),
                            default='prism')

    parser.add_argument('--black_bg',
                            action='store_true',
                            help=textwrap.dedent('Set the background to black.\nDefault: White'),
                            default=False)

    parser.add_argument('--annotate',
                            action='store_true',
                            help=textwrap.dedent('Add L and R labels to images.\nDefault: False'),
                            default=False)

    parser.add_argument('-t', '--thr',
                            type=float,
                            metavar='thr',
                            help=textwrap.dedent('Mask the input image such that all values\n  which have an absolute value less than this threshold\n  are not shown.\nIf None then no thresholding is undertaken.\nDefault: None'),
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

anat_file = arguments.anat_file
overlay_file = arguments.overlay_file
axis = arguments.axis
cmap = arguments.cmap
overlay_cmap = arguments.overlay_cmap
threshold = arguments.thr
black_bg = arguments.black_bg
annotate = arguments.annotate
dpi = arguments.dpi

# Set a couple of hard coded options
draw_cross = False

#===============================================================================
# Make a bunch of dictionaries that allow you to loop through x, y and z
#===============================================================================

# The x, y, z coord_transform dictionary contains keys that
# are either 'x', 'y', 'z' and values that are functions to
# convert that axis to alligned space.
def coord_transform_x(x, img):
    x, y, z = coord_transform(x, 0, 0, img.affine)
    return x
def coord_transform_y(y, img):
    x, y, z = coord_transform(0, y, 0, img.affine)
    return y
def coord_transform_z(z, img):
    x, y, z = coord_transform(0, 0, z, img.affine)
    return z

coord_transform_dict = { 'x' : coord_transform_x,
                         'y' : coord_transform_y,
                         'z' : coord_transform_z }

# The x, y, z dim lookup dictionary contains keys that
# are either 'x', 'y', 'z' and values that correspond to
# the axis of the image.
# For example, 'x' is the first axis of the image: 0
dim_lookup_dict = { 'x' : 0,
                    'y' : 1,
                    'z' : 2 }

# The x, y, z label lookup dictionary contains keys that
# are either 'x', 'y', 'z' and values that correspond to
# the name of the projection.
label_lookup_dict = { 'x' : 'sagittal',
                      'y' : 'coronal',
                      'z' : 'axial' }

#===============================================================================
# Get the colormap from nilearn
#===============================================================================
if hasattr(cm, cmap):
    cmap = getattr(cm, cmap)

#===============================================================================
# Set up the output directory
#===============================================================================
if '.mgz' in anat_file:
    pngs_dir = anat_file.rsplit('.nii', 1)[0] + '_PNGS'
else:
    pngs_dir = anat_file.rsplit('.mgz', 1)[0] + '_PNGS'

if not os.path.isdir(pngs_dir):
    os.makedirs(pngs_dir)

#===============================================================================
# Read in the images using nibabel
#===============================================================================
img = nib.load(anat_file)

# Convert the data to float
data = img.get_data()
data = data.astype('float')

# Reslice the image so there are no rotations in the affine.
# This step is actually included in the nilearn plot_anat command below
# but it runs faster if the image has already been resliced.
img_reslice = reorder_img(img, resample='continuous')

# Do the same if you have an overlay file too
if not overlay_file is None:
    overlay_img = nib.load(overlay_file)
    data = overlay_img.get_data()
    data = data.astype('float')
    overlay_img_reslice = reorder_img(overlayimg, resample='nearest')

#===============================================================================
# Plot each slice unless it's empty!
#===============================================================================

# Loop through all the slices in this dimension
for i in np.arange(img_reslice.shape[dim_lookup_dict['x']], dtype='float'):

    # Test to see if there is anything worth showing in the image

    # Get the co-ordinate you want
    coord = coord_transform_dict[axis](i, img_reslice)

    # Make the image
    slicer = plotting.plot_anat(img_reslice,
                                threshold=None,
                                cmap=cmap,
                                display_mode=axis,
                                black_bg=black_bg,
                                annotate=annotate,
                                draw_cross=draw_cross,
                                cut_coords=(coord,))

    # Add the overlay if given
    if not overlay_file is None:
        slicer.add_overlay(overlay_img_reslice, cmap=overlay_cmap)

    # Save the png file
    output_file = os.path.join(pngs_dir,
                       '{}_{:03.0f}.png'.format(label_lookup_dict[axis], i))

    slicer.savefig(output_file, dpi=dpi)

    slicer.close()

#===============================================================================
# Now make a gif!
#===============================================================================
gif_file = (anat_file.rsplit('.nii', 1)[0]
                + '_{}.gif'.format(label_lookup_dict[axis]))

png_list = glob(os.path.join(pngs_dir,
                    '{}*.png'.format(label_lookup_dict[axis])))


png_list.sort()

with imageio.get_writer(gif_file, mode='I') as writer:
    for fname in png_list:
        image = imageio.imread(fname)
        writer.append_data(image)

#===============================================================================
# WAY TO GO!
#===============================================================================
