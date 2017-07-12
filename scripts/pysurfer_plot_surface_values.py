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
import numpy as np

from surfer import Brain, io

import itertools as it

import matplotlib.pylab as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec


#=============================================================================
# FUNCTIONS
#=============================================================================
def setup_argparser():
    '''
    Code to read in arguments from the command line
    Aso allows you to change some settings
    '''
    # Build a basic parser.
    help_text = ('Plot values on a freesurfer surface')

    sign_off = 'Author: Kirstie Whitaker <kw401@cam.ac.uk>'

    parser = argparse.ArgumentParser(description=help_text, epilog=sign_off)

    # Now add the arguments
    parser.add_argument(dest='overlay_file',
                            type=str,
                            metavar='overlay_file',
                            help='overlay file in with the hemisphere in the first two letters')

    parser.add_argument(dest='output_dir',
                            type=str,
                            metavar='output_dir',
                            help='output directory')

    parser.add_argument('--subject_id',
                            type=str,
                            metavar='subject id',
                            help='freesurfer subject id',
                            default='fsaverage')

    parser.add_argument('-sd', '--subjects_dir',
                            type=str,
                            metavar='subjects_dir',
                            help='freesurfer subjects dir',
                            default=os.environ["SUBJECTS_DIR"])

    parser.add_argument('-c', '--cmap',
                            type=str,
                            metavar='cmap',
                            help='colormap',
                            default='RdBu_r')

    parser.add_argument('-cf', '--color_file',
                            type=str,
                            metavar='color_file',
                            help='file containing list of custom colors',
                            default=None)

    parser.add_argument('--center',
                            action='store_true',
                            help='center the color bar around 0')

    parser.add_argument('-t', '--thresh',
                            type=float,
                            metavar='thresh',
                            help='mask values below this value',
                            default=-98)

    parser.add_argument('-m', '--mask',
                            type=float,
                            metavar='mask',
                            help='mask values that are exactly this value',
                            default=0)

    parser.add_argument('-l', '--lower',
                            type=float,
                            metavar='lowerthr',
                            help='lower limit for colorbar',
                            default=None)

    parser.add_argument('-u', '--upper',
                            type=float,
                            metavar='upperthr',
                            help='upper limit for colorbar',
                            default=None)


    parser.add_argument('-s', '--surface',
                            type=str,
                            metavar='surface',
                            help='surface - one of "pial", "inflated" or "both"',
                            default='both')

    arguments = parser.parse_args()

    return arguments, parser

#------------------------------------------------------------------------------
def mask_vtx_data(overlay_fname, cortex_fname, thresh):

    vtx_data = io.read_scalar_data(overlay_fname)
    cortex_data = io.read_label(cortex_fname)

    # Create a mask of 1s where there is cortex and 0s on the medial wall
    mask = np.zeros_like(vtx_data)
    mask[cortex_data] = 1

    # Set all values that are not in cortex to thresh-1
    vtx_data[mask == 0] = thresh-1

    return vtx_data

#------------------------------------------------------------------------------
def calc_range(vtx_data_left, vtx_data_right, thresh, l, u):
    '''
    This is an important step to ensure that the colorbar is exactly
    the same for the right and left hemispheres.
    '''
    if l == None:
        # Figure out the min and max for each hemisphere
        l_l = vtx_data_left[vtx_data_left>=thresh].min()
        l_r = vtx_data_right[vtx_data_right>=thresh].min()

        # Take the smallest of these two
        l = np.min([l_l, l_r])

        # And round to a nice number
        l = np.floor(l*20)/20.0

    if u == None:
        # Figure out the min and max for each hemisphere
        u_l = vtx_data_left[vtx_data_left>=thresh].max()
        u_r = vtx_data_right[vtx_data_right>=thresh].max()

        # Take the largest of these two
        u = np.max([u_l, u_r])

        # And round to a nice number
        u = np.ceil(u*20)/20.0

    # Return the lower and upper bounds
    return l, u

#------------------------------------------------------------------------------
def plot_surface(vtx_data, subject_id, subjects_dir, hemi, surface, output_dir, prefix, l, u, cmap, center, thresh):
    # Open up a brain in pysurfer
    brain = Brain(subject_id, hemi, surface,
                  subjects_dir = subjects_dir,
                  config_opts=dict(background="white",
                                   height=665,
                                   width=800))

    if center:
        # Make sure the colorbar is centered
        if l**2 < u **2:
            l = u*-1
        else:
            u = l*-1

    # Create an empty brain if the values are all below threshold
    if np.max(vtx_data) < thresh:
        # Add your data to the brain
        brain.add_data(vtx_data*0,
                        l,
                        u,
                        thresh = thresh,
                        colormap=cmap,
                        alpha=0.0)

    # Otherwise, add the data appropriately!
    else:
        # Add your data to the brain
        brain.add_data(vtx_data,
                        l,
                        u,
                        thresh = thresh,
                        colormap=cmap,
                        alpha=.8)

    # Save the images for medial and lateral
    # putting a color bar on all of them
    brain.save_imageset(prefix = os.path.join(output_dir, prefix),
                        views = views_list,
                        colorbar = range(len(views_list)) )

#-----------------------------------------------------------------------------
def combine_pngs(surface, output_dir):
    '''
    Find four images and combine them into one nice picture
    '''
    figsize = (5,4)
    fig = plt.figure(figsize = figsize, facecolor='white')

    grid = gridspec.GridSpec(2, 2)
    grid.update(left=0, right=1, top=1, bottom = 0.08, wspace=0, hspace=0)

    f_list = [ os.path.join(output_dir, '_'.join(['lh', surface, 'lateral.png'])),
               os.path.join(output_dir, '_'.join(['rh', surface, 'lateral.png'])),
               os.path.join(output_dir, '_'.join(['lh', surface, 'medial.png'])),
               os.path.join(output_dir, '_'.join(['rh', surface, 'medial.png'])) ]

    # Plot each figure in turn
    for g_loc, f in zip(grid, f_list):
        ax = plt.Subplot(fig, g_loc)
        fig.add_subplot(ax)
        img = mpimg.imread(f)
        # Crop the figures appropriately
        # NOTE: this can change depending on which system you've made the
        # images on originally - it's a bug that needs to be sorted out!
        if 'lateral' in f:
            img_cropped = img[75:589,55:(-50),:]
        else:
            img_cropped = img[45:600,25:(-25),:]
        ax.imshow(img_cropped, interpolation='none')
        ax.set_axis_off()

    # Add the bottom of one of the images as the color bar
    # at the bottom of the combo figure
    grid_cbar = gridspec.GridSpec(1,1)
    grid_cbar.update(left=0, right=1, top=0.08, bottom=0, wspace=0, hspace=0)
    ax = plt.Subplot(fig, grid_cbar[0])
    fig.add_subplot(ax)
    img = mpimg.imread(f)
    img_cbar = img[600:,:]
    ax.imshow(img_cbar, interpolation='none')
    ax.set_axis_off()

    # Save the figure
    filename = os.path.join(output_dir, '{}_combined.png'.format(surface))
    print filename
    fig.savefig(filename, bbox_inches=0, dpi=300)


def add_four_hor_brains(grid, f_list, fig):
    '''
    Take the four pysurfer views (left lateral, left medial,
    right medial and right lateral) and arrange them in a row
    according to the grid positions given by grid

    grid    :  the gridspec list of grid placements
    f_list  :  list of four file pysurfer image files
    big_fig :  the figure to which you're adding the images

    # THIS WAS UPDATED TO INCLUDE PLOTTING IN A GRID
    # Should probably change the function name!
    '''
    for g_loc, f in zip(grid, f_list):
        img = mpimg.imread(f)
        # Crop the figures appropriately
        # NOTE: this can change depending on which system you've made the
        # images on originally - it's a bug that needs to be sorted out!
        if 'lateral' in f:
            img_cropped = img[115:564, 105:(-100),:]
        else:
            img_cropped = img[90:560, 60:(-55),:]

        # Add an axis to the figure
        ax_brain = plt.Subplot(fig, g_loc)
        fig.add_subplot(ax_brain)

        # Show the brain on this axis
        ax_brain.imshow(img_cropped, interpolation='none')
        ax_brain.set_axis_off()

    return fig

def add_colorbar(grid, big_fig, cmap_name, y_min=0, y_max=1, cbar_min=0, cbar_max=1, vert=False, label=None, show_ticks=True, pad=0):
    # Set the seaborn context and style
    sns.set(style="white")
    sns.set_context("poster", font_scale=3)
    '''
    Add a colorbar to the big_fig in the location defined by grid

    grid       :  grid spec location to add colormap
    big_fig    :  figure to which colorbar will be added
    cmap_name  :  name of the colormap
    x_min      :  the minimum value to plot this colorbar between
    x_max      :  the maximum value to plot this colorbar between
    cbar_min   :  minimum value for the colormap (default 0)
    cbar_max   :  maximum value for the colormap (default 1)
    vert       :  whether the colorbar should be vertical (default False)
    label      :  the label for the colorbar (default: None)
    ticks      :  whether to put the tick values on the colorbar (default: True)
    pad        :  how much to shift the colorbar label by (default: 0)
    '''
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap

    # Add an axis to the big_fig
    ax_cbar = plt.Subplot(big_fig, grid)
    big_fig.add_subplot(ax_cbar)

    # Normalise the colorbar so you have the correct upper and
    # lower limits and define the three ticks you want to show
    norm = mpl.colors.Normalize(vmin=cbar_min, vmax=cbar_max)

    if show_ticks:
        ticks = [y_min, np.average([y_min, y_max]), y_max]
    else:
        ticks=[]

    # Figure out the orientation
    if vert:
        orientation='vertical'
        rotation=270
    else:
        orientation='horizontal'
        rotation=0

    # Add in your colorbar:
    cb = mpl.colorbar.ColorbarBase(ax_cbar,
                                       cmap=cmap_name,
                                       norm=norm,
                                       orientation=orientation,
                                       ticks=ticks,
                                       boundaries=np.linspace(y_min, y_max, 300))

    if label:
        cb.set_label(label, rotation=rotation, labelpad=pad)

    return big_fig


def brains_in_a_row(measure, surface, output_dir, cortex_style, l, u, cmap):

    # Set up the figure
    fig, ax = plt.subplots(figsize=(20,6), facecolor='white')

    # Set up the grid
    grid = gridspec.GridSpec(1,4)
    grid.update(left=0.01, right=0.99, top=1.05, bottom=0.2, wspace=0, hspace=0)

    # Set up the file list
    f_list = [ '_'.join([os.path.join(output_dir, measure), 'lh', surface, cortex_style, 'lateral.png']),
               '_'.join([os.path.join(output_dir, measure), 'lh', surface, cortex_style, 'medial.png']),
               '_'.join([os.path.join(output_dir, measure), 'rh', surface, cortex_style, 'medial.png']),
               '_'.join([os.path.join(output_dir, measure), 'rh', surface, cortex_style, 'lateral.png']) ]

    # Add the brains
    fig = add_four_hor_brains(grid, f_list, fig)

    # Set up the colorbar grid
    cb_grid = gridspec.GridSpec(1,1)

    cb_grid.update(left=0.2,
                        right=0.8,
                        bottom=0.2,
                        top=0.25,
                        wspace=0,
                        hspace=0)

    fig = add_colorbar(cb_grid[0], fig,
                            cmap_name=cmap,
                            cbar_min=l,
                            cbar_max=u,
                            y_min=l,
                            y_max=u,
                            label='')

    # Turn off the axis
    ax.set_axis_off()

    # Save the figure
    filename = os.path.join(output_dir, '{}_{}_{}_FourHorBrains.png'.format(measure, surface, cortex_style))
    fig.savefig(filename, dpi=300)

    # Close the figure
    plt.close('all')


#=============================================================================
# SET SOME VARIABLES
#=============================================================================
# Read in the arguments from argparse
arguments, parser = setup_argparser()

overlay_file = arguments.overlay_file
output_dir = arguments.output_dir
subject_id = arguments.subject_id
subjects_dir = arguments.subjects_dir
l = arguments.lower
u = arguments.upper
cmap = arguments.cmap
color_file = arguments.color_file
center = arguments.center
surface = arguments.surface
thresh = arguments.thresh
mask = arguments.mask

if surface == 'both':
    surface_list = [ "inflated", "pial" ]
elif surface == 'inflated':
    surface_list = [ "inflated" ]
elif surface == 'pial':
    surface_list = [ "pial" ]
else:
    print "Do not recognise surface. Check {}".format(surface)
    parser.print_help()
    sys.exit()

hemi_list = [ "lh", "rh" ]
views_list = [ 'medial', 'lateral' ]

# Check that the inputs exist:
for hemi in hemi_list:
    #f = os.path.join(os.path.dirname(overlay_file), hemi + os.path.basename(overlay_file)[2:])
    f = overlay_file.replace('lh.', '{}.'.format(hemi)).replace('rh.', '{}.'.format(hemi))

    if not os.path.isfile(f):
        print "{} overlay file doesn't exist".format(hemi)
        sys.exit()

# Make the output directory if it doesn't already exist
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)


#=============================================================================
# READ IN THE VERTEX DATA
#=============================================================================

# Create a vertex data dictionary that will contain the data
# for left and right hemispheres
vtx_data_dict = {}

# Read in the left and right hemisphere surfaces
for hemi in hemi_list:

    # Define the name for the overlay surface file
    overlay_fname = overlay_file.replace('lh.', '{}.'.format(hemi)).replace('rh.', '{}.'.format(hemi))

    # Define the name for the cortex label file
    cortex_fname = os.path.join(subjects_dir, subject_id, 'label', hemi + '.cortex.label')

    # Read the data in and mask it so that non-cortex is -99
    vtx_data_dict[hemi] = mask_vtx_data(overlay_fname, cortex_fname, thresh)

#=============================================================================
# CALCULATE THE COLOR BAR RANGE
#=============================================================================
# Calculate the lower and upper values if they haven't been defined:
l, u = calc_range(vtx_data_dict['lh'], vtx_data_dict['rh'], thresh, l, u)

# Unless there's a given color file
if color_file:
    colors = [line.strip() for line in open(color_file)]
    l = 1
    u = len(colors)

#=============================================================================
# MAKE THE INDIVIDUAL PICTURES
#=============================================================================
for hemi, surface in it.product(hemi_list, surface_list):

    prefix = '_'.join([hemi, surface])

    # Show this data on a brain
    if colors:
        plot_surface(vtx_data_dict[hemi], subject_id, subjects_dir,
                     hemi, surface,
                     output_dir, prefix,
                     l, u, colors, center,
                     thresh)
    else:
        plot_surface(vtx_data_dict[hemi], subject_id, subjects_dir,
                     hemi, surface,
                     output_dir, prefix,
                     l, u, cmap, center,
                     thresh)

#=============================================================================
# COMBINE THE IMAGES
#=============================================================================
for surface in surface_list:
    combine_pngs(surface, output_dir)
    brains_in_a_row(measure, surface, output_dir, cortex_style, l, u, cmap)

# You're done :)
# Happy International Women's Day 2017
# <3 Kx
