#!/usr/bin/env python

#=============================================================================
# Created by Kirstie Whitaker
# September 2014
# Contact: kw401@cam.ac.uk
#=============================================================================

#=============================================================================
# IMPORTS
#=============================================================================
from __future__ import print_function

import os
import sys
import argparse
import numpy as np

import pandas as pd
import nibabel as nib
from surfer import Brain
import seaborn as sns

import itertools as it

import matplotlib.pylab as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec

import matplotlib.colors as mcolors
import matplotlib.cm as cm

#=============================================================================
# FUNCTIONS
#=============================================================================
def setup_argparser():
    '''
    Code to read in arguments from the command line
    Aso allows you to change some settings
    '''
    # Build a basic parser.
    help_text = ('Plot a single value for each region in the NSPN 500 parcellation of the fsaverage surface')

    sign_off = 'Author: Kirstie Whitaker <kw401@cam.ac.uk>'

    parser = argparse.ArgumentParser(description=help_text, epilog=sign_off)

    # Now add the arguments
    parser.add_argument(dest='roi_file',
                            type=str,
                            metavar='roi_file',
                            help='roi file containing list of measure values - one for each region - csv format')

    parser.add_argument('--annot_name',
                            type=str,
                            metavar='annot_name',
                            help='name of annot file that contains the parcellation',
                            default='500.aparc')

    parser.add_argument('--fsaverageid',
                            type=str,
                            metavar='fsaverage_id',
                            help='FSaverage subject id',
                            default='fsaverageSubP')

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

    parser.add_argument('-c2', '--cmap2',
                            type=str,
                            metavar='cmap2',
                            help='colormap for the second overlay',
                            default='autumn')

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

    parser.add_argument('-t2', '--thresh2',
                            type=float,
                            metavar='thresh2',
                            help='mask values below this value for the second color',
                            default=None)

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
                            default='pial')

    parser.add_argument('-cst', '--cortex_style',
                            type=str,
                            metavar='cortex_style',
                            help='cortex style - one of "classic", "bone", "high_contrast" or "low_contrast"',
                            default='classic')

    arguments = parser.parse_args()

    return arguments, parser

#------------------------------------------------------------------------------
def calc_range(roi_data, l=None, u=None, thresh=-99, center=False):
    """
    The calc_range function takes as input:
      * roi_data: a numpy array or column from a pandas data frame
      * lower threshold for color map: a float or integer value, or None
          if None, the lower threshold is calculated from the data
      * upper threshold for color map: a float or integer value, or None
          if None, the upper threshold is calculated from the data
      * threshold value (regions with values below this threhold will not be
          plotted): a float or integer value, default is -99
      * center: a boolean which if true forces 0 to be the center of the color
          bar, default is False

    If center is set to False, and l and u are both provided, then they are
      simply returned, nothing happens!

    If center is set to True and the thresholds are given then the largest one
      (according to its absolute value) is set as one limit, with the other
      calculated as that limit * -1.

    If either of the lower or upper limits are given as None, they are
      calculated by rounding to the closest 1/20th of the data in the correct
      direction.
      For example, if the minimum value is -3.236 it would be rounded to -3.20.
        If the maximum value were 10223.22 it would be rounded to 10224.
      Note, this isn't perfect, and it could be improved by rounding to 3
      decimal places before applying this rounding to the nearest 1/20th. But
      it's not too bad for now.

    If a threshold is given, and the lower limit is not given, that limit
      is calculated as being the smallest value that is greater than the
      threshold value.
      For example, if the minimum value is -183.22 but the threshold is given
        as -99, then the smallest value that is larger than -99 is used.
    """
    # Figure out the min and max for the color bar
    if l == None:
        l = roi_data[roi_data>thresh].min()
        l = np.floor(l*20)/20.0
    if u == None:
        u = roi_data[roi_data>thresh].max()
        u = np.ceil(u*20)/20.0

    if center:
        # Make sure the colorbar is centered
        if l**2 < u **2:
            l = u*-1
        else:
            u = l*-1

    return l, u

def make_colorbar(l, u, cmap, color_file):
    # If you have a custom color_file then get the lower and upper
    # values from that, and create the color map from it too.
    if color_file:
        cmap = [line.strip() for line in open(color_file)]
        l = 1
        u = len(cmap)

        # If you've passed rgb values you need to convert
        # these to tuples
        if len(cmap[0].split()) == 3:
            cmap = [ (np.float(x.split()[0]),
                      np.float(x.split()[1]),
                      np.float(x.split()[2])) for x in cmap ]

    # If you don't have a color_file then calculate the lower and
    # upper limits from the data itself
    else:
        # Set l and u so that they're the same for both hemispheres
        l, u = calc_range(df[0], l, u, thresh, center)

    return l, u, cmap

#------------------------------------------------------------------------------
def plot_surface(vtx_data,
                 subject_id,
                 hemi,
                 surface,
                 subjects_dir,
                 output_dir,
                 prefix,
                 l,
                 u,
                 cmap,
                 thresh,
                 thresh2=None,
                 cmap2='autumn',
                 cortex_style='classic'):
    """
    This function needs more documentation, but for now
    it is sufficient to know this one important fact:
	For the variable "cmap":
	    If you pass a word that defines a matplotlib
	      colormap (eg: jet, Rd_Bu etc) then the code
	      will use that for the color scheme.
	    If you pass a **list** of colors then you'll
	      just loop through those colors instead.
    """
    if cortex_style.count('_') == 2:
        cortex_style_list = cortex_style.split('_')
        cortex_name = cortex_style_list[0]
        cortex_min = np.float(cortex_style_list[1])
        cortex_max = np.float(cortex_style_list[2])

        cortex_style = ( cortex_name, cortex_min, cortex_max, False )

    # Open up a brain in pysurfer
    brain = Brain(subject_id, hemi, surface,
                      subjects_dir = subjects_dir,
                      background="white",
                      size=(800, 665),
                      cortex=cortex_style)

    # Create an empty brain if the values are all below threshold
    if np.max(vtx_data) < thresh:
        # Add your data to the brain
        brain.add_data(vtx_data*0,
                        l,
                        u,
                        thresh = thresh,
                        colormap=cmap,
                        alpha=0.0)

    # If you only have one threshold
    # then add the data!
    elif not thresh2:
        # Add your data to the brain
        brain.add_data(vtx_data,
                        l,
                        u,
                        thresh = thresh,
                        colormap=cmap,
                        alpha=.8)

    else:
        # Plot the data twice for the two
        # different settings
        vtx_data1 = np.copy(vtx_data)
        vtx_data1[vtx_data1>thresh2] = 0
        brain.add_data(vtx_data1,
                        l,
                        u,
                        thresh = thresh,
                        colormap = cmap,
                        alpha = .8)

        brain.add_data(vtx_data,
                        l,
                        u,
                        thresh = thresh2,
                        colormap = cmap2,
                        alpha = .8)

    # Save the images for medial and lateral
    # putting a color bar on all of them
    brain.save_imageset(prefix = os.path.join(output_dir, prefix),
                        views = views_list,
                        colorbar = range(len(views_list)) )

#-----------------------------------------------------------------------------
def combine_pngs(measure, surface, output_dir, cortex_style):
    '''
    Find four images and combine them into one nice picture
    '''
    figsize = (5,4)
    fig = plt.figure(figsize = figsize, facecolor='white')

    grid = gridspec.GridSpec(2, 2)
    grid.update(left=0, right=1, top=1, bottom = 0.08, wspace=0, hspace=0)

    f_list = [ '_'.join([os.path.join(output_dir, measure), 'lh', surface, cortex_style, 'lateral.png']),
               '_'.join([os.path.join(output_dir, measure), 'rh', surface, cortex_style, 'lateral.png']),
               '_'.join([os.path.join(output_dir, measure), 'lh', surface, cortex_style, 'medial.png']),
               '_'.join([os.path.join(output_dir, measure), 'rh', surface, cortex_style, 'medial.png']) ]

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
    filename = os.path.join(output_dir, '{}_{}_{}_combined.png'.format(measure, surface, cortex_style))
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

#===============================================================================

if __name__ == "__main__":

    #---------------------------------------------------------------------------
    # Read in the arguments from argparse
    #---------------------------------------------------------------------------
    arguments, parser = setup_argparser()

    annot_name = arguments.annot_name
    subject_id = arguments.fsaverageid
    subjects_dir = arguments.subjects_dir
    roi_data_file = arguments.roi_file
    l = arguments.lower
    u = arguments.upper
    cmap = arguments.cmap
    cmap2 = arguments.cmap2
    color_file = arguments.color_file
    center = arguments.center
    surface = arguments.surface
    thresh = arguments.thresh
    thresh2 = arguments.thresh2
    cortex_style = arguments.cortex_style

    # Define the output directory
    output_dir = os.path.join(os.path.dirname(roi_data_file), 'PNGS')
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Define the name of the measure you're plotting
    measure = os.path.basename(roi_data_file)
    measure = os.path.splitext(measure)[0]

    # Define the aparc names
    # Read in aparc names file
    aparc_names_file =  os.path.join(subjects_dir,
                                     subject_id, "parcellation",
                                     "{}.names.txt".format(annot_name))

    # Read in the names from the aparc names file
    aparc_names = [line.strip() for line in open(aparc_names_file)]
    # dropping the first 41            #### COMMMENTED OUT - NEED TO UPDATE THIS SO WE DON'T HAVE TO DROP THE FIRST 41...
    ## aparc_names = aparc_names[41::] #### COMMMENTED OUT - NEED TO UPDATE THIS SO WE DON'T HAVE TO DROP THE FIRST 41...

    # Figure out which surfaces you're going to use
    if surface == 'both':
        surface_list = [ "inflated", "pial" ]
    elif surface == 'inflated':
        surface_list = [ "inflated" ]
    elif surface == 'pial':
        surface_list = [ "pial" ]
    else:
        print ("Do not recognise surface. Check {}".format(surface))
        parser.print_help()
        sys.exit()

    hemi_list = [ "lh", "rh" ]
    views_list = [ 'medial', 'lateral' ]

    # Check that the inputs exist:
    if not os.path.isfile(roi_data_file):
        print ("Roi data file doesn't exist")
        sys.exit()

    if not os.path.isdir(os.path.join(subjects_dir, subject_id, "surf")):
        print ("Fsaverage directory doesn't exist")
        print ("Check subjects_dir and subject_id")
        sys.exit()

    #---------------------------------------------------------------------------
    # READ IN THE MEASURE DATA
    #---------------------------------------------------------------------------
    # Read in the data
    df = pd.read_csv(roi_data_file, index_col=False, header=None)

    # Make custom colorbar
    cmap, l, u = make_colorbar(l, u, cmap, color_file)

    # Now rearrange the data frame and match it up with
    # the aparc names
    df = df.T
    df.columns = aparc_names

    # Now make your pictures
    for hemi, surface in it.product(hemi_list, surface_list):

        prefix = '_'.join([measure, hemi, surface, cortex_style])

        # Read in aparc annot file which will be inside
        # the label folder of your fsaverage subject folder
        aparc_file = os.path.join(subjects_dir,
                              subject_id, "label",
                              "{}.{}.annot".format(hemi, annot_name))

        # Use nibabel to merge together the aparc_names and the aparc_file
        labels, ctab, names = nib.freesurfer.read_annot(aparc_file)

        # Create an empty roi_data array
        roi_data = np.ones(len(names))*(thresh-1.0)

        # Loop through the names and if they are in the data frame
        # for this hemisphere then add that value to the roi_data array
        for i, name in enumerate(names):
            roi_name = '{}_{}'.format(hemi, name)

            if roi_name in df.columns:
                roi_data[i] = df[roi_name]

        # Make a vector containing the data point at each vertex.
        vtx_data = roi_data[labels]

        # Write out the vtx_data
        #nib.freesurfer.write_annot(f_name, vtx_data, ctab, names)

        # Show this data on a brain
        plot_surface(vtx_data, subject_id, hemi,
                         surface, subjects_dir,
                         output_dir, prefix,
                         l, u, cmap,
                         thresh,
                         cmap2=cmap2,
                         thresh2=thresh2,
                         cortex_style=cortex_style)

    #=============================================================================
    # COMBINE THE IMAGES
    #=============================================================================
    for surface in surface_list:
        combine_pngs(measure, surface, output_dir, cortex_style)
        brains_in_a_row(measure, surface, output_dir, cortex_style, l, u, cmap)


# You're done :)
# Happy International Women's Day 2017
# Updated on 10th October 2017 from Nanning, China
# <3 Kx
