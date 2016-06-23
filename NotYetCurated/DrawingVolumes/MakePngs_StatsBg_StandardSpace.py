#!/usr/bin/env python

'''
This code allows you to make png images of a stats image on top of a background
image for all the slices in a 3D volume (sagittal, axial and coronal).

positional arguments:
  bg_fname              File name for background .nii.gz image
  stats_fname           File name for stats .nii.gz image
  output_dirname        Output directory for .png images

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Print verbose updates of each step to the screen
  -cm colormap, --colormap colormap
                        Colormap used to plot stats data. Default is autumn.
                        See wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
                        for all possible colormaps
  -tc1 textcolor_mni, --textcolor_mni textcolor_mni
                        Color for coordinate text information. Default is
                        black. Enter "none" to leave blank
  -tc2 textcolor_R, --textcolor_R textcolor_R
                        Color for side indicator (R) on axial slices. Default
                        is black. Enter "none" to leave blank
  -tr, --transparency   Make background transparent. Default is black
  -c crop_option, --crop_option crop_option
                        Choose to crop either the background image, or the
                        stats image. Default is background

Created on: 29th August 2013
Created by: Kirstie Whitaker
   Contact:  kw401@cam.ac.uk
'''

#==============================================================================
# IMPORT WHAT YOU NEED
import numpy as np
import matplotlib.pylab as plt
import matplotlib as mpl
import os
import sys
import nibabel as nib
import argparse

#==============================================================================
# DEFINE THE FUNCTIONS YOU NEED
def setup_argparser():
    '''
    # CODE TO READ ARGUMENTS FROM THE COMMAND LINE AND SET OPTIONS
    # ALSO INCLUDES SOME HELP TEXT
    '''
    
    # Build a basic parser.
    help_text = ('Overlay a stats nifti image on top of a standard space '
                    + 'background nifti image (eg: mean_FA, FMRIB58_FA_1mm, '
					+ 'MNI152lin_T1_2mm_brain) and create pngs of all slices '
					+ 'in all directions (axial, coronal, sagittal)')
    
    sign_off = 'Author: Kirstie Whitaker <kw401@cam.ac.uk>'
    
    parser = argparse.ArgumentParser(description=help_text, epilog=sign_off)
    
    # Now add the arguments
    # Required argument: background file
    parser.add_argument(dest='background_file', 
                            type=str,
                            metavar='bg_fname',
                            help='File name for background .nii.gz image')
    
    # Required argument: stats file
    parser.add_argument(dest='stats_file', 
                            type=str,
                            metavar='stats_fname',
                            help='File name for stats .nii.gz image')
    
    # Required argument: output dir
    parser.add_argument(dest='output_dir', 
                            type=str,
                            metavar='output_dirname',
                            help='Output directory for .png images')
    
    # Optional argument: verbosity
    parser.add_argument('-v', '--verbose',
                            dest='verbose',
                            action='store_true',
                            help='Print verbose updates of each step to the screen')
                            
    # Optional argument: use_specificcolors
    parser.add_argument('-oc', '--use_specificcolors',
                            dest='use_specificcolors',
                            action='store_true',
                            help='Ignore the colormap and generate one from a list of colors')
    
    # Optional argument: colormap
    #       default: autumn
    parser.add_argument('-cm', '--colormap',
                            dest='colormap',
                            type=str,
                            default='autumn',
                            metavar='colormap',
                            help= ('Colormap used to plot stats data. Default is autumn. '
                                    + 'See wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps '
                                    + 'for all possible colormaps') )
                                    
    # Optional argument: colorlist
    #       default: blue
    parser.add_argument('-c', '--colorlist',
                            dest='colorlist',
                            type=str,
                            nargs='+',
                            default='blue',
                            metavar='colorlist',
                            help= ('One color for every integer value in the stats file.'
                                    + 'Default is blue.') )
    
    # Optional argument: textcolor_mni
    #       default: black
    parser.add_argument('-tc1', '--textcolor_mni',
                            dest='textcolor_mni',
                            type=str,
                            default='k',
                            metavar='textcolor_mni',
                            help='Color for coordinate text information. Default is black. Enter "none" to leave blank')
    
    # Optional argument: textcolor_R
    #       default: black
    parser.add_argument('-tc2', '--textcolor_R',
                            dest='textcolor_R',
                            type=str,
                            default='k',
                            metavar='textcolor_R',
                            help='Color for side indicator (R) on axial slices. Default is black. Enter "none" to leave blank')
    
    # Optional argument: transparency
    #       default: False
    parser.add_argument('-tr', '--transparency',
                            dest='transparency',
                            action='store_true',
                            help='Make background transparent. Default is black')
    
    # Optional argument: crop_option
    #       default: background
    parser.add_argument('-cr', '--crop_option',
                            dest='crop_option',
                            type=str,
                            default='background',
                            choices=['background', 'stats'],
                            metavar='crop_option',
                            help='Choose to crop either the background image, or the stats image. Default is background')
    
    arguments = parser.parse_args()
    
    return arguments, parser

def hardcoded_variables():
    # The xyz_dict is a dictionary that allows us to loop through the dimensions
    # and label them with meaningful words
    xyz_dict = { 'name': ('sagittal', 'coronal', 'axial'),
                'letter': ('X', 'Y', 'Z'),
				'mni_const': (90, 126, 72) }
    # Name and letter are self explanatory. The mni_const is a list of values
	# that can be added to the slice number in order to convert the voxel value
	# to mni coordinates

	# The xyz_dict will actually be added to with important values from the
	# data later on
    
    # Now define the functions that convert the voxel values to mni space
    def mni_x(constant, slice_n, zoom):
        x = constant - (slice_n * zoom)
        return x
    
    def mni_y(constant, slice_n, zoom):
        y = (slice_n * zoom) - constant
        return y
    
    def mni_z(constant, slice_n, zoom):
        z = (slice_n * zoom) - constant
        return z

    mni_func_list = [ mni_x, mni_y, mni_z ]

    return xyz_dict, mni_func_list

def load_data(arguments, parser):
    '''
    READ IN THE DATA
    and check that the two files are the same size
    
    '''
    try:
        bg_img = nib.load(arguments.background_file)
        bg_data = bg_img.get_data()
    except:
        print '\n************************'
        print 'ERROR: background file can not be loaded\n'
        parser.print_help()
        sys.exit()
    
    try:
        stats_img = nib.load(arguments.stats_file)
        stats_data = stats_img.get_data()
    except:
        print '\n************************'
        print 'ERROR: stats file can not be loaded\n'
        parser.print_help()
        sys.exit()
    
    # Check that they're the same shape and have the same zoom dimensions
    if not bg_data.shape == stats_data.shape:
        print '\n************************'
        print 'ERROR: files are not the same shape\n'
        parser.print_help()
        sys.exit()
    
    if not bg_img.get_header().get_zooms() == stats_img.get_header().get_zooms():
        print '\n************************'
        print 'ERROR: files do not have the same voxel dimensions\n'
        parser.print_help()
        sys.exit()
    
    # Make sure all data is float:
    bg_data = bg_data/1.
    stats_data = stats_data/1.

    # Scale the background_data by its maximum
    bg_data = bg_data / bg_data.max()

    # And get rid of the bottom 6% of values in the background image
    bg_data[bg_data < 0.06] = 0
    stats_data[bg_data < 0.06] = 0
    
    # Now also save the pixel dimensions
    zooms = bg_img.get_header().get_zooms()

    return bg_data, stats_data, zooms

def create_colormap(arguments):
    # make a color map of fixed colors given in the arguments
    cmap = mpl.colors.ListedColormap(arguments.colorlist)
    bounds = np.linspace(0,1,len(arguments.colorlist)+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    arguments.colormap = cmap
    return arguments, cmap

def create_test_data():
    bg = np.tile(np.arange(40), (40, 40, 1)).T
    bg = bg + np.arange(40)
    bg = bg[:20,:,:30]
    
    bg[:3,:,:] = 0
    bg[(-2):,:,:] = 0
    bg[:,:11,:] = 0
    bg[:,(-17):,:] = 0
    bg[:,:6,:] = 0
    bg[:,(-10):,:] = 0
    
    bg = bg/1.
    bg = bg/np.max(bg)

    stats = np.zeros_like(bg)
    stats = stats/1.
    stats[10:15, 10:15, 10:15] = np.random.random([5,5,5])
    
    return bg, stats

def crop_data(bg, stats):
    '''
    Crop the data to get ride of large amounts of black space surrounding the
    background image.
    '''
    #---------------------------------------------------------------
    # First find all the slices that contain data you want
    slices_list_x = list(np.argwhere(np.sum(bg, (1,2)) != 0)[:,0])
    slices_list_y = list(np.argwhere(np.sum(bg, (0,2)) != 0)[:,0])
    slices_list_z = list(np.argwhere(np.sum(bg, (0,1)) != 0)[:,0])

    slices_list = [slices_list_x, slices_list_y, slices_list_z]
    
    #---------------------------------------------------------------
    # Make a copy of the data
    bg_cropped = np.copy(bg)
    stats_cropped = np.copy(stats)
    
    #---------------------------------------------------------------
    # Remove all slices that have no data in the background image
    bg_cropped = bg_cropped[ slices_list_x, :, : ]
    stats_cropped = stats_cropped[ slices_list_x, :, : ]
    
    bg_cropped = bg_cropped[ :, slices_list_y, : ]
    stats_cropped = stats_cropped[ :, slices_list_y, : ]
    
    bg_cropped = bg_cropped[ :, :, slices_list_z ]
    stats_cropped = stats_cropped[ :, :, slices_list_z ]
    
    #---------------------------------------------------------------
    # Pad with zeros for 5 slices on all sides if you're going
    # to get all the slices in the background image
    bg_cropped = np.pad(bg_cropped, 5, mode='constant')
    stats_cropped = np.pad(stats_cropped, 5, mode='constant')
    
    # Add these slices into the slices_list
    for i, sl in enumerate(slices_list):
        sl = sl + [ n - 5 for n in sl ]
        sl = sl + [ n + 5 for n in sl ]
        slices_list[i] = sorted(list(set(sl)))
    
    return bg_cropped, stats_cropped, slices_list

def stats_only(stats, slices_list):
    
    slices_list_x = list(np.argwhere(np.sum(stats, (1,2))!=0)[:,0])
    slices_list_y = list(np.argwhere(np.sum(stats, (0,2))!=0)[:,0])
    slices_list_z = list(np.argwhere(np.sum(stats, (0,1))!=0)[:,0])

    slices_list = [slices_list_x, slices_list_y, slices_list_z]
    
    return slices_list

def rotate_data(bg, stats, slices_list, axis_name, shape):
    # Rotate the data as required
    # Return the rotated data, and an updated slice list if necessary
    if axis_name == 'axial':
        # Align so that right is right
        stats = np.rot90(stats)
        stats = np.fliplr(stats)
        bg = np.rot90(bg)
        bg = np.fliplr(bg)
    
    elif axis_name == 'coronal':
        stats = np.rot90(stats)
        bg = np.rot90(bg)
        stats = np.flipud(np.swapaxes(stats, 0, 2))
        bg = np.flipud(np.swapaxes(bg, 0, 2))
        slices_list[1] = [ shape - n for n in slices_list[1] ] 
        
    elif axis_name == 'sagittal':
        stats = np.flipud(np.swapaxes(stats, 0, 2))
        bg = np.flipud(np.swapaxes(bg, 0, 2))
    
    else:
        print '\n************************'
        print 'ERROR: data could not be rotated\n'
        parser.print_help()
        sys.exit()
    
    return bg, stats, slices_list

def make_png(bg_slice, stats_slice,
                        axis_name, mni_text, 
                        png_name, arguments):
    '''
    Makes a png image from two slices overlayed, writes the 
    mni coordinate in the bottom left corner, and saves to the
    output directory 
    The transparency option makes the background transparent
    or not
    The colormap option is the color that the statistics will
    be presented in (over a gray background)
    '''

    # Set up the figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # If you don't want a transparent background, add in a
    # black background first
    if not arguments.transparency:
        black = ax.imshow(np.ones_like(bg_slice),
                                interpolation='none',
                                cmap='gray')
    # Mask the data
    m_stats_slice = np.ma.masked_where(stats_slice==0, stats_slice)
    m_bg_slice = np.ma.masked_where(bg_slice < 0.06, bg_slice)

    # First show the background slice
    im1 = ax.imshow(m_bg_slice,
                        interpolation='none',
                        cmap='gray',
                        vmin = 0,
                        vmax = 1)
    
    # Then overlay the stats_slice
    im2 = ax.imshow(m_stats_slice,
                        interpolation='none',
                        cmap=arguments.colormap,
                        alpha=0.75,
                        vmin = 0,
                        vmax = 1)
                        
    # Add a black line around the edge of the background image
    # it makes the brain look nicer :)
    CS = plt.contour(bg_slice, [0.06, 1], linewidths=3, colors='k')

    # Put a text box in the bottom center of the picture with the
    # slice number in MNI space
    text = ax.text(0.5, 0.01 , mni_text,
        horizontalalignment='center',
        verticalalignment='bottom',
        transform=ax.transAxes,
        color = arguments.textcolor_mni)
    
    # Put a little "R" in the middle right side of the image 
    # if you're making axial slices
    if axis_name == 'axial':
        text = ax.text(0.99, 0.5 , 'R',
        horizontalalignment='right',
        verticalalignment='center',
        transform=ax.transAxes,
        color = arguments.textcolor_R)
    
    # Turn off axis labels
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_frame_on(False)

    # Save the figure
    plt.savefig(os.path.join(arguments.output_dir, png_name),
                    transparent=True,
                    bbox_inches='tight',
                    edgecolor='none')
    
    plt.close()

#==============================================================================
# NOW THE FUN BEGINS

arguments, parser = setup_argparser() # Read arguments from command line

if not os.path.isdir(arguments.output_dir):     # Make the output directory if
    os.makedirs(arguments.output_dir)             # it doesn't already exist
    
xyz_dict, mni_func_list = hardcoded_variables() # Define some of the hardcoded
                                                 # variables

bg, stats, zooms = load_data(arguments, parser) # Load data

#bg, stats = create_test_data() # Use test data

xyz_dict['shape'] = bg.shape # Add shape into your xyz_dict

xyz_dict['zooms'] = zooms # Add voxel dimensions to your xyz_dict

bg_cropped, stats_cropped, slices_list = crop_data(bg, stats)
                                              # Crop data (but keep slice_ids)

if arguments.use_specificcolors:
    arguments, cmap = create_colormap(arguments) # Update the colormap if necessary    

# Loop through the three axes
for axis_id in range(3):  
    
    axis_name = xyz_dict['name'][axis_id] # Get the name of the axis
    print axis_name.capitalize()          # and print to screen

    shape = xyz_dict['shape'][axis_id] # You also need the shape of the array
    
    bg, stats, slices_list = rotate_data(bg_cropped, # Rotate the data so your
                                    stats_cropped,   # axis of interest is last 
                                    slices_list,     # and the slices look good
                                    axis_name,
                                    shape)
    
    # Loop through the slices
    for slice_id, slice_n in enumerate(slices_list[axis_id]):

        mni = mni_func_list[axis_id](xyz_dict['mni_const'][axis_id], slice_n, xyz_dict['zooms'][axis_id])

        mni_text = '{} = {}'.format(xyz_dict['letter'][axis_id], mni)
        
        png_name = '{}_slice_{:04.0f}_{:+04.0f}.png'.format(axis_name, slice_id, mni)
        
        if arguments.crop_option == 'stats':
            if np.sum(stats[:,:,slice_id]) > 0:
        
                make_png(bg[:,:,slice_id],       # Make the image ONLY from slices
                            stats[:,:,slice_id], # that have stats data and save 
                            axis_name,           # in output_directory
                            mni_text, 
                            png_name,
                            arguments)
        else:
            make_png(bg[:,:,slice_id],       # Make the image from each slice
                        stats[:,:,slice_id], # and save in output_directory 
                        axis_name,
                        mni_text, 
                        png_name,
                        arguments)

# THE END
