#!/usr/bin/env python

'''
This code allows you to make png images of a brain extracted image on
top of the whole image in the background for all the slices in a 3D volume
(sagittal, axial and coronal).

positional arguments:
  bg_fname              File name for background .nii.gz image
  overlay_fname  File name for overlay .nii.gz image
  output_dirname        Output directory for .png images

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Print verbose updates of each step to the screen
  -cm1 colormap, --colormap1 colormap
                        Colormap used to plot background data. Default is gray.
                        See wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
                        for all possible colormaps
  -cm2 colormap, --colormap2 colormap
                        Colormap used to plot overlay data. Default is autumn.
                        See wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
                        for all possible colormaps
  -tc textcolor_R, --textcolor_R textcolor_R
                        Color for side indicator (R) on axial slices. Default
                        is black. Enter "none" to leave blank
  -tr, --transparency   Make background transparent. Default is black

Created on: 11th December 2013
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
    help_text = ('Overlay a one nifti image (of the same size) on top of a '
                    + 'background nifti image (eg: highre_brain and highres) '
                    + 'and create pngs of all slices in all directions '
					+ '(axial, coronal, sagittal)')
    
    sign_off = 'Author: Kirstie Whitaker <kw401@cam.ac.uk>'
    
    parser = argparse.ArgumentParser(description=help_text, epilog=sign_off)
    
    # Now add the arguments
    # Required argument: background file
    parser.add_argument(dest='background_file', 
                            type=str,
                            metavar='bg_fname',
                            help='File name for background .nii.gz image')
    
    # Required argument: overlay file
    parser.add_argument(dest='overlay_file', 
                            type=str,
                            metavar='overlay_fname',
                            help='File name for overlay .nii.gz image')
    
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
                            
    # Optional argument: colormap1
    #       default: autumn
    parser.add_argument('-cm1', '--colormap1',
                            dest='colormap1',
                            type=str,
                            default='gray',
                            metavar='colormap1',
                            help= ('Colormap used to plot background data. Default is gray. '
                                    + 'See wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps '
                                    + 'for all possible colormaps') )
                                    
    # Optional argument: colormap2
    #       default: autumn
    parser.add_argument('-cm2', '--colormap2',
                            dest='colormap2',
                            type=str,
                            default='autumn',
                            metavar='colormap2',
                            help= ('Colormap used to plot overlay data. Default is autumn. '
                                    + 'See wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps '
                                    + 'for all possible colormaps') )
                                    
    # Optional argument: crop_option
    #       default: background
    parser.add_argument('-cr', '--crop_option',
                            dest='crop_option',
                            type=str,
                            default='background',
                            choices=['background', 'overlay'],
                            metavar='crop_option',
                            help='Choose to crop either the background image, or the overlay image. Default is background')

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
    
    # Optional argument: axial
    #       default: False
    parser.add_argument('-ax', '--axial_only',
                            dest='axial',
                            action='store_true',
                            help='Only create axial pngs. Default is false - create all 3 axes')
                           
    # Optional argument: contour
    #       default: True
    parser.add_argument('-co', '--contour',
                            dest='contour',
                            action='store_true',
                            help='Add a contour line around the overlay file. Default is false.')
                           
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
	# to mni coordinates - it isn't actually used for this code but a hangover
	# from a previous reincarnation

	# The xyz_dict will actually be added to with important values from the
	# data later on
    
    return xyz_dict

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
        overlay_img = nib.load(arguments.overlay_file)
        overlay_data = overlay_img.get_data()
    except:
        print '\n************************'
        print 'ERROR: overlay file can not be loaded\n'
        parser.print_help()
        sys.exit()
    
    # Check that they're the same shape and have the same zoom dimensions
    if not bg_data.shape == overlay_data.shape:
        print '\n************************'
        print 'ERROR: files are not the same shape\n'
        parser.print_help()
        sys.exit()
    
    if not bg_img.get_header().get_zooms() == overlay_img.get_header().get_zooms():
        print '\n************************'
        print 'ERROR: files do not have the same voxel dimensions\n'
        parser.print_help()
        sys.exit()
    
    # Make sure all data is float:
    bg_data = bg_data/1.
    overlay_data = overlay_data/1.
    
    # Scale the data by its maximum
    bg_data = bg_data / bg_data.max()
    overlay_data = overlay_data / overlay_data.max()    
    
    # Now also save the pixel dimensions
    zooms = bg_img.get_header().get_zooms()

    return bg_data, overlay_data, zooms

def crop_data(bg, overlay):
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
    overlay_cropped = np.copy(overlay)
    
    #---------------------------------------------------------------
    # Remove all slices that have no data in the background image
    bg_cropped = bg_cropped[ slices_list_x, :, : ]
    overlay_cropped = overlay_cropped[ slices_list_x, :, : ]
    
    bg_cropped = bg_cropped[ :, slices_list_y, : ]
    overlay_cropped = overlay_cropped[ :, slices_list_y, : ]
    
    bg_cropped = bg_cropped[ :, :, slices_list_z ]
    overlay_cropped = overlay_cropped[ :, :, slices_list_z ]
    
    #---------------------------------------------------------------
    # Pad with zeros for 5 slices on all sides if you're going
    # to get all the slices in the background image
    bg_cropped = np.pad(bg_cropped, 5, mode='constant')
    overlay_cropped = np.pad(overlay_cropped, 5, mode='wrap')
    
    # Add these slices into the slices_list
    for i, sl in enumerate(slices_list):
        sl = sl + [ n - 5 for n in sl ]
        sl = sl + [ n + 5 for n in sl ]
        slices_list[i] = list(set(sl))
    
    return bg_cropped, overlay_cropped, slices_list

def overlay_only(overlay, slices_list):
    
    slices_list_x = list(np.argwhere(np.sum(overlay, (1,2))!=0)[:,0])
    slices_list_y = list(np.argwhere(np.sum(overlay, (0,2))!=0)[:,0])
    slices_list_z = list(np.argwhere(np.sum(overlay, (0,1))!=0)[:,0])

    slices_list = [slices_list_x, slices_list_y, slices_list_z]
    
    return slices_list

def rotate_data(bg, overlay, slices_list, axis_name, shape):
    # Rotate the data as required
    # Return the rotated data, and an updated slice list if necessary
    if axis_name == 'axial':
        # Align so that right is right
        overlay = np.rot90(overlay)
        overlay = np.fliplr(overlay)
        bg = np.rot90(bg)
        bg = np.fliplr(bg)
    
    elif axis_name == 'coronal':
        overlay = np.rot90(overlay)
        bg = np.rot90(bg)
        overlay = np.flipud(np.swapaxes(overlay, 0, 2))
        bg = np.flipud(np.swapaxes(bg, 0, 2))
        slices_list[1] = [ shape - n - 3 for n in slices_list[1] ] 
        
    elif axis_name == 'sagittal':
        overlay = np.flipud(np.swapaxes(overlay, 0, 2))
        bg = np.flipud(np.swapaxes(bg, 0, 2))
    
    else:
        print '\n************************'
        print 'ERROR: data could not be rotated\n'
        parser.print_help()
        sys.exit()
    
    return bg, overlay, slices_list

def make_png(bg_slice, overlay_slice,
                        axis_name, 
                        png_name, arguments):
    '''
    Makes a png image from two slices overlayed and saves to the
    output directory 
    The transparency option makes the background transparent
    or not
    The colormap options (1 and 2) are the colors that the background and
    overlay files will be presented in
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
    m_overlay_slice = np.ma.masked_where(overlay_slice==0, overlay_slice)

    # First show the background slice
    im1 = ax.imshow(bg_slice,
                        interpolation='none',
                        cmap=arguments.colormap1,
                        vmin = 0,
                        vmax = 1)
    
    # Then overlay the overlay_slice
    im2 = ax.imshow(m_overlay_slice,
                        interpolation='none',
                        cmap=arguments.colormap2,
                        vmin = 0,
                        vmax = 1)
               
    if arguments.contour:
        # Add a black line around the edge of the background image
        # it makes the brain look nicer :)
        CS = plt.contour(overlay_slice, [0.01, 1], linewidths=3, colors='k')
         
    # Put a little "R" in the middle right side of the image 
    # if you're making axial slices
    if axis_name == 'axial' and not arguments.textcolor_R == 'none':
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
    
xyz_dict = hardcoded_variables() # Define some of the hardcoded
                                                 # variables

bg, overlay, zooms = load_data(arguments, parser) # Load data

xyz_dict['shape'] = bg.shape # Add shape into your xyz_dict

xyz_dict['zooms'] = zooms # Add voxel dimensions to your xyz_dict

if arguments.crop_option == 'overlay':
    overlay_cropped, bg_cropped, slices_list = crop_data(overlay, bg)
                                              # Crop data (but keep slice_ids)
else:
    bg_cropped, overlay_cropped, slices_list = crop_data(bg, overlay)
                                              # Crop data (but keep slice_ids)
    
# Figure out which axes to make pngs of
if arguments.axial:
    axes_range = range(2,3)
else:
    axes_range = range(3)
    
# Loop through the axes
for axis_id in axes_range:
    
    axis_name = xyz_dict['name'][axis_id] # Get the name of the axis
    print axis_name.capitalize()          # and print to screen

    shape = xyz_dict['shape'][axis_id] # You also need the shape of the array
    
    bg, overlay, slices_list = rotate_data(bg_cropped, # Rotate the data so your
                                    overlay_cropped,   # axis of interest is last 
                                    slices_list,       # and the slices look good
                                    axis_name,
                                    shape)
    
    # Loop through the slices
    for slice_id, slice_n in enumerate(slices_list[axis_id]):
        
        png_name = '{}_slice_{:04.0f}.png'.format(axis_name, slice_id)
                
        if arguments.crop_option == 'overlay':
            if np.sum(overlay[:,:,slice_id]) > 0:
        
                make_png(bg[:,:,slice_id],         # Make the image ONLY from slices
                            overlay[:,:,slice_id], # that have overlay data and save 
                            axis_name,             # in output_directory
                            png_name,
                            arguments)
        else:
            make_png(bg[:,:,slice_id],       # Make the image from each slice
                        overlay[:,:,slice_id], # and save in output_directory 
                        axis_name,
                        png_name,
                        arguments)

# THE END
