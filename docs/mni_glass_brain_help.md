usage: mni_glass_brain.py [-h] [--display_mode display_mode] [-c cmap] [-cb]
                          [--black_bg] [--thr_abs thr_abs] [--thr_pos thr_pos]
                          [-l lowerthr] [-u upperthr] [--dpi dpi]
                          stats_file

Plot a statistical volume in MNI space within a glass brain

positional arguments:
  stats_file            Nifti file in MNI space containing the 
                          statistical values you want to visualise

optional arguments:
  -h, --help            show this help message and exit
  --display_mode display_mode
                        The layout for the glass brain projections.
                        Options are:
                          x: sagittal
                          y: coronal
                          z: axial
                          l: sagittal left hemisphere only
                          r: sagittal right hemisphere only
                          ortho: three cuts in orthogonal directions.
                        You can combine them too, for example:
                          xz: sagittal & axial views
                          lr: left and right sagittal views
                          lyrz: left hemi sagittal, coronal,
                                right hemi sagittal & axial views
                        Default: ortho
  -c cmap, --cmap cmap  Any matplotlib colormap listed at
                          http://matplotlib.org/examples/color/colormaps_reference.html
                        Default: RdBu_r
  -cb, --cbar           Display a colorbar on the right of the plots
                        Default: False
  --black_bg            Set the background to black.
                        Default: White
  --thr_abs thr_abs     Mask the input image such that all values
                          which have an absolute value less than this threshold
                          are not shown.
                        If None then no thresholding is undertaken.
                        Default: None
  --thr_pos thr_pos     Mask the input image such that all values
                          less than this threshold are not shown.
                        If None then no thresholding is undertaken.
                        Default: None
  -l lowerthr, --lower lowerthr
                        Lower limit for colormap
  -u upperthr, --upper upperthr
                        Upper limit for colormap
  --dpi dpi             DPI of output png file
                        Default: 300

Author: Kirstie Whitaker <kw401@cam.ac.uk>
