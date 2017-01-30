```
usage: pysurfer_plot_500parcellation_surface_values.py [-h]
                                                       [--fsaverageid fsaverage_id]
                                                       [-sd subjects_dir]
                                                       [-c cmap] [-c2 cmap2]
                                                       [-cf color_file]
                                                       [--center] [-t thresh]
                                                       [-t2 thresh2]
                                                       [-l lowerthr]
                                                       [-u upperthr]
                                                       [-s surface]
                                                       [-cst cortex_style]
                                                       roi_file

Plot a single value for each region in the NSPN 500 parcellation of the 
fsaverage surface

positional arguments:
  roi_file              roi file containing list of measure values - one for
                        each region - csv format

optional arguments:
  -h, --help            show this help message and exit
  --fsaverageid fsaverage_id
                        FSaverage subject id
  -sd subjects_dir, --subjects_dir subjects_dir
                        freesurfer subjects dir
  -c cmap, --cmap cmap  colormap
  -c2 cmap2, --cmap2 cmap2
                        colormap for the second overlay
  -cf color_file, --color_file color_file
                        file containing list of custom colors
  --center              center the color bar around 0
  -t thresh, --thresh thresh
                        mask values below this value
  -t2 thresh2, --thresh2 thresh2
                        mask values below this value for the second color
  -l lowerthr, --lower lowerthr
                        lower limit for colorbar
  -u upperthr, --upper upperthr
                        upper limit for colorbar
  -s surface, --surface surface
                        surface - one of "pial", "inflated" or "both"
  -cst cortex_style, --cortex_style cortex_style
                        cortex style - one of "classic", "bone",
                        "high_contrast" or "low_contrast"

Author: Kirstie Whitaker <kw401@cam.ac.uk>
```
