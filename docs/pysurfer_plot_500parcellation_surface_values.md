# Surface values for 500 parcellation

The `pysurfer_plot_500parcellation_surface_values.py` script visualises a list of regional measures for the NSPN 500mm<sup>2</sup> freesurfer surface parcellation.


## Required Installations

1. [Freesurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)
2. [anaconda](http://continuum.io/downloads#all) with python version `2.7`
3. [VTK](http://www.vtk.org/VTK/resources/software.html#latestcand)
4. [Mayavi](http://mayavi.sourceforge.net/install.html)

Note that there can be some crazy complicated installation challenges if you try to install these three separately as there will be conflicts in different dependencies. The fix is to let anaconda do its magic by installing all at the same time:
    
    conda create -yn py27-b4p vtk mayavi python=2.7 anaconda

5. [Nibabel](http://nipy.sourceforge.net/nibabel/installation.html#installation) by typing `pip install nibabel`
6. [pysurfer](http://pysurfer.github.io/install.html) by typing `pip install pysurfer`
7. [seaborn](http://seaborn.pydata.org/index.html) `pip install seaborn`

Please note - if you try to google how to install either VTK or mayavi you will end up in a world of pain! Anaconda and pip will be your friend...I hope!

If you did this a little while ago you might need to upgrade the modules with either `conda update x` or `pip install --upgrade x`.

### Getting Kirstie's code

Download the [latest release](https://github.com/KirstieJane/BrainsForPublication/releases) of the BrainsForPublication repository.

Unzip into a location of your choice and use the command `chmod +x scripts/*.py` to ensure you can execute the commands.

You can either enter the whole (or relative) path to the command when you're running the code, or add the folder to your .bashrc file by appending:

    export PATH="<directory_containing_BrainsForPublication_code>:$PATH"

for example:

    export PATH="/home/kw401/DRIVERS/:$PATH"

### Make sure you have the appropriate parcellation files

The fsaverage folder that contains these files (for example `fsaverageSubP`) is shipped with the [BrainsForPublication](https://github.com/KirstieJane/BrainsForPublication/releases) code in a folder called `required_data`.

Specifically the folder contains all the standard fsaverage files along with:

* The aparc names file: `parcellation/500.names.txt`

* The left and right annotation files: `label/lh.500.aparc.annot` & `label/rh.500.aparc.annot`

The two directories `parcellation` and `label` must be inside the `fsaverageSubP` directory.


### Plotting parcellation values for the NSPN (area=500) parcellation

This uses the `pysurfer_plot_500parcellation_surface_values.py` command.

You can find the help information by typing:

`pysurfer_plot_500parcellation_surface_values.py --help`

Which will return:

    usage: pysurfer_plot_500parcellation_surface_values.py [-h] [--fsaverageid fsaverage_id]
                                    [-sd subjects_dir] [-c cmap] [--center]
                                    [-t thresh] [-l lowerthr] [-u upperthr]
                                    [-s surface]
                                    roi_file

    Plot values on a freesurfer surface

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
      --center              center the color bar around 0
      -t thresh, --thresh thresh
                            mask values below this value
      -l lowerthr, --lower lowerthr
                            lower limit for colorbar
      -u upperthr, --upper upperthr
                            upper limit for colorbar
      -s surface, --surface surface
                            surface - one of "pial", "inflated" or "both"

    Author: Kirstie Whitaker <kw401@cam.ac.uk>

An example command may look something like this:

`pysurfer_plot_500parcellation_surface_values.py -u 50 -c hsv -t -98 -s pial --center degree_values.txt`

Let's step through the different options:

**Inputs (positional arguments)**

A *positional argument* means that the command will not run without it, and that it has to be at the end of the optional arguments in the correct order (if there are more than 1).

* The only positional argument is the `roi_file`.

    * In the example above this is called `degree_values.txt`.
    * The roi_file must contain must contain 308 values, one on each line and each corresponding to the 308 regions in the 500 parcellation. Do not put the names in this file, nor the index numbers, but do make sure not to mix them up!

**Optional arguments**

* The first optional argument is `-h` or `--help`
    * Run the command with this flag to see the help printed above

* The `--fsaverageid` argument sets the fsaverage id name.
    * The default value is `fsaverageSubP`

* The `-sd` or `--subjects_dir` argument sets the subjects_dir.
    * This is the folder which contains the fsaverage folder.
    * The default value is `$FREESURFER_HOME/subjects` or whatever has been set as the environmental variable either previously in the terminal, or in your `.bashrc` file.
    * It can either be set here, or as an environmental variable by typing `export SUBJECTS_DIR=<insert_dirname>` in a bash shell.

* The `-c` or `--cmap` argument sets the color map for the values.
    * The default value is the `jet` color map.
    * The example above uses the hsv color map (`-c hsv`)
    * Examples of the matplotlib color maps can be found [here](http://matplotlib.org/examples/color/colormaps_reference.html)

* `-u` or `--upper` and `-l` or `--lower` set the upper and lower boundaries for the color bar
    * By default these values are calculated from the data
    * In the example above the plots will use a maximum value of 50 (`-u 50`). All regions that have a higher value than 50 will be set to 50.

* `--center` will center the color bar around 0, irregardless of the max and min of the data, or the given upper and lower bounds (see below).
    * If the absolute value of the upper limit is larger than the absolute value of the lower limit the limits will be `-upper` to `+upper`. Otherwise they will be `-lower` to `lower`.
    * If you would like to center the color bar around a different value simply set the appropriate upper and lower bounds (rather than using --center). For example, -u 10 -l 0 will center the color bar around 5.

* `-s` or `--surface` sets which surface the values are shown on.
    * Options are `pial`, `inflated` or `both`.
    * The default is `both`.
    * In the example above the plots will be made on the pial surface (`-s pial`).

* `-t` or `--thresh` sets the threshold value, and all values strictly below this limit will be masked out
    * Note that values that are equal to the threshold value will be included.
    * The default threshold is -98. This is designed such that `-99` values in the `roi_file` will be marked as missing and no color plotted for these regions.
    * You can change the threshold value to anything less than the maximum value.

**Output**

Five files will be created (ten if you've selected `both` for the surface option).

They will all be contained in a folder called `PNGS` that will be created in the same directory as the `roi_file`.

Four of the files will be named as `<roi_file_name>_<hemi>_<surface>_<view>.png`.

For example, if you give as an input an file called `degree_values.txt` you'll get as output:

* `degree_values_lh_pial_lateral.png`
* `degree_values_lh_pial_medial.png`
* `degree_values_rh_pial_lateral.png`
* `degree_values_rh_pial_medial.png`

There will be an additional file named `<roi_file_name>_<surface>_combined.png`.

For example:

* `degree_values_pial_combined.png`

### Viewing the pictures you've made

The command `gthumb` is a useful image viewer on unix systems (including the BCNI and HPC).

eg: `gthumb degree_values_pial_combined.png`
