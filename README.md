# McNaughton Lab (@Lethbridge) Two-Photon Analysis Pipeline

Collection of tools for analysing two-photon calcium imaging data in the McNaughton Lab at the Canadian Centre for Behavioural Neuroscience at the University of Lethbridge.

## Dependencies
* MATLAB R2017b or later
* Image Processing Toolbox
* Signal Processing Toolbox
* Statistics and Machine Learning Toolbox
* Parallel Computing Toolbox (optional)
* GCC or another compatible C/C++ compiler (not required for basic functionalities). All MEX functions were tested on a Linux machine; not sure whether they are cross-platform compatible.

## Installation
1. Clone repository and add the local repository's directory path to the MATLAB path.
1. (optional) Compile the MEX files under the `better_than_default` folder.

## Place cells analysis

`pc_batch_analysis()` is the primary function used for classifying and for characterizing place cells.
```
Usage:
    analysis = pc_batch_analysis(behavior, deconv, 'name', value);
Positional Arguments:
    behavior
        behavior structure constructed in "HaoRan's format". "Dun's format" must
        be converted using the function convert_behavior()
    deconv
        matrix of deconvolved time-series, where rows contain frames and columns
        represent neurons
Optional Arguments:
    'mask', maskNeurons, mimg
        required for merging two planes

    'test',
        'si'     SI shuffle test
        'ricker' new method that convolves tuning curve with a series of
                 ricker wavelets
        'mixed' (default) two tests combined

    'shuffles', 1000 (default)
        number of shuffles for tests

    'sig', 0.05 (default)
        significance threshold for p-value

    'bins', 50 (default)
        number of spatial bins

    'sd', 4 (default)
        smoothing kernel s.d. in cm

    'mad', 3 (default)
        MAD threshold for Ricker test

    'frac_trials', 1/3 (default)
        fraction of active trials in place fields (for Ricker test)

    'consecutive', true (default)
        complementary parameter to frac_trials; whether trials need to be
        consecutive

    'width', [.05 .8] (default)
        minimum and maximum width of place fields expressed in fraction

    'io_ratio', 2.5 (default)
        minimum ratio between in-field vs out-field fr

    'par', true (default)
        use parallel processing to speed up (can benefit very long
        recordings)

    'pic', false (default)
        test for alternative hypothesis of "path integrator cells"
```
### Visualising the results
The function `plot_analysis(analysis, [true true true], list)` allows for visualising your results. The first parameter is passed as the output of `pc_batch_analysis()`. The second parameter dictates the type of plots desired. The first flag gives the raster plots of individual cells. The second flag plots the mean population vector as a function of position. The third flag plots the mean running velocity. The third optional parameter defines a list of neurons to be plotted. By default, only neurons that were classified as place cells will be plotted.

The `plot_behaviour(analysis, deconv, dot)` is used to visualise the neural time-series as a function of position and movement velocity. The first parameter is the output of `pc_batch_analysis()`. The second parameter is the original `deconv` matrix. The third parameter indicates whether the time-series should be plotted as the normalized dF/F (`false`) or the "peak-per-trial" (`true`).

## The `LFP` class
If you don't like to work with scripts, but instead would like to have everything neatly organised inside one object, then the `LFP` class is for you! This class contains all the methods required to get you from loading and pre-processing your data to plotting figures for dissemination :sunglasses:

### Data organization
It is recommended for you to organise your data in the following manner:
```
.
└── Animal
    └── Date (yyyy_mm_dd)
        ├── yyyy_mm_dd_1.abf (abf file with session number appended after the date)
        └── Session number (integer)
            ├── behavior.mat
            └── Plane (plane0, plane1, plane2, etc.)
                ├── deconv.mat
                ├── masks_neurons.mat
                ├── mean_img.mat
                ├── ratio_model.mat
                └── timecourses.mat
```
If you only have a single plane, then simply store the data files under the session folder. If you're part of the Polaris Research group, than your data is likely already organised this way.

### Object instantiation
There are two ways to create an object of the LFP class: with an `abf` file or with a `behavior.mat` file. The first option is highly recommended as it performs a more robust extraction of the animal behaviour than the earlier versions. Give the full path of the behavior file during creation of an object: e.g. `obj = LFP(fullfile(pwd, '2020_01_01.abf'))` or `obj = LFP(fullfile(pwd, '2020_01_01', '1', 'behavior.mat'))`. If no file is specified (i.e. `obj = LFP()`), then the user will be prompted to select a file. Properties may be set during this stage: `obj = LFP(file, 'name', value)`.

### Object properties
Currently, it is highly advised to define all of the object properties during instantiation; changing the properties following creation may result in unpredictable behaviours.

## The `ensemble` class
TODO
