# McNaughton Lab (@Lethbridge) Two-Photon Analysis Pipeline

Collection of tools for analysing two-photon calcium imaging data in the McNaughton Lab at the Canadian Centre for Behavioural Neuroscience at the University of Lethbridge.

## Dependencies
* MATLAB R2017b or later
* Image Processing Toolbox
* Signal Processing Toolbox
* Statistics and Machine Learning Toolbox
* Parallel Computing Toolbox (optional)
* GCC or another compatible C/C++ compiler (not required for basic functionalities). All MEX functions were tested on a Linux machine; not sure whether they are cross-platform compatible.

## System Requirements
* A 64-bit Intel Processor with >=SSE4.2 instruction set is recommended to fully utilize the MEX functions
* For multi-plane analysis, a minimum of 16 GB of memory is recommended
* For general use, a minimum of 8 GB of memory is recommended
* For best performance, a discrete graphics card with hardware OpenGL acceleration is recommended

## Installation
1. Clone repository and add the local repository's directory path to the MATLAB path
1. (optional) Compile the MEX files under the `C` folder

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
If you don't like to work with scripts, but instead would like to have everything neatly organised inside one object, then the `LFP` class is for you! This class contains all the methods required to get you from loading and pre-processing your data to plotting figures for dissemination :sunglasses:.

### Data organization
It is recommended for you to organise your data in the following manner:
```
.
└── Animal
    └── Date (yyyy_mm_dd)
        ├── yyyy_mm_dd_1.abf (abf file with session number appended after the date)
        └── Session number (integer)
            ├── Experiment.xml
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
Various properties affect different stages of the analysis. In most cases, the default properties should be adequate as those are the same parameters used by Ingrid and Yours Truly. Currently, it is highly advised to define all of the object properties during instantiation; changing the properties following creation of the object may result in unpredictable behaviours. `obj = LFP(file, 'name', value)`

### Clampex Channels definition
The Clampex channels are defined inside a vector in the following order:
1. two-photon frame pulses
1. encoder channel A
1. encoder channel B
1. reward pulses
1. behaviour camera pulses
1. LFP
1. lick sensor.
Missing channels may be substituted by `NaN`. Make sure the channels are properly defined during initialization: `obj = LFP(file, 'Channels', [1 2 3 5 NaN NaN 7])`. In the latter example, the two-photon frame pulses were acquired on channel 1, the treadmill encoder signals were captures on channels 2 and 3, the reward trigger was gathered on channel 5, and the lick sensor data was obtained on channel 7. The infrared camera and LFP signals were missing. Make sure that the channels are properly defined with `obj.plot('channels')`. The method `LFP.get_channel(chan)` allows you to either retreive the channel name given the channel number (e.g. `obj.get_channel(3)`), or to obtain the channel number given the name (e.g. `obj.get_channel('lfp')`).

### Animal movement
For rest recordings, you may wish to remove all epochs during which the animal was moving, while for running data you may want to discard the frames during which the animal was idle. This can be achieved with the `LFP.remove_mvt(mode)` method, where `mode` is either `'exclude'` (rest only) or `'include'` (run only). This function __must__ be run as a preliminary step before executing the core analysis of the `ensemble` class (see below). The code __will not__ warn you if you omitted this step so please keep it in mind.

### Place cells analysis
`pc_batch_analysis()` can be invoked directly through the class method `LFP.perform_analysis()`. The parameters passed to `pc_batch_analysis()` are stored as a cell array (e.g. `ops = {'test', 'si', 'bins', 80}`) and can be set under the `'pcOps'` property (e.g. `obj.set('pcOps', ops)`).

### Topographical analysis
The FOV can be set manually through the `'FOV'` parameter. If `Experiment.xml` is available, the FOV will be detected automatically. By invoking the `LFP.topography()` method, the peak response locations of place fields are mapped to the spatial locations of neuronal ROIs. These neuronal masks are stored within the `LFP.topo.loc` field.

### Multi-plane recordings
Define the planes to be processed with the `'planes'` parameter. For example, given a dataset with 15 planes, where planes 0, 1, 12, 13 and 14 are flyback frames that should be discarded, the object should be instantiated as `obj = LFP(file, 'planes', 3:12)`. Notice that the plane indices start with 1, despite the fact that the plane folders' indices begin with 0.

The step size between light sections is defined by the `'stepsize'` parameter. The method `LFP.rm_redund()` can be used to automatically remove overlapping neurons across planes (the neuron with higher SNR is preserved). The parameter `'maxStep'` defines the maximum distance for two-neurons to be considered overlapping. The parameter `'ol'` sets the threshold for fraction of overlapping pixels. Finally, the parameter `'maxR'` sets the threshold for Pearson correlation coefficient over which the time-series of two neurons are considered to be the same. All criteria must be satisfied before an ROI is removed.

Alternatively, you can manually remove neurons with the `LFP.rm_neurons()` method. Take a look at the property `LFP.twop.planes.plane_members` to figure out which neuron corresponds with which plane.

Below is a full example:
```
session = '/home/haoran/data/Emily/2020_01_01/2020_01_01_1.abf'; % session to be processed
obj = LFP(session, 'planes', 3:12, 'stepsize', 35, 'maxstep', 100, 'ol', .8, 'maxr', .7);
% load plane2 through plane11
% define step size as 35 um between consecutive light sections (automatically set if Experiment.xml is available)
% set threshold for overlapping neurons to 100 um over the Z aspect and 80% overlapping pixels
% time-series between prospective overlapping neurons must have Pearson's Rho coefficient > .7
obj.rm_redund; % get rid of overlapping neurons
obj.perform_analysis; % run place cells analysis
```

## The `ensemble` class
The `ensemble` class is a subclass of the `LFP` class. Therefore, all methods and properties described in the previous section also apply here. The `ensemble` class contains useful methods for detecting groups of neurons that express synchronous patterns of activity and is therefore suitable for analysing reactivation-type data. Make sure to always remove movement epochs before proceeding with any analysis by running `obj.remove_mvt()`.

### Agglomerative clustering
Agglomerative clustering with average-linkage criterion is used to group neurons into synchronous ensembles. The `'thres'` parameter defines the threshold for average distance between cluster members ( distance metric defined as `1 - abs(Pearson correlation coefficient)`). The `'e_size'` parameter sets the minimum number of neurons required to form an ensemble. After defining these parameters, run agglomerative clustering with the `ensemble.hclust()` method.

### SCE detection
Synchronous calcium events (SCEs) are detected by invoking the `ensemble.detect_sce()` method. The range of durations of a SCE is set via the `'sce_dur'` parameter, where the value is given as a two-element array (e.g. `[0 .8] % from 0 to 800 msec`). Events falling short or ever this range are rejected.

### Topographical analysis
The `ensemble` class comes with its own `ensemble.topography()` method aimed at mapping synchronous ensemble clusters onto physical locations of neuronal ROIs.

### Visualising the results
The `ensemble.plot(type, param)` method is used visualize the results, where `type` dictates the plot to be shown:
* `'tree'` plots the results from agglomerative clustering as a dendrogram
* `'sce'` plots the detected SCE events
* `'clust_topo'` plots the topographical organisation of the detected ensemble, where `param` is an index array of ensembles to be plotted (e.g. `[2 5]` will only plot ensembles 2 and 5)
* `'clust_topo_stack'` plots the topographical organisation of the detected ensemble across multiple planes, where `param` is an index array of ensembles to be plotted
* `'clust_topo_stats'` plots summary statistics for "clusteredness"

## Where to go next?
The `pc_analysis` toolbox has much more to offer. Ask HaoRan about specifics.
