# Constructing a robust metric of peak quality for untargeted mass-spectrometry
## Repo information and people
This repository contains data and code processed during the 2023 eScience Winter Incubator. 
  - Project lead: [William Kumler](https://github.com/wkumler)
  - eScience advisor: [Bryna Hazelton](https://github.com/bhazelton)
  - Academic advisor: [Anitra Ingalls](https://sites.google.com/view/anitra-ingalls/home?authuser=0)

## Project summary

Mass spectrometry (MS) is a cutting-edge analysis field used to identify the molecular composition of samples taken from medical laboratories, the depths of the ocean, and even outer space. In the Ingalls Lab at UW, we use it to characterize the molecular composition of seawater and its inhabitants, a task complicated by the complex biogeochemistry of the oceans. The nascent nature of modern mass spectrometry also introduces many challenges, one of which is distinguishing biological/chemical signal from noise produced during the measurement process.

My goal in this incubator is to calibrate existing detection algorithms to a probabilistic likelihood that the signal corresponds to a real molecular feature. This will involve estimating the relative strength of various metrics used for detecting molecules, using machine learning methods to construct the probabilistic estimate, and ideally constructing packages that interface with existing software to facilitate widespread adoption.

## Project goals

First: Estimating the relative strength of various peakpicking metrics (e.g. peak shape, signal-to-noise, height, presence of isotopes)
  - Answers the question "Which peak metrics do the best job of separating real biological/chemical signal from noise?"
  - Will require a large dataset, maybe multiple, of peaks flagged as good/bad by an MS expert
  
Second: Using the above metrics to calculate the likelihood that an MS expert would flag a given peak as noise
  - Answers the question "What is the probability this peak is just noise and should be removed from the analysis?"
  - Will require the use of machine learning methods like logistic ridge regression, random forests, and others to handle highly-correlated metrics
  - This has become more of a specific request after some discussion about regression models, answering the question of "what thresholds (and metrics?) do I need to choose to ensure that no more than 5% of the peaks identified as 'good' are actually noise?" aka a strict control of the false discovery rate (ratio of false positive to true+false positives < 0.05).
  
Third: Building a package or code-sharing method that accepts input from commonly used peakpicking algorithms (e.g. [xcms](https://github.com/sneumann/xcms), [MzMine](http://mzmine.github.io/), and [MSDIAL](http://prime.psc.riken.jp/compms/msdial/main.html))

## Metrics of interest

Although part of the project is developing new metrics that will help separate good and bad peaks, I'm coming into it with some intuition about which parameters will be the most helpful. I've divided these into a few different categories I'm calling per-peak, per-file, and per-feature. Per-peak parameters are things estimated from the raw data, the mz/rt/int data points. Per-file parameters are things estimated within a given file like isotope presence/absence. Per-feature parameters are things estimated across multiple files after peak alignment and correspondence ("feature" in MS often refers to a chemical signal showing up in multiple files).

#### Per-peak metrics
Anything from XCMS directly, averaged across the different files

  - mean_mz: Average mass to charge ratio
  - sd_ppm: Average m/z deviation across multiple files in PPM space
    - Using PPM accuracy instead of absolute deviation accounts for the fact that heavier molecules tend to have larger mass errors
  - mean_rt: Average retention time
  - sd_rt: Standard deviation of retention times
  - mean_pw: Average peakwidth (in seconds)
  - sd_pw: Standard deviation of the peakwidths
  - sn: XCMS's default signal-to-noise estimate
    - Has known problems, especially for "spiky" HILIC data
  - f: Unclear parameter relating to the ROI number
  - scale: Unclear parameter relating to the Centwave peak detector
  - log_mean_height: The average of the log-scaled peak heights
  - log_sd_height: The standard deviation of the log-scaled peak heights

Custom parameters:

  - med_cor: Custom calculation estimating the correlation between the points in a perfectly smooth peak and the actual data
    - Designed to measure the "Gaussian-ness" or "peak shape"
  - med_SNR: Custom calculation estimating the Signal-to-Noise Ratio of a peak
    - Uses the residuals from the smooth peak fit in the med_cor step as an estimate of the noise *within* a peak
    - This is different from other peakpicking algorithms which require points outside of the peak to estimate the noise background - which doesn't always exist in spiky HILIC data
  - int_vs_ppm: high-intensity points tend to be more accurate, so we should expect a correlation between the intensity of a point and its deviation from the "center" of the peak
    - Currently unimplemented
  - med_missed_scans: The total number (maybe better as percentage?) of "missed" scans within a peak. Good looking peaks have data points at every retention time, while poor-quality ones often have missing data

#### Per-file metrics
Really just limited to ^13^C isotope information at this point but could eventually include anything estimated from other peaks in the same file

  - shape_cor: The correlation between the intensity values of the base peak and the intensity values of the ^13^C peak
  - area_cor: The correlation between the peak area of the base peak and the peak area of the ^13^C peak across multiple files

#### Per-feature metrics
Here's where we include some more "heuristic" metrics - usually I see these used as thresholds to remove bad peaks later on with rules like "average area has to be 3x the blank" and of course it's unlikely that noise has a significant trend with sample type.

  - smp_to_blk: Ratio of average sample area to average blank area
  - t_pval: The p-value for a statistical test measuring differences between sample types (e.g. surface samples vs deep samples)
    - Was t-tests with Falkor data, now expanded to ANOVAs for MESOSCOPE
  - smp_to_std: Ratio of average sample area to average standard area
    - Tends to be good at identifying peaks that are only found in the standard mixes which are otherwise hard to classify - a good looking peak in the standards is still worth noticing but is not a great peak in the samples
  - feat_npeaks: Total number of peaks included in the feature
  - n_found: Number of files (now % of files) in which the peak was found
  - samps_found: Number of samples (now % of samples) in which the peak was found
  - stans_found: Number of standards (% of stans) in which the peak was found
  - blank_found: Boolean. Whether or not a peak was found in the blank

#### Metric expectations
I expect a "good" peak to have

  - Large log_mean_height
  - Low sd_ppm
  - Low sd_rt
  - Middling mean_pw
  - Low sd_pw
  - High med_cor
  - High med_SNR
  - High int_vs_ppm
  - Low med_missed_scans
  - High shape_cor
  - High area_cor
  - High smp_to_blk
  - Low t_pval
  - Middling feat_npeaks
  - High n_found
  - High samps_found

Expected unimportant metrics:

  - log_sd_height
  - mean_mz
  - mean_rt
  - sn
  - f
  - scale
  - smp_to_std
  - stans_found
  - blank_found

## Repo structure

R scripts:

  - peakpicking_and_prep.R
    - Takes in the mzML files in mzMLs/ and performs peakpicking and other xcms tasks on them, then writes out peak boundaries into made_data/feature_bounds.csv
  - training_tool.R
    - Takes in the peak boundaries from made_data/feature_bounds.csv and the raw mzML files again, then enables an interactive script using plot windows that lets the user categorize features as "good", "bad", and other categories.
  - feature_extraction.R
    - Extracts information about the peaks from the raw data and output for use in later ML/classification/regression things.
  - *_attempts.R
    - Scripts that run various ML attempts - clustering, random forests, logistic regressions, etc.
    - Takes in the features_extracted.csv and outputs figures, mostly.
  - msdial_prep_script.R
    - (In progress) Runs MSDIAL on the same set of files for comparison with XCMS
    - The goal is to fiddle with peakpicking parameters to get ~the same set of peaks and see if peaks that are "good" are most likely to be found in both XCMS and MSDIAL
    - Constructs the entire msdial folder (not shared with GitHub because of large file sizes)

made_data_FT350:

  - Contains output produced by scripts during the process
  - Allows each script to be rerun independently of the one before it
  - Used to load files into memory during the "setup" chunk of each script

made_data_FT2040:

  - Same as 350 but run with slightly different parameters to include additional peaks
  - See commit 17df98e for details

figures:

  - Contains figures produced by the scripts during the process
  - Also contains the presentation 2923 given to Dave on Feb 9th, 2023

