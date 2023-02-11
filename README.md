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

