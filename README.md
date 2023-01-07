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
  
Third: Building a package or code-sharing method that accepts input from commonly used peakpicking algorithms (e.g. [xcms](https://github.com/sneumann/xcms), [MzMine](http://mzmine.github.io/), and [MSDIAL](http://prime.psc.riken.jp/compms/msdial/main.html))
