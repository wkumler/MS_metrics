# Picky with peakpicking: assessing chromatographic peak quality with simple metrics in metabolomics
Note: this README was significantly changed at the time of manuscript submission to reflect the final state of the project. For the in-progress README, see [here](https://github.com/wkumler/MS_metrics/blob/675ca3f002a96f94a6d0d5a65124ee0b2b40aa98/README.md).

## Repo information
This repository contains data and code necessary for the peakpicking manuscript submitted to BMC Bioinformatics that emerged from the 2023 eScience Winter Incubator. The manuscript has been written as a reproducible R Markdown file available in the /manuscript subdirectory but requires downloading the raw mzML data from Metabolomics Workbench as the files are too large to be hosted on GitHub directly. To recompile the manuscript, clone the repository and copy the positive-mode mzML files over from MW into the associated subdirectory's /mzML folder. The entire document can then be recreated using RStudio's `Knit` function after launching it from the .Rproj file. This will take a while (~8 hours on my laptop), as the code runs XCMS on ~500 mass-spectrometry files, performs custom feature extraction, performs the statistics and generates the figures for the manuscript. This was tested on 7/25/2023 ahead of manuscript submission.

## Authors
  - Project lead: [William Kumler](https://github.com/wkumler) (wkumler@uw.edu)
  - eScience advisor: [Bryna Hazelton](https://github.com/bhazelton)
  - Academic advisor: [Anitra Ingalls](https://sites.google.com/view/anitra-ingalls/home?authuser=0)

## Repo structure
  - manuscript
    - Contains the R Markdown document used to create the Word document submitted to BMC Bio via knit and additional necessary data. Everything other than the RMarkdown document, the /IS_integrations directory, the refs.bib, and template.docx are created via knit.
    - IS_integrations: folder containing Skyline integrations of the internal standards used to normalize the data via BMIS.
      - The is_combine.R script can be run to convert the *.sky files into clean all_IS.csv format via the SkylineRunner.exe CLI
    - refs.bib: References for the manuscript, created by exporting the Mendeley project folder
    - template.docx: Template for knitting, enabling e.g. line numbers and double-spacing and header formatting
  - MW_upload
    - Metabolomics Workbench doesn't really like it when you upload projects a little bit at a time but requires the entire project be near completion prior to providing a DOI. Since the files used here are part of both this manuscript and a followup oceanography manuscript, there's some overlap between the two when uploading to MW.
    - MW_formatting.R: Accepts as input the manually-copied CSV files and produces the Metabolomics Workbench formatted .txt files.
    - *.csv: Files are copied over from a different repository (AllMeso) where the different analyses were collated into single documents and targeted data was used as well as the untargeted approach illustrated here.
    - *.txt: Outputs from MW_formatting.R that were opened and copy-pasted into the MW upload forms.
  - made_data_*
    - Contains raw data (mzMLs) and outputs from XCMS used in the manuscript.
      - \mzMLs: Needs to be manually created within each folder by extracting the .mzML files from the Metabolomics Workbench studies.
        - made_data_CultureData/mzMLs: [Project ID](dx.doi.org/10.21228/M8QM5H) -> Study ID ST002790 -> Download raw data -> [ST002790_culturedata_HILIC_POS.zip](https://www.metabolomicsworkbench.org/studydownload/ST002790_culturedata_HILIC_POS.zip)
        - made_data_MS3000/mzMLs: [Project ID PR001738](dx.doi.org/10.21228/M82719) -> Study ID ST002789 -> Download raw data -> [ST002789_HILIC_POS.zip](https://www.metabolomicsworkbench.org/studydownload/ST002789_HILIC_POS.zip)
        - made_data_MS3000/mzMLs: [Project ID PR001738](dx.doi.org/10.21228/M82719) -> Study ID ST002788 -> Download raw data -> [ST002788_HILIC_POS.zip](https://www.metabolomicsworkbench.org/studydownload/ST002788_HILIC_POS.zip)
        - made_data_CultureData/mzMLs: [Project ID](http://dx.doi.org/10.21228/M8GH6P) -> Study ID ST002077 -> Download raw data -> [ST002790_Pttime.zip](https://www.metabolomicsworkbench.org/studydownload/ST002077_Pttime.zip)
    - classified_feats.csv: Contains manual annotations of peak quality for each feature in the dataset, named by XCMS feature number alongside the *m/z*/RT bounding box coordinates. Feature classifications can be "Good", "Bad", "Ambiguous", "Stans only", or "Unclassified". Annotations were performed manually using the `training_tool.R` script found in the root project directory after the `peakpicking_and_prep.R` script was run.
    - features_extracted.csv: Contains the extracted metrics for each mass feature. Rows are individual XCMS features and columns are the metric extracted from each mass feature. NAs have been filled and the feature classification (see above) has been added as the final column.
  - peakpicking_and_prep.R: Runs XCMS on the mzMLs in each subdirectory. The first step in the process.
  - feature_extraction.R: Performs feature extraction on the XCMS output from peakpicking_and_prep.R. Still requires access to the raw data, however.
  - training_tool.R: Interactive annotation tool built to rapidly render a mass feature to the user and capture input using the arrow keys to categorize the feature's quality, one XCMS project at a time. Currently right arrow key is bound to "Bad", left is "Good", up is "Ambiguous", and down is "Stans only".
  - *_attempts.R and *_validation.R: Various sandboxed attempts to wrap our head around the data. None of these are necessary for the final manuscript, as analyses and figures have been embedded as R chunks.
