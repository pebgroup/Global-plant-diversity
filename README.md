# Global plant diversity project

Using the WCVP taxonomy and distribution data, and Smith&Brown 2018 phylogeny to analyse environmental and evolutionary drivers of global plant diversity patterns.

## Scripts

All scripts are stored in *publish_scripts* folder. To go through the analysis step by step, open *master_analysis.R*, which sources all other scripts required to repeat analysis. The header in the master script contains system requirements and execution time information.

### Reproducibility and completeness
Please note that for the sake of reproducibility, we include *all* scripts. That includes also scripts processing large primary data files which are not included in the supplement for size reasons (i.e. CRU TS, wordclim, soil). Results from processing these files are included however. To repeat these steps, the necessary data can be obtained as free downloads as described [here](data_description.md)). 


## Data
Primary data are files obtained online / by collaborators. Processed data are files which are produced during the analysis steps and stored to shorten execution time for later analysis steps and to keep work space small. Files are listed with description here: [data description](data_description.md).

Large data files that can be freely accessed online (i.e. CRU TS, wordclim) are not part of the supplement but listed with version number and information how to obtain them.

**All steps are recorded, the final data set used for maps, figures and structural equation modeling is *sem_input_data.rds***.

## Analysis workflow
<p align="left">
<img src="flow_chart_wcvp.png"/>  
</p>  