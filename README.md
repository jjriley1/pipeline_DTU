# Pipeline DTU

## Usage

To use pipeline_DTU to investigate DTU of utrons (DTUtrons) use one of the following commands:

    "python <pipeline_location>/pipeline_DTU.py make DTUtrons -v5" <-- local interactive
    "submit_pipeline <pipeline_location>/pipeline_DTU.py make DTUtrons -v5" <-- new session

To use pipeline_DTU to investigate DTU for whole transcripts use one of the following commands:

    "python <pipeline_location>/pipeline_DTU.py make DTU -v5" <-- local session
    "submit_pipeline <pipeline_location>/pipeline_DTU.py make DAS -v5" <-- new session

## Configuration

The pipeline requires a configured 'pipeline.yml' file. An example can be found at:

    "<pipeline_dir>/pipeline_DTU/pipeline.yml"

## Requirements

The pipeline requires the results from:

* pipeline_fastqc
* pipeline_mapping


## How to use:

* Input files in '.bam' format should be placed into an "input_assemble.dir"
* Complete/modify the pipeline.yml
  * Importantly, information in the "dtu" section will need to be modified on a use-by-use basis
  * Including details on:
    * Whether data points are paired / unpaired
    * Whether analysis is on one (univariate analysis) or two (covariate analysis)
    * The design formula that the linear model will be based around
    * Significance thresholds (adjusted p-values and log2fold change thresholds)
* Design matrix must be present in the working directory
  * This allows mutliple comparisons to be made in one go
  * An example is show in "<pipeline_dir>/pipeline_DTU/design.tsv"

