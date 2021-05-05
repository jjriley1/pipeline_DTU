#Pipeline DTU

##Usage

To use pipeline_DTU to investigate DTU of utrons (DTUtrons) use one of the following commands:

    "python <pipeline_location>/pipeline_DTU.py make DTUtrons -v5" <-- in local interactive session
    "submit_pipeline <pipeline_location>/pipeline_DTU.py make DTUtrons -v5" <-- in new session (background / non-interactive)

To use pipeline_DTU to investigate DTU for whole transcripts use one of the following commands:

    "python <pipeline_location>/pipeline_DTU.py make DTU -v5" <-- in local interactive session
    "submit_pipeline <pipeline_location>/pipeline_DTU.py make DAS -v5" <-- in new session (background / non-interactive)

##Configuration

The pipeline requires a configured 'pipeline.yml' file. An example can be found at:

    "<source_dir>/pipeline_DTU/pipeline_DTU/pipeline.yml"

##Requirements

The pipeline requires the results from:

    * pipeline_fastqc
    * pipeline_mapping

Input files in '.bam' format should be placed into an "input_assemble.dir".
