##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################

"""
===========================
pipeline_DTU
===========================



:Author: Jack Riley
:Release: $Id$
:Date: |today|
:Tags: Python



Overview
========

This pipeline is designed to detect differential transcript usage (DTU) between experimental conditions
using RNA-seq data as input. The pipeline can be used for 2 purposes:
    - Primarily, (due to nature of JR's PhD Project) pipeline_DTU is designed to predict alternative splicing events in the 3'UTR of 
    transcripts between conditions. --> maybe it should be called "pipeline_DTUtrons"
        - in order to do this, code from pipeline_utrons (https://github.com/sudlab/pipeline_utrons) is used, and those transcripts
        classified as having "partnered_utrons" are filtered for. 
        - candidate genes output from this pipeline are those which show a flip (from major to minor isoform [or viceversa]) in expression
        of 3'UTR intron (utron-) splicing vs retaining transcripts between conditions 
    - Secondly, pipeline_DTU can be used to investigate differential transcript usage between conditions for whole transcripts as opposed to just 
    "utron- retaining/spliced" transcripts.
        - this will be facilitated by the classical usage of DEXSeq

There are multiple statistical packages that exist to determine DTU, and this pipeline uses 3 to maximise output:

    - DEXSeq: "DEXSeq assumes a Negative Binomial (NB) distribution for the feature counts, and considers the counts
     for each feature (originally, the exonic parts) relative to counts for all other features of the group (the gene), 
     using an interaction term in a generalized linear model (GLM). The GLM framework is an extension of the linear model
     (LM), but shares with LM the usage of a design matrix, typically represented by X, which is made up of columns of 
     covariates that are multiplied by scalar coefficients, typically represented by β. The design matrix with its multiple
     coefficients is useful for extending statistical models beyond simple group comparisons, allowing for more complex 
     situations, such as within-patient comparisons, batch correction, or testing of ratios." - Love, Soneson and Patro (2018)

    - DRIMSeq: "In contrast to the NB model, DRIMSeq assumes an Dirichlet Multinomial model (DM) for each gene, where the total
     count for the gene is considered fixed, and the quantity of interest is the proportion for the transcript within a gene for
     each sample. If the proportion for one transcript increases, it must result in a decrease for the proportions of the other
     transcripts of the gene. Genes that are detected as having statistically significant DTU are those in which the proportion
     changes observed across condition were large, considering the variation in proportions seen within condition. The variation
     in proportions across biological replicates is modeled using a single precision parameter per gene ... DRIMSeq also uses a 
     design matrix, and so can be used to analyze DTU for complex experimental designs, including within-patient comparisons, 
     batch correction, or testing of ratios." - Love, Soneson and Patro (2018)

    - swish (from the fishpond package): "We note that swish extends and builds on another method, SAMseq (Li and Tibshirani 2011), 
     implemented in the samr package, by taking into account inferential uncertainty, and allowing to control for batch effects and 
     matched samples. Additionally, swish has methods for testing changes in effect size across secondary covariates, which we refer
     to as “interactions”. swish calls functions from the qvalue (Storey and Tibshirani 2003) or samr package for calculation of 
     local FDR and q-value" - Zhu, Srivastava, Ibrahim, Patro, Love (2020)

Outputs of DTU are subset into those which are either unique to each package, or those which are present in outputs of multiple methods

N.b. for use with more complex design matrices/formulae (e.g. more than 1 condition, i.e. with 1 or more 'interaction terms'), only linear 
modelling can be used to resolve this problem. In this regard, only DEXseq and DRIMseq are suitable, and the pipeline will run accordingly
based on this. For example, "condition = day" will determine DTU between days (as factors) and is suitable input for all 3 packages, 
whereas "design = ~variant + time + variant:time" is unsuitable for swish and will only be used by DRIMSeq and DEXSeq.
                 
                        ###### SWISH MAY BE ABLE TO DO INTERACTIONS/COVARIATE ANALYSES ######



Usage
=====

To use pipeline_DTU to investigate DTU of utrons (DTUtrons) use one of the following commands:

    "python <pipeline_location>/pipeline_DTU.py make DTUtrons -v5" <-- in local interactive session
    "submit_pipeline <pipeline_location>/pipeline_DTU.py make DTUtrons -v5" <-- in new session (background / non-interactive)

To use pipeline_DTU to investigate DTU for whole transcripts use one of the following commands:

    "python <pipeline_location>/pipeline_DTU.py make DTU -v5" <-- in local interactive session
    "submit_pipeline <pipeline_location>/pipeline_DTU.py make DAS -v5" <-- in new session (background / non-interactive)



Configuration
=============

The pipeline requires a configured 'pipeline.yml' file. An example can be found at:

    "<source_dir>/pipeline_DTU/pipeline_DTU/pipeline.yml"

In order for the R scripts to work, correct naming convention must be used, and number
and names of covariates must be provided in the config file, alongside whether the data
is paired (on identifier [first input before '-']). 


Input files
===========

Files must be named using the following convention: (separated with '-')

<source/identifier/cell_line> - <variable 1> - <(optional) variable 2> - <(replicate number>.bam

(if no technical replicates then name all samples R0)

e.g. for a pooled timecourse differentiation data comparing day 0 and day 10, one may use files such as:
NPC-D00-R1.bam
NPC-D00-R2.bam
NPC-D10-R1.bam
NPC-D10-R2.bam
(i.e. technical / pooled replicates)

in the config file you can also specify whether the data is to be paired on the identifier/cell_line, example names:
NPC1-D00-R1.bam
NPC1-D10-R1.bam
NPC2-D00-R1.bam
NPC2-D10-R1.bam
(i.e. comparing different cell lines or states)

whether 1 or 2 covariates are to be tested is also specified in the config file. examples:
NPC-D00-treated-R1.bam
NPC-D10-treated-R1.bam
NPC-D00-control-R1.bam
NPC-D10-control-R1.bam

A design.tsv file must also be included within the working directory, and examples of single vs covariate design.tsv files
can be found in:

    "<source_dir>/pipeline_DTU/pipeline_DTU/example_design.tsv"

Requirements
============

The pipeline requires the results from:

    - pipeline_fastqc
    - pipeline_mapping

Input files in '.bam' format should be placed into an "input_assemble.dir".

### ULTIMATELY I WISH TO INCORPORATE fastqc + mapping INTO THE PIPELINE AT SOME STAGE ###



Pipeline output
===============

- pipeline_DTU will output quantification data from Salmon
- a list of candidates which display DTU in one/multiple/all software packages
    - alongside expression data for both retaining + spliced partners


Code
====

"""

##### imports #####

from ruffus import *
from ruffus.combinatorics import product
import sys
import os
import shutil
from gffutils import DataIterator as DataIterator
import sqlite3
import subprocess
import glob
from cgatcore import experiment as E
import cgat.Sra as Sra
from cgatcore import pipeline as P
import cgatpipelines.tasks.rnaseq as RnaSeq
import tempfile

##### main ##### 



##### run fxns ######

def DTUtrons():
    #decorate this with pipeline_utrons specific things
    pass

def DTU():
    #decorate this with DTU specific things
    pass

###### misc ######

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))