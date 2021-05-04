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
    - Secondly, pipeline_DTU can be used to investigate alternative splicing between conditions for ALL transcripts as opposed to just 
    "utron- retaining/spliced" transcripts. 

There are multiple statistical packages that exist to determine DTU, and this pipeline uses 3 to maximise output:

    - DEXSeq: <description>
    - DRIMSeq: <description>
    - swish (from the fishpond package): <description>

Outputs of DTU are subset into those which are either unique to each package, or those which are present in outputs of multiple methods

N.b. for use with more complex design matrices/formulae (e.g. more than 1 condition, i.e. with 1 or more 'interaction terms'), only linear 
modelling can be used to resolve this problem. In this regard, only DEXseq and DRIMseq are suitable, and the pipeline will run accordingly
based on this. For example, "condition = day" will determine DTU between days (as factors) and is suitable in put for all 3 packages, 
whereas "design = ~variant + time + variant:time" is unsuitable for swish and will only be used by DRIMSeq and DEXSeq.



Usage
=====

To use pipeline_DTU to investigate DTU of utrons (DTUtrons) use one of the following commands:

    "python <pipeline_location>/pipeline_DTU.py make DTUtrons -v5" <-- in local interactive session
    "submit_pipeline <pipeline_location>/pipeline_DTU.py make DTUtrons -v5" <-- in new session (background / non-interactive)

To use pipeline_DTU to investigate DTU for ALL transcripts use one of the following commands:

    "python <pipeline_location>/pipeline_DTU.py make DTU -v5" <-- in local interactive session
    "submit_pipeline <pipeline_location>/pipeline_DTU.py make DTU -v5" <-- in new session (background / non-interactive)



Configuration
=============

The pipeline requires a configured 'pipeline.yml' file. An example can be found at:

    "<source_dir>/pipeline_DTU/pipeline_DTU/pipeline.yml"



Input files
===========

### <TO COMPLETE>: need info on design matrix, file naming conventions etc ###



Requirements
============

The pipeline requires the results from:

    - pipeline_fastqc
    - pipeline_mapping

### ULTIMATELY I WISH TO INCORPORATE THESE INTO THE PIPELINE AT SOME STAGE ###



Pipeline output
===============

.. Describe output files of the pipeline here



Glossary
========

.. glossary::



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