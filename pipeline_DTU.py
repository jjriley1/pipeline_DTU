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
(i.e. looking for DTU that is common between different cell lines or states)

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
###################
##### imports #####
###################

from ruffus import *
from ruffus.combinatorics import product
import sys
import os
import shutil
from gffutils import DataIterator as DataIterator
import sqlite3
import subprocess
import glob
import csv
from cgatcore import experiment as E
import cgat.Sra as Sra
from cgatcore import pipeline as P
import cgatpipelines.tasks.rnaseq as RnaSeq
import tempfile

#############################################
##### main - substantial amount of this ##### 
#####  code comes from pipeline_utrons  #####
#############################################

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

PARAMS.update(P.peek_parameters(
    PARAMS["annotations_dir"],
    'genesets',
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))

PARAMS["project_src"]=os.path.dirname(__file__)

RnaSeq.PARAMS = PARAMS

###################
##### utility #####
###################

def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.
    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

STRINGTIE_QUANT_FILES=["i_data.ctab", "e_data.ctab", "t_data.ctab",
                       "i2t.ctab", "e2t.ctab"]

######################################
##### assemble novel transcripts #####
######################################

@follows(mkdir("assembled_transcripts.dir"), mkdir("portcullis"))
@transform(["input_assemble.dir/*.bam",
            "input_assemble.dir/*.remote"],
           formatter(),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"])),
           "assembled_transcripts.dir/{basename[0]}.gtf.gz")
def assembleWithStringTie(infiles, outfile):

    infile, reference = infiles
    basefile = os.path.basename(infile)
    job_threads = PARAMS["stringtie_threads"]
    job_memory = PARAMS["stringtie_memory"]
    tmpfile = P.get_temp_filename()
    if os.path.exists(tmpfile): 
    	os.unlink(tmpfile)
    
    statement =  '''portcullis full 
			            -t 1
                        -o portcullis/%(basefile)s/
                        -r %(portcullis_bedref)s
                        -b
                        %(portcullis_fastaref)s
                        %(infile)s &&
                    mv portcullis/%(basefile)s/portcullis.filtered.bam %(tmpfile)s &&
                    rm -r portcullis/%(basefile)s/ &&
                    stringtie %(tmpfile)s
                        -p %(stringtie_threads)s
                        -G <(zcat %(reference)s)
                        %(stringtie_options)s
                        2> %(outfile)s.log
                    | gzip > %(outfile)s &&
                    rm %(tmpfile)s'''

    P.run(statement)

###########################
##### merge/aggregate #####
###########################

@follows(mkdir("final_genesets.dir"), assembleWithStringTie)
@merge([assembleWithStringTie,
       os.path.join(
           PARAMS["annotations_dir"],
           PARAMS["annotations_interface_geneset_all_gtf"])],
       "final_genesets.dir/agg-agg-agg.gtf.gz")
def mergeAllAssemblies(infiles, outfile):

    infiles = ["<(zcat %s)" % infile for infile in infiles]
    infiles, reference = infiles[:-1], infiles[-1]

    job_threads = PARAMS["stringtie_merge_threads"]
    job_memory = PARAMS["stringtie_merge_memory"]

    infiles = " ".join(infiles)

    statement = '''stringtie --merge
                            -G %(reference)s
                            -p %(stringtie_merge_threads)s
                            %(stringtie_merge_options)s
                            %(infiles)s
                            2> %(outfile)s.log
                    | cgat gtf2gtf --method=sort
                            --sort-order=gene+transcript
                            -S %(outfile)s -L %(outfile)s.log'''

    P.run(statement) 

###############################################
##### configurable: merge by tissue/group #####
###############################################

@follows(assembleWithStringTie)
@collate([glob.glob("assembled_transcripts.dir/*%s*.gtf.gz" % PARAMS["stringtie_groups"][x]) for x in range(0, len(PARAMS["stringtie_groups"]))],
         regex("(.+)/(.+)-(.+).gtf.gz"),
         add_inputs(os.path.join(
            PARAMS["annotations_dir"],
            PARAMS["annotations_interface_geneset_all_gtf"])),
          r"final_genesets.dir/\2-agg-agg-agg.gtf.gz")
def merge_by_tissue(infiles, outfile):
    job_threads = PARAMS["stringtie_merge_threads"]
    job_memory= PARAMS["stringtie_merge_memory"]

    reference = "<(zcat %s)" % infiles[0][1]
    infiles = ["<(zcat %s)" % infile for infile in infiles[0][0]]
    infiles = " ".join(infiles)
    
    statement = '''stringtie --merge
                    -G %(reference)s
                    -p %(stringtie_merge_threads)s
                    %(stringtie_merge_options)s
                    %(infiles)s
                    2> %(outfile)s.log
                | cgat gtf2gtf --method=sort
                    --sort-order=gene+transcript
                    -S %(outfile)s -L %(outfile)s.log'''
    P.run(statement)

###################
##### utility #####
###################

@follows(mergeAllAssemblies, merge_by_tissue)
def Assembly():
    pass

#############################
##### classify/annotate #####
#############################

@transform([assembleWithStringTie, mergeAllAssemblies],
           suffix(".gtf.gz"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"])),
           ".class.gz")
def classifyTranscripts(infiles, outfile):
    '''classify transcripts.'''
    to_cluster = True
    infile, reference = infiles
    counter = PARAMS['gtf2table_classifier']
    job_memory = "16G"

    statement = ''' zcat %(infile)s
                    | cgat gtf2table
                    --counter=%(counter)s
                    --reporter=transcripts
                    --gff-file=%(reference)s
                    --log=%(outfile)s.log
                    | gzip
                    > %(outfile)s
                    '''
    P.run(statement)

#####################################
##### load .class into database #####
#####################################

@follows(mkdir("database_load"), classifyTranscripts)
@merge(classifyTranscripts,
       "database_load/transcript_class.load")
def loadTranscriptClassification(infiles, outfile):

    P.concatenate_and_load(infiles, outfile,
                         regex_filename=".+/(.+).class.gz",
                         options="-i transcript_id -i gene_id"
                         " -i match_gene_id -i match_transcript_id"
                         " -i source",
                         job_memory="64G")

############################################
##### find_utrons from pipeline_utrons #####
############################################

@follows(mkdir("utron_beds.dir"), classifyTranscripts)
@subdivide([assembleWithStringTie, mergeAllAssemblies],
           regex("(.+)/(.+).gtf.gz"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"]),
                      r"\1/\2.class.gz"),
           [r"utron_beds.dir/\2.all_utrons.bed.gz",
            r"utron_beds.dir/\2.partnered_utrons.bed.gz",
            r"utron_beds.dir/\2.novel_utrons.bed.gz"])
def find_utrons(infiles, outfiles):

    infile, reference, classfile = infiles
    job_threads=2
    job_memory="16G"

    all_out, part_out, novel_out = outfiles

    track = P.snip(all_out, ".all_utrons.bed.gz")
    current_file = __file__ 
    pipeline_path = os.path.abspath(current_file)
    pipeline_directory = os.path.dirname(pipeline_path)
    script_path = "pipeline_DTU/find_utrons.py"
    full_utron_path = os.path.join(pipeline_directory, script_path)  

    statement = ''' cgat gtf2gtf -I %(infile)s
                        --method=sort
                        --sort-order=gene+transcript
                        -L %(track)s.log
                    | python %(full_utron_path)s 
                        --reffile=%(reference)s
                        --class-file=%(classfile)s
                        --outfile %(all_out)s
                        --partfile=%(part_out)s
                        --novel-file=%(novel_out)s
                        -L %(track)s.log
                        '''
    P.run(statement)

#########################################
##### get utron ids from utron_beds #####
#########################################

@transform(find_utrons,
           suffix(".bed.gz"),
           ".ids.gz")
def getUtronIds(infile, outfile):

    statement = ''' zcat %(infile)s 
                    | cut -f 4
                    | sed 's/:/\\t/g'
                    | sort -u
                    | gzip > %(outfile)s
                    '''
    P.run(statement)

########################################
##### load utron ids into database #####
########################################

@collate(getUtronIds,
         regex(".+/(.+)\.(.+).ids.gz"),
         r"database_load/\2_ids.load")
def loadUtronIDs(infiles, outfile):
    job_threads=3

    header = "track,transcript_id"
    options = "-i track -i transcript_id"

    if not outfile == "all_utrons_ids.load":
        header += ",match_transcript_id"
        options += " -i match_transcript_id"

    P.concatenate_and_load(infiles, outfile,
                         regex_filename=".+/(.+)\..+\.ids.gz",
                         has_titles=False,
                         cat="track",
                         header=header,
                         options=options,
                         job_memory="64G")

###################
##### utility #####
###################

@follows(loadUtronIDs, loadTranscriptClassification)
def AnnotateAssemblies():
    pass

#######################################################
##### export assembled and aggregated transcripts #####
#######################################################

@follows(mkdir("export/indexed_gtfs.dir"))
@transform([assembleWithStringTie,
            mergeAllAssemblies,
            merge_by_tissue,
            "final_genesets.dir/*.gtf.gz"],
           regex("(.+)/(.+).gtf.gz"),
           r"export/indexed_gtfs.dir/\2.gtf.gz")
def exportIndexedGTFs(infile, outfile):

    statement = ''' zcat %(infile)s
                    | sort -k1,1 -k4,4n
                    | bgzip > %(outfile)s &&
                    tabix -f -p gff %(outfile)s &&
                    if [[ %(outfile)s == *"agg-agg-agg"* ]]; then 
                        cp %(outfile)s %(outfile)s.tbi export/; 
                    fi;
                    '''
    P.run(statement)

###################
##### utility #####
###################

@follows(exportIndexedGTFs)
def export():
    pass

##################################################
##### create decoy-aware transcriptome/index #####
##################################################

@follows(mkdir("salmon_index"), mkdir("salmon_index/DAT"),exportIndexedGTFs)
@transform("final_genesets.dir/agg-agg-agg.gtf.gz",
           formatter(),
           "salmon_index/{basename[0]}.salmon.index")
def makeSalmonIndex(infile,outfile):
    job_memory="64G"
    job_threads=1

    gtf_basename = P.snip(os.path.basename(infile), ".gtf.gz")
    transcript_fasta = "salmon_index/" + gtf_basename + "transcripts.fa"
    fastaref =PARAMS["portcullis_fastaref"]
    index_options=PARAMS["salmon_indexoptions"]
    tmpfile = P.get_temp_filename()

    #statement since salmon >v1.0 Selective Alignment update
    #statement now generates decoy.txt from reference genome, and concats the genome.fa and transcriptome.fa into gentrome.fa
    #gentrome.fa is then passed through salmon alongside the decoys to create decoy-aware transcriptome (DAT)
    statement = ''' gunzip -c %(infile)s > %(tmpfile)s;
                    gffread %(tmpfile)s -g %(fastaref)s -w %(transcript_fasta)s;
                    grep "^>" <%(fastaref)s | cut -d " " -f 1 > salmon_index/DAT/decoys.txt;
                    sed -i.bak -e 's/>//g' salmon_index/DAT/decoys.txt; 
                    cat %(transcript_fasta)s %(fastaref)s > salmon_index/DAT/gentrome.fa.gz;
                    salmon index
                        -p %(job_threads)s
                        %(index_options)s
                        -t salmon_index/DAT/gentrome.fa.gz
                        -d salmon_index/DAT/decoys.txt
                        -i %(outfile)s;
                    rm %(tmpfile)s
                    '''
    P.run(statement)


##########################
##### quantification #####
##########################

if not os.path.exists("temp_bams/"):
    os.makedirs("temp_bams/")

@follows(mkdir("quantification.dir"), mkdir("sorted_bams"),
         mergeAllAssemblies, makeSalmonIndex)
@product(["input_assemble.dir/*.bam",
          "input_assemble.dir/*.remote"],
         formatter(".+/(?P<TRACK>.+).([bam|remote])"),
         mergeAllAssemblies,
         formatter(".+/(?P<GENESET>.+).gtf.gz"),
         "quantification.dir/{TRACK[0][0]}_{GENESET[1][0]}")
def quantifyWithSalmon(infiles, outfile):
    '''Quantify existing samples against genesets'''
    job_threads=2
    job_memory="24G"

    infile, gtffile = infiles
    basefile = os.path.basename(infile)
    sample_name = basefile.split(os.extsep, 1)
    gtfbase = P.snip(os.path.basename(gtffile), ".gz")
    salmonIndex = "salmon_index/" + gtfbase + ".salmon.index"

    salmon_options = PARAMS["salmon_quantoptions"]
    bootstrap_options = PARAMS["num_bootstraps"]

    sorted_bam="sorted_bams/" + sample_name[0] + "_sorted.bam"
    fastq1 = P.snip(outfile, "_agg-agg-agg")+".1.fastq"
    fastq2 = P.snip(outfile, "_agg-agg-agg")+".2.fastq"
    fastq0 = P.snip(outfile, "_agg-agg-agg")+".0.fastq"

    statement = ''' samtools sort -n %(infile)s -o %(sorted_bam)s;
                    samtools fastq
                        -1 %(fastq1)s
                        -2 %(fastq2)s
                        -0 %(fastq0)s -s /dev/null -n -F 0x900
                        %(sorted_bam)s;
                    paired_end_reads=$(samtools view -c -f 1 %(sorted_bam)s);
                    if [ $paired_end_reads = 0 ]; then
                        salmon quant -i %(salmonIndex)s
                            --libType IU
                            -r %(fastq0)s
                            -o %(outfile)s
                            %(salmon_options)s
                            --numBootstraps %(bootstrap_options)s;
                    else
                        salmon quant -i %(salmonIndex)s
                            --libType IU
                            -1 %(fastq1)s
                            -2 %(fastq2)s
                            -o %(outfile)s
                            %(salmon_options)s
                            --numBootstraps %(bootstrap_options)s;
                    fi; 
                    mv %(outfile)s/quant.sf %(outfile)s.sf; 
                    rm %(fastq1)s; rm %(fastq2)s; rm %(fastq0)s; rm %(sorted_bam)s 
                    '''
    P.run(statement)

#############################################
##### load quantification into database #####
#############################################

@follows(quantifyWithSalmon)
@merge("quantification.dir/*.sf", "database_load/salmon_quant.load")
def mergeAllQuants(infiles, outfile):
    job_threads=3

    P.concatenate_and_load(infiles, outfile,
                         regex_filename="quantification.dir/(.*)_agg-agg-agg.sf",
                         options="-i Name -i Length -i EffectiveLength"
                         " -i TPM -i NumReads -i track"
                         " -i source",
                         job_memory="64G")

    if not os.path.isfile("mapping_rates.txt"):
        statement = ''' bash /shared/sudlab1/General/projects/UTRONs/MyFiles/scripts/mapping_rates_script.sh '''
        P.run(statement)
    else:
        pass

###########################################################
##### Export all_utrons, novel_utrons ids and tx2gene #####
##### text files from utrons database into .txt files #####
###########################################################

@follows(mergeAllQuants, mkdir("expression.dir", "expression.dir/csvdb_files"))
def CSVDBfiles():
    ''' utility function to connect to database. Use this method to connect to the pipeline database.
        Export all_utrons_ids, novel_utrons_ids and tx2gene data in .txt format to be used in the Rscript. 
        '''

    subprocess.call(["sqlite3", PARAMS["database_name"],
                     ".headers on", ".mode tab", ".output expression.dir/csvdb_files/tx2gene.txt",
                     "select transcript_id, match_gene_id from transcript_class where track = 'agg-agg-agg'"])

    subprocess.call(["sqlite3", PARAMS["database_name"],
                     ".headers on", ".mode tab", ".output expression.dir/csvdb_files/all_utrons_ids.txt", 
                     "select * from all_utrons_ids where track = 'agg-agg-agg'"])

    subprocess.call(["sqlite3", PARAMS["database_name"],
                     ".headers on", ".mode tab", ".output expression.dir/csvdb_files/partnered_utrons_ids.txt",
                     "select * from partnered_utrons_ids where track = 'agg-agg-agg'"])

    subprocess.call(["sqlite3", PARAMS["database_name"],
                     ".headers on", ".mode tab", ".output expression.dir/csvdb_files/novel_utrons_ids.txt",
                     "select * from novel_utrons_ids where track = 'agg-agg-agg'"])
    
    shutil.copy("/shared/sudlab1/General/projects/UTRONs/databases/60db_novel_utrons_ids.txt", "expression.dir/csvdb_files/60db_novel_utrons_ids.txt")


###########################################################################
##### Rscript for generating tpm expression values for transcript and ##### 
#####           gene level, as well as fraction expression            #####
###########################################################################

@follows(CSVDBfiles, quantifyWithSalmon)
@merge("quantification.dir/*.sf", ["database_load/utrons_expression.load", "database_load/partnered_utrons_expression.load", "database_load/novel_utrons_expression.load"])
def utrons_expression(infiles, outfiles):
    outfile_load, partnered_load, novel_load = outfiles
    job_threads = 4
    job_memory = "48G"
    
    ### all transcripts - including non-utrons ###
    outfile = "expression.dir/utrons_expression.txt"
    if not os.path.isfile(outfile):
        statement = ''' Rscript /shared/sudlab1/General/projects/UTRONs/MyFiles/scripts/utrons_Rscript.R '''
        P.run(statement)
    else:
        pass
    P.load(outfile, outfile_load, options = "-i Sample -i transcript_id -i gene_id -i tr_expr -i gene_expr -i fract_expr", job_memory="16G")
    
    ### partnered utrons only ###
    partnered = "expression.dir/partnered_utrons_expression.txt"
    if not os.path.isfile(partnered):
        subprocess.call(["sqlite3", PARAMS["database_name"],
                         ".headers on", ".mode tab", ".output expression.dir/partnered_utrons_expression.txt",
                         "SELECT * FROM utrons_expression A WHERE transcript_id in (SELECT transcript_id from partnered_utrons_ids)"])
    else:
        pass
    P.load(partnered, partnered_load, options = "-i Sample -i transcript_id -i gene_id -i tr_expr -i gene_expr -i fract_expr", job_memory="16G")

    ### novel utrons only ###
    novel = "expression.dir/novel_utrons_expression.txt"
    if not os.path.isfile(novel):
        subprocess.call(["sqlite3", PARAMS["database_name"],
                         ".headers on", ".mode tab", ".output expression.dir/novel_utrons_expression.txt",
                         "SELECT * FROM utrons_expression A WHERE transcript_id in (SELECT transcript_id from novel_utrons_ids)"])
    else:
        pass
    P.load(novel, novel_load, options = "-i Sample -i transcript_id -i gene_id -i tr_expr -i gene_expr -i fract_expr", job_memory="16G")

#################################
##### identify splice sites #####
#################################

@follows(find_utrons, mergeAllAssemblies, mergeAllQuants)
@transform("utron_beds.dir/*.bed.gz", regex("(.+)/agg-agg-agg.(.+)_utrons.bed.gz"), [r"expression.dir/\2_splice_sites.txt", r"database_load/\2_splice_sites.load"])
def identify_splice_sites(infiles, outfiles):
    infile=infiles
    outfile, outfile_load = outfiles
    
    statement = ''' python /shared/sudlab1/General/projects/stem_utrons/pipelines/pipeline_DTU/splicesites_start_end_sizes.py %(infile)s %(outfile)s;
                    sort -u %(outfile)s > %(outfile)s_2.txt; rm %(outfile)s; mv %(outfile)s_2.txt %(outfile)s;
                    sed -i $'1i transcript_id\\tstrand\\tss5\\tss3\\tcontig\\tsplice_site_start\\tsplice_site_end\\tutron_size' %(outfile)s '''
    P.run(statement)   
    P.load(outfile, outfile_load, options = "-i transcript_id -i ss5 -i ss3 -i splice_site_start -i splice_site_end -i utron_size", job_memory="16G")

####################################################################################################
#####               Extract gene_id, gene_name, strand, chr, start, end from the               #####
#####  /shared/sudlab1/General/annotations/hg38_noalt_ensembl85/ensembl.dir/geneset_all.gtf.gz #####
####################################################################################################

@follows(find_utrons, mergeAllAssemblies, mergeAllQuants)
@transform(PARAMS["annotations_interface_geneset_all_gtf"], regex("(.+).gtf.gz"), "expression.dir/gtf_stop_codons.txt")
def gtf_stop_codons(infile, gtf):
    outfile = open(gtf, "w")
    for line in DataIterator(infile):
        if line.featuretype == "stop_codon":
            if line.strand == "+":
                outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (''.join(line.attributes['gene_id']),''.join(line.attributes['transcript_id']), ''.join(line.attributes['gene_name']), line.strand, line.featuretype, line.seqid, line.start, line.end))
            elif line.strand =="-":
                outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (''.join(line.attributes['gene_id']), ''.join(line.attributes['transcript_id']), ''.join(line.attributes['gene_name']), line.strand, line.featuretype, line.seqid, line.end, line.start))
        else:
            pass
    outfile.close()

#####################################
#####    --- DTU SCRIPTS ---    #####
##### Determine how complex the #####
#####   design matrix is, and   #####
##### execute neccessary script #####
###################################################################################
###   In order for results to be interactive, this part will generate a folder  ###
### for each set of comparisons, and copy the Rmd into each of these folders so ### 
###      that they can be interogated manually/individually. The data from      ###
###          each of these will be compiled into a final output report.         ###
###################################################################################

@follows(gtf_stop_codons, mkdir("DTU.dir"))
@transform("design.tsv", suffix("design.tsv"), "DTU.dir/comparisons_to_make.txt") 
def analyseDesignMatrix(infile, outfile):
    infile = open(infile)
    design = csv.reader(infile, delimiter="\t")
    outfile = open(outfile, "w")
    
    for row in design:
        variables = len(row) - 1 
        if(variables == 1):
            folder_name = str(row[0] + "_vs_" + str(row[1]))
            formula = "~var1"
        elif(variables == 2):
            covar1, covar2 = str(row[2]).split(':', 1)
            folder_name = str(row[0] + "_vs_" + str(row[1])) + "_PLUS_" + covar1+"_vs_"+covar2
            formula = "~var1 + var2 + var1:var2"
        else:
            raise valueError("design.tsv is not configured correctly, please refer to pipeline_DTU/example_design.tsv for comparison")

        outfile.write(folder_name + "\t" + formula + "\n")

        to_cluster=False
        statement = """ mkdir DTU.dir/""" + folder_name + """ && echo \"""" + formula + """\" > DTU.dir/""" + folder_name + "/" + folder_name + """.info"""
        os.system(statement)
        
    infile.close
    outfile.close()


#####################
##### run fxns ######
#####################

@follows(Assembly, AnnotateAssemblies, export, mergeAllQuants, CSVDBfiles, utrons_expression, identify_splice_sites, gtf_stop_codons)
def DTUtrons():
    #decorate this with pipeline_utrons specific things
    pass
@follows(Assembly, AnnotateAssemblies, export, mergeAllQuants, CSVDBfiles, utrons_expression, identify_splice_sites, gtf_stop_codons)
def DTU():
    #decorate this with DTU specific things
    pass

##################
###### misc ######
##################

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))