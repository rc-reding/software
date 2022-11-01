#!/usr/bin/env python3.9

from multiprocessing import cpu_count
from helper_func import process_reads, retrieve_reads, retrieve_depth,\
                        extract_genes_seq, house_keeping, collate_gene_seq,\
                        analyse_seqs_found, plot_analysis, sanitise_options,\
                        parse_directories_init
import os, argparse
import sys  # temporal?


def main(Opts):
    """
        Main function
    """
    # Check Opts consistency
    options = sanitise_options(Opts)

    # https://github.com/ablab/spades/issues/19#issuecomment-631981051
    # SPAdes allocates fixed amount of memory per thread for some stages
    # of the assembly process. Reduce N_CPU from to circumvent errors of
    # memory limitation.
    N_CPU = cpu_count() - 2

    DIR_LIST = os.listdir(options.PATH)
    DIR_LIST = sorted([DIR for DIR in DIR_LIST if DIR.find(".") == -1 and\
                        DIR.find("reads") == -1])

    for DIR in DIR_LIST:
        PRJ_NAME, READS_PATH, DIR_DEPTH = parse_directories_init(options.PATH, DIR)

        if options.ASSEMBLE is True:
            # if options.CALC_COVERAGE is True:
            #     coverage_data_found = check_coverage_data_exists()
            #     if coverage_data_found is False:
            #         raise ValueError("MSG HERE")
            rds_count = 1
            illumina_rds = retrieve_reads(path=READS_PATH)
            for rd in illumina_rds[0:-1:2]:
                print('\nGenerating %i of %i assemblies.' % (rds_count,
                                                             len(illumina_rds[0:-1:2])))
                if options.ASSEMBLE_DENOVO is True:  # ASSEMBLE DE NOVO
                    ASSEMBLY_PATH,\
                        DEPTH_PATH,\
                            REFERENCE_PATH = process_reads(rd, options.REF_GENOME,
                                                           PRJ_NAME, N_CPU,
                                                           options.ASSEMBLE_DENOVO,
                                                           options.CALC_COVERAGE)
                else:  # ASSEMBLE USING REFERENCE GENOME
                    # FIX: adapt to reads mapped to a reference
                    sys.exit(0)
                    ASSEMBLY_PATH,\
                        DEPTH_PATH,\
                            REFERENCE_PATH = process_reads(rd, options.REF_GENOME,
                                                           PRJ_NAME, N_CPU,
                                                           options.ASSEMBLE_DENOVO,
                                                           options.CALC_COVERAGE)
                # Update count
                rds_count += 1

            # Remove alignment files to save space? (Take up more than 85% of space)
            house_keeping(ASSEMBLY_PATH, alignments=True)
        else:
            if os.path.exists(str("/").join([options.PATH, "assemblies"])) is True:
                ASSEMBLY_PATH = str("/").join([options.PATH, DIR, "de_novo/"])
                # UNMAPPED_PATH = str("/").join([options.PATH, DIR, "from_unaligned/"])
            else:
                ASSEMBLY_PATH = str("/").join([options.PATH, DIR, "assemblies/de_novo/"])
                # UNMAPPED_PATH = str("/").join([options.PATH, DIR, "assemblies/from_unaligned/"])
            house_keeping(ASSEMBLY_PATH, alignments=True)

        if options.TARGET_GENES is not None:
            # Extract information for genes of interest
            assemblies = [file for file in os.listdir(ASSEMBLY_PATH)
                            if file.find(".fa") > -1 and file.find(".fa.") == -1]
            RD_DEPTH_STATS = dict() if options.CALC_COVERAGE is True else None
            for assembly in assemblies:  # Force extracting promoter + CDS
                loci_metadata, hk_metadata, SEQ_ID, PRJ_PATH,\
                    GENES_LIST = extract_genes_seq(ASSEMBLY_PATH + assembly,
                                                   options.TARGET_GENES,
                                                   options.HK_GENES, PRJ_NAME,
                                                   options.PREFIX, DIR_DEPTH,
                                                   N_CPU+2, options.SEQ_CUTOFF,
                                                   options.ID_THRESHOLD,
                                                   DENOVO=options.ASSEMBLE_DENOVO,
                                                   CALC_COVERAGE=options.CALC_COVERAGE)
                if loci_metadata is not None and options.CALC_COVERAGE is True:
                    # Extract coverage depth for genes of interest + MLST genes
                    RD_DEPTH,\
                        STATS = retrieve_depth((ASSEMBLY_PATH + assembly,
                                                loci_metadata, SEQ_ID,
                                                options.SEQ_CUTOFF, hk_metadata),
                                                options.PREFIX, PRJ_PATH)
                    if STATS[0] is not str('NaN'):
                        RD_DEPTH_STATS[assembly] = (STATS)  # mean, std

    # Do not move inside for-loop, or will generate a phylogeny/plot per DIR
    # while updating GENES_LIST in the background if DIR_LIST has multiple DIR
    if options.PHYLOGENY is True:
        phylogeny, alignment_file = analyse_seqs_found(GENES_LIST, N_THREADS=N_CPU+2,
                                               FLT_THR=options.FILTER_THRESHOLD)

    if options.PLOT_FNAME is not None:
        # If options.PHYLOGENY is not given, discover phylogeny generated
        if options.PHYLOGENY is False:
            alignment_file, phylogeny = retrieve_phylogeny(options.PATH,
                                                           options.PREFIX)
            RD_DEPTH = None
            RD_DEPTH_STATS = None

        # RD_DEPTH needed by plot_analysis routine, correct if neccessary
        if options.CALC_COVERAGE is False:
            RD_DEPTH = None
            RD_DEPTH_STATS = None

        plot_analysis(alignment_file, phylogeny, RD_DEPTH, RD_DEPTH_STATS,
                      options.PLOT_FNAME, options.PREFIX,
                      CUTOFF=options.SEQ_CUTOFF, ROI=options.ROI_COORDS,
                      PROMOTER=options.PROMOTER, LABELS=options.XLS_DB)


if __name__ == '__main__':
    ## Retrieve options from CLI ##
    Args = argparse.ArgumentParser(description="HOUND: species-independent \
                                                gene identification system",
                                   epilog="Carlos Reding (c) Copyright 2022. \
                                   Software developed as part of FARM-SAFE \
                                   project (BBSRC grant BB/T004592/1 & Arwain \
                                   XXXX) at the University of Bristol.")

    Args.add_argument("--readsdir", metavar="DIR", type=str, nargs=1, required=True,
                      dest="PATH", help="Directory where FASTQ files can be \
                      found. It can be a directory of directories if FASTQ \
                      files are contained in a 'reads' directory. Maximum \
                      directory depth is 2.")

    Args.add_argument("--assemble", action='store_true', dest="ASSEMBLE",
                      help="Assemble reads.")

    Args.add_argument("--reference", metavar="FILE", type=str, nargs=1,
                      dest="REF_GENOME", help="Reference genome in FASTA format.\
                       Requires --assemble.")

    Args.add_argument("--de-novo", action='store_true', dest="ASSEMBLE_DENOVO",
                      help="Assemble reads de novo. Used to assemble genomes and\
                      specify assembly type for data analysis. Requires \
                      --assemble.")

    Args.add_argument("--coverage", action='store_true', dest="CALC_COVERAGE",
                      help="Compute coverage depth to estimate gene copy number.\
                       Can be used to assemble genomes or include coverage depth\
                       in the data analysis. Requires --hk-genes.")

    Args.add_argument("--hk-genes", metavar="FILE", type=str, dest="HK_GENES",
                      help="List of Multilocus sequence typing (MLST) \
                      genes, or other reference genes, in FASTA format to \
                      compute baseline coverage depth. Requires --coverage.")

    Args.add_argument("--genes", metavar="FILE", type=str, dest="TARGET_GENES",
                      help="List of genes to be found, in FASTA \
                      format. Requires --identity and --prefix.")

    Args.add_argument("--prefix", metavar="NAME", type=str, dest="PREFIX",
                      nargs=1, help="Label added to all output files. Required \
                      to do multiple searches with the same assemblies.")

    Args.add_argument("--identity", metavar="NUM", type=float, nargs=2,
                      dest="ID_THRESHOLD", help="Identity threshold required \
                      to shortlist sequences found. Two floats (min identity, \
                      max identity) between 0 and 1 are required. Requires \
                      --genes.")

    Args.add_argument("--promoter", action='store_true', dest="PROMOTER",
                      help="Isolate the promoter region of the target gene(s) \
                      sequences found and ignore coding sequences. Requires \
                      --cutoff.")

    Args.add_argument("--cutoff", metavar="NUM", type=int, default=0,
                      dest="SEQ_CUTOFF", help="Length of the promoter in \
                      nucleotides. Requires --promoter.")

    Args.add_argument("--phylo", action='store_true', dest="PHYLOGENY",
                      help="Align sequences of promoter/coding sequences found,\
                       and generate the corresponding phylogeny.")

    Args.add_argument("--phylo-thres", metavar="NUM", dest="FILTER_THRESHOLD",
                      type=float, nargs=1, default=0.5, help="Remove sequences \
                      that are a fraction of the total size of alignment. Used \
                      to improve quality alignment. Requires a number between 0\
                      and 1 (defaults to ""0.5"").")

    Args.add_argument("--plot", metavar="FILE", type=str, dest="PLOT_FNAME",
                      help="Generate plot from the multiple alignment of \
                      sequences found, and save as FILE.")

    Args.add_argument("--roi", metavar="FILE", type=str,
                      dest="ROI_COORDS", help="Sequences of interest to look \
                      for in the gene(s) found, in FASTA format. Requires \
                      --plot.")

    Args.add_argument("--labels", metavar="FILE", type=str,
                      dest="XLS_DB", help="XLS file containing assembly name \
                      (col 1), and assembly type (col 6) to label phylogeny \
                      leafs (defaults to assembly name). Requires --plot.")
    ### Run main program ###
    Opts = Args.parse_args()
    main(Opts)
# Step 1: List directory with files. [DONE]
# Step 2: Quality Control [SKIP, ALREADY DONE BY MICROBESNG]
# Step 3: Trim adaptors [SKIP, ALREADY DONE BY MICROBESNG]
# St#p 4: Align to reference. [DONE]
# Step 5: Compute read coverage to detect/visualise gene loss/amplification [DONE]
# Step 6: Filter un-aligned reads [DONE]
# Step 7: Assemble unaligned reads [DONE]
# Step 8: Gene detection using local BLAST database [DONE]
# Step 9: Parse BLAST report to detect i) genes and ii) extract their location [DONE]

# TODO promptly!
# Step 10: Detect mutations with respect to prom seq. given, and highlight it in the plot [DONE]
# Step 11: Accept user-given sequences, rather than hardcoded promoter sequences.
# Step 12: Output a CSV file with assembly ID, locus name, etc. to allow different analysis
# Step 13: Current plot is okay as demonstration, but make optional so people can do what they want.

# TODO LATER.
# Step 14: Compute INDELs and structural variations
# Step 15: Compute single-nucleotide polymorphisms (SNPs)
