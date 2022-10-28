from .unalignment_tools import process_unaligned
from .promoter_tools import analyse_seqs_found, plot_analysis
from .rtrv_reads import rtrv_reads as retrieve_reads
from .ncbi import find_amr_genes, compile_genes_detected
from .alignment_tools import align_reads, multiple_seq_alignment,\
                             collate_gene_seq, map_denovo_reads
from .assembly_utilities import retrieve_depth, calculate_assembly_depth,\
                                assemble_denovo
from argparse import Namespace
import os, warnings


def sanitise_options(options):
    """
        Check that options given by user are consistent with
        their dependencies.
    """
    import sys
    sys.tracebacklimit = 0  # No trackback needed FOR THESE errors.

    if options.PATH is None:
        raise ValueError("--readsdir DIRECTORY is required to do anything with\
                            this program. Use --help for instructions.")
    else:
        # Check if PATH has '/' and if so, remove it. Can mess with assembly.
        if options.PATH[-1] == str("/"):
            options.PATH = PATH[:-1]

    if options.PATH is not None and options.ASSEMBLE is False and\
        options.TARGET_GENES is None:
        raise ValueError("""Don't know what to do with the reads. Use --assemble
                            or --genes [FILE], or use --help for instructions.""")

    if options.ASSEMBLE is True and options.ASSEMBLE_DENOVO is False and\
        options.REF_GENOME is None:
        raise ValueError("How do I assemble the reads? Use --help for instructions.")
    # if options.ASSEMBLE_DENOVO is True and options.ASSEMBLE is False:
    #     raise ValueError("--de-novo requires --assemble. Use --help for instructions.")

    if options.CALC_COVERAGE is True and options.HK_GENES is False:
        raise ValueError("""--coverage requires reference/housekeeping genes to
                            estimate baseline coverage depth. Use --help for
                            instructions.""")
    if options.HK_GENES is True and options.CALC_COVERAGE is False:
        warnings.warn("""--hk-genes were given without --coverage. This option
                        will be ignored.""")

    if options.REF_GENOME is not None and options.ASSEMBLE is False:
        raise ValueError("--reference [FILE] requires --assemble.")

    if options.TARGET_GENES is not None and options.ID_THRESHOLD is None:
        raise ValueError("""--genes [FILE] requires --identity MIN MAX to set
                            gene selection window. Use --help for instructions.""")

    if options.TARGET_GENES is not None and options.PREFIX is None:
        raise ValueError("""--genes [FILE] requires --prefix NAME to label output
                            files. Use --help for instructions.""")

    if options.ID_THRESHOLD is not None and options.TARGET_GENES is None:
        raise ValueError("""--identity requires --genes [FILE]. Use --help for
                            instructions.""")

    if options.SEQ_CUTOFF > 0 and options.PROMOTER is False:
        raise ValueError("""--cutoff requires --promoter. Use --help for
                            instructions.""")

    if options.ROI_COORDS is not None and options.PLOT_FNAME is False:
        raise ValueError("""--roi [FILE] requires --plot [FILE]. Use --help for
                            instructions.""")

    if options.XLS_DB is not None and options.PLOT_FNAME is None:
        raise ValueError("""--labels FILE requires --plot [FILE]. Use --help for
                            instructions.""")


    if options.ID_THRESHOLD is not None:
        options.ID_THRESHOLD = tuple(options.ID_THRESHOLD)

    if options.PATH[0][-1] == str('/'):
        options.PATH = os.path.abspath(options.PATH[0][:-1])
    else:
        options.PATH = os.path.abspath(options.PATH[0])

    # TODO: CHECK THEY ARE INDEED FASTA FILES
    if options.TARGET_GENES is not None:
        options.TARGET_GENES = os.path.abspath(options.TARGET_GENES[0])

    if options.REF_GENOME is not None:
        options.REF_GENOME = os.path.abspath(options.REF_GENOME[0])

    if options.PREFIX is not None:
        options.PREFIX = options.PREFIX[0]

    # Restore Traceback for debugging
    del sys.tracebacklimit

    return options


def house_keeping(assembly_path: str, alignments: bool = False) -> str:
    """
        Remove temporary files that are not longer needed after assembly.
        BLAST index files are temporary but they will be re-used with
        every new gene search---and they can be slow to produce (2-4s each).
        However BLAST output files (.blast) can cause trouble when using the
        same assemblies to look for different genes. Fast to produce (<<1s),
        delete.
    """
    if alignments is True:
        aln_path = assembly_path.replace("/assemblies", "/alignments")
        if aln_path.split("/")[-1] != str('alignments'):
            # Is this a case of '/assemblies/de_novo' to '/alignments/de_novo'?
            # If so, remove last folder from name to retrieve abs path for 'alignments'
            aln_path = str("/").join(aln_path.split("/")[:-2])
        os.system("rm -Rf " + aln_path)
    else:
        blast_output = assembly_path.replace(".fa", ".blast")
        if os.path.exists(blast_output) is True:
            os.system("rm " + blast_output)


def parse_directories_init(PATH, DIR):
    """
        Checks the right folder structure exists, and returns
        project name, reads directory, and directory depth.
    """
    PRJ_NAME = DIR  # Keep this to make parsing results easier

    if os.path.exists(str("/").join([PATH, "reads"])) is True:
        READS_PATH = str("/").join([PATH, "reads/"])
        PRJ_NAME = PATH.split('/')[-1]  # Avoids saving files in 'reads/'
        DIR_DEPTH = 1
    elif os.path.exists(str("/").join([PATH, DIR, "reads"])) is True:
        READS_PATH = str("/").join([PATH, DIR, "reads/"])
        DIR_DEPTH = 2
    else:
        second_dir = os.listdir(str("/").join([PATH, DIR]))[0]
        if os.path.exists(str("/").join([PATH, DIR, second_dir, "reads"])) is True:
            READS_PATH = str("/").join([PATH, DIR, second_dir, "reads/"])
            DIR_DEPTH = 3
        else:
            try:
                third_dir = os.listdir(str("/").join([PATH, DIR, second_dir]))[0]
                READS_PATH = str("/").join([PATH, DIR, third_dir, "reads/"])
                DIR_DEPTH = 4
            except ValueError:
                raise ValueError("Max directory depth reached with nothing \
                                 found. Re-arrange your directories and try again.")
    return PRJ_NAME, READS_PATH, DIR_DEPTH


def process_reads(illumina_rd: str, REFERENCE: str, PRJ_NAME: str, N_CPU: int,
                  denovo_assembly: bool, calculate_coverage: bool) -> bool:
    """
        Process illumina read to map them to a reference genome, assemble de
        novo those that could not be mapped to the reference, and find genes
        of interest using BLAST.
    """
    if denovo_assembly is False:
        # Align reads to reference if given.
        if REFERENCE is None:
            raise ValueError("""Trying to align to a reference, but no reference
                                genome given. Provide a reference in FASTA
                                format or use --de-novo. Use --help for
                                instructions.""")

        assembly, REFERENCE_PATH,\
            unaligned_template = align_reads(illumina_rd, ref_chromosome=REFERENCE,
                                             project_name=PRJ_NAME, N_THREADS=N_CPU)
        # Process those reads that could not be aligned, and assemble de novo
        unaligned_reads = process_unaligned(unaligned_template, N_THREADS=N_CPU)
        unmapped_assembly = assemble_denovo(unaligned_reads, assemble_unaligned,
                                            project_name=PRJ_NAME,
                                            N_THREADS=N_CPU)
        return str("/").join(assembly.split("/")[:-1]),\
               str("/").join(unmapped_assembly.split("/")[:-1]), REFERENCE_PATH
    else:
        # Assembly de novo without aligning to reference
        denovo_assembly,\
            assembly_success = assemble_denovo(illumina_rd, N_THREADS=N_CPU,
                                               project_name=PRJ_NAME)
        if calculate_coverage is True:
            if assembly_success is True:
                mapped_assembly = map_denovo_reads(illumina_rd, denovo_assembly,
                                                   N_THREADS=N_CPU,
                                                   project_name=PRJ_NAME)
                depth_file = calculate_assembly_depth(mapped_assembly,
                                                      N_THREADS=N_CPU)
            else:
                # Assembly failed due to _very_ low quality reads
                # Generate depth file (empty) for convenience.
                depth_file = denovo_assembly.replace("de_novo/", "").replace(".fa", "_depth.txt")
                os.system("touch " + depth_file)
            return str("/").join(denovo_assembly.split("/")[:-1]) + str("/"),\
                   str("/").join(depth_file.split("/")[:-1]), None
        else:
            if assembly_success is True:
                mapped_assembly = map_denovo_reads(illumina_rd, denovo_assembly,
                                                   N_THREADS=N_CPU,
                                                   project_name=PRJ_NAME)
            return str("/").join(denovo_assembly.split("/")[:-1]) + str("/"),\
                   None, None


def extract_genes_seq(assembly: str, TARGET_GENES: str, project_name: str,
                      PREFIX: str, DIR_DEPTH: int, N_THREADS: int,
                      seq_cutoff: int, ID_THRESHOLD: tuple, DENOVO: bool = False,
                      CALC_COVERAGE: bool = False) -> str:
    """
        Process assembly folder to retrieve BLAST output files, and parse
        them to create a FASTA file containing the nucleotidic sequence of
        all genes found in the assembly. This FASTA file will later be used
        to compare the sequences across multiples assemblies.
    """
    # Wipe previous BLAST results to avoid parsing errors
    # when using the same assembly for look for different genes
    house_keeping(assembly)

    # Find matches
    find_amr_genes(assembly, TARGET_GENES, CALC_COVERAGE, project_name, N_THREADS)

    blast_output = assembly.replace(".fa", ".blast")
    gene_metadata, housekeeping_metadata, SEQ_ID, PRJ_PATH,\
        genes_list = compile_genes_detected(assembly, blast_output,
                                            TARGET_GENES, PREFIX, DIR_DEPTH,
                                            CALC_COVERAGE, seq_cutoff,
                                            ID_THRESHOLD, DENOVO)
    # Genes found?
    if os.lstat(blast_output).st_size > 0:
        return gene_metadata, housekeeping_metadata, SEQ_ID, PRJ_PATH, genes_list
    else:
        # If genes not found, create empty file to acknowledge
        # there was a search but nothing found.
        os.system("touch " + genes_list)
        return None, None, SEQ_ID, PRJ_PATH, genes_list
