import os
from Bio import Align, AlignIO
from .assembly_utilities import _calculate_median
from .sundry import _find_tool, _is_indexed, _cleanup_path,\
                    _check_project_exists


def _rtrv_reference(organism_id: str) -> str:
    """
        Retrieve reference sequence given an `organism_id', and return
        the absolute path.
    """
    REF_PATH = os.path.abspath(str("../src_references/"))  # FIX: Add as parameter?

    # Filter .fasta files (in case sequences are referenced)
    file_list = os.listdir(REF_PATH)
    ref_seqs = [file for file in file_list if file.find(organism_id) > -1]
    ref_file = [file for file in ref_seqs if file.count(".") == 1]

    # Check if sequence is indexed
    _is_indexed(str("/").join([REF_PATH, ref_file[0]]), ref_seqs)
    return str("/").join([REF_PATH, ref_file[0]])


def map_denovo_reads(illumina_rds: str, ref_assembly: str = None,
                     N_THREADS: int = 1, project_name: str = None) -> str:
    """
        Map illumina reads to de novo assembly produced with SPAdes
    """
    # Check project has a name and a directory
    rds_dir, rd_fName, project_path, ALN_dir,\
            _, ASSEMBLY_dir = _check_project_exists(illumina_rds,
                                                            project_name)
    # Tools
    samtools = _find_tool('samtools')
    samtools_convers_args = list(["view","-b", "--verbosity 1",
                                  "--threads " + str(N_THREADS)])
    samtools_sort_coord_args = list(["sort", "--verbosity 1",
                                     "--threads " + str(N_THREADS)])

    bwa = _find_tool('bwa')
    bwa_aln_args = list(["mem", "-v 1", "-t " + str(N_THREADS)])

    # Define output filenames first
    outALN = str("/").join([project_path, ALN_dir, rd_fName])  # Step 1 (DO NOT CHANGE)
    outALN = outALN.replace("_1_trimmed.fastq.gz", ".sam")  # Step 2
    outBAM = outALN.replace(".sam", ".bam")

    outFINAL = outALN.replace(".sam", "_assembly.bam")  #  Step 1
    outFINAL = outFINAL.replace(ALN_dir, ASSEMBLY_dir)  # Step 2
    outSTATS = outFINAL.replace(".bam", ".stats")

    # Index reference assembly (single-thread process)
    os.system(bwa + str(" ").join(["index", ref_assembly]))

    # Align reads to assembly (multiprocess)
    rds_File2 = illumina_rds.replace("_1_", "_2_")
    bwa_aln_args.extend(["-o " + outALN, ref_assembly, illumina_rds, rds_File2])
    os.system(bwa + str(" ").join(bwa_aln_args))

    # Convert alignment from SAM to BAM (multiprocess)
    samtools_convers_args.extend([outALN, "-o " + outBAM])
    os.system(samtools + str(" ").join(samtools_convers_args))

    # Sort reads now by coordinate for further processing (multiprocess)
    # (needed to calculate coverage)
    samtools_sort_coord_args.extend(["-o " + outFINAL, outBAM])
    os.system(samtools + str(" ").join(samtools_sort_coord_args))

    # Remove temporary files
    tmp_files = list([outBAM, outALN])
    [os.remove(file) for file in tmp_files]

    idx_dir = str("/").join(ref_assembly.split("/")[:-1])
    idx_template = ref_assembly.split("/")[-1] + str(".")
    idx_tmp_files = [f for f in os.listdir(idx_dir) if f.find(idx_template) > -1]
    [os.remove(idx_dir + '/' + file) for file in idx_tmp_files]
    return outFINAL


def align_reads(illumina_rds: str, ref_chromosome: str = None,
                N_THREADS: int = 1, project_name: str = None) -> str:
    """
        Align illumina reads to reference genome. This genome must be
        previously indexed.
    """
    if ref_chromosome is None:
        raise ValueError('Reference genome cannot be "NoneType".')
    else:
        # Check project has a name and a directory
        rds_dir, rd_fName, project_path, ALN_dir,\
                UNALN_dir, ASSEMBLY_dir = _check_project_exists(illumina_rds,
                                                                project_name)
        # Tools
        samtools = _find_tool('samtools')
        samtools_convers_args = list(["view","-b", "--verbosity 1",
                                      "--threads " + str(N_THREADS)])
        samtools_sort_coord_args = list(["sort", "--verbosity 1",
                                         "--threads " + str(N_THREADS)])
        samtools_sort_name_args = list(["sort", "--verbosity 1", "-n",
                                        "--threads " + str(N_THREADS)])
        samtools_fixmate_args = list(["fixmate", "-m", "--verbosity 1",
                                      "--threads " + str(N_THREADS)])
        samtools_calmd_args = list(["calmd", "-b", "--verbosity 1",
                                    "--threads " + str(N_THREADS)])

        bwa = _find_tool('bwa')
        bwa_aln_args = list(["mem", "-v 1", "-t " + str(N_THREADS)])

        # Process...
        REF_SEQ = _rtrv_reference(ref_chromosome)

        # Define output filenames first
        outALN = str("/").join([project_path, ALN_dir, rd_fName])  # Step 1 (DO NOT CHANGE)
        outALN = outALN.replace("1_trimmed.fastq.gz", "Aligned.sam")  # Step 2

        outBAM = outALN.replace(".sam", ".bam")
        outCoordSorted = outALN.replace(".sam", ".sorted")
        outNameSorted = outALN.replace(".sam", ".sorted.name")
        outFixmate = outALN.replace(".sam", ".fixmate")
        outCALM = outALN.replace(".sam", ".calmd.bam")
        outUNMAPPED = outCALM.replace(ALN_dir, UNALN_dir)  # Use this fName to process unaligned reads

        outFINAL = outALN.replace("Aligned.sam", "Assembly.bam")  #  Step 1
        outFINAL = outFINAL.replace(ALN_dir, ASSEMBLY_dir)  # Step 2
        outSTATS = outFINAL.replace(".bam", ".stats")

        # Align reads to reference (multiprocess)
        rds_File2 = illumina_rds.replace("_1_", "_2_")
        bwa_aln_args.extend(["-o " + outALN, REF_SEQ, illumina_rds, rds_File2])
        os.system(bwa + str(" ").join(bwa_aln_args))

        # Convert alignment from SAM to BAM (multiprocess)
        samtools_convers_args.extend(["-T " + REF_SEQ, outALN, "-o " + outBAM])
        os.system(samtools + str(" ").join(samtools_convers_args))

        # Sort reads by name to allow removing duplicates (multiprocess)
        # This extra step is required since rmdup is DEPRECATED
        samtools_sort_name_args.extend(["-o " + outNameSorted, outBAM])
        os.system(samtools + str(" ").join(samtools_sort_name_args))

        # Add ms and MC score to remove duplicates with markdup
        # This extra step is required since rmdup is DEPRECATED
        samtools_fixmate_args.extend([outNameSorted, outFixmate])
        os.system(samtools + str(" ").join(samtools_fixmate_args))

        # Sort reads now by coordinate for further processing (multiprocess)
        # (needed to calculate coverage)
        samtools_sort_coord_args.extend(["-o " + outCoordSorted, outFixmate])
        os.system(samtools + str(" ").join(samtools_sort_coord_args))

        # Mark mismatches and insertions
        samtools_calmd_args.extend([outCoordSorted, REF_SEQ, "> " + outCALM])
        os.system(samtools + str(" ").join(samtools_calmd_args))

        # Remove any PCR duplicates left (rmdup is DEPRECATED, use markdup)
        samtools_markdup_args = list(["markdup", "--threads " + str(N_THREADS),
                                      "-r", outCALM, outFINAL])
        os.system(samtools + str(" ").join(samtools_markdup_args))

        # Index assembly and generate stats
        samtools_index_args = list(["index", "-@ " + str(N_THREADS), outFINAL])
        os.system(samtools + str(" ").join(samtools_index_args))

        samtools_stats_args = list(["flagstat", "--threads " + str(N_THREADS),
                                    outFINAL, ">", outSTATS])
        os.system(samtools + str(" ").join(samtools_stats_args))

        # Temporary files are removed at a later stage
    return outFINAL, REF_SEQ, outUNMAPPED


def collate_gene_seq(DIR_LIST: list, PATH: str, DE_NOVO: bool = False) -> str:
    """
        Scan all sub-folders and collate all sequences (FASTA) filtered
        into one file, easier to handle for downstream analysis
        than having to continuously screen through all sub-folder
        structure.
    """
    if DE_NOVO is True:
        input_file = str("genes_found_denovo_assembly.fa")
        coverage_file = str("coverage_denovo_genes.txt")
    else:
        input_file = str("genes_found_assembly.fa")
        coverage_file = str("coverage_genes.txt")

    target_file = str("genes_found_compilation.fa")
    coverage_target_file = str("coverage_compilation.txt")

    for DIR in DIR_LIST:
        if os.path.exists(PATH + str("/") + DIR + str("/") + input_file):
            # Collate sequences
            fIn = PATH + str("/") + DIR + str("/") + input_file
            fOut = PATH + str("/") + target_file
            os.system('touch ' + fOut)
            os.system('cat ' + fIn + ' >> ' + fOut)
            # Collate coverage
            fIn_cov = PATH + str("/") + DIR + str("/") + coverage_file
            fOut_cov = PATH + str("/") + coverage_target_file
            os.system('touch ' + fOut_cov)
            os.system('cat ' + fIn_cov + ' >> ' + fOut_cov)
    return


def _mark_duplicates(sequence_list: str) -> bool:
    """
        Screens through entry names in an alignment to detect duplicated names,
        and add 'LOCX' (locus + number) to avoid errors downstream.
    """
    sequences = AlignIO.read(sequence_list, 'fasta')
    renamed_sequences = Align.MultipleSeqAlignment([])

    n_copies = dict()
    seq_seen = set()
    # with open(sequence_list, 'w') as fOut:
    for seq in sequences:
        # Has this sequence been seen?
        if seq.id.split("__")[0] in seq_seen:
            # Mark number of copies if needed
            if seq.id.split("__")[0] in n_copies.keys():
                n_copies[seq.id.split("__")[0]] += 1
            else:
                n_copies[seq.id.split("__")[0]] = 2
            # Re-name duplicated entry to avoid issues with 'phyml'
            new_id = seq.id.split("__")[0] + str("_L") +\
                     str(n_copies[seq.id.split("__")[0]]) +\
                     str("__") + seq.id.split("__")[-1]
            seq.id = new_id
        # Annotate sequence has been seen
        seq_seen.add(seq.id.split("__")[0])  # Here first, or ERR
        renamed_sequences.append(seq)

    # Write sequence
    AlignIO.write(renamed_sequences, sequence_list, 'fasta')
    return


def multiple_seq_alignment(seqs_file: str, FILTER_THRESHOLD: float,
                           N_THREADS: int) -> str:
    """
        Filter sequences that are too short to be meaninful, and
        use clustalO/muscle to align the remainder.
    """
    from Bio import SeqIO
    from Bio.Align.Applications import MuscleCommandline

    # First: filter sequences, if it hasn't been done already
    filtered_seqs_file = seqs_file.replace(".fa", "_FILTERED.fa")
    if os.path.exists(filtered_seqs_file) is False and\
        os.lstat(seqs_file).st_size > 0:
        # Calculate average sequence length (not worth importing numpy)
        seqs = SeqIO.parse(seqs_file, 'fasta')  # Iterator!
        seq_len = [len(entry.seq) for entry in seqs]
        median_seq_len = _calculate_median(seq_len)

        # Remove sequences that are too short to be meaninful
        seqs = SeqIO.parse(seqs_file, 'fasta')  # Re-generate iterator
        filtered_seqs = [entry for entry in seqs\
                            if len(entry.seq) >= median_seq_len * FILTER_THRESHOLD]
        # Write filtered sequences to disk
        SeqIO.write(filtered_seqs, filtered_seqs_file, 'fasta')

    if os.lstat(seqs_file).st_size > 0:
        # Now: align sequences. TODO: adapt for clustalO for multiprocessing as muscle5 has no gapopen penalty.
        muscle3_bin = _find_tool('muscle')[:-1]  # Remove final space, ERR otherwise
        alignment_file = filtered_seqs_file.replace(".fa", ".aln")
        muscle3_cmd = MuscleCommandline(cmd=muscle3_bin, input=filtered_seqs_file,
                                        out=alignment_file, gapopen=-1950.0)
        muscle3_cmd()
        return alignment_file
    else:
        print('*** NO HITS FOUND/ADDED, NOTHING TO ALIGN. ***\n')
        return None


def retrieve_phylogeny(PATH: str, PREFIX: str):
    """
        Discover alignment and phylogeny files in PATH.
    """
    FILE_LIST = os.listdir(PATH)
    ALN_EXT = str(".aln")  # Safe, autogenerated by alignment routine
    TREE_EXT = str("tree.txt")  # Safe, autogenerated by phylogeny routine

    # Filter only those files with PREFIX
    FILTERED_LIST = [FILE for FILE in FILE_LIST if PREFIX in FILE]

    alignment_file = [FILE for FILE in FILTERED_LIST if ALN_EXT in FILE][0]  # TODO: check there is onle 1
    tree_file = [FILE for FILE in FILTERED_LIST if TREE_EXT in FILE][0]  # TODO: check there is onle 1

    # Import phylogeny if it exists
    if len(tree_file) == 1:
        tree = PhyloTree(str("/").join([PATH, tree_file]))
        putative_outgroup = tree.get_farthest_leaf()[0].name
        tree.set_outgroup(putative_outgroup)
        tree.ladderize(direction=1)
        return str("/").join([PATH, alignment_file]), tree
    else:
        print('*** NO HITS FOUND/ADDED, NOTHING TO ALIGN. ***\n')
        return None, None
