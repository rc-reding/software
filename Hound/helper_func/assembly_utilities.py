import os, sys, subprocess, math
from .sundry import _find_tool, _check_directory, _cleanup_path
from .ncbi import _sanitise_gene_metadata

# Define constants
SEQS_FNAME = "/coverage_mapped_genes.txt"
SEQS_DN_FNAME = "/coverage_denovo_genes.txt"
samtools = _find_tool('samtools')

def _bam_to_fastq(input_file: str, N_THREADS: int) -> str:
    """
        Convert BAM file containing unaligned reads into FastQ for
        downstream processing.
    """
    samtools_bam2fastq_args = list(["fastq", "--threads " + str(N_THREADS)])

    output_file = input_file.replace(".bam", ".fastq")
    samtools_bam2fastq_args.extend(["-o " + output_file, input_file])
    os.system(samtools + str(" ").join(samtools_bam2fastq_args))

    return output_file


def _assemble(fastq_file: str, output_dir: str, denovo_assembly: bool,
              N_THREADS: int, KMER_SIZE: int = 131) -> str:
    """
        Assemble reads using SPAdes. Output filename is fixed to scaffold.fasta
        and cannot be changed, enforcing a serial pipeline whereby one file
        can be processed at a time.

        NOTE: The default K-Mer value (131) was checked by hand as optimal for
              Escherichia coli.
    """
    # Prepare directory/file names and structure
    target_dir = fastq_file.split("/")[-2]  # [-1] ==
    output_dir = _check_directory(fastq_file,
                                  str("/").join([target_dir, output_dir]),
                                  target_dir, export_filename=False)
    # Step 1: run spades
    spades = _find_tool('spades.py')
    spades_args = list(["--isolate", "--threads " + str(N_THREADS),
                        "-k " + str(KMER_SIZE), "-m 1024"])
    if denovo_assembly is False:
        spades_args.extend(["--pe-12 1 " + fastq_file, "-o " + output_dir])
    else:
        if fastq_file.find("_U1_") > -1:
            fastq2_file = fastq_file.replace("_U1_trimmed", "_U2_trimmed")
        else:
            fastq2_file = fastq_file.replace("_1_trimmed", "_2_trimmed")
        spades_args.extend(["--pe-1 1 " + fastq_file, "--pe-2 1 " + fastq2_file,
                            "-o " + output_dir])
    os.system(spades + str(" ").join(spades_args))
    return output_dir, "scaffolds.fasta"


def assemble_denovo(fastq_file: str, N_THREADS: int = 1,
                    project_name: str = None) -> str:
    """
        Assemble de novo reads contained in reads_file. Useful to analyse
        reads that could not be mapped to a reference (merged into one file),
        or create a scaffold that can be used to detect species or calculate
        coverage depth.
    """
    scaffold_dir, scaffold_file = _assemble(fastq_file, "tmp", True,
                                            N_THREADS, KMER_SIZE=127)

    # Step 2: Use spades (-> dustmask -> prinseq not yet implemented)
    relative_path, output_file = _cleanup_path(fastq_file, project_name)
    if output_file.find("_U1_") > -1:
        # Is it a dataset of unpaired reads?
        output_file = output_file.replace("_U1_trimmed.fastq.gz",
                                          "_U_assembly.fa")
    else:
        output_file = output_file.replace("_1_trimmed.fastq.gz", "_assembly.fa")
    output_dir = _check_directory(fastq_file, "assemblies/de_novo",
                                  relative_path , export_filename=False)
    scaffold_exists = True
    # Move assembled scaffold to corresponding folder
    scaffold_path = str("/").join([scaffold_dir, scaffold_file])
    if os.path.exists(scaffold_path) is False:
        # Assemble failed due to _very_ low quality reads? Create empty file
        # or the rest of the program will fail.
        os.system("touch " + scaffold_path)
        scaffold_exists = False

    output_path = str("/").join([output_dir, output_file])
    os.system("mv " + str(" ").join([scaffold_path, output_path]))

    # SPAdes will NOT overwrite a folder containing assembly files: delete
    os.system("rm -Rf " + scaffold_dir)
    return output_path, scaffold_exists


def calculate_assembly_depth(assembly: str, N_THREADS: int = 1) -> str:
    """
        Given a SPAdes assembly, calculate read depth at every position of
        each scaffold. NOTE the assembly provided (.bam) was generated to
        allow this step, and can safely be deleted after read depth is
        calculated to free disk space.
    """
    samtools_depth_args = list(['depth', '--threads ' + str(N_THREADS)])

    # Calculate readp depth across the whole dataset
    depth_file = assembly.replace(".bam", "_depth.txt")
    samtools_depth_args.extend(['-o ' + depth_file, assembly])
    os.system(samtools + str(" ").join(samtools_depth_args))

    # Delete .bam file
    os.system("rm " + assembly)
    return depth_file


def _extract_contig_coords(depth_filename: str, contig_id: str) -> int:
    """
        Screen through depth data to find first and last lines where
        'contig_id' is found.
    """
    depth = open(depth_filename, 'r')

    INIT_FOUND = False
    END_FOUND = False

    line_n = 0
    for line in depth:
        if contig_id in line and INIT_FOUND is False:
            INIT_FOUND = True
            INIT = line_n
        elif contig_id not in line and INIT_FOUND is True:
            END = line_n-1  # Previous line was the last with contig_id
            END_FOUND = True
            break
        line_n += 1

    del depth

    if INIT_FOUND and END_FOUND:
        return INIT, END
    else:
        return False, False


def _extract_contig_rd_depth(depth_filename: str, contig_id: str) -> tuple:
    """
        Retrieve coverage depth for a specific contig
    """
    init_pos, end_pos = _extract_contig_coords(depth_filename, contig_id)

    # If gene is in complementary strand...
    if init_pos > end_pos:
        tmp = end_pos
        end_pos = init_pos
        init_pos = tmp

    depth_raw = subprocess.check_output(['head -n ' + str(end_pos+1) + ' ' +\
                                         depth_filename + ' | tail -n ' +\
                                         str(end_pos-init_pos+1)],
                                         shell=True).decode()

    return tuple(int(i.split('\t')[-1]) for i in depth_raw.split('\n') if len(i) > 0)


def _depth_little_helper(assembly: str, GENE_ID: str, metadata: dict,
                         SEQ_ID: str, COORDS: tuple, SEQ_CUTOFF: int) -> tuple:
    """
        Parse coverage file to retrieve coverage depth of a genomic region,
        modifying the coordinates to include start codon (ATG) and termination
        codon (TAA)
    """
    # File containing coverage depth per position (here, to load once per asmbl)
    depth_fName = assembly.replace('/de_novo', '').replace('.fa', '_depth.txt')

    # Retrieve INIT and STOP position to narrow region
    END_POS, START_POS, TARGET_SEQID = COORDS
    GENE_START_POS = metadata[GENE_ID][START_POS]
    GENE_END_POS = metadata[GENE_ID][END_POS]
    TARGET_SEQ = metadata[GENE_ID][TARGET_SEQID]

    contig_rd_depth = _extract_contig_rd_depth(depth_fName, TARGET_SEQ)

    if GENE_START_POS < GENE_END_POS:
        INIT_POS = max([0, GENE_START_POS-SEQ_CUTOFF-1])  # INIT_POS always >= 0
        STOP_POS = GENE_END_POS+3
    elif GENE_START_POS > GENE_END_POS:
        INIT_POS = max([0, GENE_END_POS-4])  # INIT_POS always >= 0
        STOP_POS = GENE_START_POS+SEQ_CUTOFF

    return contig_rd_depth[INIT_POS:STOP_POS]


def _calculate_median(read_depth: list) -> int:
    """
        Calculate median coverage (more robust than mean to outliers)
    """
    n = len(read_depth)
    read_depth = sorted(read_depth)

    if n % 2 == 0:
        median1 = read_depth[n//2]
        median2 = read_depth[n//2 - 1]
        median = (median1 + median2)/2
    else:
        median = read_depth[n//2]
    return int(median)


def retrieve_depth(seq_data: tuple, PREFIX: str, PRJ_PATH: str):
    """
        Retrieve how many times the reads cover `seq_data' in the new assembly.
    """
    # Include DB_TYPE in output files
    global SEQS_DN_FNAME
    if PREFIX not in SEQS_DN_FNAME:
        SEQS_DN_FNAME = SEQS_DN_FNAME.split('.')[0] + str('_') + PREFIX + \
                        str('.') + SEQS_DN_FNAME.split('.')[1]

    END_POS = 3
    START_POS = 2
    TARGET_SEQID = 5

    assembly, loci_metadata, SEQ_ID, SEQ_CUTOFF, hk_metadata = seq_data

    # Sanitise metadata list
    loci_metadata = _sanitise_gene_metadata(loci_metadata)
    hk_metadata = _sanitise_gene_metadata(hk_metadata)

    # Calculate coverage of assembly (one baseline per assembly, out of loop)
    if hk_metadata is not False:
        b_line = list()
        for GENE_ID in hk_metadata:
            rd_depth = _depth_little_helper(assembly, GENE_ID, hk_metadata, SEQ_ID,
                                            (END_POS, START_POS, TARGET_SEQID), 0)
            b_line.append(_calculate_median(rd_depth))

        mean_b_line = sum(b_line) / len(hk_metadata)  # mean depth
        var = sum(math.pow(i - mean_b_line, 2) for i in b_line) / len(b_line)
        std_b_line = math.sqrt(var)  # standard deviation of mean depth
        print("Baseline depth for this assembly is", str(mean_b_line), "±",
              round(std_b_line, 3))
    else:
        mean_b_line = str("NaN")
        std_b_line = str("NaN")
        print("Assembly '"+ assembly +"' has not enough data to calculate baseline depth.")

    if loci_metadata is not False:
        for GENE_ID in loci_metadata:
            rd_depth = _depth_little_helper(assembly, GENE_ID, loci_metadata, SEQ_ID,
                                            (END_POS, START_POS, TARGET_SEQID),
                                            SEQ_CUTOFF)
            # Write to file
            with open(PRJ_PATH + SEQS_DN_FNAME, 'a') as fOut:
                fOut.write(">" + SEQ_ID + str("\t") + GENE_ID + str("\t") +
                           loci_metadata[GENE_ID][TARGET_SEQID] + str("\t") +
                           str(sum(rd_depth) / len(rd_depth)) + str("\t") +
                           str(mean_b_line) + str("\n"))
    return PRJ_PATH + SEQS_DN_FNAME, tuple([mean_b_line, std_b_line])
