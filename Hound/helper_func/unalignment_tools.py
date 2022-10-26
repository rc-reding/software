import os
from .sundry import _find_tool, _check_directory


def process_unaligned(output_template: str, N_THREADS: int = 1) -> str:
    """
        Process unaligned reads (and mates) that could not be mapped to a
        reference genome, and proceed to assemble de novo.
    """
    # Prepare variables and directories (CHECK: https://www.biostars.org/p/56246/)
    # Tentative pipeline from https://github.com/jalwillcox/Analyze-Unmapped-Reads
    reads_file = output_template.replace("/unaligned_reads", "")
    output_file = _check_directory(output_template, [], [])  # Changes extension to .txt as byproduct...

    outUNALN = output_file.replace("Aligned.calmd.bam", "_UNALIGNED.calmd.bam")
    outTMP = output_file.replace("Aligned.calmd.bam", "_UNALIGNED_TMP.calmd.bam")
    outMATES = output_file.replace("Aligned.calmd.bam", "_UNALIGNED_MATES.calmd.bam")
    outMERGED = output_file.replace("Aligned.calmd.bam", "_MERGED.calmd.bam")
    outSORTED = output_file.replace("Aligned.calmd.bam", "_UNALIGNED_FINAL.calmd.bam")

    samtools = _find_tool('samtools')
    samtools_merge_args = list(["merge", "--verbosity 1",
                                "--threads " + str(N_THREADS)])
    samtools_sort_args = list(["sort", "--verbosity 1",
                               "--threads " + str(N_THREADS)])
    # samtools_unaligned_args = list(["view", "--verbosity 1", "-bh",
    #                                 "-f 4 -F 2048",
    #                                 "--threads " + str(N_THREADS)])
    # samtools_mates_args = list(["view", "--verbosity 1", "-bh", "-f 8 -F 2048",
    #                             "--threads " + str(N_THREADS)])
    # samtools_rmdup_args = list(["view", "--verbosity 1", "-bh", "-F 8",
    #                             "--threads " + str(N_THREADS)])
    samtools_unaligned_args = list(["view", "--verbosity 1", "-bh",
                                    "-f 4",
                                    "--threads " + str(N_THREADS)])
    samtools_mates_args = list(["view", "--verbosity 1", "-bh", "-f 8",
                                "--threads " + str(N_THREADS)])
    samtools_rmdup_args = list(["view", "--verbosity 1", "-bh", "-f 12",  # 12 = 8 + 4 = read1 AND read2 unmapped
                                "--threads " + str(N_THREADS)])
    # Extract unaligned reads for further processing
    samtools_unaligned_args.extend(["-o " + outUNALN, reads_file])
    os.system(samtools + str(" ").join(samtools_unaligned_args))

    # Extract unmapped mates of mapped reads for further processing
    samtools_mates_args.extend(["-o " + outTMP, reads_file])
    os.system(samtools + str(" ").join(samtools_mates_args))

    # Remove duplicates (unaligned and mates are both unaligned seqs)
    samtools_rmdup_args.extend(["-o " + outMATES, outTMP])
    os.system(samtools + str(" ").join(samtools_rmdup_args))

    # Merge unaligned and mates reads
    samtools_merge_args.extend(["-o " + outMERGED, outUNALN, outMATES,])
    os.system(samtools + str(" ").join(samtools_merge_args))

    # Sort unaligned and mates reads by NAME
    samtools_sort_args.extend(["-o " + outSORTED, "-n", outMERGED])
    os.system(samtools + str(" ").join(samtools_sort_args))

    # Remove temporary files
    temp_files = list([outTMP, outMATES, outMERGED, outUNALN])
    # [os.remove(file) for file in temp_files]
    # os.system("rm " +str(" ").join(temp_files))
    return outSORTED
