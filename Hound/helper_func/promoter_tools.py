from .alignment_tools import multiple_seq_alignment, _mark_duplicates
from genoplot.plot_align import plot_alignment, plot_depth_stats
from .sundry import _find_tool
from Bio import Align, AlignIO, SeqIO, Seq
from ete3 import PhyloTree


def _convert_phylip(fasta_alignment: str, phylip_alignment: str) -> str:
    """
        Convert alignment file from FASTA to ***Relaxed** Phylip to use
        with PhymlCommandline. WARNING: phyml does NOT support standard
        phylip format.
    """
    SeqIO.convert(fasta_alignment, 'fasta', phylip_alignment, 'phylip-relaxed')
    return phylip_alignment


def _generate_phylogeny(alignment_file: str) -> str:
    """
        Generate phylogeny for the alignment given.
    """
    from Bio.Phylo.Applications import PhymlCommandline
    from Bio import Phylo

    # Convert .aln file (FASTA) to .phy (relaxed phylip) to use phyml
    alignment_file = _convert_phylip(alignment_file,
                                     alignment_file.replace(".aln", ".phy"))
    # phyml does not permit custom outputfiles, use defaults
    tree_file = alignment_file.replace(".phy", ".phy_phyml_tree.txt")
    phyml_cmd = _find_tool('phyml')[:-1]  # remove final space, not needed
    phylogeny = PhymlCommandline(cmd='"'+ phyml_cmd +'"', input=alignment_file,
                                 model='GTR', alpha='e', bootstrap=-4,
                                 prop_invar='e', optimize='tlr', search='NNI',
                                 pars=True, frequencies='e', r_seed=100100,
                                 quiet=True)
    # Generate phylogeny
    print("Generating phylogeny...")
    phylogeny()
    tree = PhyloTree(tree_file)
    putative_outgroup = tree.get_farthest_leaf()[0].name
    tree.set_outgroup(putative_outgroup)
    tree.ladderize(direction=1)
    return tree_file, tree


def _sort_alignment(alignment_file: str, tree: str) -> Seq:
    """
        Re-arrange alignment based on phylogenetic tree.
    """
    # Extracts all clades following the phylogeny
    clades_list = [i.name for i in tree.traverse("postorder") if len(i.name)>0]
    # Sort alignment
    alignment = AlignIO.read(alignment_file, 'fasta')
    alignment_sorted = Align.MultipleSeqAlignment([])
    for entry in clades_list[::-1]:
        for seq in alignment:
            if entry == seq.id.replace("[", "").replace("]", ""):
                seq.id = seq.id.split("__")[0]
                seq.name = seq.name.split("__")[0]
                alignment_sorted.append(seq)
    # Overwrite old alignment
    AlignIO.write(alignment_sorted, alignment_file, 'fasta')
    return


def _estimate_consensus(aln_file: str, CONS_THRESHOLD: float = 0.8) -> Seq:
    """
        Estimate a consensus sequence given an alignment.
    """
    from Bio.Align import AlignInfo
    alignment = AlignIO.read(aln_file, 'fasta')
    aln_summary = Align.AlignInfo.SummaryInfo(alignment)
    consensus = aln_summary.gap_consensus(threshold=CONS_THRESHOLD,
                                          ambiguous='X', require_multiple=True)
    # mostly for plotting purposes, estimate conserved sequence
    conserved = aln_summary.gap_consensus(threshold=0.95, ambiguous='X',
                                          require_multiple=True)
    return consensus, conserved


def _extract_reference_locations(consensus_seq: Seq,
                                 seqs_of_interest: list) -> list:
    """
        Given a series of sequences of interest, find their start and end
        position in the consensus sequence to aid plotting.
    """
    reference_locations = list()
    for seq in seqs_of_interest:
        START_POS = consensus_seq.find(seq)
        END_POS = consensus_seq.find(seq) + len(seq) if START_POS != -1 else -1
        reference_locations.append([START_POS, END_POS])
    return reference_locations


def _normalise_coverage(loci_coverage: float, baseline_coverage: float) -> float:
    """
        Use coverage depth baseline to normalise coverage depth of detected
        loci.
    """
    return loci_coverage / baseline_coverage


def _sort_coverage_data(loci_coverage_dict: dict, hk_coverage_dict: dict,
                        phylogeny: Align) -> list:
    """
        Given an imported coverage dataset, re-arrange it to match
        the order of the corresponding phylogeny for plotting.

        WARNING!!! SPAdes PROVIDES K-MER COVERAGE NOT COVERAGE DEPTH!!!!!!
    """
    loci_cov_sorted = list()
    for sequence in phylogeny:
        seq_id = sequence.name.split('__')[0]
        if seq_id in loci_coverage_dict.keys():
            norm_cov = _normalise_coverage(float(loci_coverage_dict[seq_id]),
                                           float(hk_coverage_dict[seq_id]))
            loci_cov_sorted.append(norm_cov)
    return loci_cov_sorted


def _extract_coverage(coverage_file: str, phylogeny: Align) -> tuple:
    """
        Parse 'coverage_file' to import coverage data *in the same order*
        as the phylogeny, and normalised.
    """
    cov_raw = open(coverage_file, 'r').read().split("\n")
    cov_loci = dict()
    cov_hk = dict()

    for entry in cov_raw:
        if len(entry) > 0:
            data = entry.split('\t')
            seq_id = data[0].removeprefix('>').removesuffix('_assembly.fa')
            loci_cov = data[-2]
            hk_cov = data[-1]

            # Isolate locus number, if more than 1
            if len(data[1].split('_locus')) == 1:
                locus_n = 0
            else:
               locus_n = data[1].split('_locus')[1]

            if locus_n == 0:
                cov_loci[seq_id] = loci_cov
                cov_hk[seq_id] = hk_cov
            else:
                cov_loci[seq_id + str('_L') + locus_n] = loci_cov
                cov_hk[seq_id + str('_L') + locus_n] = hk_cov
    ## More elegant... but excludes genes with multiple copies
    # cov_loci = dict([itemgetter(*[0, -2])(i.split("\t")) for i in cov_raw if len(i) > 0])
    # cov_hk = dict([itemgetter(*[0, -1])(i.split("\t")) for i in cov_raw if len(i) > 0])
    return _sort_coverage_data(cov_loci, cov_hk, phylogeny)


def _find_mutations(sample_sequence: Seq, reference_sequence: Seq,
                    location: list) -> list:
    """
        Return the exact location of a mutation, using the coordinates from
        the general alignment, for downstream plotting.
    """
    from difflib import ndiff

    sequence_difference = ndiff(reference_sequence, sample_sequence)

    init_location = location[0]
    coordinates = list()

    correction_coeff = 0  # Account for presence of added nt (+C)
    for pos, change in enumerate(sequence_difference):
        if change[0] == '-':
            coordinates.append(init_location + pos - correction_coeff)
        elif change[0] == '+':
            correction_coeff +=1
    return coordinates


def _detect_mutations(alignment_file: str, ref_seqs_locations: list,
                      ref_seqs: list) -> dict:
    """
        Retrieve promoter sequences and compare them to a given reference.
        If there are differences, highlight them for plotting.
    """
    # Load alignment sequences
    alignment = SeqIO.parse(alignment_file, 'fasta')

    mutations = dict()
    for seq in alignment:
        mutations[seq.id] = list()  # Create an entry for each assembly
        for id, location in enumerate(ref_seqs_locations):
            mut_coordinates = _find_mutations(seq[location[0]:location[1]],
                                              ref_seqs[id], location)
            if len(mut_coordinates) > 0:
                # Add to dictionary only if mutations are found
                mutations[seq.id].append(mut_coordinates)
    return mutations


def analyse_seqs_found(seqs_file: str, FLT_THR: float = 0.5,
                       N_THREADS: int = 1) -> str:
    """
        Align the sequences found, calculate consensus/conserved sequences,
        and generate the corresponding phylogeny.
    """
    import os
    if os.path.exists(seqs_file) is True:
        # Align sequences
        alignment_file = multiple_seq_alignment(seqs_file, FLT_THR, N_THREADS)
        # Mark duplicated sequence names
        _mark_duplicates(alignment_file)
        # Generate phylogeny to re-arrange alignment sequences
        # This step is slow (phyml!), do it only if the phylogeny doesn't exist
        tree_file, phylogeny = _generate_phylogeny(alignment_file)
        _sort_alignment(alignment_file, phylogeny)
        # Consensus is 80% ID, conserved is 95%. Conserved is for plotting *
        consensus_seq, conserved_seq = _estimate_consensus(alignment_file,
        CONS_THRESHOLD=0.7)
        return phylogeny, alignment_file, consensus_seq, conserved_seq


def plot_analysis(alignment_file: str, cov_file: str, cov_stats: dict,
                  PREFIX: str, CUTOFF: int = 0, ROI: bool = False,
                  PROMOTER: bool = True, LABELS: str = None) -> str:
    """
    """
    # Extract metadata for TEM promoter boxes
    # TODO: Avoid hard-coding sequences, take a user-given parameter
    if ROI is not False:
        # Target sequences
        # TODO: Write FASTA parser with outputs ref_seqs, and ref_seqs_label
        Pb_10_SEQ = str("TAATGT")
        Pa_35_SEQ = str("TTGAAG")
        Pb_10_35_SEQ = str("GTGATACGCC")
        P3_35 = str("TTCAAA")
        P3_4_10 = str("GACAAT")
        ref_seqs = tuple([P3_35, P3_4_10, Pb_10_SEQ, Pa_35_SEQ,
                          Pb_10_35_SEQ])
        ref_seqs_label = tuple(["P3 -35", "P3/4-10", "Pb -10",
                                "Pa -35", "Pb -35/\nPa -10"])
        ref_seqs_locations = _extract_reference_locations(consensus_seq,
                                                          ref_seqs)
        mutations_found = _detect_mutations(alignment_file,
                                            ref_seqs_locations, ref_seqs)
    else:
        mutations_found = None
        ref_seqs_locations = None
        ref_seqs_label = None

    alignment = AlignIO.read(alignment_file, 'fasta')
    cov_data = _extract_coverage(cov_file, phylogeny)
    # Plot alignment
    plot_alignment(alignment, PREFIX, cov_data, phylogeny,
                   consensus_seq, conserved_seq, mutations_found,
                   cov_file, CUTOFF, (ref_seqs_label,
                   ref_seqs_locations), PROMOTER_ONLY=PROMOTER,
                   show_plot=False)

    # Plot [baseline] coverage depth stats for completeness
    plot_depth_stats(cov_stats)
