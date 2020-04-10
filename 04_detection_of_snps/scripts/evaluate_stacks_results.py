# -*- coding: utf-8 -*-
"""Evaluate a stack mapping versus a simulated rage dataset.
"""
import argparse

import sys

from Bio.pairwise2 import align
from Bio.pairwise2 import format_alignment

from collections import Counter
from functools import partial

from yaml import dump
from yaml import CDumper as Dumper

import sketching as sk
import file_parser


def join_seqs(seq_p5, seq_p7, join_seq):
    """Join two sequences like the ones written by stacks."""
    return "".join((seq_p5, join_seq, seq_p7))


def find_matching_loci(gt_data, stacks_data, similarity, verbose=True):
    """Compare clusterings with minhash sketching.
            similar sequences. Default: 0.3

    Returns:
        dict: Mapping a tuple of (RAGE locus name, RAGE locus reference
        sequence) to the RAGE locus record and a list of associated Stacks
        locus records.
    """
    # sketching parameters
    k = 7
    s = 50
    sketch = partial(sk.bottom_sketch, k=k, s=s)

    # initialize assembly
    assembly = {
        (gt_record.name, gt_record.seq_p5):
        (gt_record, []) for gt_record in gt_data}
    # print("Sketching stacks data")
    # pre-sketch the data to avoid quadratic (re-) sketching
    sketched_stacks_data = [(sketch(stacks_record.seq.decode()), stacks_record)
                            for stacks_record in stacks_data]
    
    print(assembly)
    input("waiting for you")
    # for s, record in sketched_stacks_data:
    #     print(record, s)

    # print("Searching")
    for index, (gt_record) in enumerate(gt_data):
        # print user output
        if verbose and index % 100 == 0:
            print(index, file=sys.stderr)
        # print(joined_gt_seq)
        sketch_gt = sketch(gt_record.seq_p5)

        # compare all stacks locus sketches against the active RAGE loc. sketch
        # Append the loci of all similar sequences to the assembly list.

        for sketch_stacks, stacks_record in sketched_stacks_data:
            if sk.compare_sketches(sketch_gt, sketch_stacks, s) > similarity:
                assembly[(gt_record.name, gt_record.seq_p5)][1].append(stacks_record)
                record.found = True

    for _, record in sketched_stacks_data:
        if record.found is False:
            print("Not found:", record.seq)
    print(assembly)
    input("waiting for you")
    return assembly


def evaluate_assembly(assembly, gt_data, stacks_data, gt_stats, args):
    """Analyze an assembly of RAGE vs Stacks data.

    Arguments:
        assembly (dict): Mapping of RAGE sequences to RAGE records
            and associated Stacks records.
        gt_data (list of GTRecords): From a RAGE _gt.yaml file.
        stacks_data (list of TSVRecords): From a Stacks _export.tsv file.
        args (argparse.NameSpace): User parameters.
    """
    nr_undiscovered_gt_loci = 0
    nr_loci_with_undiscovered_mutations = 0
    nr_loci_with_discovered_mutations = 0
    nr_evaluated_loci = 0
    nr_successfully_aligned_loci = 0
    # write mapping to file

    with open(args.output, "w") as outfile:
        outdata = {
            "Loci": {},
        }
        for (gt_name, gt_seq), (gt_locus, stacks_loci) in assembly.items():
            outdata["Loci"][gt_name] = {
                "ground_truth_seq": gt_seq,
                "ground_truth_alleles": gt_locus.alleles,
                "stacks_loci": []
            }

            # if len(stacks_loci) > 1:
            #     # This means that more than one stacks locus associated
            #     # with the active ground truth locus
            #     print("HIT")
            #     print(gt_seq)
            #     print(stacks_loci)

            successfully_detected, successfully_aligned = False, False
            undetected, no_mutations = False, False

            # compute semiglobal alignments of the loci to verify that
            # they actually match
            for stacks_locus in stacks_loci:

                stacks_locus_info = {
                    "seq": stacks_locus.seq.decode(),
                }
                all_alns = align.globalms(
                    gt_seq,
                    stacks_locus.seq.decode(),
                    1,   # match score
                    0,   # mismatch panalty
                    -5,  # gap open penalty
                    -3,  # gap extend penalty
                    penalize_end_gaps=(False, False),
                    one_alignment_only=True,
                )
                # pick the first reported alignments
                # these are either unique or good enough
                aln = all_alns[0]

                if aln[2] >= 80:
                    successfully_aligned = True
                    stacks_locus_info["SNPs"] = [
                        {
                            # TODO: this is not 100% accurate. The precise
                            # position can variate with different spacer length
                            # etc. This is a conservative estimate.
                            "orientation": "p7" if entry.pos > 98 else "p5",
                            "pos": entry.pos,
                            "ref": entry.ref,
                            "alt": "".join(entry.alts),
                        } for entry in stacks_locus.data
                    ]
                    if stacks_locus.data:
                        successfully_detected = True
                        # NOTE: Currently this only checks, if the stacks locus
                        # was successfully aligned and also detected a mutation
                        # the kind of the mutation is not checked. Neither, if
                        # the correct mutation (e.g. SNP C>G at position 42 in
                        # the forward read) was detected.
                        #
                        # TODO: Evaluate if right SNPs were found
                        #       (mind that Stacks might not call the allele
                        #       RAGE simulated as root as main allele
                        #         -> consider x>y == y>x)
                        #
                        # WARNING how to manage split up loci?
                        # consider: stacks split up a heterozygous locus
                        # to two loci. How do we count that? Two misses?
                        # Check if the union of all variants covers the
                        # set of simulated variants?
                        #
                        # TODO: Evaluate if the right allele frequencies were
                        #       detected by stacks
                    elif gt_locus.mutations:
                        undetected = True
                    else:
                        no_mutations = True

                else:
                    print(f"MISMATCH with {stacks_locus.data[0].chrom}")
                    print(format_alignment(*aln))
                outdata["Loci"][gt_name]["stacks_loci"].append(
                    stacks_locus_info)

            if not stacks_loci:
                # print("No matching stack locus found", file=outfile)
                outdata["Loci"][gt_name]["stacks_loci"] = "No matches found"
                nr_undiscovered_gt_loci += 1

            if successfully_detected:
                nr_loci_with_discovered_mutations += 1
            elif undetected:
                nr_loci_with_undiscovered_mutations += 1
                print(f"The variants\n{gt_locus.alleles}\nin RAGE locus:\n{gt_name:<10} {gt_seq} were not detected in stacks loci\n{stacks_loci}", file=sys.stderr)
            elif no_mutations:
                ...

            nr_evaluated_loci += 1
            if successfully_aligned:
                nr_successfully_aligned_loci += 1

        # find stacks loci that are not in the RAGE loci (singletons and HRLs?)
        # create one entry for each locus assembled by stacks.
        # then find out in the assembly which ones were never assigned to a
        # rage locus:
        #
        # 1. create a dictionary with zero for all stacks locus names
        # 2. iterate through the gt_locus -> stacks-locus mapping
        # 3. increment the counter of each stacks locus, every time it occurs
        #    in the mapping
        stacks_locus_occurence_count = Counter(
            {rec.name: 0 for rec in stacks_data})
        for (gt_name, gt_seq), (gt_locus, stacks_loci) in assembly.items():
            for locus in stacks_loci:
                # count every gt locus that has an assigned stacks locus
                stacks_locus_occurence_count[locus.name] += 1

        # Report stacks loci that have no counterpart in the ground truth
        # (RAGE) loci
        stacks_only_loci = [
            (l, c) for l, c in stacks_locus_occurence_count.items() if c == 0]
        if stacks_only_loci:
            print(f"The following {len(stacks_only_loci)} loci were not "
                  "simulated by rage, but identified by stacks. These might "
                  "include incompletely digested reads, Null Alleles, "
                  "Singletons, and HRLs/ Lumberjack stacks.", file=sys.stderr)
            print([name for name, _ in stacks_only_loci], file=sys.stderr)

        # Check to how many RAGE loci each Stacks loci was assigned to.
        # This should always be 1.
        # A value >1 would suggest that a stacks locus is similar to two or
        # more RAGE loci, which is highly unlikely
        most_assigned_gt_loci = stacks_locus_occurence_count.most_common(5)
        stacks_locus_name, assigned_gt_loci = most_assigned_gt_loci[0]
        if assigned_gt_loci > 1:
            print(f"The 5 most assigned stacks locus id. This should always "
                  f"be 1:\n  {most_assigned_gt_loci}\n", file=sys.stderr)
        else:
            print(f"Each stacks locus was assigned to at most one ground "
                  f"truth (RAGE) loci", file=sys.stderr)

        outdata["metadata"] = {
            "Loci with mutations that were successfully discovered": nr_loci_with_discovered_mutations,
            "Total simulated mutation loci": gt_stats.nr_loci_with_muts,
            "Loci with mutations that were not discovered by stacks": nr_loci_with_undiscovered_mutations,
            "SNP discovery ratio": nr_loci_with_discovered_mutations / gt_stats.nr_loci_with_muts,
            }

        outfile.write(dump(outdata, default_flow_style=False, Dumper=Dumper,
                           explicit_start=True, sort_keys=False))


def main(args):
    """Compare ground truth with stacks assembly."""
    print(f"Loading gt data", file=sys.stderr)
    gt_data, gt_stats = file_parser.parse_rage_gt_file(args)

    print(f"Loading stacks data", file=sys.stderr)
    stacks_data = file_parser.get_stacks_data(args)

    print("Analyzing:", file=sys.stderr)
    assembly = find_matching_loci(gt_data, stacks_data,
                                  similarity=args.similarity_threshold,
                                  )

    print("\n\nLocus Analysis:\n", file=sys.stderr)
    evaluate_assembly(assembly, gt_data, stacks_data, gt_stats, args)

    # print("\n\nSNPs Analysis:\n", file=sys.stderr)
    # evaluate_snps(assembly, gt_data, stacks_data, args)


def get_argparser():
    """Manage user parameters"""
    parser = argparse.ArgumentParser()
    # input
    parser.add_argument(
        "-s", "--stacks-snps-file",
        help="Path to a stacks snps vcf file",
        dest="stacks_haplo",
        )
    parser.add_argument(
        "-f", "--stacks-fasta-file",
        help="Path to a stacks catalog fasta file",
        dest="stacks_fa",
        )
    parser.add_argument(
        "-g", "-y", "--ground-truth", "--yaml",
        help="Path the a YAML gt file",
        default="RAGEdataset_ATCACG_gt.yaml",
        dest="yaml",
        )
    # input
    parser.add_argument(
        "-o", "--output",
        help="If an output file should be written",
        default=False,
        )
    parser.add_argument(
        "-v", "--verbose",
        help="Print supplementary information",
        dest="verbose",
        action="store_true",
        default=False,
        )
    # analysis parameters
    parser.add_argument(
        "-t", "--similarity-threshold",
        help="The minimal estimated Jaccard similarity for two loci to be"
             "considered similar.",
        dest="similarity_threshold",
        type=float,
        default=0.2,
    )
    return parser


if __name__ == '__main__':
    parser = get_argparser()
    args = parser.parse_args()
    main(args)
