"""Parse the output of aggregate_snp_detection_results rule.

This output contains blocks of results.
"""
import yaml
import copy
import sys

def snp_to_text(snp):
    return f"{snp['orientation']}@{snp['pos']}:{snp['ref']}>{snp['alt']}"

def snps_to_set(snps):
    return set([snp_to_text(snp) for snp in snps])


def offset_snp(snp, offset, invert=False):
    orientation, tail = snp.split("@")
    pos, base_change = tail.split(":")
    base_from, base_to = base_change.split(">")
    pos = int(pos) + offset
    if invert:
        tmp = base_from
        base_from = base_to
        base_to = tmp
    
    return f"{orientation}@{pos}:{base_from}>{base_to}"


def evaluate_file(path, stats_file):
    nr_of_misidentified_snps = 0
    nr_of_missed_snps = 0
    nr_of_correct_snps = 0
    nr_loci_missed_by_stacks = 0
    nr_loci_split_by_stacks = 0
    correctly_classified = 0
    total_simulated_mutations = set()
    with open(path, "r") as yf:
        yaml_data = yaml.safe_load(yf)

    for name, l in yaml_data["Loci"].items():
        if not l["stacks_loci"]:
            nr_loci_missed_by_stacks += 1
        elif len(l["stacks_loci"]) > 1:
            nr_loci_split_by_stacks += 1

        all_stacks_mutations = set()
        for sl in l["stacks_loci"]:
            try:
                all_stacks_mutations |= snps_to_set(sl["SNPs"])
            except:
                all_stacks_mutations = set()

        all_simulated_mutations = set()
        for _, gta in l["ground_truth_alleles"].items():
            if gta:
                all_simulated_mutations |= snps_to_set(gta["mutations"])
                total_simulated_mutations |= set([f"name_{mut}" for mut in snps_to_set(gta["mutations"])])

        # after removing all detected SNPs from the simulated set
        # check again using simulated position +1 and remove
        # then simulated position +2 etc.
        # This should correctly identify problems through locus alignment
        
        remaining_simulated_mutations = all_simulated_mutations - all_stacks_mutations
        remaining_stacks_mutations = all_stacks_mutations - (all_simulated_mutations & all_stacks_mutations)
        
        corrected_remaining_simulated_mutations = copy.deepcopy(remaining_simulated_mutations)
        for sim_mut in remaining_simulated_mutations:
            offset_1 = offset_snp(sim_mut, 1)
            offset_2 = offset_snp(sim_mut, 2)
            offset_0_switched = offset_snp(sim_mut, 0, True)
            offset_1_switched = offset_snp(sim_mut, 1, True)
            offset_2_switched = offset_snp(sim_mut, 2, True)            
            if offset_1 in remaining_stacks_mutations:
                remaining_stacks_mutations.remove(offset_1)
                corrected_remaining_simulated_mutations.remove(sim_mut)
            if offset_2 in remaining_stacks_mutations:
                remaining_stacks_mutations.remove(offset_2)
                corrected_remaining_simulated_mutations.remove(sim_mut)

            if offset_0_switched in remaining_stacks_mutations:
                remaining_stacks_mutations.remove(offset_0_switched)
                corrected_remaining_simulated_mutations.remove(sim_mut)
            if offset_1_switched in remaining_stacks_mutations:
                remaining_stacks_mutations.remove(offset_1_switched)
                corrected_remaining_simulated_mutations.remove(sim_mut)
            if offset_2_switched in remaining_stacks_mutations:
                remaining_stacks_mutations.remove(offset_2_switched)
                corrected_remaining_simulated_mutations.remove(sim_mut)

        # relabel for easier use
        missed_snps = corrected_remaining_simulated_mutations
        bad_snps = remaining_stacks_mutations
        
        if missed_snps:
            nr_of_missed_snps += len(missed_snps)

        if bad_snps:
            nr_of_misidentified_snps += len(bad_snps)
            
        if not bad_snps and not missed_snps:
            nr_of_correct_snps += len(all_simulated_mutations)
        elif bad_snps:
            print(f"Simulated: {all_simulated_mutations}", file=sys.stderr)
            print(f"Stacks:    {all_stacks_mutations}", file=sys.stderr)
            print(f"NOT Corrected:\nundetected {missed_snps}\n{bad_snps}", file=sys.stderr)
            print(f"Position Corrected:\nundetected {corrected_remaining_simulated_mutations}\n{remaining_stacks_mutations}", file=sys.stderr)
            print("\n\n", file=sys.stderr)
        if all_simulated_mutations == all_stacks_mutations:
            correctly_classified += 1
        print(f"All simulated mutations: {all_simulated_mutations} (len(all_simulated_mutations)) = {(len(all_simulated_mutations))}")
    print(f"Total simulated mutations: {total_simulated_mutations} (len(total_simulated_mutations)) = {(len(total_simulated_mutations))}")
    output = [
        f"Analyzed {path}",
        f"Misidentified SNPs (FP):   {nr_of_misidentified_snps:>5}",
        f"Missed SNPs (FN):          {nr_of_missed_snps:>5}",
        f"Corectly called SNPs (TP): {nr_of_correct_snps:>5}",
        f"Total SNPs:                {len(total_simulated_mutations):>5}",
        f"Split loci:                {nr_loci_split_by_stacks:>5}",
        f"Missed loci:               {nr_loci_missed_by_stacks:>5}",
        f"Correctly classified loci: {correctly_classified:>5}",
        "\n",
    ]
    print("\n".join(output), file=stats_file)


if __name__ == "__main__":
    with open(snakemake.output.stats_file, "w") as stats_file:
        for f in sorted(snakemake.input.validation_files):
            evaluate_file(f, stats_file)
