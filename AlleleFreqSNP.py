#! /usr/bin/env python

from Preppers.Population import genotypes


def parse_popfile(filename):

    f = open(filename)
    pop_dict = {}

    for line in f:
        line = line.strip().split()
        pop_dict[line[0]] = int(line[1])

    return pop_dict


def get_most_frequent_alleles(records):

    seq_length = len(records[0].seq)
    alleles = []
    for x in range(seq_length):
        site = []

        for rec in records:
            gt = genotypes[rec.seq[x]]
            site += gt

        bases = list(set(site))

        if site.count(bases[0]) >= site.count(bases[1]):
            alleles.append(bases[0])
        else:
            alleles.append(bases[1])

    return alleles


def assign_groups(population_dict):
    pops = sorted(set(list(population_dict.values())))

    groups = []
    for val in pops:
        curr_group = []
        for key in population_dict:
            if population_dict[key] == val:
                curr_group.append(key)
        groups.append(curr_group)

    groups.append(list(population_dict.keys()))

    return groups


def main():
    from Parsers.SeqIO import PhylipReader
    phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/sback_snps.phy")
    pop_assign = parse_popfile("/home/jamc/Data/GitHub/PopPyGen/Data/sback_popfile.txt")
    out = open("/home/jamc/Data/GitHub/PopPyGen/Data/sback_allele_counts.txt", "w")
    pops = assign_groups(pop_assign)

    records = list(phy)
    seq_length = len(records[0].seq)

    mfa = get_most_frequent_alleles(records)

    for x in range(seq_length):
        major_allele = mfa[x]
        allele_count = []

        for population in pops:
            groups_records = [y for y in records if y.name in population]
            alleles = []

            for rec in groups_records:
                alleles += genotypes[rec.seq[x]]

            allele_count.append(alleles.count(major_allele))

        counts = "\t".join([str(k) for k in allele_count])
        print("%s\t%s" % (x, counts), file=out)

    out.close()

    return 0


if __name__ == '__main__':
    main()
