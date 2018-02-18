#! /usr/bin/env python

from Parsers.Fasta import FastaReader


def transpose_list_of_lists(matrix):
    return list(map(list, zip(*matrix)))


def extract_segregating_sites(seqs, outgroup=0):
    """ outgroup variable is the row number of the
        seq of the outgroup in the alignment   """

    seq_length = len(seqs[0])

    outgroup_seq = seqs[outgroup]
    # allow outgroup to be placed anywhere in alignment
    seqs = seqs[0:outgroup] + seqs[outgroup+1:]

    site_matrix = []
    for x in range(seq_length):
        position = []
        for seq in seqs:
            position.append(seq[x])
        # len(set) should be 1 or 2 (Infinite Sites Model)
        unique_bases = set(position)
        if len(unique_bases) == 2:
            if outgroup_seq[x] in unique_bases:
                position.insert(0, outgroup_seq[x])
                site_matrix.append(position)

    return site_matrix


def convert_to_binary_matrix(site_mat):

    binmat = []
    for site in site_mat:
        outgroup = site[0]
        position = site[1:]
        binrow = [0 if base == outgroup else 1 for base in position]
        binmat.append(binrow)

    return binmat


def plot_binary_matrix(binmat):

    # binmat = transpose_list_of_lists(binmat)
    for row in binmat:
        out = "".join([str(x) for x in row])
        print(out)


def expected_sfs_coal(sample_size):

    total = sum([1/k for k in range(1,sample_size)])
    sfs = []
    for j in range(1, sample_size):
        sfs.append((1/j)/total)

    return sfs


def observed_sfs_coal(binmat):

    sfs = [0 for x in range(len(binmat))]
    for row in binmat:
        idx = row.count(1)

    return 0


def count_singletons(binmat):
    singletons = 0
    for site in binmat:
        if site.count(1) == 1:
            singletons += 1

    return singletons


def coalest_proportion_of_singletons(sample_size):
    total = sum([1/k for k in range(1,sample_size)])

    return 1/total


def main():
    fasta_file = "/home/jamc/Data/GitHub/PopPyGen/Data/TNFSF5_All_aligned.fas"
    fasta = FastaReader(fasta_file)

    seqs = []
    for record in fasta:
        seqs.append(record.seq)

    site_matrix = extract_segregating_sites(seqs)
    binary_matrix = convert_to_binary_matrix(site_matrix)
    plot_binary_matrix(binary_matrix)

    obs_singletons = count_singletons(binary_matrix)/len(binary_matrix)
    exp_singletons = coalest_proportion_of_singletons(len(seqs))
    print("Obs Single: %s" % obs_singletons)
    print("Exp Single: %s" % exp_singletons)
    exp_sfs = expected_sfs_coal(len(seqs))
    print(exp_sfs)


if __name__ == '__main__':
    main()
