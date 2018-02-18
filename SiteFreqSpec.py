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


def convert_to_binary_matrix(data_matrix):

    binmat = []
    for site in data_matrix:
        outgroup = site[0]
        position = site[1:]
        binrow = ["0" if base == outgroup else "1" for base in position]
        binmat.append(binrow)

    return binmat


def main():
    fasta_file = "/home/jamc/Data/GitHub/PopPyGen/Data/TNFSF5_All_aligned.fas"
    fasta = FastaReader(fasta_file)

    seqs = []
    for record in fasta:
        seqs.append(record.seq)

    site_matrix = extract_segregating_sites(seqs)
    binary_matrix = convert_to_binary_matrix(site_matrix)
    print(site_matrix)
    print(binary_matrix)
    print("SM: nrows = %s :: ncols = %s" % (len(site_matrix), len(site_matrix[0])))
    print("BM: nrows = %s :: ncols = %s" % (len(binary_matrix), len(binary_matrix[0])))


if __name__ == '__main__':
    main()
