#! /usr/bin/env python

from Parsers.Fasta import FastaReader


def transpose_list_of_lists(matrix):
    return list(map(list, zip(*matrix)))


def extract_segregating_sites(seqs, outgroup=0):
    """ outgroup variable is the row number of the 
        seq of the outgroup in the alignment   """

    n = len(seqs)
    S = 0
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
                position.insert(0, outgroup[x])
                site_matrix.append(position)

    return site_matrix



def main():
    fasta_file = "/mnt/Data/Files/OpenCourseWare/PopGenData/TNFSF5_aligned.fas"
    fasta = FastaReader(fasta_file)

    seqs = []
    for record in fasta:
        seqs.append(record.seq)

    theta_d = compute_tajima_estimator(seqs)
    print("Td = %s" % theta_d)
    theta_w = compute_watterson_estimator(seqs)
    print("Ts = %s" % theta_w)


if __name__ == '__main__':
    main()
