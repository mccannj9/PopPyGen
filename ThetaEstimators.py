#! /usr/bin/env python

from Parsers.Fasta import FastaReader


def compute_pairwise_differences(seq1, seq2):
    zipper = zip(seq1,seq2)

    diffs = 0
    for p,q in zipper:
        if p != q:
            diffs += 1

    return diffs


def compute_tajima_estimator(seqlist):
    import itertools

    seqs_pairwise = itertools.combinations(seqlist, 2)

    number_comparisons = len(seqlist)*(len(seqlist) - 1)/2
    diffs_pairwise = 0

    for seq1, seq2 in seqs_pairwise:
        diffs_pairwise += compute_pairwise_differences(seq1, seq2)

    # this is averaged over sequence length and number of comparisons
    # exp number mutations per site between two lineages
    # how DNASP calculates it (nucleotide diversity)
    return diffs_pairwise/(len(seqlist[0])*number_comparisons)*100


def compute_watterson_estimator(seqs):
    """ S = number of segregating sites
        n = the sample size """

    n = len(seqs)
    S = 0

    # seqs should obviously be all the same length
    sequence_length = len(seqs[0])

    for x in range(sequence_length):
        position = []
        for seq in seqs:
            position.append(seq[x])
        # len(set) should be 1 or 2 (Infinite Sites Model)
        unique_bases = set(position)

        if len(unique_bases) == 2:
            S += 1

        position = []

    print("S = %s" % S)

    # this is averaged over sequence length
    # how DNASP calculates it (Watterson's Estimator)
    return S/(sum([1/i for i in range(1,n)])*sequence_length)*100


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
