#! /usr/bin/env python

from Parsers.SeqIO import *


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


def compute_pairwise_differences_radseq(seq1, seq2):
    zipper = zip(seq1,seq2)
    bases = set(["A", "T", "G", "C", "R", "Y", "S", "W", "K", "M"])

    diffs = 0
    for p,q in zipper:
        position = set([p,q])
        check = bases.intersection(position)
        # print(position, check)
        if check == position:
            if p != q:
                diffs += 1
        # print(diffs)

    return diffs


def compute_watterson_estimator_radseq(seqdata):
    sequence_length = seqdata.nloci
    n = seqdata.nind
    S = 0

    seqdata = list(seqdata)
    bases = set(["A", "T", "G", "C", "R", "Y", "S", "W", "K", "M"])

    for x in range(sequence_length):
        position = []
        for seq in seqdata:
            position.append(seq.seq[x])
        # len(set) should be 1 or 2 (Infinite Sites Model)
        unique_bases = set(position)
        if unique_bases.intersection(bases) == unique_bases:
            if len(unique_bases) == 2:
                S += 1

    return S/(sum([1/i for i in range(1,n)])*sequence_length)*100


def compute_tajima_estimator_radseq(seqdata):
    import itertools

    sequence_length = seqdata.nloci
    # print(seqdata)
    seqlist = list(seqdata)
    seqs_pairwise = itertools.combinations(seqlist, 2)

    number_comparisons = len(seqlist)*(len(seqlist) - 1)/2
    diffs_pairwise = 0

    print(sequence_length, number_comparisons)

    for seq1, seq2 in seqs_pairwise:
        diffs_pairwise += compute_pairwise_differences_radseq(seq1.seq, seq2.seq)
        # print(diffs_pairwise)

    # this is averaged over sequence length and number of comparisons
    # exp number mutations per site between two lineages
    # how DNASP calculates it (nucleotide diversity)
    return diffs_pairwise/(sequence_length*number_comparisons)*100


def main():
    # fasta_file = "/mnt/Data/Files/OpenCourseWare/PopGenData/TNFSF5_aligned.fas"
    # fasta = FastaReader(fasta_file)

    # seqs = []
    # for record in fasta:
    #     seqs.append(record.seq)

    # theta_d = compute_tajima_estimator(seqs)
    # print("Td = %s" % theta_d)
    # theta_w = compute_watterson_estimator(seqs)
    # print("Ts = %s" % theta_w)

    # phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/sback.phy")

    # # theta_w = compute_watterson_estimator_radseq(phy)
    # # print("Tw = %s" % theta_w)
    # # print("Hew = %s" % (theta_w/(theta_w + 1)))
    # theta_d = compute_tajima_estimator_radseq(phy)
    # print("Td = %s" % theta_d)
    # print("Hed = %s" % (theta_d/(theta_d + 1)))

    print(" ")
    phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/sback_snps_with_constants_pop1.phy")
    print(phy)
    theta_w = compute_watterson_estimator_radseq(phy)
    print("Tw = %s" % theta_w)
    print("Hew = %s" % (theta_w/(theta_w + 1)))

    print(" ")
    phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/sback_snps_with_constants_pop1.phy")
    theta_d = compute_tajima_estimator_radseq(phy)
    print("Td = %s" % theta_d)
    print("Hed = %s" % (theta_d/(theta_d + 1)))

    print(" ")
    phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/sback_snps_with_constants_pop2.phy")
    print(phy)
    theta_w = compute_watterson_estimator_radseq(phy)
    print("Tw = %s" % theta_w)
    print("Hew = %s" % (theta_w/(theta_w + 1)))
    print(" ")

    phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/sback_snps_with_constants_pop2.phy")
    theta_d = compute_tajima_estimator_radseq(phy)
    print("Td = %s" % theta_d)
    print("Hed = %s" % (theta_d/(theta_d + 1)))
    print(" ")

    phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/sback_snps_with_constants_pop3.phy")
    print(phy)
    theta_w = compute_watterson_estimator_radseq(phy)
    print("Tw = %s" % theta_w)
    print("Hew = %s" % (theta_w/(theta_w + 1)))
    print(" ")

    phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/sback_snps_with_constants_pop3.phy")
    theta_d = compute_tajima_estimator_radseq(phy)
    print("Td = %s" % theta_d)
    print("Hed = %s" % (theta_d/(theta_d + 1)))

    phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/sback_snps_with_constants.phy")
    print(phy)
    theta_w = compute_watterson_estimator_radseq(phy)
    print("Tw = %s" % theta_w)
    print("Hew = %s" % (theta_w/(theta_w + 1)))
    print(" ")

    phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/sback_snps_with_constants.phy")
    theta_d = compute_tajima_estimator_radseq(phy)
    print("Td = %s" % theta_d)
    print("Hed = %s" % (theta_d/(theta_d + 1)))


if __name__ == '__main__':
    main()
