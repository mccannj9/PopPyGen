#! /usr/bin/env python3


class Sequence:

    """Sequence object that can be used for generating sequence objects
       from fasta files, things can be added as needed for other applications"""

    def __init__(self, name="", seq="", qual="", seqform=""):
        self.name = name
        self.seq = seq
        self.qual = qual
        self.seqform = seqform
        self.is_revcomp = False

    def check_duplicate(self, other):
        return self.seq == other.seq

    def reverse_complement(self):
        self.seq = reverseCompSequence(self.seq)
        self.qual = reverseSequence(self.qual)
        self.is_revcomp = not(self.is_revcomp)

    def fasta_oneline(self):
        return ">%s\n%s" % (self.name, self.seq)

    def fasta_wrap(self, n=80):
        seq_out = ">%s\n" % self.name
        for x in range(0,len(self.seq),n):
            seq_out += self.seq[x:x+n] + "\n"
        seq_out = seq_out.rstrip('\n')
        return seq_out

    def fastq_out(self):
        return "%s\n%s\n+\n%s" % (self.name,self.seq,self.qual)

    def __str__(self):
        return "%s\t%s" % (self.name, self.seq[0:50] + " ...")

    def __repr__(self):
        return "%s\t%s" % (self.name, self.seq[0:50] + " ...")

    def __len__(self):
        return len(self.seq)

    def __eq__(self, other):
        return self.name == other.name


rclookup = {
    'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
    'w': 'w', 's': 's', 'm': 'k', 'k': 'm', 'r': 'y', 'y': 'r',
    'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b', 'n': 'n',
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
    'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N', "-": "-"}


def reverseCompSequence(sequence):
    """
    Defines a generic method for reverse complementing a sequence of
    nucleotide codes.  This method fully supports all of the IUPAC
    ambiguity codes.
    """
    tmp = ""

    for base in sequence[::-1]:
        tmp += rclookup[base]

    return tmp

