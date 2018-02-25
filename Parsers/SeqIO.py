#! /usr/bin/env python3


class Sequence:

    """
        Sequence object that can be used for generating sequence objects
        from fasta files, things can be added as needed for other applications
    """

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
        self.qual = self.qual[::-1]
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


class FastaReader(object):
    """
    Reader for FASTA files.
    """
    def __init__(self, file, wholefile=False, keep_linebreaks=False, sequence_class=Sequence):
        """
        file is a filename or a file-like object.
        If file is a filename, then .gz files are supported.
        If wholefile is True, then it is ok to read the entire file
        into memory. This is faster when there are many newlines in
        the file, but may obviously need a lot of memory.
        keep_linebreaks -- whether to keep the newline characters in the sequence
        """
        #if isinstance(file, basestring):
        #   file = xopen(file)
        file = open(file)
        self.fp = file
        self.sequence_class = sequence_class
        self.delivers_qualities = False

    def __iter__(self):
        """
        Read next entry from the file (single entry at a time).

        # TODO this can be quadratic since += is used for the sequence string
        """
        name = None
        seq = ''
        for line in self.fp:
            # strip() should also take care of DOS line breaks
            line = line.strip()
            if line and line[0] == '>':
                if name is not None:
                    assert seq.find('\n') == -1
                    yield self.sequence_class(name, seq, None, "fasta")
                name = line[1:]
                seq = ''
            else:
                seq += line
        if name is not None:
            assert seq.find('\n') == -1
            yield self.sequence_class(name, seq, None, "fasta")

    def __enter__(self):
        if self.fp is None:
            raise ValueError("I/O operation on closed FastaReader")
        return self

    def __exit__(self, *args):
        self.fp.close()


def trim_read(fq_obj, trim=0):
    if trim > 0:
        fq_obj.seq = fq_obj.seq[trim:]
    elif trim < 0:
        fq_obj.seq = fq_obj.seq[::-1][trim:][::-1]
    return fq_obj


class PhylipReader(object):
    """
    Reader for Phylip files.
    """
    def __init__(self, file, keep_linebreaks=False, sequence_class=Sequence):
        """
        file is a filename or a file-like object.
        If file is a filename, then .gz files are supported.
        If wholefile is True, then it is ok to read the entire file
        into memory. This is faster when there are many newlines in
        the file, but may obviously need a lot of memory.
        keep_linebreaks -- whether to keep the newline characters in the sequence
        """

        file = open(file)
        self.fp = file
        self.sequence_class = sequence_class
        header = self.fp.readline().strip().split()
        self.nind = int(header[0])
        self.nloci = int(header[1])

    def __iter__(self):
        """
        Read next entry from the file (single entry at a time).

        # TODO this can be quadratic since += is used for the sequence string
        """

        for line in self.fp:
            # strip() should also take care of DOS line breaks
            line = line.strip().split()
            name = line[0]
            seq = line[1]
            yield self.sequence_class(name, seq, None, "fasta")

    def __enter__(self):
        if self.fp is None:
            raise ValueError("I/O operation on closed FastaReader")
        return self

    def __exit__(self, *args):
        self.fp.close()

    def __str__(self):
        return "Phylip Filename: %s\nInd: %s -- Loci: %s" % (
            self.fp.name, self.nind, self.nloci
        )
