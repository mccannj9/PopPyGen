#! /usr/bin/env python3


from SeqIO import Sequence


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
        #if isinstance(file, basestring):
        #   file = xopen(file)
        file = open(file)
        self.fp = file
        self.sequence_class = sequence_class
        header = self.fp.readline().strip().split()
        self.nindividuals = int(header[0])
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


if __name__ == '__main__':
    phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/sback.phy")
    for seq in phy:
        print(seq)
