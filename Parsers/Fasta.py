#! /usr/bin/env python3

from SeqIO import Sequence


def trim_read(fq_obj, trim=0):
    if trim > 0:
        fq_obj.seq = fq_obj.seq[trim:]
    elif trim < 0:
        fq_obj.seq = fq_obj.seq[::-1][trim:][::-1]
    return fq_obj


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
