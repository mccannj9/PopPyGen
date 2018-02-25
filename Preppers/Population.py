#! /usr/bin/env python3


genotypes = {
    "A": ["A", "A"],
    "T": ["T", "T"],
    "G": ["G", "G"],
    "C": ["C", "C"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
}

excluded_bases = set(["B","D","H","V","N", "U"])


def transpose_list_of_lists(matrix):
    return list(map(list, zip(*matrix)))


def ex_seg_sites(seqs):

    seq_length = len(seqs[0].seq)
    site_matrix = []

    for x in range(seq_length):
        position = []
        bases = []
        for seq in seqs:
            # if seq.seq[x] in genotypes:
            position.append(seq.seq[x])
            try:
                bases += genotypes[seq.seq[x]]
            except KeyError:
                bases.append(seq.seq[x])
        # len(set) should be 1 or 2 (Infinite Sites Model)
        unique_bases = set(bases)
        # check if set contains any unwanted bases, i.e. non-empty
        if unique_bases.intersection(excluded_bases):
            continue
        if len(unique_bases) == 2:
            # print("".join(position))
            site_matrix.append(position)
            # print(site_matrix)

    return transpose_list_of_lists(site_matrix)


def ex_seg_sites_with_constant(seqs):

    seq_length = len(seqs[0].seq)
    site_matrix = []

    for x in range(seq_length):
        position = []
        bases = []
        for seq in seqs:
            # if seq.seq[x] in genotypes:
            position.append(seq.seq[x])
            try:
                bases += genotypes[seq.seq[x]]
            except KeyError:
                bases.append(seq.seq[x])
        # len(set) should be 1 or 2 (Infinite Sites Model)
        unique_bases = set(bases)
        # check if set contains any unwanted bases, i.e. non-empty
        if unique_bases.intersection(excluded_bases):
            continue
        if len(unique_bases) <= 2:
            # print("".join(position))
            site_matrix.append(position)
            # print(site_matrix)

    return transpose_list_of_lists(site_matrix)


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


def extract_segregating_sites_no_outgroup(seqs):
    """ outgroup variable is the row number of the
        seq of the outgroup in the alignment   """

    seq_length = len(seqs[0])

    outgroup_seq = seqs[0]
    # no outgroup >> take first sequence as outgroup
    seqs = seqs[1:]

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
