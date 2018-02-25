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

    sfs = [0 for x in range(len(binmat[0]))]
    for row in binmat:
        idx = row.count(1) - 1
        sfs[idx] += 1

    return [x/sum(sfs) for x in sfs]


def plot_sfs_bokeh(obs_sfs, exp_sfs, graph_output_file, gene_name):

    from bokeh.core.properties import value
    from bokeh.io import show, output_file
    from bokeh.models import ColumnDataSource
    from bokeh.plotting import figure
    from bokeh.transform import dodge

    output_file(graph_output_file)

    ncats = len(obs_sfs)
    cats = ["%s" % (x,) for x in range(1,ncats)]
    dtypes = ["Obs", "Exp"]

    data = {
        "cats": cats,
        "Obs": obs_sfs,
        "Exp": exp_sfs
    }

    x = [(cat, dtype) for cat in cats for dtype in dtypes]
    # counts = sum(zip(data['Obs'], data['Exp']), ())

    source = ColumnDataSource(data=data)

    p = figure(x_range=cats, y_range=(0, 0.6), plot_height=500, plot_width=1000, title="%s SFS" % gene_name,
               toolbar_location=None, tools="")

    p.vbar(x=dodge('cats', -0.20, range=p.x_range), top='Obs', width=0.3, source=source,
           color="#e84d60", legend=value("Obs"))

    p.vbar(x=dodge('cats', 0.20, range=p.x_range), top='Exp', width=0.3, source=source,
           color="#718dbf", legend=value("Exp"))

    p.x_range.range_padding = 0.1
    p.xgrid.grid_line_color = None
    p.legend.location = "top_left"
    p.legend.orientation = "horizontal"

    show(p)


def main():
    # data from Hammer et al 2004 (Human X-linked genes)
    fasta_file = "/home/jamc/Data/GitHub/PopPyGen/Data/TNFSF5_All_aligned.fas"
    fasta = FastaReader(fasta_file)

    seqs = []
    for record in fasta:
        seqs.append(record.seq)

    site_matrix = extract_segregating_sites(seqs)
    binary_matrix = convert_to_binary_matrix(site_matrix)
    plot_binary_matrix(binary_matrix)

    exp_sfs = expected_sfs_coal(len(binary_matrix[0]))
    obs_sfs = observed_sfs_coal(binary_matrix)
    print(exp_sfs)
    print(obs_sfs)

    print(len(exp_sfs), len(obs_sfs))
    plot_sfs_bokeh(
        obs_sfs, exp_sfs, "./Data/dodged_bars_TNFSF5.html", "TNFSF5"
    )

    fasta_file = "/home/jamc/Data/GitHub/PopPyGen/Data/AMELX_All_aligned.fas"
    fasta = FastaReader(fasta_file)

    seqs = []
    for record in fasta:
        seqs.append(record.seq)

    site_matrix = extract_segregating_sites(seqs)
    binary_matrix = convert_to_binary_matrix(site_matrix)
    plot_binary_matrix(binary_matrix)

    exp_sfs = expected_sfs_coal(len(binary_matrix[0]))
    obs_sfs = observed_sfs_coal(binary_matrix)
    print(exp_sfs)
    print(obs_sfs)

    print(len(exp_sfs), len(obs_sfs))
    plot_sfs_bokeh(obs_sfs, exp_sfs, "./Data/dodged_bars_AMELX.html", "AMELX")

if __name__ == '__main__':
    main()
