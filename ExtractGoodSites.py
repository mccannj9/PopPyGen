#! /usr/bin/env python


def main():
    from Preppers.Population import ex_seg_sites_with_constant
    from Parsers.SeqIO import PhylipReader
    phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/sback.phy")
    out = open("/home/jamc/Data/GitHub/PopPyGen/Data/sback_snps_with_constants.phy", "w")

    # phy = PhylipReader("/home/jamc/Data/GitHub/PopPyGen/Data/test.phy")
    # out = open("/home/jamc/Data/GitHub/PopPyGen/Data/test_snps.phy", "w")

    seqs = list(phy)
    print(len(seqs[0].seq))
    segsites = ex_seg_sites_with_constant(seqs)
    # print(segsites)

    taxa = len(segsites)
    sites = len(segsites[0])
    print("%s %s" % (taxa, sites))
    print("%s %s" % (taxa, sites), file=out)

    for x in range(taxa):
        taxa_id = seqs[x].name + "\t"
        site = "".join(segsites[x])
        print(taxa_id + site, file=out)

    out.close()
    return 0


if __name__ == '__main__':
    main()
