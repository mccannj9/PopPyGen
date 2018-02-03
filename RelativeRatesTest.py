#! /usr/bin/env python

from Parsers.Fasta import FastaReader


def n_choose_k(n, k):

    partials = 1
    for val in range(1, k+1):
        partials *= (n-k+val)/(val)

    return partials


def compute_binomial_probability(n, k, p):

    combos = n_choose_k(n, k)
    return combos*(p**k)*((1-p)**(n-k))


def compute_binomial_cdf(n_, p_, from_, to_):

    prob = 0
    for val in range(from_, to_+1):
        prob += compute_binomial_probability(n_, val, p_)

    return prob


def compute_significance(n, k, p):
    if k > (n - k):
        k = (n - k)
    p1 = compute_binomial_cdf(n, p, 0, k)
    p2 = compute_binomial_cdf(n, p, n-k, n)

    return p1 + p2




def main():
    from math import factorial
    n = 485
    k = 286

    # computes bin prob for exercise 1.5
    # performing something like two-tailed binomial test
    sig = compute_significance(n, k, 1/2)
    print("PROB = %s" % (sig))

    n = 8
    k = 4

    # computes bin prob for exercise 1.6
    # performing something like two-tailed binomial test
    sig = compute_significance(n, k, 1/2)
    print("PROB = %s" % (sig))

if __name__ == '__main__':
    main()
