#! /usr/bin/env python

import random

def population_simulation(pop_size, num_gens, init_freq=0.5):

    trials = 2*pop_size

    allele_freqs = []

    for x in range(num_gens):
        success = 0
        for y in range(trials):
            if random.uniform(0,1) < init_freq:
                success += 1

        init_freq = success/trials
        allele_freqs.append(init_freq)

    return allele_freqs


def main():

    from Plotters.BokehWrappers import plot_single_drift_simulation

    nsims = 15

    # simulations = []
    # for x in range(nsim):
    #     simulations.append()

    freqs = population_simulation(50, 200, init_freq=0.5)
    for x in range(len(freqs)):
        print("Gen %s: %s" % (x+1, freqs[x]))

    plot_single_drift_simulation(list(range(1,201)), freqs, "./Data/test.html")

    return 0


if __name__ == '__main__':
    main()
