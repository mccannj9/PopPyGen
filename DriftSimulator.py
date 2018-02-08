#! /usr/bin/env python

import random

def population_simulation(pop_size, num_gens, init_freq=0.5):

    trials = 2*pop_size

    # adding initial freq so all sims start at same value in plot
    allele_freqs = [init_freq]

    for x in range(num_gens):
        success = 0
        for y in range(trials):
            if random.uniform(0,1) < init_freq:
                success += 1

        init_freq = success/trials
        allele_freqs.append(init_freq)

    return allele_freqs


def main():

    from Plotters.BokehWrappers import plot_multiple_drift_simulations

    nsims = 10
    popsize = 500
    gens = 500

    simulations = []
    for x in range(nsims):
        freqs = population_simulation(popsize, gens, init_freq=0.5)
        simulations.append(freqs)

    # xdata for bokeh
    xgenerations = [list(range(0,gens+1)) for y in range(nsims)]

    plot_multiple_drift_simulations(
        xgenerations, simulations, "./Data/test.html"
    )

    return 0


if __name__ == '__main__':
    main()
