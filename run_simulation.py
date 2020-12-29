import math
import pandas as pd
import numpy as np
import random


def simulate_dat(curr_pen, prevalence, n_iterations, bottleneck_gen = 30, N_bottleneck = 500, fecundity_ratio_females = 0.5, fecundity_ratio_males = 0.25):
    
    #Values are set based on empirical observations 
    bottleneck_N_chr =  N_bottleneck*2
    max_pop_N_chr = 20000000

    total_growth = max_pop_N_chr / bottleneck_N_chr

    growth_rate = math.exp(np.log(total_growth)/bottleneck_gen)

    mean_off_per_woman = 2 * growth_rate
    mean_off_per_man = 2 * growth_rate
    mean_off = growth_rate 

    offspring_cases = (mean_off*(fecundity_ratio_females + fecundity_ratio_males))/2

    mean_off_carriers = (curr_pen * offspring_cases) + ((1 - curr_pen) * mean_off)

    carriers = np.zeros([bottleneck_gen, n_iterations], 'int')
    noncarriers = np.zeros([bottleneck_gen, n_iterations], 'int')

    noncarriers[0, :] = N_bottleneck - 1
    carriers[0, :] = 1

    #Poisson process
    for gen in range(1, bottleneck_gen): 
        noncarriers[gen, :] = np.random.poisson(mean_off  * noncarriers[gen - 1, :])
        carriers[gen, :] = np.random.poisson(mean_off_carriers  * carriers[gen - 1, :])
    return noncarriers[gen, :], carriers[gen, :]




n_iterations = 10000
penetrance = np.arange(0.05, 0.3, 0.05)
prevalence = 0.01
output_dat = np.zeros([n_iterations * 5, 5])

for index in range(0, len(penetrance)):
    curr_pen = penetrance[index]
    noncars, cars = simulate_dat(curr_pen, prevalence, n_iterations, 25, 300, 0.5, 0.25)
    noncar_cases = np.floor(noncars * prevalence)
    noncar_controls = noncars - noncar_cases
    car_cases = np.floor(cars * curr_pen)
    car_controls = cars - car_cases
    first = n_iterations * index
    last = n_iterations * (index + 1)
    pens = np.ones([n_iterations,]) * curr_pen
    output_dat[range(first, last), :] = curr_pen
    output_dat[range(first, last), 1] = car_cases
    output_dat[range(first, last), 2] = car_controls
    output_dat[range(first, last), 3] = car_controls + noncar_controls
    output_dat[range(first, last), 4] = car_cases + noncar_cases
df = pd.DataFrame(output_dat)
df.to_csv('scenario_1.csv', columns = ['Penetrance', 'CaseCars', 'ControlCars', 'TotControls', 'TotCases'])




