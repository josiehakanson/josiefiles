# Uploading data, cleaning, and defining constants

import numpy as np
import math as mt
import pandas as pd
import matplotlib.pyplot as plt

# %matplotlib inline

# skipping separator rows
skip_rows = [999999, 1999999, 2999999, 3999999, 4999999, 5999999, 6999999, 7999999]

data_func = pd.read_csv('MCSTAS.txt', delim_whitespace = True, header = None, skiprows = skip_rows)

data_func = data_func.apply(pd.to_numeric, errors='coerce')

# Drop any rows with NaN values
data_func = data_func.dropna()

data = data_func.to_numpy()

weight_values = data[:, 0]


x_values = data[:, 1]
y_values = data[:, 2]
# z_values = data[:, 3]
v_x = data[:, 4]
v_y = data[:, 5]
v_z = data[:, 6]

test_values = [weight_values[0], x_values[0], y_values[0], v_x[0], v_y[0], v_z[0]]

weight_sum = np.sum(weight_values)
events = len(weight_values)
avg_vel = np.average(v_z, weights = weight_values)




# Differing magnetic field for each neutron event

# events_lol = 100000           # test amount of events to look at

# Constants
hbar = 6.58e-16 # eV*s
kappa = 1e-5 # relative nTmm (goes to 1e-6?)
    # unknown value - fraction of the nTmm, using 1e-5  for testing
eps_1 = 10e-17 # epsilon in eV
mu = 6.3e-12 # normal dipole in eV/G
length = 1.727         # lenght of one magnet in meters


# constants that will change
t = 0
i = 0
B = 24 # Magnetic field in gauss (Will be varied from 24-26 in steps)

# step sizes
dt = length/avg_vel          # weighted average velocity used for step size of magnetic field
    #1E-7 # 1 millisec
   # = 2m/v
dB = (0.04/3600) * dt # step size for B_gauss
    # delta_Bgauss = 0.04 (gauss/hr)
    # steps of 0.04 gauss, going from 24-26 gauss in one hour
dz = 0.001 # step size im meters
ms_hr = 3.6e6               # milliseconds per hour
time_steps = 100000         # test amount of step sizes to look at

# stuff i'm keeping track of
Prob_array = [0] * events # assign an amount of values to this
event_prob_avg = [0] * events
eps_array = [0] * int(2/dB)         # amount of steps is 2/dB bc B varies across 2 gauss (from 24-26)
B_gauss = [0] * int(2/dB)
time = [0] * events
Prob_sum = [0] * int(2/dB)
P_avg = [0] * int(2/dB)

# equations needed
V_f = 1.51e-10 # Optical potential in eV (determined by B gauss)
    # above number is 6.03e-12*25 for the 25 gauss determined for the test
# V_f = (6.03e-12) * B_gauss
delta_u = V_f - (mu*B)
eps = eps_1 + (kappa*mu*B)
omega = np.sqrt(((delta_u/2)**2) + (eps**2))
# Prob_nn = ((eps**2)/(omega**2)) * ((np.sin((omega*t)/hbar))**2)
# eps_bv = eps_1 + (kappa*mu*B)

# print("Probability at t=0:", Prob_nn)

k = 0
j = 0
sum = 0
angle = 0

while B <= 26:           # loop that varies magnetic field
    B_gauss[i] = B                  # keeping track of magnetic field
    eps = eps_1 + (kappa*mu*B)      # calculating epsilon for each magnetic field
    V_f = (6.03e-12) * B            # optical potential for this magnetic field

    k = 0
    t_tot = 0               # resetting total time to cross magnet for each event
    time_steps = 0
    while k < events:       # cycling through each event (out of 8 million)
        vel = v_z[k]        # velocity of event being looked at
        dt = dz/vel         # determining step size based on velocity
        t_tot = length/vel                     # total time to cross the first magnet
        time_steps = int(t_tot/dt) + 1      # amount of steps needed to account for time it takes to cross the magnet

        sum = 0
        j = 0
        t = 0               # setting everything to zero before going in the loop
        while j < time_steps:               # loop calculating probability of this neutron with a specific magnetic field (varying through time)
            B = B_gauss[i]                  # ensuring we're looking at the right magnetic field                                                           # calculating optical potential and everything else for this event

            delta_u = V_f - (mu*B)
            omega = np.sqrt(((delta_u/2)**2) + (eps**2))

            angle = (omega * t) / hbar  # Confirm units result in radians

            Prob_array[j] = ((eps**2)/(omega**2)) * ((np.sin(angle))**2)           # is this in degrees or radians? how do i make sure?

            sum += Prob_array[j]            # accumulated sum from this event
            #Prob_sum[i] += Prob_array[j]
            #eps_array[k] = eps_bv
            time[j] = t               # keeping track of time just in case

            t += dt
            #print(Prob_sum[k])
            j += 1

        event_prob_avg[k] = sum / time_steps        # average probability of the event we're looking at
        Prob_sum[i] += event_prob_avg[k]            # sum of the all the average probabilities for this specific magnetic field
        k += 1

    P_avg[i] = Prob_sum[i]/weight_sum   # array of averge probabilities
    eps_array[i] = eps      # keeping track of epsilon for each magnetic field just in case
    B += dB                 # updating values for the next run of the loop
    i += 1


# printing valuable information to a .txt file
f = open("output_dec.txt", "a")
print("Probabilies:", P_avg, file=f)
print("Sum of probabilities:", Prob_sum, file=f)
print("B in Gauss:", B_gauss, file=f)
print("Epsilon:", eps_array, file=f)
f.close()


# plotting Probability as a function of magnetic field, and saving
plt.plot(B_gauss, P_avg)
plt.ylabel('Probability n to n')
plt.xlabel('Magnetic Field (gauss)')

plt.savefig('simulation_dec.png')

# if this plot is linear again im going to kms
# i am hoping that factoring in all the events will fix this