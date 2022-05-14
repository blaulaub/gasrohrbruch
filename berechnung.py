# -*- coding: utf-8 -*-

from math import pi, sqrt

#%% Problem Description
#
# We look at a section of a gas pipeline where we want to detect damages
# to the pipe that lead to gas leakages, by means of a sensor that measures
# pressure.


#%% Model Description
#
# For this analysis, there is one "sensor position" where nonimal conditions
# (temperature, pressure, velocity) are defined, and where pressure is
# measured.
#
# We considere three consecutive pipe segments. The first pipe segment is
# the observed segment, the leakage can occur anywhere within this segment.
# The sensor position is at the start of the observed segment.
#
# To approximate a realistic system, we assume there are two gas tanks
# upstream and downstream of the observed segment. Both are assumed to be of
# infinite volume, in order to provide stationary boundary conditions for the
# analysis. The tanks are assumed to be connected to the observed segments by
# an upstream and a downstream pipe segment, respectively.
#
# The combination of tank and pipe segments on both sides, upstream and
# downstream, acts as a model for a real gas source or sink, respectively;
# that should be a sufficient approximation for moderate system variations,
# but might be incorrect for extreme conditions (e.g., when the mass flow
# deviates significantly from the nominal mass flow, or when the pipe segments
# are not significantly longer than the observed pipe segment).
#
# We use index 1 for the upstream tank and the beginning of the upstream
# segment, index s for the sensor position, which coincides with the end of
# the upstream segment and the beginning of the observed segment, index 2 for
# the leakage position along the observed segment, and index 4 for the
# downstream tank an the end of the downstream segment.
#
# We reserve but don't use index 3 for the connection from observed to
# downstream; but presumably there will be the next sensor, so this could be
# interesting, too.


#%% Given: gas properties of methane
R     =  518.28     # Gas constant in [J/(kg·K)]
gamma =    1.31     # specific heat ratio in [-]
cp    = 2232        # specific heat in [J/(kg·K)]


#%% Given: pipe parameters
D              = 0.3              # duct diameter in [m]
cross_section  = pi/4 * D**2      # duct cross section
L_observed     =    700           # length of the observed pipe section in [m]


#%% Given: nominal conditions at the upstream sensor
T =        10  + 273.15   # fluid temperature (~10 °C) in [K]
p = 6_000_000             # duct pressure in [Pa] (1 bar = 1e5 Pa)
u =        10             # flow velocity in [m/s]


#%% Given: pipe friction parameter

# Note: I have no idea what the correct value is for the given pipes, so I
# picked some value from https://de.wikipedia.org/wiki/Moody-Diagramm
rohrreibungsbeiwert = 0.01       # a.k.a. "Rohrreibungszahl"
# see https://de.wikipedia.org/wiki/Rohrreibungszahl for how to estimate


#%% Chosen: model parameters
L_upstream   = 7_000   # length of the upstream pipe section in [m]
L_downstream = 7_000   # length of the downstream pipe section in [m]


#%% Dependent: loss coefficients
zeta_upstream   = L_upstream   / D * rohrreibungsbeiwert
zeta_observed   = L_observed   / D * rohrreibungsbeiwert
zeta_downstream = L_downstream / D * rohrreibungsbeiwert


#%% Helper function: ratio of total to static temperature
# This is a frequently used function of Ma, hence I name it M.
M = lambda Ma : 1 + 0.5*(gamma-1) * Ma**2


#%% System Identification preliminaries: sensor site
speed_of_sound = sqrt(gamma * R * T)   # local speed of sound in [m/s]
Ma             = u / speed_of_sound    # local Mach number


#%% System Identification: upstream tank
pt1 = p * (M(Ma)**(gamma/(gamma-1)) + zeta_upstream*gamma*Ma**2)
Tt1 = T *  M(Ma)


#%% System Identification: downstream tank
zeta_total = zeta_upstream + zeta_observed + zeta_downstream
p4 = pt1 / (M(Ma)**(gamma/(gamma-1)) + zeta_total*gamma*Ma**2)
# Note on accuracy: Here we used the Mach number computed for the nominal
# conditions, in particular for the nominal pressure, and neglected the
# Mach number increase due to the pressure drop due to the pipe friction
# along the flow path. This should be ok for the current analysis goal.






#%% other values at the upstream sensor position
density        = p / (R*T)             # density in [kg/m³]


#%% system identification

# pressure loss along the upstream section due to pipe friction
# (see https://de.wikipedia.org/wiki/Rohrreibungszahl)
p_loss_upstream = zeta_upstream * density/2 * u**2

# upstream reservoir conditions including pressure loss (due to pipe friction,
# pressure reduces along the way) and assuming adiabatic flow (no heat loss,
# since the gas is roughly at ambient temperature)
p_upstream  = p + p_loss_upstream
pt_upstream = p_upstream * (1 + (gamma-1)/2 * Ma**2)**(gamma/(gamma-1))
Tt_upstream = T * (pt_upstream/p_upstream)**((gamma-1)/gamma)

# This is simply the upstream pressure minus all pressure losses along the way
zeta_total = zeta_upstream + zeta_observed + zeta_downstream
p_downstream = p_upstream - zeta_total * density/2 * u**2


#%% sanity check: compute real flow at sensor position, from upstream

def real_massflow_density(gamma, R, pt, Tt, zeta, Ma):
    M = 1 + 0.5*(gamma-1) * Ma**2
    c1 = M**(0.5*(1+gamma)/(1-gamma))
    c2 = zeta * 0.5*Ma**2 * M**(0.5*(3*gamma-1)/(gamma-1))
    assert c2 < c1  # otherwise, value of Ma is not possible
    return pt/sqrt(Tt) * sqrt(gamma/R) * Ma * (c1 - c2)

real_massflow = cross_section * real_massflow_density(gamma, R, pt1, Tt1, zeta_upstream, Ma)

#%% LEAKAGE -- from undamaged to damaged system

# So far we had an upstream reservoir, three pipe sections (upstream, observed,
# downstream), and a downstream reservoir.

# When a leakage occurs in the observed pipe section, a third reservoir
# (ambient) is connected. Because the pressure ratio from pipe to ambient is
# so high (60 bar duct pressure to 1 bar ambient pressure), this results in
# high-Mach-number flow (instead of 10 meter per second along the pipe, we get
# velocities in the order of to the speed of sound).

# The leakage occurs within the observed section
# (the observed section starts with the pressure sensor)
L_leakage = 300        # distance downstream from pressure sensor
assert L_leakage > 0 and L_leakage <= L_observed

# The leakage is assumed to have a circular cross section/diameter
# (the maximum value possible is the pipe diameter)
D_leakage = 0.1 * D    # leakage hole diameter in [m]
assert D_leakage > 0 and D_leakage <= D

# The leakage opens to ambient, with ambient pressure
# (with respect to a scale of 60 bar, this may be considered constant)
p_amb  = 101_300              # nominal ambient pressure in [Pa]
