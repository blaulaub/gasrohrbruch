# -*- coding: utf-8 -*-

from math import pi, sqrt


#%% methane gas properties

R     =  518.28     # Gas constant in [J/(kg·K)]
gamma =    1.31     # specific heat ratio in [-]
cp    = 2232        # specific heat in [J/(kg·K)]


#%% pipe parameters

D              = 0.3              # duct diameter in [m]
cross_section  = pi/4 * D**2      # duct cross section


#%% conditions at the upstream sensor

T     =        10  + 273.15   # fluid temperature (~10 °C) in [K]
p     = 6_000_000             # duct pressure in [Pa] (1 bar = 1e5 Pa)
u     =        10             # flow velocity in [m/s]


#%% other values at the upstream sensor position

speed_of_sound = sqrt(gamma * R * T)   # local speed of sound in [m/s]
Ma             = u / speed_of_sound    # local Mach number
density        = p / (R*T)             # density in [kg/m³]


#%% ideal upstream reservoir conditions (neglecting pipe friction)
#   (these are only estimates, may be used for plausibility checks)

# Note: "upstream reservoir" means an infinite volume of gas at rest;
# this is similar to an "ideal voltage source" in electrical engineering

Tt   = T * (1 + (gamma-1)/2 * Ma**2)
pt   = p * (1 + (gamma-1)/2 * Ma**2)**(gamma/(gamma-1))


#%% ideal flow
#   (neglecting friction, hence the pipe length does not matter)

def ideal_massflow_density(gamma, R, pt, Tt, Ma):
    M = 1 + 0.5*(gamma-1) * Ma**2
    c1 = M**(0.5*(1+gamma)/(1-gamma))
    return pt/sqrt(Tt) * sqrt(gamma/R) * Ma * c1

massflow = cross_section * ideal_massflow_density(gamma, R, pt, Tt, Ma)
# computes to 28.9 kg/s for the given conditions


#%% pipe friction parameter

# Note: I have no idea what the correct value is for the given pipes, so I
# picked some value from https://de.wikipedia.org/wiki/Moody-Diagramm
rohrreibungsbeiwert = 0.01       # a.k.a. "Rohrreibungszahl"
# see https://de.wikipedia.org/wiki/Rohrreibungszahl for how to estimate


#%% real-world model parameters

# Note: just like a real-world voltage source is modeled by combining
# an ideal voltage source with an inner resistance, the real-world upstream
# reservoir is modelled by combining an ideal upstream reservoir (an infinite
# volume of gas at rest) with an outlet pipe of finite length; here, the
# outlet pipe is called "upstream" pipe.

# Note: the "upstream length" here is very similar to the "inner resistance"
# in electrical engineering; it needs to BE TUNED ACCORDING TO MEASUREMENTS or
# BE ESTIMATED ON SOUND ASSUMPTIONS, otherwise model predictions with regard
# to leakage flows (pipe damage in fluid systems, similar to leak currents or
# short circuits in electrical circuits) will be far off.

# Also note: in reality, there probably isn't a pressurized reservoir, but some
# compressor driven by some motor; I am not sure how adequate the "reservoir +
# outlet pipe" modeling is for such a system.

# Note: I have no idea what the correct value is for the given system, so I
# picked some numbers...

L_upstream   = 10_000   # length of the upstream pipe section in [m]
L_observed   =    500   # length of the observed pipe section in [m]
L_downstream = 10_000   # length of the downstream pipe section in [m]

# ...and the corresponding loss coefficients
zeta_upstream = L_upstream / D * rohrreibungsbeiwert
zeta_observed = L_observed / D * rohrreibungsbeiwert
zeta_downstream = L_downstream / D * rohrreibungsbeiwert


#%% real-world upstream reservoir parameters

# pressure loss due to pipe friction
# (see https://de.wikipedia.org/wiki/Rohrreibungszahl)
p_upstream    = p + zeta_upstream * density/2 * u**2

# upstream reservoir conditions including pressure loss (due to pipe friction,
# pressure reduces along the way) and assuming adiabatic flow (no heat loss,
# since the gas is roughly at ambient temperature)
pt = p_upstream * (1 + (gamma-1)/2 * Ma**2)**(gamma/(gamma-1))
Tt = T * (pt/p)**((gamma-1)/gamma)




#%% real-world downstream exit conditions

# This is simply the upstream pressure minus all pressure losses along the way
zeta_total = zeta_upstream + zeta_observed + zeta_downstream
p_downstream = p_upstream - zeta_total * density/2 * u**2

#%% NOMINALly operating system so far

assert 'Tt' in vars()              # total temperature upstream
assert 'pt' in vars()              # total pressure upstream
assert 'zeta_upstream' in vars()   # friction losses while comming from upstream
assert 'zeta_observed' in vars()   # friction losses while going through observed
assert 'zeta_downstream' in vars() # friction losses while going downstream
assert 'p_downstream' in vars()    # back pressure downstream

#%% sanity check: compute real flow at sensor position, from upstream

def real_massflow_density(gamma, R, pt, Tt, zeta, Ma):
    M = 1 + 0.5*(gamma-1) * Ma**2
    c1 = M**(0.5*(1+gamma)/(1-gamma))
    c2 = zeta * 0.5*Ma**2 * M**(0.5*(3*gamma-1)/(gamma-1))
    assert c2 < c1  # otherwise, value of Ma is not possible
    return pt/sqrt(Tt) * sqrt(gamma/R) * Ma * (c1 - c2)

real_massflow = cross_section * real_massflow_density(gamma, R, pt, Tt, zeta_upstream, Ma)

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





#%% remainder

L     =     1_000             # nominal duct length in [m]

# ambient conditions

# derived values
p              = P * p_amb                     # duct fluid pressure in [Pa]
mass_flow      = density * cross_section * u   # mass flow in kg/s

# SO FAR
# - velocity is small compared to sound (may be neglected)
# - sound takes ~2 seconds to travel 1 km
# - mass flow is ~29 kilogram per second

# QUESTIONS
# - how to model the system ?
# - how to model a damage event (sudden duct breakage) ?
# - how to model the leakage ?
# - what is the detectable wave (amplitude, dispersion) ?
# - what are the detectable characteristics (e.g., permanent pressure drop) ?
# - what is the impact of the leakage area (from 0 to 1 duct cross section) ?

# KEYWORDS
# - acoustic waves in a duct
# - compression waves in a duct




# UPFRONT: static analysis of pipe an leakage flow

# SIMPLE ESTIMATE: Mach number for given duct-to-ambient pressure ratio
# (maximum possible Mach number for discharge to ambient)
max_discharge_Mach_number = sqrt(2/(gamma-1)*(P**((gamma-1)/gamma)-1))
# - computes to Mach 3.2
#   (a higher value means more energy potential in the fluid in the duct)
#   (a Mach number >1 means that the duct conditions are not the limiting factor for leakage flow)

# SIMPLE ESTIMATE: discharge mass flow for given duct conditions and
# given a leakage cross section equal to the duct cross section, at Ma=1
coeff1 = cross_section * P * p_amb/sqrt(T)
coeff3 = 0.5 * (1+gamma)/(1-gamma)
max_discharge_mass_flow = coeff1 * sqrt(gamma/R) * (1 + (gamma-1)/2)**coeff3
# - computes to 750 kg/s for leakage = 100% duct cross section
#   (under optimum conditions, e.g., breakage right after the pipe start)

# SIMPLE ESTIMATE: fraction of leakage-to-duct cross section, where at Ma=1
# the leakage mass flow is equal to the pipe mass flow
f = mass_flow/cross_section / coeff1 * sqrt(R/gamma) * (1+(gamma-1)/2)**(-coeff3)
# - computes to 0.55
#   (means that given a hole of 55% of the duct cross section and sufficient
#   downstream pipe resistance, all gas would leak)



#%% methane hydraulic properties
#   (required to compute friction effects)

# Note: I did not find precises numbers for that on Google Search; there are
# equations to approximate this based on temperature and pressure, but that
# takes more time.

# Here, I just use some approximate numbers
# (from Google Search, but not for the correct pressure and/or temperature):
dynamic_viscosity_cp = 0.01107                      # in Centipoises [cP]
dynamic_viscosity    = dynamic_viscosity_cp * 0.001 # in Pascal seconds [P·s]

