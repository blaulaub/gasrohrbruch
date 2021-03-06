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
R     =  518.28     # Gas constant in [J/(kg??K)]
gamma =    1.31     # specific heat ratio in [-]
cp    = 2232        # specific heat in [J/(kg??K)]


#%% Given: pipe parameters
D              = 0.3              # duct diameter in [m]
cross_section  = pi/4 * D**2      # duct cross section
L_observed     =    700           # length of the observed pipe section in [m]


#%% Given: nominal conditions at the upstream sensor and at ambient
T     =        10  + 273.15   # fluid temperature (~10 ??C) in [K]
p     = 6_000_000             # duct pressure in [Pa] (1 bar = 1e5 Pa)
u     =        10             # flow velocity in [m/s]

p_amb = 101_300               # ambient pressure in [Pa]


#%% Given: pipe friction parameter

# Note: I have no idea what the correct value is for the given pipes, so I
# picked some value from https://de.wikipedia.org/wiki/Moody-Diagramm
rohrreibungsbeiwert = 0.01       # a.k.a. "Rohrreibungszahl"
# see https://de.wikipedia.org/wiki/Rohrreibungszahl for how to estimate


#%% Chosen: model parameters
L_upstream   = 7_000   # length of the upstream pipe section in [m]
L_downstream = 7_000   # length of the downstream pipe section in [m]


#%% Helper function: ratio of total to static temperature
# This is a frequently used function of Ma, hence I name it M.
def M(Ma): return 1 + 0.5*(gamma-1) * Ma**2


#%% Helper function: lossy massflow
def compute_massflow(A, pt, Tt, Ma, zeta = 0):
    """Compute the massflow, with pressure losses due to friction being
    subtracted from the total pressure upfront.

    A    -- the cross section
    pt   -- total pressure
    Tt   -- total temperature
    Ma   -- Mach number
    zeta -- loss coefficient
    """
    M_ = M(Ma)
    c_effectiv = M_**(0.5*(gamma+1)/(gamma-1))
    c_lost    = 0 if zeta == 0 else zeta*gamma*Ma**2/sqrt(M_)
    return A * pt/sqrt(Tt) * sqrt(gamma/R) * Ma / (c_effectiv + c_lost)


#%% Helper function: loss coefficients

zeta_1s   = L_upstream / D * rohrreibungsbeiwert
zeta_s4   = (L_observed + L_downstream) / D * rohrreibungsbeiwert
zeta_t    = (L_upstream + L_observed + L_downstream) / D * rohrreibungsbeiwert

def zeta_12(L_leak):
    return (L_upstream + L_observed) / D * rohrreibungsbeiwert

def zeta_24(L_leak):
    return (L_observed-L_leak + L_downstream) / D * rohrreibungsbeiwert


#%% Helper function: Newton method

def seek(f, x_guess, x_step, y_goal, y_accuracy, tries = 10):
    """Implements the Newton method. Seeks an x for the given y_goal,
    in the limit of the given y_accuracy and the given tries.

    Keyword arguments:
    f          -- some function that takes an x and computes a y
    x_guess    -- the first x to try
    x_step     -- the increment of x for successive tries
    y_goal     -- the desired outcome
    y_accuracy -- the maximum allowable delta for termination
    tries      -- the maximum number of refinements before giving up
    """
    if tries <= 0:
        raise RecursionError("did not converge")
    y1 = f(x_guess)
    if (abs(y1-y_goal) <= y_accuracy):
        return x_guess
    else:
        y2 = f(x_guess + x_step)
        x_next_guess = x_guess + x_step * (y1-y_goal)/(y1-y2)
        return seek(f, x_next_guess, x_step, y_goal, y_accuracy, tries-1)


#%% System Identification preliminaries: sensor site
speed_of_sound = sqrt(gamma * R * T)   # local speed of sound in [m/s]
Ma             = u / speed_of_sound    # local Mach number


#%% System Identification: upstream tank
pt1 = p * (M(Ma)**(gamma/(gamma-1)) + zeta_1s*gamma*Ma**2)
Tt1 = T *  M(Ma)


#%% System Identification: downstream tank
p4 = pt1 / (M(Ma)**(gamma/(gamma-1)) + zeta_t*gamma*Ma**2)
# Note on accuracy: Here we used the Mach number computed for the nominal
# conditions, in particular for the nominal pressure, and neglected the
# Mach number increase due to the pressure drop due to the pipe friction
# along the flow path. This should be ok for the current analysis goal.

# Better accuracy could be achieved by observing the massflow balance.
# Given is the massflow from the upstream tank to the sensor:
massflow_1s = compute_massflow(cross_section, pt1, Tt1, Ma, zeta_1s)
# ...what we need is the p4 that achieves the same downstream massflow

pts = p * M(Ma)**(gamma/(gamma-1))
Tts = T * M(Ma)

Ma_s4 = seek(lambda Ma_s4 : compute_massflow(cross_section, pts, Tts, Ma_s4, zeta_s4), Ma, Ma*0.01, massflow_1s, 0.001)
massflow_s4 = compute_massflow(cross_section, pts, Tts, Ma_s4, zeta_s4)
p4 = pts / (M(Ma_s4)**(gamma/(gamma-1)) + zeta_s4*gamma*Ma_s4**2)


#%% Helper functions

def compute_p2(Ma2, L_leak):
    """Compute the static pressure at the leak site for a given Mach number.

    Keyword arguments:
    Ma2    -- the Mach number just before the leak site
    L_leak -- the position of the leak site, on the observed segment
    """
    M2 = M(Ma2)
    return pt1 / (M2**(gamma/(gamma-1)) + zeta_12(L_leak)*gamma*Ma2**2)

def compute_p4(pt2, Ma3, L_leak):
    """Compute the static pressure at the leak site for a given Mach number.

    Keyword arguments:
    pt2    -- the total pressure available at the leak site
    Ma3    -- the Mach number just after the leak site
    L_leak -- the position of the leak site, on the observed segment
    """
    M3 = M(Ma3)
    return pt2 / (M3**(gamma/(gamma-1)) + zeta_24(L_leak)*gamma*Ma3**2)


#%% Solve for leakage of 0% of pipe diameter at 50% of observed segment

A_leak = 0.0 * cross_section
L_leak = 0.5 * L_observed

def massflow_error_for_p2(p2):
    seek_Ma2_for_p2 = lambda p2 : seek(lambda Ma2: compute_p2(Ma2, L_leak), Ma, 0.01*Ma, p2, 1.)
    Ma2 = seek_Ma2_for_p2(p2)
    pt2 = p2 * M(Ma2)**(gamma/(gamma-1))
    seek_Ma3_for_p4 = lambda p4: seek(lambda Ma3: compute_p4(pt2, Ma3, L_leak), Ma2, 0.01*Ma2, p4, 1.)
    Ma3 = seek_Ma3_for_p4(p4)
    m_in  = compute_massflow(cross_section, pt1, Tt1, Ma2, zeta_12(L_leak))
    m_out = compute_massflow(cross_section, pt2, Tt1, Ma3, zeta_24(L_leak))
    return m_in - m_out

p2 = seek(massflow_error_for_p2, p, 100, 0, 0.1)


#%% Solve for leakage of 1% of pipe diameter at 50% of observed segment

def solve(L_leak, A_leak, Ma2_guess, Ma3_guess):

    def compute_before_2(p2, Ma2_guess):
        seek_Ma2_for_p2 = lambda p2 : seek(lambda Ma2: compute_p2(Ma2, L_leak), Ma2_guess, 0.01*Ma2_guess, p2, 1)
        Ma2 = seek_Ma2_for_p2(p2)
        pt2 = p2 * M(Ma2)**(gamma/(gamma-1))
        m_in  = compute_massflow(cross_section, pt1, Tt1, Ma2, zeta_12(L_leak))
        return (Ma2, pt2, m_in)

    def compute_after_2(pt2, Ma3_guess):
        seek_Ma3_for_p4 = lambda p4: seek(lambda Ma3: compute_p4(pt2, Ma3, L_leak), Ma3_guess, 0.01*Ma3_guess, p4, 1)
        Ma3 = seek_Ma3_for_p4(p4)
        m_out = compute_massflow(cross_section, pt2, Tt1, Ma3, zeta_24(L_leak))
        return (Ma3, m_out)

    def compute_leakage(pt2):
        Ma_leak = min(1., sqrt(((pt2/p_amb)**((gamma-1)/gamma)-1)*2/(gamma-1)))
        m_leak = compute_massflow(A_leak, pt2, Tt1, Ma_leak)
        return m_leak

    def massflow_error_for_p2(p2):
        (Ma2, pt2, m_in) = compute_before_2(p2, Ma2_guess=Ma)
        (Ma3, m_out)     = compute_after_2(pt2, Ma3_guess=Ma2)
        m_leak           = compute_leakage(pt2)
        return m_in - m_out - m_leak

    p2 = seek(massflow_error_for_p2, p, 100, 0, 0.01)
    (Ma2, pt2, m_in) = compute_before_2(p2, Ma2_guess)
    (Ma3, m_out)     = compute_after_2(pt2, Ma3_guess)
    m_leak           = compute_leakage(pt2)

    return (p2, Ma2, pt2, m_in, Ma3, m_out, m_leak)

A_leak = 0.01 * cross_section
L_leak = 0.5 * L_observed

(p2, Ma2, pt2, m_in, Ma3, m_out, m_leak) = solve(L_leak, A_leak, Ma2_guess=Ma, Ma3_guess=Ma)

#%% Solve for range of leakages at 50% of observed segment

def solve_all(L_leak, A_leaks, Ma_guess):
    def s(L_leak, A_leaks, Ma2_guess, Ma3_guess, result):
        if not A_leaks: return result
        A_leak = A_leaks.pop(0)
        (p2, Ma2, pt2, m_in, Ma3, m_out, m_leak) = solve(L_leak, A_leak, Ma2_guess=Ma, Ma3_guess=Ma)
        result.append((
            A_leak/cross_section,    # leakage area fraction
            p2,                      # gas pressure at leakage site
            p2/p,                    # fraction of pressure relative to sensor
            m_in,                    # massflow from upstream
            m_out,                   # massflow going downstream
            m_leak                   # massflow lost to environment
        ))
        s(L_leak, A_leaks, Ma2, Ma3, result)
        return s(L_leak, A_leaks, Ma2, Ma3, result)
    return s(L_leak, A_leaks, Ma_guess, Ma_guess, [])

#A_leaks = [0., 0.01*cross_section]
A_leaks = [ cross_section * i/100. for i in range(6) ]
result = solve_all(L_leak, A_leaks, Ma)

# Note: for the given model (where rohrreibungsbeiwert, L_upstream and
# L_downstream are just guesses), at 6% leakage cross section, the pressure
# already drops so low that fluid comes back from downstream.
