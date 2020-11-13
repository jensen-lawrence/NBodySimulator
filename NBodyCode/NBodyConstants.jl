# -----------------------------------------------------------------------------------------------------------------------------------------

#= Global Constant Definitions =#

# Physical constants
const M_SUN = 1.98847e30
const R_SUN = 6.957e8
const M_JUPITER = 1.898e27
const R_JUPITER = 6.9911e7
const M_EARTH = 5.9722e24
const R_EARTH = 6.371e6

# Reference scales
const MINUTE = 60
const HOUR = 60*MINUTE
const DAY = 24*HOUR
const WEEK = 7*DAY
const MONTH = DAY*(365/12)
const YEAR = DAY*365
const AU = 149597870700
const LY = 9454254955488000
const KMS = 1000
const KMH = 5/18

# -----------------------------------------------------------------------------------------------------------------------------------------

#= Body Structure Definition =#

"""
    mutable struct Body

Body : mutable struct
    | mutable struct that represents a celestial body.

name : String
    | The name of the celestial body. Accessed by body.name.

mass : Float64
    | The mass of the celestial body. Accessed by body.mass.
    Measured in kg.

radius : Float64
    | The radius of the celestial body. Accessed by body.radius.
    Measured in m.

x : Vector{Float64}
    | The position of the celestial body in [x, y, z]. Accessed by body.x.
    Measured in m.

v : Vector{Float64}
    | The velocity of the celestial body in [u, v, w]. Accessed by body.v.
    Measured in m/s.

h : Float64
    | The step size of the celestial body. Determines how often the body is
    updated in simulations.
    Measured in s.
"""
mutable struct Body
    name::String
    mass::Float64
    radius::Float64
    x::Vector{Float64}
    v::Vector{Float64}
    h::Float64
end

# -----------------------------------------------------------------------------------------------------------------------------------------

#= Object Group Definitions =#

# Earth-Moon System
Earth_EM = Body("Earth", M_EARTH, R_EARTH, [-4.67066e6, 0.0, 0.0], [0.0, -12.57, 0.0], 150)
Moon_EM = Body("Moon", 0.0123*M_EARTH, 0.2727*R_EARTH, [3.79728e8, 0.0, 0.0], [0.0, 1022.0, 0.0], 150)
EarthMoon = [Earth_EM, Moon_EM]

# Galilean Moons System
Jupiter_JM = Body("Jupiter", M_JUPITER, R_JUPITER, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 150)
Io_JM = Body("Io", 8.931938e22, 1.8216e6, [4.217e8, 0.0, 0.0], [0.0, 17334.0, 0.0], 150)
Europa_JM = Body("Europa", 4.799844e22, 1.5608e6, [6.709e8, 0.0, 0.0], [0.0, 13740.0, 0.0], 150)
Ganymede_JM = Body("Ganymede", 1.4819e23, 2.6341e6, [1.0704e9, 0.0, 0.0], [0.0, 10880.0, 0.0], 150)
Callisto_JM = Body("Callisto", 1.075938e23, 2.4103e6, [1.8827e9, 0.0, 0.0], [0.0, 8204.0, 0.0], 150)
GalileanMoons = [Jupiter_JM, Io_JM, Europa_JM, Ganymede_JM, Callisto_JM]

# Saturn Moons System
Saturn_SM = Body("Saturn", 0.29944*M_JUPITER, 0.83294*R_JUPITER, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 150)
Enceladus_SM = Body("Enceladus", 1.08022e20, 2.521e5, [2.37948e8, 0.0, 0.0], [0.0, 12622.0, 0.0], 150)
Dione_SM = Body("Dione", 1.095452e21, 5.614e5, [3.77396e8, 0.0, 0.0], [0.0, 10022.0, 0.0], 150)
Rhea_SM = Body("Rhea", 2.306518e21, 7.638e5, [5.27108e8, 0.0, 0.0], [0.0, 8480.0, 0.0], 150)
Titan_SM = Body("Titan", 1.3452e23, 2.57473e6, [1.22187e9, 0.0, 0.0], [0.0, 5570.0, 0.0], 150)
Iapetus_SM = Body("Iapetus", 1.805635e21, 7.345e5, [3.56082e9, 0.0, 0.0], [0.0, 3260.0, 0.0], 150)
SaturnMoons = [Saturn_SM, Enceladus_SM, Dione_SM, Rhea_SM, Titan_SM, Iapetus_SM]

# Solar System
Sun_SS = Body("Sun", 1.0*M_SUN, 1.0*R_SUN, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.01*DAY)
Mercury_SS = Body("Mercury", 0.055*M_EARTH, 0.3829*R_EARTH, [0.387098*AU, 0.0, 0.0], [0.0, 47362.0, 0.0], 0.01*DAY)
Venus_SS = Body("Venus", 0.815*M_EARTH, 0.9499*R_EARTH, [0.723332*AU, 0.0, 0.0], [0.0, 35020.0, 0.0], 0.04*DAY)
Earth_SS = Body("Earth", 1.0*M_EARTH, 1.0*R_EARTH, [1*AU, 0.0, 0.0], [0.0, 29780.0, 0.0], 0.08*DAY)
Mars_SS = Body("Mars", 0.107*M_EARTH, 0.532*R_EARTH, [1.523679*AU, 0.0, 0.0], [0.0, 24007.0, 0.0], 0.15*DAY)
Jupiter_SS = Body("Jupiter", M_JUPITER, R_JUPITER, [5.2044*AU, 0.0, 0.0], [0.0, 13070.0, 0.0], 1.0*DAY)
Saturn_SS = Body("Saturn", 0.29944*M_JUPITER, 0.83294*R_JUPITER, [9.5826*AU, 0.0, 0.0], [0.0, 9680.0, 0.0], 2.0*DAY)
Uranus_SS = Body("Uranus", 14.536*M_EARTH, 3.981*R_EARTH, [19.2184*AU, 0.0, 0.0], [0.0, 6800.0, 0.0], 4.0*DAY)
Neptune_SS = Body("Neptune", 17.147*M_EARTH, 3.865*R_EARTH, [30.07*AU, 0.0, 0.0], [0.0, 5430.0, 0.0], 10.0*DAY)
SolarSystem = [Sun_SS, Mercury_SS, Venus_SS, Earth_SS, Mars_SS, Jupiter_SS, Saturn_SS, Uranus_SS, Neptune_SS]

# Inner Solar System
InnerSolarSystem = [Sun_SS, Mercury_SS, Venus_SS, Earth_SS, Mars_SS]

# TRAPPIST-1 System
TRAPPIST1 = Body("TRAPPIST-1", 0.0898*M_SUN, 0.1192*R_SUN, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 250)
TRAPPIST1B = Body("TRAPPIST-1b", 1.374*M_EARTH, 1.116*R_EARTH, [0.01154*AU, 0.0, 0.0], [0.0, 83096.0, 0.0], 250)
TRAPPIST1C = Body("TRAPPIST-1c", 1.308*M_EARTH, 1.097*R_EARTH, [0.01580*AU, 0.0, 0.0], [0.0, 71016.0, 0.0], 250)
TRAPPIST1D = Body("TRAPPIST-1d", 0.388*M_EARTH, 0.778*R_EARTH, [0.02227*AU, 0.0, 0.0], [0.0, 59817.0, 0.0], 250)
TRAPPIST1E = Body("TRAPPIST-1e", 0.692*M_EARTH, 0.920*R_EARTH, [0.02925*AU, 0.0, 0.0], [0.0, 52194.0, 0.0], 250)
TRAPPIST1F = Body("TRAPPIST-1f", 1.039*M_EARTH, 1.045*R_EARTH, [0.03849*AU, 0.0, 0.0], [0.0, 45450.0, 0.0], 250)
TRAPPIST1G = Body("TRAPPIST-1g", 1.321*M_EARTH, 1.129*R_EARTH, [0.04683*AU, 0.0, 0.0], [0.0, 41250.0, 0.0], 250)
TRAPPIST1H = Body("TRAPPIST-1h", 0.326*M_EARTH, 0.775*R_EARTH, [0.06189*AU, 0.0, 0.0], [0.0, 35882.0, 0.0], 250)
Trappist1 = [TRAPPIST1, TRAPPIST1B, TRAPPIST1C, TRAPPIST1D, TRAPPIST1E, TRAPPIST1F, TRAPPIST1G, TRAPPIST1H]

# Stars
Star1 = Body("Star 1", M_SUN, 1, [0.0, 0.0, 0.0], [1e4, 1e4, 1e4], 10)
Star2 = Body("Star 2", M_SUN, 1, [1e10, 0.0, 0.0], [-3e4, 2e4, 0.0], 10)
Star3 = Body("Star 3", M_SUN, 1, [0.0, 1e10, 0.0], [-1e3, 1e4, 4e4], 10)
Star4 = Body("Star 4", M_SUN, 1, [0.0, 0.0, 1e10], [0.0, 0.0, 1e4], 10)
Star5 = Body("Star 5", M_SUN, 1, [1e10, 1e10, 0.0], [-3e4, -3e4, 2e4], 10)
Star6 = Body("Star 6", M_SUN, 1, [0.0, 1e10, 1e10], [0.0, 0.0, 0.0], 10)
Star7 = Body("Star 7", M_SUN, 1, [1e10, 0.0, 1e10], [3e4, 2e4, -4e4], 10)
Star8 = Body("Star 8", M_SUN, 1, [1e10, 1e10, 1e10], [-2e4, 2e4, -1e4], 10)
Stars = [Star1, Star2, Star3, Star4, Star5, Star6, Star7, Star8]

# -----------------------------------------------------------------------------------------------------------------------------------------