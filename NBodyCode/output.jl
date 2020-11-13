# ---------------------------------------------------------------------------------------------------------------------

#= NBodySimulation Imports =#
include("NBodyConstants.jl")
include("NBodySolver.jl")
include("NBodyGraphing.jl")

# ---------------------------------------------------------------------------------------------------------------------

#= Simulation Parameters and Execution =#

# Pick one of the pre-defined object groups, or define your own
bodies = Stars

# Pick the time over which the simulation will be run
# Default time is in seconds. You can can also use other time constants defined in NBodyConstants.jl
T = 5.0*DAY

# Running the simulation and acquiring the data
data = runsimulation(bodies, T)

# Change this to the absolute path to where the plots should be stored
graphdir = "C:\\Users\\Jensen\\Documents\\Projects\\NBodySimulator\\"

# ---------------------------------------------------------------------------------------------------------------------

#= Plotting Simulation Results =#

# Uncomment then run code to generate system momentum vs. time plot
# momentumplot(
#     data,
#     "Weeks",
#     "All",
#     graphdir * "momentumexample.png",
#     showlegend = true,
#     dpi = 400
# )

# Uncomment then run code to generate system energy vs. time plot
# energyplot(
#     data,
#     "Weeks",
#     graphdir * "energyexample.png",
#     showlegend = true,
#     dpi = 400
# )

# Uncomment then run code to generate object position vs. time plot
# positionplot(
#     data,
#     ["x", "y"],
#     "AU",
#     "Years",
#     graphdir * "positionexample.png",
#     showlegend = true,
#     dpi = 400
# )

# Uncomment then run code to generate object velocity vs. time plot
# velocityplot(
#     data,
#     ["x", "y"],
#     "m/s",
#     "Years",
#     graphdir * "velocityexample.png",
#     showlegend = true,
#     dpi = 400
# )

# Uncomment then run code to generate object space plot
# staticspaceplot(
#     data,
#     ["x", "y", "z"],
#     "AU",
#     "Days",
#     graphdir * "3Dstaticexample.png",
#     showlegend = true,
#     squareplot = false,
#     dpi = 400,
#     angle3D = (30, 30)
# )

# Uncomment then run code to generate object space animation
animatedspaceplot(
    data,
    ["x", "y", "z"],
    "AU",
    "Days",
    graphdir * "3Danimatedexample.gif",
    showlegend = true,
    showpaths = true,
    squareplot = false,
    dpi = 400,
    displayevery = 200,
    fps = 30,
    angle3D = (30, 30)
)

# ---------------------------------------------------------------------------------------------------------------------