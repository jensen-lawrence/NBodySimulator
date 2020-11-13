# NBodySimulator
This is an n-body simulator written in Julia. The solver uses the Position Extended Forest-Ruth Like algorithm to solve a system of Newton's equations and then update the positions and velocities of each of the bodies in the system at each timestep.

## NBodyConstants.jl
This file contains many relevant physical constants and reference scales, along with the Body mutable struct, and some pre-defined object groups
#### Physical Constants
The available physical constants are: ```M_SUN```, ```R_SUN```, ```M_JUPITER```, ```R_JUPITER```, ```M_EARTH```, and ```R_EARTH```. These are primarily useful when defining new bodies, as astronomical literature often lists the masses and radii of stars and planets in terms of the masses and radii of the Sun, Jupiter, and Earth.
#### Reference Scales
The available reference scales can be broken down into time scales, position scales, and velocity scales. Time is defined in seconds in the simulator, so the available time scales are ```MINUTE```, ```HOUR```, ```DAY```, ```WEEK```, ```MONTH```, and ```YEAR```. These are useful for defining the a body's step size and for defining the period of time over which the simulation takes place. Position is defined in metres in the simulator, so the available position scales are ```AU``` and ```LY```. These are useful for defining a body's initial position. Velocity is defined in m/s in the simulator, so the available velocity scales are ```KMS``` and ```KMH```. These are useful for defining a body's initial velocity.
#### The Body Mutable Struct
The ```Body``` mutable struct is the core object of the simulation. Any instance of the ```Body``` struct is defined as follows:

```julia
body = Body(name, mass, radius, x, v, h)
```

```name``` is a string representing the name of the body, ```mass``` is a float representing the mass of the body in kg, ```radius``` is a float representing the radius of the body in m, ```x``` is a three-element vector of floats representing the position of the body in m, ```v``` is a three-element vector of floats representing the velocity of the body in m/s, and ```h``` is the step size of the body in s.

To run a simulation on a system of interacting bodies, the simulator must be provided with a vector that contains all of the relevant bodies. Each body in the system should be instantiated with the desired initial position and initial velocity. Bodies in the same system can have different time steps (useful if one is simulating the orbits of Mercury and Neptune, for example), but a given time step in a system must be an integer multiple of the system's smallest time step.
#### Pre-defined Groups
There are a number of pre-defined groups of interacting bodies: ```EarthMoon```, ```GalileanMoons```, ```SaturnMoons```, ```SolarSystem```, ```InnerSolarSystem```, ```Trappist1```, and ```Stars```. More are likely to be added in the future. The user can easily define their own groups of bodies simply defining the desired instances of the ```Body``` struct, and then grouping them all together in a vector.

## NBodySolver.jl
This file contains all of the functions that contribute to solving the n-body problem. Included are functions capable of computing momentum, kinetic energy, gravitational potential energy, total energy, and acceleration, along with the function that numerically solves the n-body problem, and the function that iterates the solution function over time.
#### Numerically Solving the N-Body Problem
Suppose we want to solve the n-body problem numericall from ```t = 0``` to ```t = T```. To do so, we determine the force (given by Newton's law of gravitation) on each body due to every other body at a given time, update the positions and velocities of each body based on the force they experience, and then increase ```t``` by the value of the smallest time step. In pseudocode, this reads as
```
t = 0
while t <= T do
  for body in bodies do
    F = get_force(body)
    body_position = new_position(F)
    body_velocity = new_velocity(F)
  end
  t += smallest_h
 end
```
The results are typically saved somewhere as the process iterates over values of ```t```. This simulator uses the Position Extended Forest-Ruth Like algorithm. More information on the algorithm can be found [here](https://arxiv.org/pdf/cond-mat/0110585.pdf).
#### Running a Simulation
To run an n-body simulation, first define a vector of ```Body``` structs (called ```bodies``` for this example). Then, define the time over which the simulation is run (called ```T``` for this example). Then, execute
```julia
data = runsimulation(bodies, T)
```
All of the relevant results from the simulation will be stored in ```data```, which can then be used for graphing (detailed in next section).
#### Advanced Usage
The functions ```syspotential``` and ```acceleration``` are defined with the optional arguments ```G``` and ```soften```. The functions ```sysenergy``` and ```totalacceleration``` call ```syspotential``` and ```acceleration```, respectively. By default, ```G = 6.67408e-11``` and ```soften = 0.0```. However, these values can be changed in the definitions of ```sysenergy``` and ```totalacceleration```, resulting in a simulation that follows different physical laws. Changing ```G``` will change the strength of gravity, and changing ```soften``` will reduce the strength of gravity. The ```soften``` parameter is typically changed to reduce unrealistic (approaching infinite) accelerations caused by integration errors when bodies get very close to each other. The value of ```soften``` must be chosen carefully for the given scenario so that the unrealistic accelerations are reduced but the simulation results are not rendered wholly inaccurate.

## NBodyGraphing.jl
This file contains functions that graph all of the noteworthy results from a simulation.
#### Graphing Change in Momentum vs. Time
Given ```data``` (as defined in Running a Simulation), a plot of change in momentum vs. time can be produced as follows:
```julia
momentumplot(data, timescale, component, saveas, showlegend = true, dpi = 300)
```
The ```timescale``` argument is a string chosen from ```Seconds```, ```Minutes```, ```Hours```, ```Days```, ```Weeks```, ```Months```, and ```Years```, and scales the time values displayed on the graph appropriately. The ```component``` argument is a string chosen from ```x```, ```y```, ```z```, ```Total```, and ```All```, and determines the which component of change in momentum will be displayed on the graph. The first three options have obvious effects, ```Total``` plots the change in the norm of momentum, and ```All``` produces a 2x2 plot of all four previous options. The ```saveas``` argument is the absolute path to where the graph will be saved. The ```showlegend``` argument is an optional Boolean argument, set to ```true``` by default, and determines if a legend is displayed on the graph. The ```dpi``` argument is an optional integer argument, set to ```300``` by default, and determines the dpi (resolution of the resulting graph).
#### Graphing Change in Energy vs. Time
Given ```data``` (as defined in Running a Simulation), a plot of change in energy vs. time can be produced as follows:
```julia
energyplot(data, timescale, saveas, showlegend = true, dpi = 300)
```
The ```timescale``` argument is a string chosen from ```Seconds```, ```Minutes```, ```Hours```, ```Days```, ```Weeks```, ```Months```, and ```Years```, and scales the time values displayed on the graph appropriately. The ```saveas``` argument is the absolute path to where the graph will be saved. The ```showlegend``` argument is an optional Boolean argument, set to ```true``` by default, and determines if a legend is displayed on the graph. The ```dpi``` argument is an optional integer argument, set to ```300``` by default, and determines the dpi (resolution of the resulting graph).
#### Graphing Position or Velocity vs. Time
Given ```data``` (as defined in Running a Simulation), a plot of position vs. time can be produced as follows:
```julia
positionplot(data, coordinates, positionscale, timescale, saveas, showlegend = true, dpi = 300)
```
Similarly, a plot of velocity vs. time can be produced as follows:
```julia
velocityplot(data, coordinates, velocityscale, timescale, saveas, showlegend = true, dpi = 300)
```
The ```coordinates``` argument is an array of strings chosen from (or any permutation of) ```["x"]```, ```["y"]```, ```["z"]```, ```["x", "y"]```, ```["x", "z"]```, ```["y", "z"]```, and ```["x", "y", "z"]```. The plots of the position/velocity value associated with each coordinate vs. time will appear left to right in the order presented in ```coordinates```. The ```positionscale``` argument is a string chosen from ```m```, ```AU```, and ```LY```; the ```velocityscale``` argument is a string chosen from ```m/s```, ```km/s```, and ```km/h```. They scale the position/velocity values displayed on the graph appropriately. The ```timescale``` argument is a string chosen from ```Seconds```, ```Minutes```, ```Hours```, ```Days```, ```Weeks```, ```Months```, and ```Years```, and scales the time values displayed on the graph appropriately. The ```saveas``` argument is the absolute path to where the graph will be saved. The ```showlegend``` argument is an optional Boolean argument, set to ```true``` by default, and determines if a legend is displayed on the graph. The ```dpi``` argument is an optional integer argument, set to ```300``` by default, and determines the dpi (resolution of the resulting graph).
#### Graphing Trajectories Through Space
Given ```data``` (as defined in Running a Simulation), a plot of the trajectories through space can be produced as follows:
```julia
staticspaceplot(data, coordinates, positionscale, timescale, saveas, showlegend = true, squareplot = false, dpi = 300, angle3D = (30, 30))
```
The ```coordinates``` argument is an array of string chosen from (or any permutation of) ```["x", "y"]```, ```["x", "z"]```, ```["y", "z"]```, and ```["x", "y", "z"]```. A two-element array results in a 2D plot; a three-element array results in a 3D plot. The x, y, and z axes of the graph plot the first, second, and third coordinates in the array, respectively. The ```positionscale``` argument is a string chosen from ```m```, ```AU```, and ```LY```, and scales the values displayed on each axis appropriately. The ```timescale``` argument is a string chosen from ```Seconds```, ```Minutes```, ```Hours```, ```Days```, ```Weeks```, ```Months```, and ```Years```, and scales the time displayed in the plot title appropriately. The ```saveas``` argument is the absolute path to where the graph will be saved. The ```showlegend``` argument is an optional Boolean argument, set to ```true``` by default, and determines if a legend is displayed on the graph. The ```squareplot``` argument is an optional Boolean argument, set to ```false``` by default, and determines if each axis has the same minimum and maximum values. The ```dpi``` argument is an optional integer argument, set to ```300``` by default, and determines the dpi (resolution of the resulting graph). The ```angle3D``` is an optional tuple argument, set to ```(30, 30)``` by default, and determines the (horizontal, vertical) angle at which the plot is generated if the plot is 3D.
#### Animating Trajectories Through Space
Given ```data``` (as defined in Running a Simulation), a plot of the trajectories through space can be produced as follows:
```julia
animatedspaceplot(data, coordinates, positionscale, timescale, saveas, showlegend = true, showpaths = false, squareplot = false, dpi = 300, displayevery = 50, fps = 30, angle3D = (30, 30))
```
The ```coordinates``` argument is an array of string chosen from (or any permutation of) ```["x", "y"]```, ```["x", "z"]```, ```["y", "z"]```, and ```["x", "y", "z"]```. A two-element array results in a 2D plot; a three-element array results in a 3D plot. The x, y, and z axes of the graph plot the first, second, and third coordinates in the array, respectively. The ```positionscale``` argument is a string chosen from ```m```, ```AU```, and ```LY```, and scales the values displayed on each axis appropriately. The ```timescale``` argument is a string chosen from ```Seconds```, ```Minutes```, ```Hours```, ```Days```, ```Weeks```, ```Months```, and ```Years```, and scales the time displayed in the plot title appropriately. The ```saveas``` argument is the absolute path to where the graph will be saved. The ```showlegend``` argument is an optional Boolean argument, set to ```true``` by default, and determines if a legend is displayed on the graph. The ```squareplot``` argument is an optional Boolean argument, set to ```false``` by default, and determines if each axis has the same minimum and maximum values. The ```dpi``` argument is an optional integer argument, set to ```300``` by default, and determines the dpi (resolution of the resulting graph). The ```displayevery``` argument is an optional integer argument, set to 50 by default, and determines the number of time steps between plots being generated and added to the animation. The ```fps``` argument is an optional integer argument, set to 30 by default, and determines the frames-per-second value of the animation. The ```angle3D``` is an optional tuple argument, set to ```(30, 30)``` by default, and determines the (horizontal, vertical) angle at which the plot is generated if the plot is 3D.
**Note:** generating animations with high ```dpi``` and low ```displayevery``` values can take a very long time.
#### Examples
For examples of each graphing output, see the ```ExampleGraphs``` folder.

## Errors
If you spot any errors in the code or find any bugs while using the code, please contact me at jptlawre@uwaterloo.ca
