# -----------------------------------------------------------------------------------------------------------------------------------------

#= Required Imports =#
using Plots
using ColorSchemes

# -----------------------------------------------------------------------------------------------------------------------------------------

#= Function Definitions =#

"""
    timevalues(data::Vector, timescale::String)

timevalues : Vector, String -> Vector{Float64}
    | Converts the time values in s obtained from a simulation to a different
    time unit. Returns the new time values as a vector.

data : Vector
    | An output from runsimulation.

timescale : String 
    | The units to which the time values in s will be converted. Possible values
    are "Seconds", "Minutes", "Hours", "Days", "Weeks", "Months", and "Years".
"""
function timevalues(data::Vector, timescale::String)
    # Assert statements
    @assert timescale ∈ ["Seconds", "Minutes", "Hours", "Days", "Weeks", "Months", "Years"]

    # Constants
    minute = 60; hour = 60*minute; day = 24*hour
    week = 7*day; month = (365/12)*day; year = 365*day

    # Function body
    t = data[2]
    if timescale == "Minutes"
        t /= minute
    elseif timescale == "Hours"
        t /= hour
    elseif timescale == "Days"
        t /= day
    elseif timescale == "Weeks"
        t /= week
    elseif timescale == "Months"
        t /= month
    elseif timescale == "Years"
        t /= year
    end
    return t
end


"""
    momentumplot(data::Vector, timescale::String, component::String, saveas::String;
                 showlegend = true, dpi = 300)

momentumplot : Vector, String, String, String; Int -> Plot
    | Using data from a simulation, generates a plot of a given component of
    momentum against time. The plot is then saved locally.

data : Vector
    | An output from runsimulation.

timescale : String
    | The unit of time that will be used to display time on the plot. Possible
    values are "Seconds", "Minutes", "Hours", "Days", "Weeks", "Months", and
    "Years".

component : String
    | The component of the momentum that will be plotted against time. Possible
    values are "x", "y", "z", "Total" (norm of momentum vector), and "All"
    (displays "x", "y", "z", and "Total" in a 2x2 plot).

saveas : String
    | The absolute path to where the plot will be saved.

showlegend : Bool
    | Optional argument. Determines whether or not the legend is shown on the
    plot. Default value is true.

dpi : Int
    | Optional argument. Sets the dpi (resolution) of the generated plot. Default 
    value is 300.
"""
function momentumplot(data::Vector, timescale::String, component::String, saveas::String;
                      showlegend = true, dpi = 300)
    # Assert statements
    @assert component ∈ ["x", "y", "z", "Total", "All"]
    @assert saveas != ""

    # Constants
    t = timevalues(data, timescale)
    pvals = data[5]
    L = length(pvals)

    px = [pvals[i][1] for i ∈ 1:L]
    py = [pvals[i][2] for i ∈ 1:L]
    pz = [pvals[i][3] for i ∈ 1:L]
    P = [norm(p) for p ∈ pvals]

    x₀ = px[1]; y₀ = py[1]; z₀ = pz[1]; P₀ = P[1]

    Δpx = [x - x₀ for x ∈ px]
    Δpy = [y - y₀ for y ∈ py]
    Δpz = [z - z₀ for z ∈ pz]
    ΔP = [p - P₀ for p ∈ P]

    # Function body
    if component == "x"
        plt = plot()
        plot!(t, Δpx, legend = false, dpi = dpi)
        xlabel!("t ($(timescale))")
        ylabel!("p (kg m/s)")
        title!("Change in System x Momentum Over Time")

    elseif component == "y"
        plt = plot()
        plot!(t, Δpy, legend = false, dpi = dpi)
        xlabel!("t ($(timescale))")
        ylabel!("p (kg m/s)")
        title!("Change in System y Momentum Over Time")

    elseif component == "z"
        plt = plot()
        plot!(t, Δpz, legend = false, dpi = dpi)
        xlabel!("t ($(timescale))")
        ylabel!("p (kg m/s)")
        title!("Change in System z Momentum Over Time")

    elseif component == "Total"
        plt = plot()
        plot!(t, ΔP, legend = false, dpi = dpi)
        xlabel!("t ($(timescale))")
        ylabel!("p (kg m/s)")
        title!("Change in Total System Momentum Over Time")

    elseif component == "All"
        plotx = plot(t, Δpx, color = :green, label = "x", dpi = dpi)
        ploty = plot(t, Δpy, color = :blue, label = "y", dpi = dpi)
        plotz = plot(t, Δpz, color = :red, label = "z", dpi = dpi)
        plotP = plot(t, ΔP, color = :black, label = "Total", dpi = dpi)
        plt = plot(plotx, ploty, plotz, plotP, layout = (2, 2))
        xlabel!("t ($(timescale))")
        ylabel!("p (kg m/s)")
    end

    if !showlegend
        plot!(legend = false)
    end
    savefig(saveas)
end


"""
    energyplot(data::Vector, timescale::String, saveas::String;
               showlegend = true, dpi = 300)

energyplot : Vector, String, String; Int -> Plot
    | Using data from a simulation, generates a plot of a the kinetic, potential,
    and total energy against time. The plot is then saved locally.

data : Vector
    | An output from runsimulation.

timescale : String
    | The unit of time that will be used to display time on the plot. Possible
    values are "Seconds", "Minutes", "Hours", "Days", "Weeks", "Months", and
    "Years".

saveas : String
    | The absolute path to where the plot will be saved.

showlegend : Bool
    | Optional argument. Determines whether or not the legend is shown on the
    plot. Default value is true.

dpi : Int
    | Optional argument. Sets the dpi (resolution) of the generated plot. Default 
    value is 300.
"""
function energyplot(data::Vector, timescale::String, saveas::String;
                    showlegend = true, dpi = 300)
    # Assert statements
    @assert saveas != ""

    # Constants
    t = timevalues(data, timescale)
    Tvals = data[6]; Vvals = data[7]; Evals = data[8]
    T₀ = Tvals[1]; V₀ = Vvals[1]; E₀ = Evals[1]

    ΔT = [T - T₀ for T ∈ Tvals]
    ΔV = [V - V₀ for V ∈ Vvals]
    ΔE = [E - E₀ for E ∈ Evals]

    # Function body
    plt = plot()
    plot!(t, ΔT, color = :red, label = "Kinetic", legend = :best, dpi = dpi)
    plot!(t, ΔV, color = :blue, label = "Potential", dpi = dpi)
    plot!(t, ΔE, color = :black, label = "Total", dpi = dpi)
    xlabel!("t ($(timescale))")
    ylabel!("E (J)")
    title!("Change in System Energy Over Time")

    if !showlegend
        plot!(legend = false)
    end
    savefig(saveas)
end


"""
    posvelplot(data::Vector, posvel::Vector{String}, coordinate::String, timescale::String)

posvelplot : Vector, Vector{String}, String, String -> Plot
    | Using data from a simulation, generates a plot of one component of either
    position or velocity against time. This function is a reference for other
    functions and not meant to be directly used.

data : Vector
    | An output from runsimulation.

posvel : Vector{String}
    | A vector of strings of length 2. The first element determines whether the
    plot is of position or velocity; possible values are "Position" and 
    "Velocity". The second element determines the units that will be used to
    plot the first component; possible values are "m", "AU", and "LY" for 
    "Position", and "m/s", "km/s", and "km/h" for "Velocity".

coordinate : String
    | The component of the given quantity that will be plotted. Possible values
    are "x", "y", and "z".

timescale : String
    | The unit of time that will be used to display time on the plot. Possible
    values are "Seconds", "Minutes", "Hours", "Days", "Weeks", "Months", and
    "Years".
"""
function posvelplot(data::Vector, posvel::Vector{String}, coordinate::String, timescale::String)
    # Assert statements
    poscondition = (posvel[1] == "Position" && posvel[2] ∈ ["m", "AU", "LY"])
    velcondition = (posvel[1] == "Velocity" && posvel[2] ∈ ["m/s", "km/s", "km/h"])
    @assert poscondition || velcondition
    @assert coordinate ∈ ["x", "y", "z"]

    # Constants
    t = timevalues(data, timescale)
    bodies = data[1]
    L₁ = length(bodies)
    au = 149597870700; ly = 9454254955488000
    kms = 1000; kmh = 3.6

    # Function body
    if posvel[1] == "Position"
        valindex = 3
    elseif posvel[1] == "Velocity"
        valindex = 4
    end

    if coordinate == "x"
        dimindex = 1
    elseif coordinate == "y"
        dimindex = 2
    elseif coordinate == "z"
        dimindex = 3
    end

    L₂ = length(data[valindex])

    plot(legend = :best)
    for i ∈ 1:L₁
        posveldata = [data[valindex][j][i][dimindex] for j ∈ 1:L₂]
        if posvel[2] == "AU"
            posveldata /= au
        elseif posvel[2] == "LY"
            posveldata /= ly
        elseif posvel[2] == "km/s"
            posveldata /= kms
        elseif posvel[2] == "km/h"
            posveldata *= kmh
        end
        plot!(t, posveldata, label = bodies[i])
    end
    xlabel!("t ($(timescale))")
    ylabel!("$(coordinate) $(posvel[1]) ($(posvel[2]))")
end


"""
    positionplot(data::Vector, coordinates::Vector{String}, positionscale::String, timescale::String, saveas::String;
                 showlegend = true, dpi = 300)

positionplot : Vector, Vector{String}, String, String, String; Bool, Int -> Plot
    | Using data from a simulation, generates a plot of one or more components
    of the position values against time. The plot is then saved locally.

data : Vector
    An output from runsimulation.

coordinates : Vector{String}
    | An array of strings that determines which components of position will be
    plotted. Possible elements are "x", "y", and "z".

positionscale : String
    | The unit of position that will be used to display position on the plot.
    Possible values are "m", "AU", and "LY".

timescale : String
    | The unit of time that will be used to display time on the plot. Possible
    values are "Seconds", "Minutes", "Hours", "Days", "Weeks", "Months", and
    "Years".

saveas : String
    | The absolute path to where the plot will be saved.

showlegend : Bool
    | Optional argument. Determines whether or not the legend is shown on the
    plot. Default value is true.

dpi : Int
    | Optional argument. Sets the dpi (resolution) of the generated plot. Default
    value is 300.
"""
function positionplot(data::Vector, coordinates::Vector{String}, positionscale::String, timescale::String, saveas::String;
                      showlegend = true, dpi = 300)
    # Assert statements
    @assert 1 ≤ length(coordinates) ≤ 3
    @assert saveas != ""

    # Function body
    plotsarray = []
    if "x" ∈ coordinates
        xplot = posvelplot(data, ["Position", positionscale], "x", timescale)
        push!(plotsarray, xplot)
    end
    if "y" ∈ coordinates
        yplot = posvelplot(data, ["Position", positionscale], "y", timescale)
        push!(plotsarray, yplot)
    end
    if "z" ∈ coordinates
        zplot = posvelplot(data, ["Position", positionscale], "z", timescale)
        push!(plotsarray, zplot)
    end

    L = length(plotsarray)
    if L == 1
        plt = plot(plotsarray[1], dpi = dpi)
    elseif L == 2
        plt = plot(plotsarray[1], plotsarray[2], layout = (1, 2), dpi = dpi)
    elseif L == 3
        plt = plot(plotsarray[1], plotsarray[2], plotsarray[3], layout = (1, 3), dpi = dpi)
    end

    if !showlegend
        plot!(legend = false)
    end
    savefig(saveas)
end


"""
    velocityplot(data::Vector, coordinates::Vector{String}, velocityscale::String, timescale::String, saveas::String;
                 showlegend = true, dpi = 300)

velocityplot : Vector, Vector{String}, String, String, String; Bool, Int -> Plot
    | Using data from a simulation, generates a plot of one or more components
    of the velocity values against time. The plot is then saved locally.

data : Vector
    | An output from runsimulation.

coordinates : Vector{String}
    | An array of strings that determines which components of velocity will be
    plotted. Possible elements are "x", "y", and "z".

velocityscale : String
    | The unit of velocity that will be used to display velocity on the plot.
    Possible values are "m/s", "km/s", and "km/h".

timescale : String
    | The unit of time that will be used to display time on the plot. Possible
    values are "Seconds", "Minutes", "Hours", "Days", "Weeks", "Months", and
    "Years".

saveas : String
    | The absolute path to where the plot will be saved.

showlegend : Bool
    | Optional argument. Determines whether or not the legend is shown on the
    plot. Default value is true.

dpi : Int
    | Optional argument. Sets the dpi (resolution) of the generated plot. Default
    value is 300.
"""
function velocityplot(data::Vector, coordinates::Vector{String}, velocityscale::String, timescale::String, saveas::String;
                      showlegend = true, dpi = 300)
    # Assert statements
    @assert 1 ≤ length(coordinates) ≤ 3
    @assert saveas != ""

    # Function body
    plotsarray = []
    if "x" ∈ coordinates
        xplot = posvelplot(data, ["Velocity", velocityscale], "x", timescale)
        push!(plotsarray, xplot)
    end
    if "y" ∈ coordinates
        yplot = posvelplot(data, ["Velocity", velocityscale], "y", timescale)
        push!(plotsarray, yplot)
    end
    if "z" ∈ coordinates
        zplot = posvelplot(data, ["Velocity", velocityscale], "z", timescale)
        push!(plotsarray, zplot)
    end

    L = length(plotsarray)
    if L == 1
        plt = plot(plotsarray[1], dpi = dpi)
    elseif L == 2
        plt = plot(plotsarray[1], plotsarray[2], layout = (1, 2), dpi = dpi)
    elseif L == 3
        plt = plot(plotsarray[1], plotsarray[2], plotsarray[3], layout = (1, 3), dpi = dpi)
    end

    if !showlegend
        plot!(legend = false)
    end
    savefig(saveas)
end


"""
    staticspaceplot(data::Vector, coordinates::Vector{String}, positionscale::String, timescale::String, saveas::String;
                    showlegend = true, squareplot = false, dpi = 300, angle3D = (30, 30))

staticspaceplot : Vector, Vector{String}, String, String; Bool, Int, Tuple{Int} -> Plot
    | Using data from a simulation, generates a plot of the paths of all of the
    bodies in the simulation through space in 2D or 3D. The plot is then saved
    locally.

data : Vector
    | An output from runsimulation.

coordinates : Vector{String}
    | An array of strings that determines which components of velocity will be
    plotted. Possible values are (permutations of) ["x", "y"], ["x", "z"],
    ["y", "z"], and ["x", "y", "z"]. The first element in the list is plotted
    on the x-axis, the second on the y-axis, and the third on the z-axis.

positionscale : String
    | The unit of position that will be used to display position on the plot.
    Possible values are "m", "AU", and "LY".

timescale : String
    | The unit of time that will be used in the plot title. Possible values are
    "Seconds", "Minutes", "Hours", "Days", "Weeks", "Months", and "Years".

saveas : String
    | The absolute path to where the plot will be saved.

showlegend : Bool
    | Optional argument. Determines whether or not the legend is shown on the
    plot. Default value is true.

squareplot : Bool
    | Optional argument. Determines whether all the axes have the same ranges.
    Default value is false.

dpi : Int
    | Optional argument. Sets the dpi (resolution) of the generated plot. Default
    value is 300.

angle3D : Tuple{Int}
    | Optional argument. For a 3D plot, adjusts the angle at which the plot is
    generated. Default value is (30, 30).
"""
function staticspaceplot(data::Vector, coordinates::Vector{String}, positionscale::String, timescale::String, saveas::String;
                         showlegend = true, squareplot = false, dpi = 300, angle3D = (30, 30))
    # Assert statements
    @assert length(coordinates) ∈ [2, 3]
    for coord ∈ coordinates
        @assert coord ∈ ["x", "y", "z"]
    end
    @assert positionscale ∈ ["m", "AU", "LY"]
    @assert saveas != ""

    # Constants
    possiblecoords = ["x", "y", "z"]
    au = 149597870700; ly = 9454254955488000
    bodies = data[1]; t = timevalues(data, timescale); xvals = data[3]
    L₁ = length(bodies); L₂ = length(t)

    if positionscale == "AU"
        xvals /= au
    elseif positionscale == "LY"
        xvals /= ly
    end

    if squareplot
        xmax = 0.0; ymax = 0.0; zmax = 0.0
        for i ∈ 1:L₂
            for j ∈ 1:L₁
                if abs(xvals[i][j][1]) > xmax
                    xmax = abs(xvals[i][j][1])
                end
                if abs(xvals[i][j][2]) > ymax
                    ymax = abs(xvals[i][j][2])
                end
                if abs(xvals[i][j][3]) > zmax
                    zmax = abs(xvals[i][j][3])
                end
            end
        end
        xmax *= 1.1; ymax *= 1.1; zmax *= 1.1

        plotmax = max(xmax, ymax, zmax)
    end

    index1 = findfirst(isequal(coordinates[1]), possiblecoords)
    index2 = findfirst(isequal(coordinates[2]), possiblecoords)

    # Function body
    plt = plot(title = "System Evolution Over $(maximum(t)) $(timescale)")
    if !showlegend
        plot!(legend = false)
    end

    if length(coordinates) == 2
        for i ∈ 1:L₁
            axis1data = [xvals[j][i][index1] for j ∈ 1:L₂]
            axis2data = [xvals[j][i][index2] for j ∈ 1:L₂]
            plotcolour = get(ColorSchemes.rainbow1, i./L₁)
            plot!(axis1data, axis2data, color = plotcolour, label = "", dpi = dpi)
            scatter!([axis1data[end]], [axis2data[end]], color = plotcolour, label = bodies[i], dpi = dpi)
            xlabel!("$(coordinates[1]) ($(positionscale))")
            ylabel!("$(coordinates[2]) ($(positionscale))")
        end

    elseif length(coordinates) == 3
        index3 = findfirst(isequal(coordinates[3]), possiblecoords)
        for i ∈ 1:L₁
            axis1data = [xvals[j][i][index1] for j ∈ 1:L₂]
            axis2data = [xvals[j][i][index2] for j ∈ 1:L₂]
            axis3data = [xvals[j][i][index3] for j ∈ 1:L₂]
            plotcolour = get(ColorSchemes.rainbow1, i./L₁)
            plot!(axis1data, axis2data, axis3data, color = plotcolour, label = "", zlabel = "$(coordinates[3]) ($(positionscale))", dpi = dpi)
            scatter!([axis1data[end]], [axis2data[end]], [axis3data[end]], color = plotcolour, label = bodies[i], dpi = dpi)
            xlabel!("$(coordinates[1]) ($(positionscale))")
            ylabel!("$(coordinates[2]) ($(positionscale))")
        end
        if squareplot
            zlims!(-plotmax, plotmax)
        end
        plot!(camera = angle3D)
    end
    if squareplot
        xlims!(-plotmax, plotmax)
        ylims!(-plotmax, plotmax)
    end
    savefig(saveas)
end


"""
    animatedspaceplot(data::Vector, coordinates::Vector{String}, positionscale::String, timescale::String, saveas::String;
                      showlegend = true, showpaths = false, squareplot = false, dpi = 300, displayevery = 50, fps = 30, angle3D = (30, 30))

animatedspaceplot : Vector, Vector{String}, String, String, String; Bool, Int, Int, Int, Tuple{Int} -> Plot
    | Using data from a simulation, generates a gif that shows how the positions
    of all the bodies in the simulation move through space in 2D or 3D. The 
    plot is then saved locally.

data : Vector
    | An output from runsimulation.

coordinates : Vector{String}
    | An array of strings that determines which components of velocity will be
    plotted. Possible values are (permutations of) ["x", "y"], ["x", "z"],
    ["y", "z"], and ["x", "y", "z"]. The first element in the list is plotted
    on the x-axis, the second on the y-axis, and the third on the z-axis.

positionscale : String
    | The unit of position that will be used to display position on the plot.
    Possible values are "m", "AU", and "LY".

timescale : String
    | The unit of time that will be used to display time on the plot. Possible
    values are "Seconds", "Minutes", "Hours", "Days", "Weeks", "Months", and
    "Years". Currently only works for 2D plots.

saveas : String
    | The absolute path to where the plot will be saved.

showlegend : Bool
    | Optional argument. Determines whether or not the legend is shown on the
    plot. Default value is true.

showpaths : Bool
    | Optional argument. Determines whether the bodies leave a trail through
    space. Default value is false.

squareplot : Bool
    | Optional argument. Determines whether all the axes have the same ranges.
    Default value is false.

dpi : Int
    | Optional argument. Sets the dpi (resolution) of the generated plot. Default
    value is 300.

displayevery : Int
    | Optional argument. Sets the frequency at which a frame is added to the gif.
    Default value is 10.

fps : Int
    | Optional argument. Sets the frames-per-second value of the gif. Default
    value is 30.

angle3D : Tuple{Int}
    | Optional argument. For a 3D plot, adjusts the angle at which the plot is
    generated. Default value is (30, 30).
"""
function animatedspaceplot(data::Vector, coordinates::Vector{String}, positionscale::String, timescale::String, saveas::String;
                           showlegend = true, showpaths = false, squareplot = false, dpi = 300, displayevery = 50, fps = 30, angle3D = (30, 30))
    # Assert statements
    @assert length(coordinates) ∈ [2, 3]
    for i ∈ coordinates
        @assert i ∈ ["x", "y", "z"]
    end
    @assert positionscale ∈ ["m", "AU", "LY"]
    @assert saveas != ""

    # Constants
    possiblecoords = ["x", "y", "z"]
    au = 149597870700; ly = 9454254955488000
    bodies = data[1]; t = timevalues(data, timescale); xvals = data[3]
    L₁ = length(bodies); L₂ = length(t)

    if positionscale == "AU"
        xvals /= au
    elseif positionscale == "LY"
        xvals /= ly
    end

    xmax = 0.0; ymax = 0.0; zmax = 0.0
    for i ∈ 1:L₂
        for j ∈ 1:L₁
            if abs(xvals[i][j][1]) > xmax
                xmax = abs(xvals[i][j][1])
            end
            if abs(xvals[i][j][2]) > ymax
                ymax = abs(xvals[i][j][2])
            end
            if abs(xvals[i][j][3]) > zmax
                zmax = abs(xvals[i][j][3])
            end
        end
    end
    xmax *= 1.1; ymax *= 1.1; zmax *= 1.1

    if squareplot
        xmax = max(xmax, ymax, zmax)
        ymax = max(xmax, ymax, zmax)
        zmax = max(xmax, ymax, zmax)
    end

    axislims = [xmax, ymax, zmax]

    index1 = findfirst(isequal(coordinates[1]), possiblecoords)
    index2 = findfirst(isequal(coordinates[2]), possiblecoords)

    # Function body
    if length(coordinates) == 2
        animation = @animate for i ∈ 1:L₂
            plot()
            for j ∈ 1:L₁
                axis1vals = [xvals[k][j][index1] for k ∈ 1:i]; axis1point = [axis1vals[end]]
                axis2vals = [xvals[k][j][index2] for k ∈ 1:i]; axis2point = [axis2vals[end]]
                plotcolour = get(ColorSchemes.rainbow1, j./L₁)
                scatter!(axis1point, axis2point, color = plotcolour, label = bodies[j], dpi = dpi)
                if showpaths
                    plot!(axis1vals, axis2vals, color = plotcolour, label = "", dpi = dpi)
                end
                xlabel!("$(coordinates[1]) ($(positionscale))")
                ylabel!("$(coordinates[2]) ($(positionscale))")
                xlims!(-axislims[index1], axislims[index1])
                ylims!(-axislims[index2], axislims[index2])
            end
            if !showlegend
                plot!(legend = false)
            end
            title!("System Evolution Over $(round(t[i], digits = 1)) $(timescale)")
        end every displayevery

    elseif length(coordinates) == 3
        index3 = findfirst(isequal(coordinates[3]), possiblecoords)
        animation = @animate for i ∈ 1:L₂
            plot()
            for j ∈ 1:L₁
                axis1vals = [xvals[k][j][index1] for k ∈ 1:i]; axis1point = [axis1vals[end]]
                axis2vals = [xvals[k][j][index2] for k ∈ 1:i]; axis2point = [axis2vals[end]]
                axis3vals = [xvals[k][j][index3] for k ∈ 1:i]; axis3point = [axis3vals[end]]
                plotcolour = get(ColorSchemes.rainbow1, j./L₁)
                scatter!(axis1point, axis2point, axis3point, color = plotcolour, label = bodies[j], zlabel = "$(coordinates[3]) ($(positionscale))", dpi = dpi)
                if showpaths
                    plot!(axis1vals, axis2vals, axis3vals, color = plotcolour, label = "", dpi = dpi)
                end
                xlabel!("$(coordinates[1]) ($(positionscale))")
                ylabel!("$(coordinates[2]) ($(positionscale))")
                xlims!(-axislims[index1], axislims[index1])
                ylims!(-axislims[index2], axislims[index2])
                zlims!(-axislims[index3], axislims[index3])
            end
            if !showlegend
                plot!(legend = false)
            end
            plot!(camera = angle3D)
            title!("System Evolution Over $(round(t[i], digits = 1)) $(timescale)")
        end every displayevery
    end
    gif(animation, saveas, fps = fps)
end

# -----------------------------------------------------------------------------------------------------------------------------------------