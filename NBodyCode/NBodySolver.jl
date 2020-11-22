# -----------------------------------------------------------------------------------------------------------------------------------------

#= Required Imports =#
using LinearAlgebra
include("NBodyConstants.jl")

# -----------------------------------------------------------------------------------------------------------------------------------------

#= Function Definitions =#

"""
    distance(x₁::Vector{Float64}, x₂::Vector{Float64})

distance : Vector{Float64}, Vector{Float64} -> Float64
    | Returns the Euclidean distance between two vectors.

x₁ : Vector{Float64}
    | The first of two vectors involved in the calculation.

x₂ : Vector{Float64}
    | The second of two vectors involved in the calculation.
"""
function distance(x₁::Vector{Float64}, x₂::Vector{Float64})
    d = norm(x₂ - x₁)
    return d
end


"""
    pairing(items::Vector)

pairing : Vector -> Vector{Vector}
    | Returns a vector of all unique pairs of elements in the input vector.

items : Vector
    | Vector of elements that will be sorted into unique pairs.
"""
function pairing(items::Vector)
    pairs = []
    for item1 ∈ items, item2 ∈ items
        if (item1 != item2) && ([item2, item1] ∉ pairs)
            push!(pairs, [item1, item2])
        end
    end
    return pairs
end


"""
    sysmomentum(bodies::Vector{Body})

sysmomentum : Vector{Body} -> Vector{Float64}
    | Computes and returns the total 3-momentum of a given system of Body structs.

bodies : Vector{Body}
    | The system of Body structs whose total 3-momentum will be computed.
"""
function sysmomentum(bodies::Vector{Body})
    p = [0.0, 0.0, 0.0]
    for body ∈ bodies
        m = body.mass; v = body.v 
        p += m*v
    end
    return p
end


"""
    syskinetic(bodies::Vector{Body})

syskinetic : Vector{Body} -> Float64
    | Computes and returns the total kinetic energy of a given system of Body 
    structs.

bodies : Vector{Body}
    | The system of Body structs whose total kinetic energy will be computed.
"""
function syskinetic(bodies::Vector{Body})
    T = 0.0
    for body ∈ bodies 
        m = body.mass; v = body.v 
        T += 1/2 * m * norm(v)^2
    end
    return T
end


"""
    syspotential(bodies::Vector{Body};
                 G::Float64 = 6.67408e-11, soften::Float64 = 0.0)

syspotential : Vector{Body}; Float64, Float64 -> Float64
    | Computes and returns the total gravitational potential energy of a given
    system of Body structs.

bodies : Vector{Body}
    | The system of Body structs whose total gravitational potential energy 
    will be computed.

G : Float64
    | Optional argument. Allows the strength of the gravitational constant to
    be changed. Set to 6.67408e-11 by default.
"""
function syspotential(bodies::Vector{Body};
                      G::Float64 = 6.67408e-11, soften::Float64 = 0.0)
    L = length(bodies)
    if L == 0
        V = 0.0
    else
        V = 0.0
        bodypairs = pairing(bodies)
        for pair ∈ bodypairs
            m₁ = pair[1].mass; x₁ = pair[1].x 
            m₂ = pair[2].mass; x₂ = pair[2].x
            r = distance(x₁, x₂)
            V += -(G*m₁*m₂)/(sqrt(r^2 + soften^2))
        end
    end
    return V
end


"""
    sysenergy(bodies::Vector{Body})

sysenergy : Vector{Body} -> Float64
    | Computes and returns the total energy of a given system of Body structs.

bodies : Vector{Body}
    | The system of Body structs whose total energy will be computed.
"""
function sysenergy(bodies::Vector{Body})
    E = syskinetic(bodies) + syspotential(bodies, G = 6.67408e-11, soften = 0.0)
    return E
end


"""
    acceleration(x₁::Vector{Float64}, M::Float64, x₂::Vector{Float64};
                 G::Float64 = 6.67408e-11, soften::Float64 = 0.0)

acceleration : Float64, Vector{Float64}, Vector{Float64}; Float64, Float64 -> Vector{Float64}
    | Computes and returns the acceleration of a body at x₁ due to the force
    from a body of mass M at x₂.

x₁ : Vector{Float64}
    | The position of the body whose acceleration is being computed
    Measured in m.

M : Float64
    | The mass of the body causing the body at x₁ to accelerate
    Measured in kg.

x₂ : Vector{Float64}
    | The position of the body causing the body at x₁ to accelerate.
    Measured in m.

G : Float64
    | Optional argument. Allows the strength of the gravitational constant to.
    be changed. Set to 6.67408e-11 by default.

soften : Float64
    | Optional argument. Softening parameter to reduce unrealistic acceleration.
    between close bodies. Set to 0.0 by default.
"""
function acceleration(x₁::Vector{Float64}, M::Float64, x₂::Vector{Float64};
                      G::Float64 = 6.67408e-11, soften::Float64 = 0.0)
    r = distance(x₁, x₂)
    a = (G*M)/(r^2 + soften^2) * (1/r) * (x₂ - x₁)
    return a
end


"""
    totalacceleration(centralx::Vector{Float64}, othermass::Vector{Float64}, otherx::Vector{Vector{Float64}})

totalacceleration : Vector{Float64}, Vector{Float64}, Vector{Vector{Float64}} -> Vector{Float64}
    | Computes and returns the total acceleration of a body at x₁ due to the
    forces of bodies with masses in othermass at positions in otherx.

centralx : Vector{Float64}
    | The position of the body whose acceleration is being computed.
    Measured in m.

othermass : Vector{Float64}
    | The masses of the other bodies causing the body at x₁ to accelerate.
    Measured in kg.

otherx : Vector{Vector{Float64}}
    | The positions of the other bodies causing the body at x₁ to accelerate.
    Measured in m.
"""
function totalacceleration(centralx::Vector{Float64}, othermass::Vector{Float64}, otherx::Vector{Vector{Float64}})
    # Assert statements
    @assert length(othermass) == length(otherx)

    # Constants
    L = length(othermass)
    x₁ = centralx

    # Function body
    a = [0.0, 0.0, 0.0]
    for i ∈ 1:L
        M = othermass[i]; x₂ = otherx[i]
        a += acceleration(x₁, M, x₂, G = 6.67408e-11, soften = 0.0)
    end
    return a
end


"""
    evolution!(bodies::Vector{Body}, method::String, stepcount::Int64)

evolution! : Vector{Body}, String, Int -> nothing
    | Computes the evolution in space and time over one time step for every Body 
    struct in bodies using the Position Extended Forest-Ruth Like algorithm. 
    Updates the x and v values of each Body struct accordingly.
        
bodies : Vector{Body}
    | The system of Body structs whose evolution in space and time is being
    computed.

method : String
    | Determines the method used to numerically solve the system of Newton's 
    laws. Possible values are "Euler", "Euler-Cromer", "Position Verlet",
    "Velocity Verlet", "Forest-Ruth", and "PEFRL". The recommended method is
    "PEFRL".

stepcount : Int
    | Determines whether a given Body struct in bodies is updated based on 
    that body's h (stepsize) value.
"""
function evolution!(bodies::Vector{Body}, method::String, stepcount::Int)
    # Assert statements
    @assert method ∈ ["Euler", "Euler-Cromer", "Position Verlet", "Velocity Verlet", "Forest-Ruth", "PEFRL"]
    
    # Constants
    L = length(bodies)

    # Function body
    if L == 1
        central = bodies[1]
        v₀ = central.v; h = central.h 
        central.x += v₀*h

    else
        update = [false for i ∈ 1:L]
        hmin = minimum([body.h for body ∈ bodies])
        updateon = [body.h/hmin for body ∈ bodies]
        xvals = [[] for i ∈ 1:L]; vvals = [[] for i ∈ 1:L]

        for i ∈ 1:L
            if stepcount % updateon[i] == 0
                # Constants
                central = bodies[i]
                h = central.h; x₀ = central.x; v₀ = central.v 
                othermass = [other.mass for other ∈ bodies if other != central]
                otherx = [other.x for other ∈ bodies if other != central]
                a(x) = totalacceleration(x, othermass, otherx)

                if method == "Euler"
                    x₁ = x₀ + (h*v₀)
                    v₁ = v₀ + (h*a(x₀))
                    xnew = x₁; vnew = v₁

                elseif method == "Euler-Cromer"
                    x₁ = x₀ + (h*v₀)
                    v₁ = v₀ + (h*a(x₁))
                    xnew = x₁; vnew = v₁

                elseif method == "Position Verlet"
                    x₁ = x₀ + (h/2 * v₀)
                    v₁ = v₀ + (h*a(x₁))
                    x₂ = x₁ + (h/2 * v₁)
                    xnew = x₂; vnew = v₁

                elseif method == "Velocity Verlet"
                    v₁ = v₀ + (h/2 * a(x₀))
                    x₁ = x₀ + (h*v₁)
                    v₂ = v₁ + (h/2 * a(x₁))
                    xnew = x₁; vnew = v₂

                elseif method == "Forest-Ruth"
                    ζ = 1/(2 - 2^(1/3))
                    x₁ = x₀ + (ζ * h/2 * v₀)
                    v₁ = v₀ + (ζ*h*a(x₁))
                    x₂ = x₁ + ((1 - ζ) * h/2 * v₁)
                    v₂ = v₁ + ((1 - 2*ζ) * h * a(x₂))
                    x₃ = x₂ + ((1 - ζ) * h/2 * v₂)
                    v₃ = v₂ + (ζ*h*a(x₃))
                    x₄ = x₃ + (ζ * h/2 * v₃)
                    xnew = x₄; vnew = v₃

                elseif method == "PEFRL"
                    α = 0.1786178958448091; β = -0.2123418310626054; γ = -0.06626458266981849
                    x₁ = x₀ + (α*h*v₀)
                    v₁ = v₀ + ((1 - 2*β) * h/2 * a(x₁))
                    x₂ = x₁ + (γ*h*v₁)
                    v₂ = v₁ + (β*h*a(x₂))
                    x₃ = x₂ + ((1 - 2*(γ + α)) * h * v₂)
                    v₃ = v₂ + (β*h*a(x₃))
                    x₄ = x₃ + (γ*h*v₃)
                    v₄ = v₃ + ((1 - 2*β) * h/2 * a(x₄))
                    x₅ = x₄ + (α*h*v₄)
                    xnew = x₅; vnew = v₄
                end
                
                xvals[i] = xnew; vvals[i] = vnew
                update[i] = true
            else
                update[i] = false
            end
        end

        for i ∈ 1:L
            if update[i]
                central = bodies[i]
                central.x = xvals[i]; central.v = vvals[i]
            end
        end
    end
end


"""
    runsimulation(bodies::Vector{Body}, method::String, T::Real)

runsimulation : Vector{Body}, String, Real - > Vector{Vector{String}, Vector{Float64}}
    | Computes the evolution in space and time evolution for every Body struct in
    bodies from t = 0 to t = T. Returns a vector containing the names, time
    values, position values, velocity values, momentum values, kinetic energy
    values, gravitational potential energy values, and total energy values for
    each body across the whole simulation.

bodies : Vector{Body}
    | The system of Body structs whose evolution in space and time is being
    computed.

method : String
    | Determines the method used to numerically solve the system of Newton's 
    laws. Possible values are "Euler", "Euler-Cromer", "Position Verlet",
    "Velocity Verlet", "Forest-Ruth", "PEFRL", "RK4", and "RKNFD". The
    recommended method is "PEFRL".

T : Real
    | The time over which the simulation is run.
    Measured in s.
"""
function runsimulation(bodies::Vector{Body}, method::String, T::Real)
    names = [body.name for body in bodies]
    tvals = [0.0]
    xvals = [[body.x for body ∈ bodies]]
    vvals = [[body.v for body ∈ bodies]]
    pvals = [sysmomentum(bodies)]
    Tvals = [syskinetic(bodies)]
    Vvals = [syspotential(bodies)]
    Evals = [sysenergy(bodies)]

    stepcount = 1
    hmin = minimum([body.h for body ∈ bodies])
    t = hmin

    while t ≤ T
        evolution!(bodies, method, stepcount)
        push!(tvals, t)
        push!(xvals, [body.x for body ∈ bodies])
        push!(vvals, [body.v for body ∈ bodies])
        push!(pvals, sysmomentum(bodies))
        push!(Tvals, syskinetic(bodies))
        push!(Vvals, syspotential(bodies))
        push!(Evals, sysenergy(bodies))

        stepcount += 1
        t += hmin
    end
    return [names, tvals, xvals, vvals, pvals, Tvals, Vvals, Evals]
end

# -----------------------------------------------------------------------------------------------------------------------------------------