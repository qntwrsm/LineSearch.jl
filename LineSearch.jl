#=
LineSearch.jl

    Provides line search methods, such as backtracking, to determine the optimal 
    step size in optimization and root-finding algorithms.

@author: Quint Wiersma <q.wiersma@vu.nl>

@date: 2022/02/02
=#

module LineSearch

using LinearAlgebra

export
    # Structs
    BackTrack, 
    # Backtracking
    backtrack!

# Structs
Base.@kwdef mutable struct BackTrack{Tf, Gf, Ti}
    α::Tf	            # stepsize
	α_init::Tf 	        # initial stepsize
    ρ::Gf = .5		    # contraction factor
    c_1::Tf = 1e-4      # sufficient decrease constant
    order::Ti = 0       # interpolation order
    ρ_hi::Gf = .5       # safeguard upper bound
    ρ_lo::Gf = .1       # safeguard lower bound
    max_iter::Ti = 1000 # max iterations
end

# Include programs
include("backtrack.jl")

end