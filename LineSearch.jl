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

# Include programs
include("backtrack.jl")

end