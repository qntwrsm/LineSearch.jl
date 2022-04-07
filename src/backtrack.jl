#=
backtrack.jl

    Backtracking line search algorithms based on the Armijo or sufficient 
    decrease condition to determine the optimal step length, using either no, 
    quadratic, or cubic safeguarded interpolation.

@author: Quint Wiersma <q.wiersma@vu.nl>

@date: 2022/02/02
=#

# Backtracking struct
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

"""
    fixed!(ls, x, x_prev, p, f_prev, ∇f_prev, f)

Backtracking line search wtih a fixed contraction parameter ``ρ`` to find ``α``,
storing the result in `x`.

#### Arguments
  - `ls::BackTrack`             : line search parameters
  - `x_prev::AbstractVector`    : previous state
  - `p::AbstractVector`         : search direction
  - `f_prev::Real`              : previous objective function value
  - `∇f_prev::AbstractVector`   : previous gradient
  - `f::Function`               : ``f(x)``

#### Returns
  - `x::AbstractVector` : updated state
  - `f_new::Real`       : objective function value
"""
function fixed!(ls::BackTrack, x::AbstractVector, x_prev::AbstractVector, 
                p::AbstractVector, f_prev::Real, ∇f_prev::AbstractVector, 
                f::Function)
    # Unpack
    α= ls.α_init
    c_1= ls.c_1
    ρ= ls.ρ
    max_iter= ls.max_iter
    
    # Objective function
    x.= x_prev .+ α .* p
    f_new= f(x)

    # Directional derivative
    d= dot(∇f_prev, p)
    
    # Iteration counter
    iter= 1
    # Backtracking
    while f_new > f_prev + c_1 * α * d && iter < max_iter
        # Contract
        α*= ρ

        # Update objective function
        x.= x_prev .+ α .* p
        f_new= f(x)

        # Update iteration counter
        iter+= 1
    end

    # Store stepsize
    ls.α= α

    return f_new
end

"""
    quad_interp!(ls, x, x_prev, p, f_prev, ∇f_prev, f)

Backtracking line search wtih  quadratic interpolation to find ``α``, storing
the result in `x`.

#### Arguments
  - `ls::BackTrack`             : line search parameters
  - `x_prev::AbstractVector`    : previous state
  - `p::AbstractVector`         : search direction
  - `f_prev::Real`              : previous objective function value
  - `∇f_prev::AbstractVector`   : previous gradient
  - `f::Function`               : ``f(x)``

#### Returns
  - `x::AbstractVector` : updated state
  - `f_new::Real`       : objective function value
"""
function quad_interp!(ls::BackTrack, x::AbstractVector, x_prev::AbstractVector, 
                        p::AbstractVector, f_prev::Real, ∇f_prev::AbstractVector, 
                        f::Function)
    # Unpack
    α_1= ls.α_init
    c_1= ls.c_1
    ρ_hi, ρ_lo= ls.ρ_hi, ls.ρ_lo
    max_iter= ls.max_iter

    # Objective function
    x.= x_prev .+ α_1 .* p
    f_new= f(x)

    # Directional derivative
    d= dot(∇f_prev, p)

    # Iteration counter
    iter= 1
    # Backtracking
    while f_new > f_prev + c_1 * α_1 * d && iter < max_iter
        # Quadratic interpolation using f(x), ∇f(x)'p, and f(x + α_1 * p)
        α_tmp= -d * α_1^2 * inv(2 * (f_new - f_prev - d * α_1))

        # Safeguard
        α_tmp= min(α_tmp, ρ_hi * α_1)   # avoid too small reductions
        α_1= max(α_tmp, ρ_lo * α_1) # avoid too big reductions

        # Update objective function
        x.= x_prev .+ α_1 .* p
        f_new= f(x)

        # Update iteration counter
        iter+= 1
    end

    # Store stepsize
    ls.α= α_1
    
    return f_new
end

"""
    cubic_interp!(ls, x, x_prev, p, f_prev, ∇f_prev, f)

Backtracking line search wtih cubic interpolation to find ``α``, storing the
result in `x`.

#### Arguments
  - `ls::BackTrack`             : line search parameters
  - `x_prev::AbstractVector`    : previous state
  - `p::AbstractVector`         : search direction
  - `f_prev::Real`              : previous objective function value
  - `∇f_prev::AbstractVector`   : previous gradient
  - `f::Function`               : ``f(x)``

#### Returns
  - `x::AbstractVector` : updated state
  - `f_new::Real`       : objective function value
"""
function cubic_interp!(ls::BackTrack, x::AbstractVector, x_prev::AbstractVector, 
                        p::AbstractVector, f_prev::Real, ∇f_prev::AbstractVector, 
                        f::Function)
    # Infer type
    Tα= typeof(ls.α_init)

    # Unpack
    α_2, α_1= ls.α_init, ls.α_init
    c_1= ls.c_1
    ρ_hi, ρ_lo= ls.ρ_hi, ls.ρ_lo
    max_iter= ls.max_iter

    # Objective function values
    x.= x_prev .+ α_2 .* p
    f_new= f(x)
    f_tmp= f_new

    # Directional derivative
    d= dot(∇f_prev, p)

    # Iteration counter
    iter= 1
    # Backtracking
    while f_new > f_prev + c_1 * α_2 * d && iter < max_iter
        if iter == 1
            # Quadratic interpolation using f(x), ∇f(x)'p, and f(x + α_2 * p)
            α_tmp= -d * α_2^2 * inv(2 * (f_new - f_prev - d * α_2))
        else
            # Cubic interpolation using f(x), ∇f(x)'p, f(x + α_1 * p), and f(x + α_2 * p)
            div= inv(α_1^2 * α_2^2 * (α_2 - α_1))
            tmp_2= f_new - f_prev - d * α_2
            tmp_1= f_tmp - f_prev - d * α_1
            a= (α_1^2 * tmp_2  - α_2^2 * tmp_1) * div
            b= (-α_1^3 * tmp_2  + α_2^3 * tmp_1) * div

            if isapprox(a, zero(a), atol=eps(Tα))
                # approximate quadratic interpolation
                α_tmp= -d * inv(2 * b)
            else
                # discriminant
                D= max(b^2 - 3 * a * d, Tα(0))
                # quadratic equation root
                α_tmp= (-b + sqrt(D)) * inv(3 * a)
            end
        end

        # Store stepsize
        α_1= α_2

        # Safeguard
        α_tmp= min(α_tmp, ρ_hi * α_2)   # avoid too small reductions
        α_2= max(α_tmp, ρ_lo * α_2) # avoid too big reductions

        # Update objective function values
        x.= x_prev .+ α_2 .* p
        f_tmp, f_new= f_new, f(x)

        # Update iteration counter
        iter+= 1
    end

    # Store stepsize
    ls.α= α_2

    return f_new
end

"""
    backtrack!(ls, x, x_prev, p, f_prev, ∇f_prev, f)

Backtracking line search to find optimal step size ``α``, storing the result in
`x`.

#### Arguments
  - `ls::BackTrack`             : line search parameters
  - `x_prev::AbstractVector`    : previous state
  - `p::AbstractVector`         : search direction
  - `f_prev::Real`              : previous objective function value
  - `∇f_prev::AbstractVector`   : previous gradient
  - `f::Function`               : ``f(x)``

#### Returns
  - `x::AbstractVector` : updated state
  - `f_new::Real`       : objective function value
"""
function backtrack!(ls::BackTrack, x::AbstractVector, x_prev::AbstractVector, 
                    p::AbstractVector, f_prev::Real, ∇f_prev::AbstractVector, 
                    f::Function)
    # Check interpolation
    if ls.order == 0
        # Fixed backtracking
        f_new= fixed!(ls, x, x_prev, p, f_prev, ∇f_prev, f)
    elseif ls.order == 2
        # Quadratic interpolation
        f_new= quad_interp!(ls, x, x_prev, p, f_prev, ∇f_prev, f)
    elseif ls.order == 3
        # Cubic interpolation
        f_new= cubic_interp!(ls, x, x_prev, p, f_prev, ∇f_prev, f)
    end

    return f_new
end