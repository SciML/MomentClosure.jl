abstract type MomentEquations end

"""
$(TYPEDEF)

Raw moment equations generated for the given system plus a number of
helper parameters (used internally).

# Fields
$(FIELDS)
"""
struct RawMomentEquations <: MomentEquations
    """[`ModelingToolkit.System`](https://docs.sciml.ai/ModelingToolkit/stable/API/System/#ModelingToolkitBase.System)
    consisting of the time-evolution equations of raw moments."""
    odes::System
    """Dictionary mapping species of the associated `ReactionSystem` to their index in the moment equations."""
    smap::Dict
    """Symbolic variables defining the raw moments."""
    μ::OrderedDict
    """Number of species within the system."""
    N::Int
    """Order of moment equations."""
    m_order::Int
    """Expansion order."""
    q_order::Int
    """Vector of all index combinations up to `q_order`."""
    iter_all::Vector
    """Vector of all index combinations up to `m_order`."""
    iter_m::Vector
    """Vector of all index combinations of order greater than `m_order`
    up to `q_order`."""
    iter_q::Vector
    """Vector of index combinations of order 1."""
    iter_1::Vector
end

"""
$(TYPEDEF)

Central moment equations generated for the given system plus a number of
helper parameters (used internally).

# Fields
$(FIELDS)
"""
struct CentralMomentEquations <: MomentEquations
    """[`ModelingToolkit.System`](https://docs.sciml.ai/ModelingToolkit/stable/API/System/#ModelingToolkitBase.System)
    consisting of the time-evolution equations of central moments."""
    odes::System
    """Dictionary mapping species of the associated `ReactionSystem` to their index in the moment equations."""
    smap::Dict
    """Symbolic variables defining the means."""
    μ::OrderedDict
    """Symbolic variables defining the central moments."""
    M::Dict
    """Number of species within the system."""
    N::Int
    """Order of moment equations."""
    m_order::Int
    """Expansion order."""
    q_order::Int
    """Vector of all index combinations up to `q_order`."""
    iter_all::Vector
    """Vector of all index combinations up to `m_order`."""
    iter_m::Vector
    """Vector of all index combinations of order greater than `m_order`
    up to `q_order`."""
    iter_q::Vector
    """Vector of index combinations of order 1."""
    iter_1::Vector
end

"""
$(TYPEDEF)

Closed moment equations and the corresponding closure functions.

# Fields
$(FIELDS)
"""
struct ClosedMomentEquations <: MomentEquations
    """[`ModelingToolkit.System`](https://docs.sciml.ai/ModelingToolkit/stable/API/System/#ModelingToolkitBase.System)
    consisting of the time-evolution equations of *closed* moments."""
    odes::System
    """Dictionary of moment closure functions for each higher order moment."""
    closure::OrderedDict
    """Original raw or central moment equations (before closure was applied)."""
    open_eqs::MomentEquations
end

"""
```julia
SciMLBase.ODEProblem(sys::MomentEquations, u0, tspan, p=NullParameters();
                     use_deterministic_IC::Bool=true, kwargs...)
```

Generates a [`SciMLBase.ODEProblem`](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEProblem) from `MomentEquations`.

Note that the initial condition `u0` can be defined as follows:
- If `u0` contains the initial molecule numbers (either as a vector or a mapping) and `use_deterministic_IC` is set to `true`,
  `u0` will be mapped to the corresponding raw/central moments under the assumption of deterministic
  initial conditions (using the function [`deterministic_IC`](@ref)). This is the default use case.
- To handle a more complicated initial condition, `u0` can be defined with the initial values of each moment directly,
  but `use_deterministic_IC` must be set to `false` to ensure correctness.
- If `u0` is a plain vector, its elements are read in the order of the species in
  [`Catalyst.speciesmap`](https://docs.sciml.ai/Catalyst/stable/api/core_api/#Catalyst.speciesmap)`(sys)`
  (or of the moments in `unknowns(sys)` if `use_deterministic_IC = false`).
"""
function SciMLBase.ODEProblem(
        sys::MomentEquations, u0, tspan, p = NullParameters();
        use_deterministic_IC::Bool = true, kwargs...
    )
    odes = complete(get_odes(sys))
    u0map = use_deterministic_IC ? deterministic_IC(u0, sys) : u0
    op = merge(_u0_opmap(u0map, unknowns(odes)), _p_opmap(p))
    return ODEProblem(odes, op, tspan; kwargs...)
end

# Convert `u0` and `p` to a single map
# Plain `u0` vector follows the ordering set by MomentClosure
_u0_opmap(x::AbstractDict, _) = x
_u0_opmap(x::Union{Tuple, AbstractArray{<:Pair}}, _) = Dict(x)
_u0_opmap(x::AbstractArray, syms) = isempty(x) ? Dict() : Dict(syms .=> x)

# Plain `p` vector is rejected (ordering is internal to ModelingToolkit)
_p_opmap(::Union{Nothing, NullParameters}) = Dict()
_p_opmap(x::AbstractDict) = x
_p_opmap(x::Union{Tuple, AbstractArray{<:Pair}}) = Dict(x)
function _p_opmap(x::AbstractArray)
    isempty(x) && return Dict()
    throw(
        ArgumentError(
            "parameters must be given as a map, e.g. `[k₁ => 1.0, …]` or a `Dict`; " *
                "a plain vector is ambiguous because the parameter order is not guaranteed."
        )
    )
end

function Base.nameof(eqs::MomentEquations)
    return nameof(get_odes(eqs))
end

# Basic `MomentEquations`-specific accessors
"""
    get_odes(sys::MomentEquations)

Return the `ModelingToolkit.System` of moment equations.
"""
get_odes(sys::MomentEquations) = getfield(sys, :odes)

"""
    get_closure(sys::ClosedMomentEquations)

Return the dictionary of moment closure functions for each higher-order moment.
"""
get_closure(sys::ClosedMomentEquations) = getfield(sys, :closure)

Catalyst.speciesmap(sys::MomentEquations) = getfield(sys, :smap)
Catalyst.speciesmap(sys::ClosedMomentEquations) = speciesmap(getfield(sys, :open_eqs))
ModelingToolkitBase.get_iv(sys::MomentEquations) = get_iv(get_odes(sys))
ModelingToolkitBase.get_eqs(sys::MomentEquations) = get_eqs(get_odes(sys))
ModelingToolkitBase.unknowns(sys::MomentEquations) = unknowns(get_odes(sys))
ModelingToolkitBase.get_ps(sys::MomentEquations) = get_ps(get_odes(sys))
