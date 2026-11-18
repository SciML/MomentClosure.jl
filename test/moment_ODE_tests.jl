using MomentClosure
using Test
using Catalyst
using OrdinaryDiffEqTsit5
using ModelingToolkit: t_nounits as t
using SciMLBase: ODEProblem
using ModelingToolkitBase: complete, unknowns

@parameters k₁ k₂ k₃ k₄ k₅ k₆ k₇
@species x(t) y(t) y₀(t)

rn = @reaction_network begin
    @parameters k₁ k₂ k₃ k₄ k₅ k₆ k₇
    (k₁), 0 → x
    (k₂), x → 0
    (k₃ * x * y / (x + k₇)), x ⇒ 0
    (k₄ * x), 0 ⇒ y₀
    (k₅), y₀ → y
    (k₆), y → 0
end

pmap = [:k₁ => 90, :k₂ => 0.002, :k₃ => 1.7, :k₄ => 1.1, :k₅ => 0.93, :k₆ => 0.96, :k₇ => 0.01]
tspan = (0.0, 50.0)
x0 = Dict(x => 70, y => 60, y₀ => 30)
# Reference mean of x is the converged solution (stable to ~10 digits for reltol/abstol ≤ 1e-10);
# integrate at a tight tolerance so the comparison is independent of solver defaults.
res = [35.156249675, 60.255177096, 29.237221341, 65.59590462, 24.900095364]
solve_tight(oprob) = solve(oprob, Tsit5(), saveat = 10:10:50, reltol = 1.0e-10, abstol = 1.0e-10)

# the default species ordering and a cyclic shift of it
smap_default = speciesmap(rn)
smap_shifted = Dict(s => mod1(i + 1, numspecies(rn)) for (s, i) in smap_default)
@test !isequal(smap_shifted, smap_default)

@testset "species ordering: $label" for (label, smap) in
    [("default", smap_default), ("shifted", smap_shifted)]
    eqs = generate_central_moment_eqs(rn, 2, 3; combinatoric_ratelaws = false, smap)
    closed_eqs = moment_closure(eqs, "normal")
    @test isequal(speciesmap(closed_eqs), smap)

    # all forms of u0 (Symbol pairs, Dict, symbolic pairs, plain vector in smap order)
    # give the same initial moments
    u0vec = [x0[s] for s in first.(sort(collect(smap), by = last))]
    u0forms = [(:y₀ => 30, :x => 70, :y => 60), x0, collect(x0), u0vec]
    @test allequal(deterministic_IC(u0, closed_eqs) for u0 in u0forms)
    oprobs = [ODEProblem(closed_eqs, u0, tspan, pmap) for u0 in u0forms]
    @test allequal(oprob.u0 for oprob in oprobs)

    # the mean of x is the unknown at index smap[x]
    @test isapprox(solve_tight(oprobs[1])[smap[x], :], res, atol = 1.0e-6)

    # moment-level initial conditions: a plain vector relies on the order of the unknowns
    @test isequal(unknowns(complete(get_odes(closed_eqs))), unknowns(closed_eqs))
    oprob = ODEProblem(closed_eqs, oprobs[1].u0, tspan, pmap, use_deterministic_IC = false)
    @test isapprox(solve_tight(oprob)[smap[x], :], res, atol = 1.0e-6)

    # a plain parameter vector is ambiguous and must be rejected
    @test_throws ArgumentError ODEProblem(closed_eqs, x0, tspan, last.(pmap))
end
