using MomentClosure, BenchmarkTools
using Catalyst, Symbolics

const SUITE = BenchmarkGroup()

rn = @reaction_network begin
    (c₁ / Ω^2), 2X + Y → 3X
    (c₂), X → Y
    (Ω * c₃, c₄), 0 ↔ X
end

# =============================================================================
# Moment equation generation
# =============================================================================

SUITE["generate"] = BenchmarkGroup()

SUITE["generate"]["central_2"] = @benchmarkable generate_central_moment_eqs(
    $rn, 2; combinatoric_ratelaws = false
)
SUITE["generate"]["central_2_4"] = @benchmarkable generate_central_moment_eqs(
    $rn, 2, 4; combinatoric_ratelaws = false
)
SUITE["generate"]["raw_2"] = @benchmarkable generate_raw_moment_eqs(
    $rn, 2; combinatoric_ratelaws = false
)
SUITE["generate"]["raw_3"] = @benchmarkable generate_raw_moment_eqs(
    $rn, 3; combinatoric_ratelaws = false
)

# =============================================================================
# Moment closure
# =============================================================================

sys = generate_central_moment_eqs(rn, 2, 4; combinatoric_ratelaws = false)
raw_sys = generate_raw_moment_eqs(rn, 2; combinatoric_ratelaws = false)

SUITE["closure"] = BenchmarkGroup()

SUITE["closure"]["zero_central"] = @benchmarkable moment_closure($sys, "zero")
SUITE["closure"]["normal_raw"] = @benchmarkable moment_closure($raw_sys, "normal")
