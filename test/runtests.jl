using Test
using TimerOutputs


include("../src/Yggdrasil.jl")
using .Yggdrasil

include("solvers/PIC/test_particle_pusher.jl")
include("vector_ops/test_cpu_vector_ops.jl")