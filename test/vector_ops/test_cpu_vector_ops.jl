using Test
include("../../src/vector_ops/VectorOps.jl")


using .VectorOps: CPU, GPU, ddx_up

@testset "ddx_up Tests for CPUVectorOps" begin
    dx = 0.1

    # Prefactors (you can adapt these from your `VectorOps` module)
#    const pf2 = VectorOps.prefactors_2nd
    # Test 1D
    backend = CPU()
    field_1D = [1.0, 2.0, 4.0, 7.0, 11.0]  # Linear test field
    expected_1D = [0.0, 20.0, 30.0, 40.0, 0.0]  # Analytical derivative (zeros at boundaries)
    result_1D = ddx_up(backend,field_1D, dx, 2)
    @test size(result_1D) == size(field_1D)
    @test isapprox(result_1D, expected_1D; atol=1e-5)

end
