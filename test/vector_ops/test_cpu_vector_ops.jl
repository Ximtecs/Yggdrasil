using Test
include("../../src/vector_ops/VectorOps.jl")


using .VectorOps: CPU, GPU, ddx_up, ddy_up, ddz_up, ddx_dn, ddy_dn, ddz_dn



@testset "ddx_up Tests for CPUVectorOps" begin
    dx = 0.01  # Spacing in x-direction
    backend = CPU()
    #---------------- Define a 1D test function and its derivative --------------------
    f_1D(x) = sin(x)                # Function: sin(x)
    df_1D(x) = cos(x)               # Derivative: cos(x)
    x_1D = 0:dx:(2π)                # 1D grid
    field_1D = f_1D.(x_1D)          # Evaluate the function on the grid
    expected_1D = df_1D.(x_1D .+ 0.5 * dx)      # Analytical derivative
    expected_1D[end] = 0.0          # Set boundary derivative to 0
    result_1D = ddx_up(backend, field_1D, dx, 2)
    @test size(result_1D) == size(field_1D)
    @test isapprox(result_1D, expected_1D; rtol=1e-5) #atol=1e-5, rtol=0)

    #------- now 4th order test
    result_1D = ddx_up(backend, field_1D, dx, 4)
    expected_1D[end-1] = 0.0          # Set boundary derivative to 0
    expected_1D[1] = 0.0          # Set boundary derivative to 0
    @test isapprox(result_1D, expected_1D; rtol=1e-5) #atol=1e-5, rtol=0)

    #------- now 6th order test
    result_1D = ddx_up(backend, field_1D, dx, 6)
    expected_1D[end-2] = 0.0          # Set boundary derivative to 0
    expected_1D[2] = 0.0          # Set boundary derivative to 0
    @test isapprox(result_1D, expected_1D; rtol=1e-5) #atol=1e-5, rtol=0)
    #-----------------------------------------------------------------------------------
    #------------- Define a 2D test function and its derivative along x -----------------
    f_2D(x, y) = sin(x) + cos(y)    # Function: sin(x) + cos(y)
    df_2D(x, y) = cos(x)            # Derivative w.r.t x: cos(x)
    x_2D = 0:dx:(2π)                # x-grid
    y_2D = 0:dx:(2π)                # y-grid
    field_2D = [f_2D(x, y) for x in x_2D, y in y_2D]
    expected_2D = [df_2D(x .+ 0.5 * dx, y) for x in x_2D, y in y_2D]
    expected_2D[end, :] .= 0.0      # Set boundary derivatives to 0
    result_2D = ddx_up(backend, field_2D, dx, 2)
    @test size(result_2D) == size(field_2D)
    @test isapprox(result_2D, expected_2D; rtol=1e-5)
    #------- now 4th order test
    result_2D = ddx_up(backend, field_2D, dx, 4)
    expected_2D[end-1, :] .= 0.0      # Set boundary derivatives to 0
    expected_2D[1, :] .= 0.0      # Set boundary derivatives to 0
    @test isapprox(result_2D, expected_2D; rtol=1e-5)
    #------- now 6th order test
    result_2D = ddx_up(backend, field_2D, dx, 6)
    expected_2D[end-2, :] .= 0.0      # Set boundary derivatives to 0
    expected_2D[2, :] .= 0.0      # Set boundary derivatives to 0
    @test isapprox(result_2D, expected_2D; rtol=1e-5)
    #--------------------------------------------------------------------------------------
    #------------- Define a 3D test function and its derivative along x -------------------
    f_3D(x, y, z) = sin(x) + cos(y) + z^2  # Function: sin(x) + cos(y) + z^2
    df_3D(x, y, z) = cos(x)                # Derivative w.r.t x: cos(x)
    x_3D = 0:dx:(2π)                       # x-grid
    y_3D = 0:dx:(2π)                       # y-grid
    z_3D = 0:dx:0.2                        # z-grid
    field_3D = [f_3D(x, y, z) for x in x_3D, y in y_3D, z in z_3D]
    expected_3D = [df_3D(x .+ 0.5 * dx, y, z) for x in x_3D, y in y_3D, z in z_3D]
    expected_3D[end, :, :] .= 0.0          # Set boundary derivatives to 0
    result_3D = ddx_up(backend, field_3D, dx, 2)
    @test size(result_3D) == size(field_3D)
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #------- now 4th order test
    result_3D = ddx_up(backend, field_3D, dx, 4)
    expected_3D[end-1, :, :] .= 0.0          # Set boundary derivatives to 0
    expected_3D[1, :, :] .= 0.0          # Set boundary derivatives to 0
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #------- now 6th order test
    result_3D = ddx_up(backend, field_3D, dx, 6)
    expected_3D[end-2, :, :] .= 0.0          # Set boundary derivatives to 0
    expected_3D[2, :, :] .= 0.0          # Set boundary derivatives to 0
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #--------------------------------------------------------------------------------------
end

@testset "ddy_up Tests for CPUVectorOps" begin
    dy = 0.01  # Spacing in y-direction
    backend = CPU()
    #---------------- Define a 1D test function and its derivative --------------------
    f_1D(y) = sin(y)                # Function: sin(y)
    df_1D(y) = cos(y)               # Derivative: cos(y)
    y_1D = 0:dy:(2π)                # 1D grid
    field_1D = f_1D.(y_1D)          # Evaluate the function on the grid
    expected_1D = df_1D.(y_1D .+ 0.5 * dy)      # Analytical derivative
    expected_1D[end] = 0.0          # Set boundary derivative to 0
    result_1D = ddy_up(backend, field_1D, dy, 2)
    @test size(result_1D) == size(field_1D)
    @test isapprox(result_1D, expected_1D; rtol=1e-5)

    #------ no need to test further 1D cases for ddy as it falls back to ddx 

    #-----------------------------------------------------------------------------------
    #------------- Define a 2D test function and its derivative along y -----------------
    f_2D(x, y) = sin(x) + cos(y)    # Function: sin(x) + cos(y)
    df_2D(x, y) = -sin(y)           # Derivative w.r.t y: -sin(y)
    x_2D = 0:dy:(2π)                # x-grid
    y_2D = 0:dy:(2π)                # y-grid
    field_2D = [f_2D(x, y) for x in x_2D, y in y_2D]
    expected_2D = [df_2D(x, y .+ 0.5 * dy) for x in x_2D, y in y_2D]
    expected_2D[:, end] .= 0.0      # Set boundary derivatives to 0
    result_2D = ddy_up(backend, field_2D, dy, 2)
    @test size(result_2D) == size(field_2D)
    @test isapprox(result_2D, expected_2D; rtol=1e-5)
    #------- now 4th order test
    result_2D = ddy_up(backend, field_2D, dy, 4)
    expected_2D[:, end-1] .= 0.0      # Set boundary derivatives to 0
    expected_2D[:, 1] .= 0.0      # Set boundary derivatives to 0
    @test isapprox(result_2D, expected_2D; rtol=1e-5)
    #------- now 6th order test
    result_2D = ddy_up(backend, field_2D, dy, 6)
    expected_2D[:, end-2] .= 0.0      # Set boundary derivatives to 0
    expected_2D[:, 2] .= 0.0      # Set boundary derivatives to 0
    @test isapprox(result_2D, expected_2D; rtol=1e-5)
    #--------------------------------------------------------------------------------------
    #------------- Define a 3D test function and its derivative along y -------------------
    f_3D(x, y, z) = sin(x) + cos(y) + z^2  # Function: sin(x) + cos(y) + z^2
    df_3D(x, y, z) = -sin(y)               # Derivative w.r.t y: -sin(y)
    x_3D = 0:dy:(2π)                       # x-grid
    y_3D = 0:dy:(2π)                       # y-grid
    z_3D = 0:dy:0.2                        # z-grid
    field_3D = [f_3D(x, y, z) for x in x_3D, y in y_3D, z in z_3D]
    expected_3D = [df_3D(x, y .+ 0.5 * dy, z) for x in x_3D, y in y_3D, z in z_3D]
    expected_3D[:, end, :] .= 0.0          # Set boundary derivatives to 0
    result_3D = ddy_up(backend, field_3D, dy, 2)
    @test size(result_3D) == size(field_3D)
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #------- now 4th order test
    result_3D = ddy_up(backend, field_3D, dy, 4)
    expected_3D[:, end-1, :] .= 0.0          # Set boundary derivatives to 0
    expected_3D[:, 1, :] .= 0.0          # Set boundary derivatives to 0
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #------- now 6th order test
    result_3D = ddy_up(backend, field_3D, dy, 6)
    expected_3D[:, end-2, :] .= 0.0          # Set boundary derivatives to 0
    expected_3D[:, 2, :] .= 0.0          # Set boundary derivatives to 0
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #--------------------------------------------------------------------------------------
end


@testset "ddz_up Tests for CPUVectorOps" begin
    dz = 0.01  # Spacing in z-direction
    backend = CPU()
    #---------------- Define a 1D test function and its derivative --------------------
    f_1D(z) = sin(z)                # Function: sin(z)
    df_1D(z) = cos(z)               # Derivative: cos(z)
    z_1D = 0:dz:(2π)                # 1D grid
    field_1D = f_1D.(z_1D)          # Evaluate the function on the grid
    expected_1D = df_1D.(z_1D .+ 0.5 * dz)      # Analytical derivative
    expected_1D[end] = 0.0          # Set boundary derivative to 0
    result_1D = ddz_up(backend, field_1D, dz, 2)
    @test size(result_1D) == size(field_1D)
    @test isapprox(result_1D, expected_1D; rtol=1e-5)

    #------ no need to test further 1D cases for ddz as it falls back to ddx 

    #-----------------------------------------------------------------------------------
    #------------- Define a 2D test function and its derivative along z -----------------
    f_2D(x, z) = sin(x) + cos(z)    # Function: sin(x) + cos(z)
    df_2D(x, z) = -sin(z)           # Derivative w.r.t z: -sin(z)
    x_2D = 0:dz:(2π)                # x-grid
    z_2D = 0:dz:(2π)                # z-grid
    field_2D = [f_2D(x, z) for x in x_2D, z in z_2D]
    expected_2D = [df_2D(x, z .+ 0.5 * dz) for x in x_2D, z in z_2D]
    expected_2D[:, end] .= 0.0      # Set boundary derivatives to 0
    result_2D = ddz_up(backend, field_2D, dz, 2)
    @test size(result_2D) == size(field_2D)
    @test isapprox(result_2D, expected_2D; rtol=1e-5)

    #------ no need to test further 2D cases for ddz as it falls back to ddy 
    #--------------------------------------------------------------------------------------
    #------------- Define a 3D test function and its derivative along z -------------------
    f_3D(x, y, z) = sin(x) + cos(y) + z^2  # Function: sin(x) + cos(y) + z^2
    df_3D(x, y, z) = 2z                    # Derivative w.r.t z: 2z
    x_3D = 0:0.1:(2π)                       # x-grid
    y_3D = 0:0.1:(2π)                       # y-grid
    z_3D = 0:dz:(2)                       # z-grid
    field_3D = [f_3D(x, y, z) for x in x_3D, y in y_3D, z in z_3D]
    expected_3D = [df_3D(x, y, z .+ 0.5 * dz) for x in x_3D, y in y_3D, z in z_3D]
    expected_3D[:, :, end] .= 0.0          # Set boundary derivatives to 0
    result_3D = ddz_up(backend, field_3D, dz, 2)
    @test size(result_3D) == size(field_3D)
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #------- now 4th order test
    result_3D = ddz_up(backend, field_3D, dz, 4)
    expected_3D[:, :, end-1] .= 0.0          # Set boundary derivatives to 0
    expected_3D[:, :, 1] .= 0.0          # Set boundary derivatives to 0
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #------- now 6th order test
    result_3D = ddz_up(backend, field_3D, dz, 6)
    expected_3D[:, :, end-2] .= 0.0          # Set boundary derivatives to 0
    expected_3D[:, :, 2] .= 0.0          # Set boundary derivatives to 0
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #--------------------------------------------------------------------------------------
end

@testset "ddx_dn Tests for CPUVectorOps" begin
    dx = 0.01  # Spacing in x-direction
    backend = CPU()
    #---------------- Define a 1D test function and its derivative --------------------
    f_1D(x) = sin(x)                # Function: sin(x)
    df_1D(x) = cos(x)               # Derivative: cos(x)
    x_1D = 0:dx:(2π)                # 1D grid
    field_1D = f_1D.(x_1D)          # Evaluate the function on the grid
    expected_1D = df_1D.(x_1D .- 0.5 * dx)      # Analytical derivative
    expected_1D[1] = 0.0            # Set boundary derivative to 0
    result_1D = ddx_dn(backend, field_1D, dx, 2)
    @test size(result_1D) == size(field_1D)
    @test isapprox(result_1D, expected_1D; rtol=1e-5)

    #------- now 4th order test
    result_1D = ddx_dn(backend, field_1D, dx, 4)
    expected_1D[2] = 0.0            # Set boundary derivative to 0
    expected_1D[end] = 0.0          # Set boundary derivative to 0
    @test isapprox(result_1D, expected_1D; rtol=1e-5)

    #------- now 6th order test
    result_1D = ddx_dn(backend, field_1D, dx, 6)
    expected_1D[3] = 0.0            # Set boundary derivative to 0
    expected_1D[end-1] = 0.0        # Set boundary derivative to 0
    @test isapprox(result_1D, expected_1D; rtol=1e-5)
    #-----------------------------------------------------------------------------------
    #------------- Define a 2D test function and its derivative along x -----------------
    f_2D(x, y) = sin(x) + cos(y)    # Function: sin(x) + cos(y)
    df_2D(x, y) = cos(x)            # Derivative w.r.t x: cos(x)
    x_2D = 0:dx:(2π)                # x-grid
    y_2D = 0:dx:(2π)                # y-grid
    field_2D = [f_2D(x, y) for x in x_2D, y in y_2D]
    expected_2D = [df_2D(x .- 0.5 * dx, y) for x in x_2D, y in y_2D]
    expected_2D[1, :] .= 0.0        # Set boundary derivatives to 0
    result_2D = ddx_dn(backend, field_2D, dx, 2)
    @test size(result_2D) == size(field_2D)
    @test isapprox(result_2D, expected_2D; rtol=1e-5)
    #------- now 4th order test
    result_2D = ddx_dn(backend, field_2D, dx, 4)
    expected_2D[2, :] .= 0.0        # Set boundary derivatives to 0
    expected_2D[end, :] .= 0.0      # Set boundary derivatives to 0
    @test isapprox(result_2D, expected_2D; rtol=1e-5)
    #------- now 6th order test
    result_2D = ddx_dn(backend, field_2D, dx, 6)
    expected_2D[3, :] .= 0.0        # Set boundary derivatives to 0
    expected_2D[end-1, :] .= 0.0    # Set boundary derivatives to 0
    @test isapprox(result_2D, expected_2D; rtol=1e-5)
    #--------------------------------------------------------------------------------------
    #------------- Define a 3D test function and its derivative along x -------------------
    f_3D(x, y, z) = sin(x) + cos(y) + z^2  # Function: sin(x) + cos(y) + z^2
    df_3D(x, y, z) = cos(x)                # Derivative w.r.t x: cos(x)
    x_3D = 0:dx:(2π)                       # x-grid
    y_3D = 0:dx:(2π)                       # y-grid
    z_3D = 0:dx:0.2                        # z-grid
    field_3D = [f_3D(x, y, z) for x in x_3D, y in y_3D, z in z_3D]
    expected_3D = [df_3D(x .- 0.5 * dx, y, z) for x in x_3D, y in y_3D, z in z_3D]
    expected_3D[1, :, :] .= 0.0            # Set boundary derivatives to 0
    result_3D = ddx_dn(backend, field_3D, dx, 2)
    @test size(result_3D) == size(field_3D)
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #------- now 4th order test
    result_3D = ddx_dn(backend, field_3D, dx, 4)
    expected_3D[2, :, :] .= 0.0            # Set boundary derivatives to 0
    expected_3D[end, :, :] .= 0.0          # Set boundary derivatives to 0
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #------- now 6th order test
    result_3D = ddx_dn(backend, field_3D, dx, 6)
    expected_3D[3, :, :] .= 0.0            # Set boundary derivatives to 0
    expected_3D[end-1, :, :] .= 0.0        # Set boundary derivatives to 0
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #--------------------------------------------------------------------------------------
end

@testset "ddy_dn Tests for CPUVectorOps" begin
    dy = 0.01  # Spacing in y-direction
    backend = CPU()
    #---------------- Define a 1D test function and its derivative --------------------
    f_1D(y) = sin(y)                # Function: sin(y)
    df_1D(y) = cos(y)               # Derivative: cos(y)
    y_1D = 0:dy:(2π)                # 1D grid
    field_1D = f_1D.(y_1D)          # Evaluate the function on the grid
    expected_1D = df_1D.(y_1D .- 0.5 * dy)      # Analytical derivative
    expected_1D[1] = 0.0            # Set boundary derivative to 0
    result_1D = ddy_dn(backend, field_1D, dy, 2)
    @test size(result_1D) == size(field_1D)
    @test isapprox(result_1D, expected_1D; rtol=1e-5)

    #------ no need to test further 1D cases for ddy as it falls back to ddx 

    #-----------------------------------------------------------------------------------
    #------------- Define a 2D test function and its derivative along y -----------------
    f_2D(x, y) = sin(x) + cos(y)    # Function: sin(x) + cos(y)
    df_2D(x, y) = -sin(y)           # Derivative w.r.t y: -sin(y)
    x_2D = 0:dy:(2π)                # x-grid
    y_2D = 0:dy:(2π)                # y-grid
    field_2D = [f_2D(x, y) for x in x_2D, y in y_2D]
    expected_2D = [df_2D(x, y .- 0.5 * dy) for x in x_2D, y in y_2D]
    expected_2D[:, 1] .= 0.0        # Set boundary derivatives to 0
    result_2D = ddy_dn(backend, field_2D, dy, 2)
    @test size(result_2D) == size(field_2D)
    @test isapprox(result_2D, expected_2D; rtol=1e-5)
    #------- now 4th order test
    result_2D = ddy_dn(backend, field_2D, dy, 4)
    expected_2D[:, 2] .= 0.0        # Set boundary derivatives to 0
    expected_2D[:, end] .= 0.0      # Set boundary derivatives to 0
    @test isapprox(result_2D, expected_2D; rtol=1e-5)
    #------- now 6th order test
    result_2D = ddy_dn(backend, field_2D, dy, 6)
    expected_2D[:, 3] .= 0.0        # Set boundary derivatives to 0
    expected_2D[:, end-1] .= 0.0    # Set boundary derivatives to 0
    @test isapprox(result_2D, expected_2D; rtol=1e-5)
    #--------------------------------------------------------------------------------------
    #------------- Define a 3D test function and its derivative along y -------------------
    f_3D(x, y, z) = sin(x) + cos(y) + z^2  # Function: sin(x) + cos(y) + z^2
    df_3D(x, y, z) = -sin(y)               # Derivative w.r.t y: -sin(y)
    x_3D = 0:dy:(2π)                       # x-grid
    y_3D = 0:dy:(2π)                       # y-grid
    z_3D = 0:dy:0.2                        # z-grid
    field_3D = [f_3D(x, y, z) for x in x_3D, y in y_3D, z in z_3D]
    expected_3D = [df_3D(x, y .- 0.5 * dy, z) for x in x_3D, y in y_3D, z in z_3D]
    expected_3D[:, 1, :] .= 0.0            # Set boundary derivatives to 0
    result_3D = ddy_dn(backend, field_3D, dy, 2)
    @test size(result_3D) == size(field_3D)
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #------- now 4th order test
    result_3D = ddy_dn(backend, field_3D, dy, 4)
    expected_3D[:, 2, :] .= 0.0            # Set boundary derivatives to 0
    expected_3D[:, end, :] .= 0.0          # Set boundary derivatives to 0
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #------- now 6th order test
    result_3D = ddy_dn(backend, field_3D, dy, 6)
    expected_3D[:, 3, :] .= 0.0            # Set boundary derivatives to 0
    expected_3D[:, end-1, :] .= 0.0        # Set boundary derivatives to 0
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #--------------------------------------------------------------------------------------
end

@testset "ddz_dn Tests for CPUVectorOps" begin
    dz = 0.01  # Spacing in z-direction
    backend = CPU()
    #---------------- Define a 1D test function and its derivative --------------------
    f_1D(z) = sin(z)                # Function: sin(z)
    df_1D(z) = cos(z)               # Derivative: cos(z)
    z_1D = 0:dz:(2π)                # 1D grid
    field_1D = f_1D.(z_1D)          # Evaluate the function on the grid
    expected_1D = df_1D.(z_1D .- 0.5 * dz)      # Analytical derivative
    expected_1D[1] = 0.0            # Set boundary derivative to 0
    result_1D = ddz_dn(backend, field_1D, dz, 2)
    @test size(result_1D) == size(field_1D)
    @test isapprox(result_1D, expected_1D; rtol=1e-5)

    #------ no need to test further 1D cases for ddz as it falls back to ddx 

    #-----------------------------------------------------------------------------------
    #------------- Define a 2D test function and its derivative along z -----------------
    f_2D(x, z) = sin(x) + cos(z)    # Function: sin(x) + cos(z)
    df_2D(x, z) = -sin(z)           # Derivative w.r.t z: -sin(z)
    x_2D = 0:dz:(2π)                # x-grid
    z_2D = 0:dz:(2π)                # z-grid
    field_2D = [f_2D(x, z) for x in x_2D, z in z_2D]
    expected_2D = [df_2D(x, z .- 0.5 * dz) for x in x_2D, z in z_2D]
    expected_2D[:, 1] .= 0.0        # Set boundary derivatives to 0
    result_2D = ddz_dn(backend, field_2D, dz, 2)
    @test size(result_2D) == size(field_2D)
    @test isapprox(result_2D, expected_2D; rtol=1e-5)

    #------ no need to test further 2D cases for ddz as it falls back to ddy 
    #--------------------------------------------------------------------------------------
    #------------- Define a 3D test function and its derivative along z -------------------
    f_3D(x, y, z) = sin(x) + cos(y) + z^2  # Function: sin(x) + cos(y) + z^2
    df_3D(x, y, z) = 2z                    # Derivative w.r.t z: 2z
    x_3D = 0:0.1:(2π)                       # x-grid
    y_3D = 0:0.1:(2π)                       # y-grid
    z_3D = 0:dz:(2)                       # z-grid
    field_3D = [f_3D(x, y, z) for x in x_3D, y in y_3D, z in z_3D]
    expected_3D = [df_3D(x, y, z .- 0.5 * dz) for x in x_3D, y in y_3D, z in z_3D]
    expected_3D[:, :, 1] .= 0.0            # Set boundary derivatives to 0
    result_3D = ddz_dn(backend, field_3D, dz, 2)
    @test size(result_3D) == size(field_3D)
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #------- now 4th order test
    result_3D = ddz_dn(backend, field_3D, dz, 4)
    expected_3D[:, :, 2] .= 0.0            # Set boundary derivatives to 0
    expected_3D[:, :, end] .= 0.0          # Set boundary derivatives to 0
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #------- now 6th order test
    result_3D = ddz_dn(backend, field_3D, dz, 6)
    expected_3D[:, :, 3] .= 0.0            # Set boundary derivatives to 0
    expected_3D[:, :, end-1] .= 0.0        # Set boundary derivatives to 0
    @test isapprox(result_3D, expected_3D; rtol=1e-5)
    #--------------------------------------------------------------------------------------
end