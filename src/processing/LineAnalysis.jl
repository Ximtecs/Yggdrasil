"""
    get_line(midpoint::Vector{Float64}, direction::Vector{Float64}, half_length::Float64, num_points::Int)

Construct a straight line in 3D space by sampling `num_points` along a line of length `2 * half_length`,
centred at `midpoint` and oriented along the unit `direction` vector.

# Arguments
- `midpoint`: A 3-element vector `[x, y, z]` representing the centre of the line.
- `direction`: A 3-element unit vector `[dx, dy, dz]` giving the line's direction.
- `half_length`: Half the total length of the line.
- `num_points`: Number of points to sample along the line.

# Returns
- A vector of 3-element arrays representing the points along the line.
"""
function get_line(midpoint::Vector{Float64},
                  direction::Vector{Float64},
                  half_length::Float64,
                  num_points::Int)

    # Validate input lengths
    if length(midpoint) != 3 || length(direction) != 3
        error("Both `midpoint` and `direction` must be 3-element vectors.")
    end

    # Compute endpoints of the line
    p1 = midpoint .- half_length .* direction
    p2 = midpoint .+ half_length .* direction

    # Direction vector from p1 to p2
    dvec = p2 .- p1

    # Linearly interpolate points along the line
    ts = range(0.0, stop=1.0, length=num_points)
    line_points = [p1 .+ t .* dvec for t in ts]

    return line_points
end

"""
    generate_line_plane(midpoint, direction, normal, length, width, n_perp; n_par=nothing)

Generate a 2D sampling plane aligned with a given line and normal vector.

# Arguments
- `midpoint::Vector{Float64}`: Centre of the line and plane (2D or 3D).
- `direction::Vector{Float64}`: Vector along the main line direction.
- `normal::Vector{Float64}`: Vector normal to the 2D sampling plane.
- `length::Float64`: Half-length of the line along the direction vector.
- `width::Float64`: Half-width of the plane in the perpendicular direction.
- `n_perp::Int`: Number of grid points in the perpendicular direction.
- `n_par::Int` (optional): Number of points along the line. If not given, chosen to match spacing in `n_perp`.

# Returns
- `line`: Points along the central line.
- `perp_vals`: Array of perpendicular offsets.
- `par_vals`: Array of line-aligned (parallel) offsets.
- `Δpar`, `Δperp`: Grid spacing in each direction.
- `grid`: 2D array of global coordinates [parallel, perpendicular].
"""
function generate_line_plane(midpoint::Vector{Float64},
                             direction::Vector{Float64},
                             normal::Vector{Float64},
                             half_length::Float64,
                             width::Float64,
                             n_perp::Int; 
                             n_par::Union{Nothing, Int}=nothing)

    if length(direction) != 3 || length(normal) != 3 || length(midpoint) != 3
        error("All input vectors must have the same length (2D or 3D).")
    end

    # Unit vectors
    û = direction ./ norm(direction)
    n̂_raw = cross(normal, û)
    if norm(n̂_raw) ≈ 0
        error("Direction and normal vectors must not be parallel.")
    end
    n̂ = n̂_raw ./ norm(n̂_raw)

    # Check orthogonality
    if abs(dot(û, n̂)) > 1e-10
        error("Direction and normal vectors must be orthogonal.")
    end

    # Build axis values
    perp_vals = range(-width, width, length=n_perp)
    Δperp = step(perp_vals)

    if isnothing(n_par)
        len_width_ratio = half_length / width
        n_par = Int(round(n_perp * len_width_ratio))
    end

    par_vals = range(-half_length, half_length, length=n_par)
    Δpar = step(par_vals)

    # Main line points
    line = [midpoint .+ s_par .* û for s_par in par_vals]

    # 2D grid in global coordinates
    grid = [midpoint .+ s_par .* û .+ s_perp .* n̂ for s_par in par_vals, s_perp in perp_vals]

    return line, collect(perp_vals), collect(par_vals), Δpar, Δperp, grid
end



"""
    interpolate_to_line(data, x, y, z, line_pts)

Interpolate scalar 3D data [nx, ny, nz] to a list of 3D line points.

Returns a vector of length `n_line`.
"""
function interpolate_to_line(
    data::Array{Float32,3},
    x::Vector{Float64},
    y::Vector{Float64},
    z::Vector{Float64},
    line_pts::Vector{Vector{Float64}}
)
    has_x = length(x) > 1
    has_y = length(y) > 1
    has_z = length(z) > 1

    dx = has_x ? x[2] - x[1] : 0.0
    dy = has_y ? y[2] - y[1] : 0.0
    dz = has_z ? z[2] - z[1] : 0.0

    x_vals = has_x ? x : [x[1]]
    y_vals = has_y ? y : [y[1]]
    z_vals = has_z ? z : [z[1]]

    itp = if has_x && has_y && has_z
        interpolate((x_vals, y_vals, z_vals), data, Gridded(Linear()))
    elseif has_x && has_y
        interpolate((x_vals, y_vals), data[:, :, 1], Gridded(Linear()))
    elseif has_x && has_z
        interpolate((x_vals, z_vals), data[:, 1, :], Gridded(Linear()))
    elseif has_y && has_z
        interpolate((y_vals, z_vals), data[1, :, :], Gridded(Linear()))
    elseif has_x
        @warn "Only x-dimension is active — line interpolation is likely redundant"
        interpolate((x_vals,), data[:, 1, 1], Gridded(Linear()))
    elseif has_y
        @warn "Only y-dimension is active — line interpolation is likely redundant"
        interpolate((y_vals,), data[1, :, 1], Gridded(Linear()))
    elseif has_z
        @warn "Only z-dimension is active — line interpolation is likely redundant"
        interpolate((z_vals,), data[1, 1, :], Gridded(Linear()))
    else
        error("No active spatial dimensions.")
    end

    n_line = length(line_pts)
    result = Vector{Float64}(undef, n_line)

    for i in 1:n_line
        xq, yq, zq = line_pts[i]

        in_bounds = (!has_x || (xq ≥ first(x_vals) && xq ≤ last(x_vals))) &&
                    (!has_y || (yq ≥ first(y_vals) && yq ≤ last(y_vals))) &&
                    (!has_z || (zq ≥ first(z_vals) && zq ≤ last(z_vals)))

        if in_bounds
            if has_x && has_y && has_z
                result[i] = itp(xq, yq, zq)
            elseif has_x && has_y && !has_z
                result[i] = itp(xq, yq)
            elseif has_x && !has_y && has_z
                result[i] = itp(xq, zq)
            elseif !has_x && has_y && has_z
                result[i] = itp(yq, zq)
            elseif has_x && !has_y && !has_z
                result[i] = itp(xq)
            elseif !has_x && has_y && !has_z
                result[i] = itp(yq)
            elseif !has_x && !has_y && has_z
                result[i] = itp(zq)
            else
                error("No active spatial dimensions in the data.")
            end

        else
            result[i] = NaN
        end
    end

    return result
end

"""
    interpolate_to_line(data, x, y, z, line_pts)

Interpolate 4D data [nx, ny, nz, n_var] to line points.

Returns array of size [n_line, n_var].
"""
function interpolate_to_line(
    data::Array{Float32,4},
    x::Vector{Float64},
    y::Vector{Float64},
    z::Vector{Float64},
    line_pts::Vector{Vector{Float64}}
)
    nx, ny, nz, n_var = size(data)
    n_line = length(line_pts)
    result = Array{Float64}(undef, n_line, n_var)

    for var in 1:n_var
        result[:, var] = interpolate_to_line(data[:, :, :, var], x, y, z, line_pts)
    end

    return result
end

"""
    interpolate_to_line(data, x, y, z, line_pts)

Interpolate 5D data [nx, ny, nz, n_var, n_step] to line points.

Returns array of size [n_line, n_var, n_step].
"""
function interpolate_to_line(
    data::Array{Float32,5},
    x::Vector{Float64},
    y::Vector{Float64},
    z::Vector{Float64},
    line_pts::Vector{Vector{Float64}}
)
    nx, ny, nz, n_var, n_step = size(data)
    n_line = length(line_pts)
    result = Array{Float64}(undef, n_line, n_var, n_step)

    for step in 1:n_step
        result[:, :, step] = interpolate_to_line(data[:, :, :, :, step], x, y, z, line_pts)
    end

    return result
end


"""
    interpolate_to_plane(data, x, y, z, grid)

Interpolate scalar field data [nx, ny, nz] to a 2D grid of sampling points.

# Arguments
- `data::Array{Float32,3}`: Scalar field on Cartesian grid.
- `x`, `y`, `z`: Coordinate vectors.
- `grid::Array{Vector{Float64},2}`: 2D array of global coordinates, shape [n_par, n_perp].

# Returns
- `Array{Float64,2}`: Interpolated values at each grid point.
"""
function interpolate_to_plane(
    data::Array{Float32,3},
    x::Vector{Float64},
    y::Vector{Float64},
    z::Vector{Float64},
    grid::Array{Vector{Float64},2}
)
    has_x, has_y, has_z = length(x) > 1, length(y) > 1, length(z) > 1
    x_vals = has_x ? x : [x[1]]
    y_vals = has_y ? y : [y[1]]
    z_vals = has_z ? z : [z[1]]

    if has_x && has_y && has_z
        itp = interpolate((x_vals, y_vals, z_vals), data, Gridded(Linear()))
    elseif has_x && has_y
        itp = interpolate((x_vals, y_vals), data[:, :, 1], Gridded(Linear()))
    elseif has_x && has_z
       itp =  interpolate((x_vals, z_vals), data[:, 1, :], Gridded(Linear()))
    elseif has_y && has_z
       itp =  interpolate((y_vals, z_vals), data[1, :, :], Gridded(Linear()))
    elseif has_x
        @warn "Only x-dimension is active — plane interpolation is likely degenerate"
        itp = interpolate((x_vals,), data[:, 1, 1], Gridded(Linear()))
    elseif has_y
        @warn "Only y-dimension is active — plane interpolation is likely degenerate"
        itp = interpolate((y_vals,), data[1, :, 1], Gridded(Linear()))
    elseif has_z
        @warn "Only z-dimension is active — plane interpolation is likely degenerate"
        itp = interpolate((z_vals,), data[1, 1, :], Gridded(Linear()))
    else
        error("No active spatial dimensions.")
    end

    n_par, n_perp = size(grid)
    result = Array{Float64}(undef, n_par, n_perp)

    for i in 1:n_par, j in 1:n_perp
        xq, yq, zq = grid[i, j]

        in_bounds = (!has_x || (xq ≥ first(x_vals) && xq ≤ last(x_vals))) &&
                    (!has_y || (yq ≥ first(y_vals) && yq ≤ last(y_vals))) &&
                    (!has_z || (zq ≥ first(z_vals) && zq ≤ last(z_vals)))

        if in_bounds
            if has_x && has_y && has_z
                result[i, j] = itp(xq, yq, zq)
            elseif has_x && has_y
                result[i, j] = itp(xq, yq)
            elseif has_x && has_z
                result[i, j] = itp(xq, zq)
            elseif has_y && has_z
                result[i, j] = itp(yq, zq)
            elseif has_x
                result[i, j] = itp(xq)
            elseif has_y
                result[i, j] = itp(yq)
            elseif has_z
                result[i, j] = itp(zq)
            else
                error("No active spatial dimensions.")
            end
        else
            result[i, j] = NaN
        end
    end

    return result
end

function interpolate_to_plane(
    data::Array{Float32,4},
    x::Vector{Float64},
    y::Vector{Float64},
    z::Vector{Float64},
    grid::Array{Vector{Float64},2}
)
    nx, ny, nz, n_var = size(data)
    n_par, n_perp = size(grid)
    result = Array{Float64}(undef, n_par, n_perp, n_var)

    for var in 1:n_var
        result[:, :, var] = interpolate_to_plane(data[:, :, :, var], x, y, z, grid)
    end

    return result
end

function interpolate_to_plane(
    data::Array{Float32,5},
    x::Vector{Float64},
    y::Vector{Float64},
    z::Vector{Float64},
    grid::Array{Vector{Float64},2}
)
    nx, ny, nz, n_var, n_step = size(data)
    n_par, n_perp = size(grid)
    result = Array{Float64}(undef, n_par, n_perp, n_var, n_step)

    for step in 1:n_step
        result[:, :, :, step] = interpolate_to_plane(data[:, :, :, :, step], x, y, z, grid)
    end

    return result
end
