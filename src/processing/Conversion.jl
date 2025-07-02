
function mom_to_vel!(data::AbstractArray, idx_nml::IDX_NML; var_idx::Union{Nothing, Int}=nothing)
    nd = ndims(data)
    vidx = isnothing(var_idx) ? (nd == 5 ? nd - 1 : nd) : var_idx

    # Variable indices from name list
    ρ_idx  = get_idx_value(idx_nml, "rho")
    px_idx = get_idx_value(idx_nml, "px")
    py_idx = get_idx_value(idx_nml, "py")
    pz_idx = get_idx_value(idx_nml, "pz")

    # Construct slicing tuples for each variable
    slicer_ρ  = ntuple(i -> i == vidx ? ρ_idx  : Colon(), nd)
    slicer_px = ntuple(i -> i == vidx ? px_idx : Colon(), nd)
    slicer_py = ntuple(i -> i == vidx ? py_idx : Colon(), nd)
    slicer_pz = ntuple(i -> i == vidx ? pz_idx : Colon(), nd)

    ρ = view(data, slicer_ρ...)
    @views data[slicer_px...] ./= ρ
    @views data[slicer_py...] ./= ρ
    @views data[slicer_pz...] ./= ρ

    return data
end

function vel_to_mom!(data::AbstractArray, idx_nml::IDX_NML; var_idx::Union{Nothing, Int}=nothing)
    nd = ndims(data)
    vidx = isnothing(var_idx) ? (nd == 5 ? nd - 1 : nd) : var_idx

    ρ_idx  = get_idx_value(idx_nml, "rho")
    px_idx = get_idx_value(idx_nml, "px")
    py_idx = get_idx_value(idx_nml, "py")
    pz_idx = get_idx_value(idx_nml, "pz")

    slicer_ρ  = ntuple(i -> i == vidx ? ρ_idx  : Colon(), nd)
    slicer_px = ntuple(i -> i == vidx ? px_idx : Colon(), nd)
    slicer_py = ntuple(i -> i == vidx ? py_idx : Colon(), nd)
    slicer_pz = ntuple(i -> i == vidx ? pz_idx : Colon(), nd)

    ρ = view(data, slicer_ρ...)
    @views data[slicer_px...] .*= ρ
    @views data[slicer_py...] .*= ρ
    @views data[slicer_pz...] .*= ρ

    return data
end


"""
    rotate_aligned_vectors!(aligned_data, û, n̂, idx_nml; component_order = (:par, :perp, :norm), var_idx = nothing)

Rotate vector components in a data array from global Cartesian coordinates to a line-aligned frame defined by:
- `û`: unit vector along the line ("parallel")
- `n̂`: unit vector normal to the plane ("normal")
- The third direction is the perpendicular direction within the plane ("perpendicular")

# Arguments
- `aligned_data::AbstractArray`: Data array (can be 3D, 4D, or 5D) containing vector components.
- `û::Vector{Float64}`: 3D unit vector in the parallel direction.
- `n̂::Vector{Float64}`: 3D unit vector normal to the plane.
- `idx_nml::IDX_NML`: Object that maps variable names to array indices.

# Optional
- `component_order::NTuple{3,Symbol}`: Permutation of (:par, :perp, :norm) to define how rotated components are assigned to (x, y, z) respectively. Default: `(:par, :perp, :norm)`.
- `var_idx::Union{Nothing, Int}`: Dimension along which variables are stored. Defaults to last dimension for ≤3D, second-to-last for ≥4D arrays.

# Notes
- The function modifies `aligned_data` in-place.
- If a variable group (e.g. `("bx", "by", "bz")`) is not found in `idx_nml`, it is skipped with a warning.
- Variable groups will also be skipped if their indices exceed the number of variables in `aligned_data`.

"""
function rotate_aligned_vectors!(
    aligned_data::AbstractArray,
    û::Vector{Float64},
    n̂::Vector{Float64},
    idx_nml::IDX_NML;
    component_order::NTuple{3,Symbol} = (:par, :perp, :norm),
    var_idx::Union{Nothing, Int} = nothing
)
    # --- Validate input vectors ---
    if length(û) != 3 || length(n̂) != 3
        error("û and n̂ must be 3D vectors")
    end

    if Set(component_order) != Set((:par, :perp, :norm))
        error("component_order must be a permutation of (:par, :perp, :norm)")
    end

    # --- Determine variable index ---
    nd = ndims(aligned_data)
    vidx = isnothing(var_idx) ? (nd >= 4 ? nd - 1 : nd) : var_idx

    n_var = size(aligned_data, vidx)

    # --- Normalise basis vectors ---
    b̂_raw = cross(n̂, û)
    if norm(b̂_raw) ≈ 0
        error("û and n̂ must not be parallel")
    end

    b̂ = b̂_raw / norm(b̂_raw)
    û = û / norm(û)
    n̂ = n̂ / norm(n̂)

    # Rotation matrix: rows = new basis vectors
    R = transpose(hcat(û, b̂, n̂))  # 3x3 matrix

    # Variable groups to rotate
    vector_groups = [
        ("px", "py", "pz"),   # momentum / velocity
        ("bx", "by", "bz"),   # magnetic field
        ("ex", "ey", "ez"),   # electric field
        ("jx", "jy", "jz"),   # current density
        ("v_e_x", "v_e_y", "v_e_z"),  # electron velocity
        ("v_p_x", "v_p_y", "v_p_z")   # proton velocity
    ]

    for (xkey, ykey, zkey) in vector_groups
        try
            ix = get_idx_value(idx_nml, xkey)
            iy = get_idx_value(idx_nml, ykey)
            iz = get_idx_value(idx_nml, zkey)

            if ix > n_var || iy > n_var || iz > n_var
                continue
            end

            slicer_x = ntuple(i -> i == vidx ? ix : Colon(), nd)
            slicer_y = ntuple(i -> i == vidx ? iy : Colon(), nd)
            slicer_z = ntuple(i -> i == vidx ? iz : Colon(), nd)

            vx = @view aligned_data[slicer_x...]
            vy = @view aligned_data[slicer_y...]
            vz = @view aligned_data[slicer_z...]

            v_par  = R[1,1]*vx .+ R[1,2]*vy .+ R[1,3]*vz
            v_perp = R[2,1]*vx .+ R[2,2]*vy .+ R[2,3]*vz
            v_norm = R[3,1]*vx .+ R[3,2]*vy .+ R[3,3]*vz

            component_map = Dict(:par => v_par, :perp => v_perp, :norm => v_norm)

            aligned_data[slicer_x...] .= component_map[component_order[1]]
            aligned_data[slicer_y...] .= component_map[component_order[2]]
            aligned_data[slicer_z...] .= component_map[component_order[3]]

        catch e
            @warn "Skipping rotation of ($xkey, $ykey, $zkey): $(e.msg)"
        end
    end
end