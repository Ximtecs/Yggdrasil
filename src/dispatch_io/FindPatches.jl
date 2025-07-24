# Functions for locating patches intersected by squares, cubes, or lines

"""
    square_intersects_patch(points, patch)

Return `true` if the axis-aligned bounding box formed by `points` intersects the
bounding box of `patch`.
"""
function square_intersects_patch(
    points::Union{AbstractMatrix{<:AbstractFloat},Vector{<:AbstractVector{<:AbstractFloat}}},
    patch::Patch_NML,
)
    mat = _convert_points(points, 4)
    mins, maxs = _bounding_box(mat)
    return _boxes_intersect(mins, maxs, patch)
end

"""
    cube_intersects_patch(points, patch)

Return `true` if the axis-aligned bounding box formed by `points` intersects the
bounding box of `patch`.
"""
function cube_intersects_patch(
    points::Union{AbstractMatrix{<:AbstractFloat},Vector{<:AbstractVector{<:AbstractFloat}}},
    patch::Patch_NML,
)
    mat = _convert_points(points, 8)
    mins, maxs = _bounding_box(mat)
    return _boxes_intersect(mins, maxs, patch)
end

function _convert_points(points, expected_rows)
    if points isa AbstractMatrix
        @assert size(points) == (expected_rows, 3)
        mat = Array{Float64}(points)
    elseif points isa AbstractVector
        @assert length(points) == expected_rows
        mat = hcat(points...)'
        @assert size(mat) == (expected_rows, 3)
        mat = Array{Float64}(mat)
    else
        error("points must be a matrix or vector of vectors")
    end
    return mat
end

function _bounding_box(mat)
    mins = vec(mapreduce(min, dims=1, mat))
    maxs = vec(mapreduce(max, dims=1, mat))
    return mins, maxs
end

function _boxes_intersect(mins, maxs, patch::Patch_NML)
    pmin = patch.LLC_NAT
    pmax = patch.LLC_NAT .+ patch.SIZE
    for i in 1:3
        if maxs[i] < pmin[i] || mins[i] > pmax[i]
            return false
        end
    end
    return true
end

"""
    find_patches_square(Snapshot_meta, points)

Return all patches for which the square defined by `points` intersects the patch.
`points` can be a 4×3 matrix or a vector of four 3-component vectors.
"""
function find_patches_square(
    Snapshot_meta::Snapshot_metadata,
    points::Union{AbstractMatrix{<:AbstractFloat},Vector{<:AbstractVector{<:AbstractFloat}}},
)
    mat = _convert_points(points, 4)
    patches = Patch_NML[]
    for patch in Snapshot_meta.PATCHES
        if square_intersects_patch(mat, patch)
            push!(patches, patch)
        end
    end
    return patches
end

"""
    find_patches_square(Snapshot_meta, points, level)

As above but restricts the search to patches at a given refinement `level`.
"""
function find_patches_square(
    Snapshot_meta::Snapshot_metadata,
    points::Union{AbstractMatrix{<:AbstractFloat},Vector{<:AbstractVector{<:AbstractFloat}}},
    level::Int,
)
    mat = _convert_points(points, 4)
    patches = Patch_NML[]
    for patch in Snapshot_meta.PATCHES
        if patch.LEVEL == level && square_intersects_patch(mat, patch)
            push!(patches, patch)
        end
    end
    return patches
end

"""
    find_patches_cube(Snapshot_meta, points)

Return all patches for which the cube defined by `points` intersects the patch.
`points` can be a 8×3 matrix or a vector of eight 3-component vectors.
"""
function find_patches_cube(
    Snapshot_meta::Snapshot_metadata,
    points::Union{AbstractMatrix{<:AbstractFloat},Vector{<:AbstractVector{<:AbstractFloat}}},
)
    mat = _convert_points(points, 8)
    patches = Patch_NML[]
    for patch in Snapshot_meta.PATCHES
        if cube_intersects_patch(mat, patch)
            push!(patches, patch)
        end
    end
    return patches
end

"""
    find_patches_cube(Snapshot_meta, points, level)

As above but restricts the search to patches at a given refinement `level`.
"""
function find_patches_cube(
    Snapshot_meta::Snapshot_metadata,
    points::Union{AbstractMatrix{<:AbstractFloat},Vector{<:AbstractVector{<:AbstractFloat}}},
    level::Int,
)
    mat = _convert_points(points, 8)
    patches = Patch_NML[]
    for patch in Snapshot_meta.PATCHES
        if patch.LEVEL == level && cube_intersects_patch(mat, patch)
            push!(patches, patch)
        end
    end
    return patches
end

# ---------------------------------------------------------------------
# Line intersection utilities

"""
    line_intersects_patch(start_point, end_point, patch)

Return `true` if the line segment defined by `start_point` and `end_point`
intersects the axis-aligned bounding box of `patch`.
"""
function line_intersects_patch(
    start_point::AbstractVector{<:AbstractFloat},
    end_point::AbstractVector{<:AbstractFloat},
    patch::Patch_NML,
)
    pmin = patch.LLC_NAT
    pmax = patch.LLC_NAT .+ patch.SIZE

    dir = end_point .- start_point
    tmin = 0.0
    tmax = 1.0

    for i in 1:3
        if abs(dir[i]) < eps(Float64)
            if start_point[i] < pmin[i] || start_point[i] > pmax[i]
                return false
            end
        else
            invd = 1.0 / dir[i]
            t1 = (pmin[i] - start_point[i]) * invd
            t2 = (pmax[i] - start_point[i]) * invd
            tlow = min(t1, t2)
            thigh = max(t1, t2)
            tmin = max(tmin, tlow)
            tmax = min(tmax, thigh)
            if tmin > tmax
                return false
            end
        end
    end

    return true
end

"""
    find_patches_line(Snapshot_meta, start_point, end_point)

Return all patches in `Snapshot_meta` for which the line segment
`start_point` → `end_point` intersects the patch.
"""
function find_patches_line(
    Snapshot_meta::Snapshot_metadata,
    start_point::AbstractVector{<:AbstractFloat},
    end_point::AbstractVector{<:AbstractFloat},
)
    @assert length(start_point) == 3 "start_point must have 3 values (x, y, z)"
    @assert length(end_point) == 3 "end_point must have 3 values (x, y, z)"
    patches = Patch_NML[]
    for patch in Snapshot_meta.PATCHES
        if line_intersects_patch(start_point, end_point, patch)
            push!(patches, patch)
        end
    end
    return patches
end

"""
    find_patches_line(Snapshot_meta, start_point, end_point, level)

As above but restricts the search to patches at a given refinement `level`.
"""
function find_patches_line(
    Snapshot_meta::Snapshot_metadata,
    start_point::AbstractVector{<:AbstractFloat},
    end_point::AbstractVector{<:AbstractFloat},
    level::Int,
)
    @assert length(start_point) == 3 "start_point must have 3 values (x, y, z)"
    @assert length(end_point) == 3 "end_point must have 3 values (x, y, z)"
    patches = Patch_NML[]
    for patch in Snapshot_meta.PATCHES
        if patch.LEVEL == level && line_intersects_patch(start_point, end_point, patch)
            push!(patches, patch)
        end
    end
    return patches
end

