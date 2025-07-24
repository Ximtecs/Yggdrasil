# Functions for locating patches intersected by squares and cubes

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

