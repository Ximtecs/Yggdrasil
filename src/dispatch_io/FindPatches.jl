# Functions for locating patches intersected by squares, cubes, or lines

"""
    square_contains_patch(points, patch)

Return `true` if the axis-aligned square defined by `points` fully contains the
patch in the two in-plane directions and intersects the square's plane.
"""
function square_contains_patch(
    points::Union{AbstractMatrix{<:AbstractFloat},Vector{<:AbstractVector{<:AbstractFloat}}},
    patch::Patch_NML,
)
    mat = _convert_points(points, 4)
    mins, maxs = _bounding_box(mat)
    ranges = maxs .- mins
    naxis = argmin(ranges)
    plane = mat[1, naxis]
    if patch.LLC_NAT[naxis] > plane || patch.LLC_NAT[naxis] + patch.SIZE[naxis] < plane
        return false
    end
    for ax in setdiff(1:3, naxis)
        pmin = patch.LLC_NAT[ax]
        pmax = patch.LLC_NAT[ax] + patch.SIZE[ax]
        if pmin < mins[ax] || pmax > maxs[ax]
            return false
        end
    end
    return true
end

"""
    cube_contains_patch(points, patch)

Return `true` if the axis-aligned cube defined by `points` fully contains the
patch.
"""
function cube_contains_patch(
    points::Union{AbstractMatrix{<:AbstractFloat},Vector{<:AbstractVector{<:AbstractFloat}}},
    patch::Patch_NML,
)
    mat = _convert_points(points, 8)
    mins, maxs = _bounding_box(mat)
    pmin = patch.LLC_NAT
    pmax = patch.LLC_NAT .+ patch.SIZE
    for i in 1:3
        if pmin[i] < mins[i] || pmax[i] > maxs[i]
            return false
        end
    end
    return true
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
    mins = [minimum(col) for col in eachcol(mat)]
    maxs = [maximum(col) for col in eachcol(mat)]
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
    find_patches_square(Snapshot_meta, points; level, all_levels=false)

Return patches for which the square defined by `points` fully contains the patch
in-plane and intersects the square's plane. If `all_levels` is `true` patches from
all refinement levels are considered. Otherwise, the search is restricted to the
specified `level` (defaults to `Snapshot_meta.LEVELMIN`). `points` can be a 4×3
matrix or a vector of four 3-component vectors.
"""
function find_patches_square(
    Snapshot_meta::Snapshot_metadata,
    points::Union{AbstractMatrix{<:AbstractFloat},Vector{<:AbstractVector{<:AbstractFloat}}};
    level::Union{Int,Nothing}=nothing,
    all_levels::Bool=false,
)
    mat = _convert_points(points, 4)
    lvl = level === nothing ? Snapshot_meta.LEVELMIN : level
    patches = Patch_NML[]
    for patch in Snapshot_meta.PATCHES
        if (all_levels || patch.LEVEL == lvl) && square_contains_patch(mat, patch)
            push!(patches, patch)
        end
    end
    return patches
end

"""
    find_patches_cube(Snapshot_meta, points; level, all_levels=false)

Return patches for which the cube defined by `points` fully contains the patch. If
`all_levels` is `true` patches from all refinement levels are considered. Otherwise,
the search is restricted to the specified `level` (defaults to
`Snapshot_meta.LEVELMIN`). `points` can be an 8×3 matrix or a vector of eight
3-component vectors.
"""
function find_patches_cube(
    Snapshot_meta::Snapshot_metadata,
    points::Union{AbstractMatrix{<:AbstractFloat},Vector{<:AbstractVector{<:AbstractFloat}}};
    level::Union{Int,Nothing}=nothing,
    all_levels::Bool=false,
)
    mat = _convert_points(points, 8)
    lvl = level === nothing ? Snapshot_meta.LEVELMIN : level
    patches = Patch_NML[]
    for patch in Snapshot_meta.PATCHES
        if (all_levels || patch.LEVEL == lvl) && cube_contains_patch(mat, patch)
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
    find_patches_line(Snapshot_meta, start_point, end_point; level, all_levels=false)

Return patches in `Snapshot_meta` for which the line segment `start_point` →
`end_point` intersects the patch. If `all_levels` is `true` patches from all
refinement levels are considered. Otherwise, the search is restricted to the
specified `level` (defaults to `Snapshot_meta.LEVELMIN`).
"""
function find_patches_line(
    Snapshot_meta::Snapshot_metadata,
    start_point::AbstractVector{<:AbstractFloat},
    end_point::AbstractVector{<:AbstractFloat};
    level::Union{Int,Nothing}=nothing,
    all_levels::Bool=false,
)
    @assert length(start_point) == 3 "start_point must have 3 values (x, y, z)"
    @assert length(end_point) == 3 "end_point must have 3 values (x, y, z)"
    lvl = level === nothing ? Snapshot_meta.LEVELMIN : level
    patches = Patch_NML[]
    for patch in Snapshot_meta.PATCHES
        if (all_levels || patch.LEVEL == lvl) && line_intersects_patch(start_point, end_point, patch)
            push!(patches, patch)
        end
    end
    return patches
end

