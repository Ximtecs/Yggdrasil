# New functions to locate patches intersected by a line

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
    patches = Patch_NML[]
    for patch in Snapshot_meta.PATCHES
        if patch.LEVEL == level && line_intersects_patch(start_point, end_point, patch)
            push!(patches, patch)
        end
    end
    return patches
end

