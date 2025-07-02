function get_xyz(Snapshot_meta :: Snapshot_metadata, drop_dims :: Bool)
    BOX = Snapshot_meta.SNAPSHOT.BOX
    ORIGIN = Snapshot_meta.SNAPSHOT.ORIGIN


    #------------ gets coordinates for lowest level patch
    level = Snapshot_meta.LEVELMIN
    found_level = false
    i = 1
    while (! found_level)
        patch_level = Snapshot_meta.PATCHES[i].LEVEL
        if patch_level == level
            found_level = true
        else
            i += 1
        end
    end
    #----------------------------
    ds = Snapshot_meta.PATCHES[i].DS






    
    mem_size = get_mem_size(Snapshot_meta)[1:3]

    x = zeros(mem_size...)
    y = zeros(mem_size...)
    z = zeros(mem_size...)


    #-------------- Gets the x, y, z coordinates of the cell centers
    x_range = collect(range(ORIGIN[1] + 0.5 * ds[1], stop=ORIGIN[1] + BOX[1] - 0.5*ds[1], length=mem_size[1]))
    y_range = collect(range(ORIGIN[2] + 0.5 * ds[2], stop=ORIGIN[2] + BOX[2] - 0.5*ds[2], length=mem_size[2]))
    z_range = collect(range(ORIGIN[3] + 0.5 * ds[3], stop=ORIGIN[3] + BOX[3] - 0.5*ds[3], length=mem_size[3]))

    for i in 1:mem_size[1]
        for j in 1:mem_size[2]
            for k in 1:mem_size[3]
                x[i,j,k] = x_range[i]
                y[i,j,k] = y_range[j]
                z[i,j,k] = z_range[k]
            end
        end
    end
            
    if drop_dims
        x = drop_unit_dims(x)
        y = drop_unit_dims(y)
        z = drop_unit_dims(z)
    end


    return x, y, z, ds
end 


function get_xyz(Snapshot_meta :: Snapshot_metadata, drop_dims :: Bool, level :: Int)
    BOX = Snapshot_meta.SNAPSHOT.BOX
    ORIGIN = Snapshot_meta.SNAPSHOT.ORIGIN

    #-------- assumes no AMR 

    found_level = false
    i = 1
    while (! found_level)
        patch_level = Snapshot_meta.PATCHES[i].LEVEL
        if patch_level == level
            found_level = true
        else
            i += 1
        end
    end
    ds = Snapshot_meta.PATCHES[i].DS

    mem_size = get_mem_size(Snapshot_meta, level)[1:3]

    x = zeros(mem_size...)
    y = zeros(mem_size...)
    z = zeros(mem_size...)


    #-------------- Gets the x, y, z coordinates of the cell centers
    if mem_size[1] > 1
    x_range = collect(range(ORIGIN[1] + 0.5 * ds[1], stop=ORIGIN[1] + BOX[1] - 0.5*ds[1], length=mem_size[1]))
    else 
        x_range = [ORIGIN[1] + 0.5 * ds[1]]
    end
    if mem_size[2] > 1
    y_range = collect(range(ORIGIN[2] + 0.5 * ds[2], stop=ORIGIN[2] + BOX[2] - 0.5*ds[2], length=mem_size[2]))
    else
        y_range = [ORIGIN[2] + 0.5 * ds[2]]
    end
    if mem_size[3] > 1
    z_range = collect(range(ORIGIN[3] + 0.5 * ds[3], stop=ORIGIN[3] + BOX[3] - 0.5*ds[3], length=mem_size[3]))
    else
        z_range = [ORIGIN[3] + 0.5 * ds[3]]
    end


    for i in 1:mem_size[1]
        for j in 1:mem_size[2]
            for k in 1:mem_size[3]
                x[i,j,k] = x_range[i]
                y[i,j,k] = y_range[j]
                z[i,j,k] = z_range[k]
            end
        end
    end
            
    if drop_dims
        x = drop_unit_dims(x)
        y = drop_unit_dims(y)
        z = drop_unit_dims(z)
    end


    return x, y, z, ds
end 


function get_xyz(
    Snapshot_meta::Snapshot_metadata,
    drop_dims::Bool,
    level::Int;
    llc::Union{Nothing,Vector{Int}}=nothing,
    urc::Union{Nothing,Vector{Int}}=nothing
)
    BOX = Snapshot_meta.SNAPSHOT.BOX
    ORIGIN = Snapshot_meta.SNAPSHOT.ORIGIN

    #-------- assumes no AMR 
    found_level = false
    i = 1
    while !found_level
        patch_level = Snapshot_meta.PATCHES[i].LEVEL
        if patch_level == level
            found_level = true
        else
            i += 1
        end
    end
    ds = Snapshot_meta.PATCHES[i].DS

    mem_size = get_mem_size(Snapshot_meta, level)[1:3]

    # Default to full domain
    llc = llc === nothing ? [1, 1, 1] : llc
    urc = urc === nothing ? mem_size : urc

    # Construct coordinate ranges
    x_range = mem_size[1] > 1 ? collect(range(ORIGIN[1] + 0.5 * ds[1], stop=ORIGIN[1] + BOX[1] - 0.5 * ds[1], length=mem_size[1])) :
                                [ORIGIN[1] + 0.5 * ds[1]]
    y_range = mem_size[2] > 1 ? collect(range(ORIGIN[2] + 0.5 * ds[2], stop=ORIGIN[2] + BOX[2] - 0.5 * ds[2], length=mem_size[2])) :
                                [ORIGIN[2] + 0.5 * ds[2]]
    z_range = mem_size[3] > 1 ? collect(range(ORIGIN[3] + 0.5 * ds[3], stop=ORIGIN[3] + BOX[3] - 0.5 * ds[3], length=mem_size[3])) :
                                [ORIGIN[3] + 0.5 * ds[3]]

    # Slice subregions
    x_sub = x_range[llc[1]:urc[1]]
    y_sub = y_range[llc[2]:urc[2]]
    z_sub = z_range[llc[3]:urc[3]]

    # Allocate arrays
    nx, ny, nz = length(x_sub), length(y_sub), length(z_sub)
    x = zeros(nx, ny, nz)
    y = zeros(nx, ny, nz)
    z = zeros(nx, ny, nz)

    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                x[i,j,k] = x_sub[i]
                y[i,j,k] = y_sub[j]
                z[i,j,k] = z_sub[k]
            end
        end
    end

    if drop_dims
        x = drop_unit_dims(x)
        y = drop_unit_dims(y)
        z = drop_unit_dims(z)
    end

    return x, y, z, ds
end


function get_ds(
    Snapshot_meta::Snapshot_metadata,
    level::Int
)
    found_level = false
    i = 1
    while !found_level
        patch_level = Snapshot_meta.PATCHES[i].LEVEL
        if patch_level == level
            found_level = true
        else
            i += 1
        end
    end

    if i > length(Snapshot_meta.PATCHES)
        error("Level $level not found in Snapshot_metadata.")
    end
    ds = Snapshot_meta.PATCHES[i].DS
    return ds
end 

function get_ds(
    Snapshot_meta::Snapshot_metadata
)
    found_level = false
    level = Snapshot_meta.LEVELMIN
    i = 1
    while !found_level
        patch_level = Snapshot_meta.PATCHES[i].LEVEL
        if patch_level == level
            found_level = true
        else
            i += 1
        end
    end
    if i > length(Snapshot_meta.PATCHES)
        error("Level $level not found in Snapshot_metadata.")
    end
    ds = Snapshot_meta.PATCHES[i].DS
    return ds
end 

function get_patch_float_size(
    Snapshot_meta::Snapshot_metadata,
    level::Int
)
    ds = get_ds(Snapshot_meta, level)
    patch_size = Snapshot_meta.SNAPSHOT.N .* ds
    return patch_size
end

function get_patch_float_size(
    Snapshot_meta::Snapshot_metadata
)
    ds = get_ds(Snapshot_meta)
    patch_size = Snapshot_meta.SNAPSHOT.N .* ds
    return patch_size
end


"""
# Function to get the area of a patch in a snapshot
# Area will be aligned with the patches at the given level
# Arguments:
# - `Snapshot_meta`: Metadata of the snapshot
# - `midpoint`: Midpoint of the area in the format [x, y, z]
# - `patches`: Number of patches to extend in each direction
#             fx. 4 meants 8x8 patches in total (Maybe be +-1 as it rounds to nearest patch)
# - `level`: Level of the snapshot
"""
function get_area( Snapshot_meta::Snapshot_metadata,
    midpoint:: Vector{Float64},
    patches :: Int,
    level ::Int)
    

    #------------------ check input arguments ------------------
    @assert length(midpoint) == 3 "midpoint must have 3 values (x, y, z)"
    @assert patches > 0 "patches must be a positive integer"
    #------------------------------------------------------------

    
    ds = get_ds(Snapshot_meta, level)
    patch_size = get_patch_float_size(Snapshot_meta, level)

    area_llc = midpoint .- (patch_size * patches)
    area_urc = midpoint .+ (patch_size * patches)


    area_llc_fixed = floor.(area_llc ./ patch_size) .* patch_size .+ ds
    area_urc_fixed = ceil.(area_urc ./ patch_size) .* patch_size

    llc_int = [Int(area_llc_fixed[1] / ds[1]),
               Int(area_llc_fixed[2] / ds[2]),
               Int(area_llc_fixed[3] / ds[3])]

    urc_int = [Int(area_urc_fixed[1] / ds[1]),
               Int(area_urc_fixed[2] / ds[2]),
               Int(area_urc_fixed[3] / ds[3])]


    #------------- if negative integer index, set to 1 -------------
    for i in 1:3
        if llc_int[i] <= 0 || urc_int[i] <= 0
            llc_int[i] = 1
            urc_int[i] = 1

            area_llc_fixed[i] = 0.0
            area_urc_fixed[i] = 0.0
        end
    end 
    #---------------------------------------------------------------

    return llc_int, urc_int, area_llc_fixed, area_urc_fixed
end


"""
    get_area_from_size(Snapshot_meta, midpoint, size, level)

Compute grid-aligned area bounds given a physical size (in code units).

# Arguments:
- `Snapshot_meta::Snapshot_metadata`: Metadata of the snapshot.
- `midpoint::Vector{Float64}`: Centre of the area, e.g., `[x, y, z]`.
- `size::Vector{Float64}`: Total size of the area in each dimension (in code units).
- `level::Int`: Refinement level to use.

# Returns:
- `llc_int::Vector{Int}`: Lower-left corner as integer grid indices.
- `urc_int::Vector{Int}`: Upper-right corner as integer grid indices.
- `area_llc_fixed::Vector{Float64}`: Adjusted lower-left corner (floored to grid).
- `area_urc_fixed::Vector{Float64}`: Adjusted upper-right corner (ceiled to grid).
"""
function get_area(
    Snapshot_meta::Snapshot_metadata,
    midpoint::Vector{Float64},
    size::Vector{Float64},
    level::Int
)
    @assert length(midpoint) == 3 "midpoint must have 3 values (x, y, z)"
    @assert length(size) == 3 "size must have 3 values (dx, dy, dz)"
    @assert all(size .>= 0) "size must contain only positive values"

    ds = get_ds(Snapshot_meta, level)
    patch_size = get_patch_float_size(Snapshot_meta, level)

    half_size = 0.5 .* size

    area_llc = midpoint .- half_size
    area_urc = midpoint .+ half_size

    area_llc_fixed = floor.(area_llc ./ patch_size) .* patch_size .+ ds
    area_urc_fixed = ceil.(area_urc ./ patch_size) .* patch_size

    llc_int = [Int(area_llc_fixed[1] / ds[1]),
               Int(area_llc_fixed[2] / ds[2]),
               Int(area_llc_fixed[3] / ds[3])]

    urc_int = [Int(area_urc_fixed[1] / ds[1]),
               Int(area_urc_fixed[2] / ds[2]),
               Int(area_urc_fixed[3] / ds[3])]

    for i in 1:3
        if llc_int[i] <= 0 || urc_int[i] <= 0
            llc_int[i] = 1
            urc_int[i] = 1
            area_llc_fixed[i] = 0.0
            area_urc_fixed[i] = 0.0
        end
    end

    return llc_int, urc_int, area_llc_fixed, area_urc_fixed
end


"""
    get_area_from_size(Snapshot_meta, midpoint, size, level)

Compute grid-aligned area bounds given a physical size (in code units).

# Arguments:
- `Snapshot_meta::Snapshot_metadata`: Metadata of the snapshot.
- `midpoint::Vector{Float64}`: Centre of the area, e.g., `[x, y, z]`.
- `size::Vector{Float64}`: Total size of the area in each dimension (in code units).
- `level::Int`: Refinement level to use.

# Returns:
- `llc_int::Vector{Int}`: Lower-left corner as integer grid indices.
- `urc_int::Vector{Int}`: Upper-right corner as integer grid indices.
- `area_llc_fixed::Vector{Float64}`: Adjusted lower-left corner (floored to grid).
- `area_urc_fixed::Vector{Float64}`: Adjusted upper-right corner (ceiled to grid).
"""
function get_area(
    Snapshot_meta::Snapshot_metadata,
    midpoint::Vector{Float64},
    size::Float64,
    level::Int
)
    @assert length(midpoint) == 3 "midpoint must have 3 values (x, y, z)"
    @assert size .>= 0 "size must contain only positive values"

    ds = get_ds(Snapshot_meta, level)
    patch_size = get_patch_float_size(Snapshot_meta, level)

    half_size = 0.5 * size

    area_llc = midpoint .- half_size
    area_urc = midpoint .+ half_size

    area_llc_fixed = floor.(area_llc ./ patch_size) .* patch_size .+ ds
    area_urc_fixed = ceil.(area_urc ./ patch_size) .* patch_size

    llc_int = [Int(area_llc_fixed[1] / ds[1]),
               Int(area_llc_fixed[2] / ds[2]),
               Int(area_llc_fixed[3] / ds[3])]

    urc_int = [Int(area_urc_fixed[1] / ds[1]),
               Int(area_urc_fixed[2] / ds[2]),
               Int(area_urc_fixed[3] / ds[3])]

    for i in 1:3
        if llc_int[i] <= 0 || urc_int[i] <= 0
            llc_int[i] = 1
            urc_int[i] = 1
            area_llc_fixed[i] = 0.0
            area_urc_fixed[i] = 0.0
        end
    end

    return llc_int, urc_int, area_llc_fixed, area_urc_fixed
end