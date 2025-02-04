#module VectorOps

# Import backend-specific implementations
include("Prefactors.jl")
include("VectorCPU.jl")
include("VectorGPU.jl")


#export Backend, CPU, GPU, 
#      ddx_up, ddy_up, ddz_up, ddx_dn, ddy_dn, ddz_dn,
#      curl_up, curl_dn

# Define backend types
abstract type Backend end
struct CPU <: Backend end
struct GPU <: Backend end


function curl_up(backend::Backend, field::AbstractArray{<:AbstractFloat, 4}, 
    dx::AbstractFloat,dy::AbstractFloat,dz::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.curl_up(field, dx, dy, dz, order)
    elseif backend isa GPU
        return GPUVectorOps.curl_up(field, dx, dy, dz, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function curl_dn(backend::Backend, field::AbstractArray{<:AbstractFloat, 4}, 
    dx::AbstractFloat,dy::AbstractFloat,dz::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.curl_dn(field, dx, dy, dz, order)
    elseif backend isa GPU
        return GPUVectorOps.curl_dn(field, dx, dy, dz, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function curl_up(backend::Backend, field::AbstractArray{<:AbstractFloat, 4}, result::AbstractArray{<:AbstractFloat, 4}, 
    dx::AbstractFloat,dy::AbstractFloat,dz::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.curl_up(field, result, dx, dy, dz, order)
    elseif backend isa GPU
        return GPUVectorOps.curl_up(field, result, dx, dy, dz, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function curl_dn(backend::Backend, field::AbstractArray{<:AbstractFloat, 4}, result::AbstractArray{<:AbstractFloat, 4}, 
    dx::AbstractFloat,dy::AbstractFloat,dz::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.curl_dn(field, result, dx, dy, dz, order)
    elseif backend isa GPU
        return GPUVectorOps.curl_dn(field, result, dx, dy, dz, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end



# Exported functions
function ddx_up(backend::Backend, field::AbstractArray, dx::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddx_up(field, dx, order)
    elseif backend isa GPU
        return GPUVectorOps.ddx_up(field, dx, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function ddy_up(backend::Backend, field::AbstractArray, dy::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddy_up(field, dy, order)
    elseif backend isa GPU
        return GPUVectorOps.ddy_up(field, dy, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end 
end

function ddz_up(backend::Backend, field::AbstractArray, dz::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddz_up(field, dz, order)
    elseif backend isa GPU
        return GPUVectorOps.ddz_up(field, dz, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

# Exported functions
function ddx_dn(backend::Backend, field::AbstractArray, dx::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddx_dn(field, dx, order)
    elseif backend isa GPU
        return GPUVectorOps.ddx_dn(field, dx, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function ddy_dn(backend::Backend, field::AbstractArray, dy::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddy_dn(field, dy, order)
    elseif backend isa GPU
        return GPUVectorOps.ddy_dn(field, dy, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function ddz_dn(backend::Backend, field::AbstractArray, dz::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddz_dn(field, dz, order)
    elseif backend isa GPU
        return GPUVectorOps.ddz_dn(field, dz, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end
function ddx_up(backend::Backend, field::AbstractArray, result::AbstractArray{<:AbstractFloat, 4}, 
    dx::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddx_up(field, result, dx, order)
    elseif backend isa GPU
        return GPUVectorOps.ddx_up(field, result, dx, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function ddy_up(backend::Backend, field::AbstractArray, result::AbstractArray{<:AbstractFloat, 4}, 
    dy::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddy_up(field, result, dy, order)
    elseif backend isa GPU
        return GPUVectorOps.ddy_up(field, result, dy, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function ddz_up(backend::Backend, field::AbstractArray, result::AbstractArray{<:AbstractFloat, 4}, 
    dz::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddz_up(field, result, dz, order)
    elseif backend isa GPU
        return GPUVectorOps.ddz_up(field, result, dz, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function ddx_dn(backend::Backend, field::AbstractArray, result::AbstractArray{<:AbstractFloat, 4}, 
    dx::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddx_dn(field, result, dx, order)
    elseif backend isa GPU
        return GPUVectorOps.ddx_dn(field, result, dx, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function ddy_dn(backend::Backend, field::AbstractArray, result::AbstractArray{<:AbstractFloat, 4}, 
    dy::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddy_dn(field, result, dy, order)
    elseif backend isa GPU
        return GPUVectorOps.ddy_dn(field, result, dy, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function ddz_dn(backend::Backend, field::AbstractArray, result::AbstractArray{<:AbstractFloat, 4},
     dz::AbstractFloat, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddz_dn(field, result, dz, order)
    elseif backend isa GPU
        return GPUVectorOps.ddz_dn(field, result, dz, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

#function curl(::Type{<:Backend}, field, spacing) end
#function divergence(::Type{<:Backend}, field, spacing) end

#end