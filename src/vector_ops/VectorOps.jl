module VectorOps

# Import backend-specific implementations
include("VectorCPU.jl")
include("VectorGPU.jl")


export Backend, CPU, GPU, 
      ddx_up, ddy_up, ddz_up, ddx_dn, ddy_dn, ddz_dn

# Define backend types
abstract type Backend end
struct CPU <: Backend end
struct GPU <: Backend end



# Exported functions
function ddx_up(backend::Backend, field, dx, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddx_up(field, dx, order)
    elseif backend isa GPU
        return GPUVectorOps.ddx_up(field, dx, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function ddy_up(backend::Backend, field, dx, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddy_up(field, dx, order)
    elseif backend isa GPU
        return GPUVectorOps.ddy_up(field, dx, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function ddz_up(backend::Backend, field, dx, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddz_up(field, dx, order)
    elseif backend isa GPU
        return GPUVectorOps.ddz_up(field, dx, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

# Exported functions
function ddx_dn(backend::Backend, field, dx, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddx_dn(field, dx, order)
    elseif backend isa GPU
        return GPUVectorOps.ddx_dn(field, dx, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function ddy_dn(backend::Backend, field, dx, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddy_dn(field, dx, order)
    elseif backend isa GPU
        return GPUVectorOps.ddy_dn(field, dx, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end

function ddz_dn(backend::Backend, field, dx, order::Int)
    if backend isa CPU
        return CPUVectorOps.ddz_dn(field, dx, order)
    elseif backend isa GPU
        return GPUVectorOps.ddz_dn(field, dx, order)
    else
        throw(ArgumentError("Unsupported backend: $backend"))
    end
end


#function curl(::Type{<:Backend}, field, spacing) end
#function divergence(::Type{<:Backend}, field, spacing) end

end