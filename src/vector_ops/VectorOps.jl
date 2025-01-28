module VectorOps

# Import backend-specific implementations
include("VectorCPU.jl")
include("VectorGPU.jl")


export Backend, CPU, GPU, ddx_up  #curl, divergence

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


#function curl(::Type{<:Backend}, field, spacing) end
#function divergence(::Type{<:Backend}, field, spacing) end

end