module GPUVectorOps

include("Prefactors.jl")
using .Prefactors: pf2, pf4, pf6


export ddx_up   


    function ddx_up(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YER"))
    end


end