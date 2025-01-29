module GPUVectorOps

include("Prefactors.jl")
using .Prefactors: pf2, pf4, pf6


export ddx_up, ddy_up, ddz_upd, dx_dn, ddy_dn, ddz_dn



    function ddx_up(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YER"))
    end
    function ddy_up(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YER"))
    end
    function ddz_up(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YER"))
    end

    function ddx_dn(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YER"))
    end
    function ddy_dn(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YER"))
    end
    function ddz_dn(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YER"))
    end


end