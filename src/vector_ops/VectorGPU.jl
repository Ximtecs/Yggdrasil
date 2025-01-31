module GPUVectorOps_mod

include("Prefactors.jl")
using .Prefactors: pf2, pf4, pf6


export ddx_up, ddy_up, ddz_upd, dx_dn, ddy_dn, ddz_dn,
     curl_up, curl_dn



    function curl_up(field::AbstractArray{<:AbstractFloat, 4}, dx::AbstractFloat, dy::AbstractFloat, dz::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YER"))
    end    
    function curl_up(field::AbstractArray{<:AbstractFloat, 4},result::AbstractArray{<:AbstractFloat, 4}, dx::AbstractFloat, dy::AbstractFloat, dz::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YER"))
    end
    function curl_dn(field::AbstractArray{<:AbstractFloat, 4}, dx::AbstractFloat, dy::AbstractFloat, dz::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YER"))
    end
    function curl_dn(field::AbstractArray{<:AbstractFloat, 4},result::AbstractArray{<:AbstractFloat, 4}, dx::AbstractFloat, dy::AbstractFloat, dz::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YER"))
    end

    
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

    function ddx_up!(result::AbstractArray{<:AbstractFloat}, field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YET"))
    end

    function ddy_up!(result::AbstractArray{<:AbstractFloat}, field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YET"))
    end

    function ddz_up!(result::AbstractArray{<:AbstractFloat}, field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YET"))
    end

    function ddx_dn!(result::AbstractArray{<:AbstractFloat}, field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YET"))
    end

    function ddy_dn!(result::AbstractArray{<:AbstractFloat}, field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YET"))
    end

    function ddz_dn!(result::AbstractArray{<:AbstractFloat}, field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
        throw(ArgumentError("NOT IMPLEMENTED YET"))
    end


end