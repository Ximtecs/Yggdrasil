module CPUVectorOps
#---------------- Module containing CPU implementations of vector operations ---------------------------------------------
# expects the following layout:
#                                      1st dimension: x-direction   
#                                      2nd dimension: y-direction
#                                      3rd dimension: z-direction
#       1D arrays which call ddy or ddz will fall back to ddx
#       2D arrays which call ddz will fall back to ddy 
# ------------------------------------------------------------------------------------------------------------------------



include("Prefactors.jl")
using .Prefactors: pf2, pf4, pf6
using LoopVectorization

export ddx_up, ddy_up, ddz_up, ddx_dn, ddy_dn, ddz_dn  # Export only what's necessary


#---------------------- ddx_up --------------------------------------------------------------------------------------------
function ddx_up(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
    if order == 2
        return ddx_up_2nd(field, dx)
    elseif order == 4
        return ddx_up_4th(field, dx)
    elseif order == 6
        return ddx_up_6th(field, dx)
    else
        throw(ArgumentError("Unsupported order: $order"))
    end
end

function ddx_up_2nd(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf2.a / dx
    if dims[1] < 3
        return result  # Return zeros for insufficient data
    end
    if ndims == 1
        @turbo for i in 1:dims[1]-1
            result[i] = pfa * (field[i+1] - field[i])
        end
    elseif ndims == 2
        @turbo for i in 1:dims[1]-1, j in 1:dims[2]
            result[i, j] = pfa * (field[i+1, j] - field[i, j])
        end
        
    elseif ndims == 3
        @turbo for i in 1:dims[1]-1, j in 1:dims[2], k in 1:dims[3]
            result[i, j, k] = pfa * (field[i+1, j, k] - field[i, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end


function ddx_up_4th(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf4.a / dx
    pfb = pf4.b / dx
    
    if dims[1] < 4
        return result  # Return zeros for insufficient data
    end
    
    if ndims == 1
        @turbo for i in 2:dims[1]-2
            result[i] = pfa * (field[i+1] - field[i]) +
                        pfb * (field[i+2] - field[i-1])
        end
    elseif ndims == 2
        @turbo for i in 2:dims[1]-2, j in 1:dims[2]
            result[i, j] = pfa * (field[i+1, j] - field[i, j]) +
                           pfb * (field[i+2, j] - field[i-1, j])
        end
    elseif ndims == 3
        @turbo for i in 2:dims[1]-2, j in 1:dims[2], k in 1:dims[3]
            result[i, j, k] = pfa * (field[i+1, j, k] - field[i, j, k]) +
                              pfb * (field[i+2, j, k] - field[i-1, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end


function ddx_up_6th(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf6.a / dx
    pfb = pf6.b / dx
    pfc = pf6.c / dx

    if dims[1] < 6
        return result  # Return zeros for insufficient data
    end

    if ndims == 1
        @turbo for i in 3:dims[1]-3
            result[i] = pfa * (field[i+1] - field[i]) +
                        pfb * (field[i+2] - field[i-1]) +
                        pfc * (field[i+3] - field[i-2])
        end
    elseif ndims == 2
        @turbo for i in 3:dims[1]-3, j in 1:dims[2]
            result[i, j] = pfa * (field[i+1, j] - field[i, j]) +
                           pfb * (field[i+2, j] - field[i-1, j]) +
                           pfc * (field[i+3, j] - field[i-2, j])
        end
    elseif ndims == 3
        @turbo for i in 3:dims[1]-3, j in 1:dims[2], k in 1:dims[3]
            result[i, j, k] = pfa * (field[i+1, j, k] - field[i, j, k]) +
                              pfb * (field[i+2, j, k] - field[i-1, j, k]) +
                              pfc * (field[i+3, j, k] - field[i-2, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end

#----------------------------------------------------------------------------------------------------------------------


#---------------------- ddy_up --------------------------------------------------------------------------------------------

function ddy_up(field::AbstractArray{<:AbstractFloat}, dy::AbstractFloat, order::Int)


    #------------- for 1D arrays use ddx 
    if length(size(field)) == 1
        return ddx_up(field, dy, order)
    end 


    if order == 2
        return ddy_up_2nd(field, dy)
    elseif order == 4
        return ddy_up_4th(field, dy)
    elseif order == 6
        return ddy_up_6th(field, dy)
    else
        throw(ArgumentError("Unsupported order: $order"))
    end
end

function ddy_up_2nd(field::AbstractArray{<:AbstractFloat}, dy::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf2.a / dy

    if dims[2] < 3
        return result  # Return zeros for insufficient data
    end

    if ndims == 2
        @turbo for i in 1:dims[1], j in 1:dims[2]-1
            result[i, j] = pfa * (field[i, j+1] - field[i, j])
        end
    elseif ndims == 3
        @turbo for i in 1:dims[1], j in 1:dims[2]-1, k in 1:dims[3]
            result[i, j, k] = pfa * (field[i, j+1, k] - field[i, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end


function ddy_up_4th(field::AbstractArray{<:AbstractFloat}, dy::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf4.a / dy
    pfb = pf4.b / dy

    if dims[2] < 4
        return result  # Return zeros for insufficient data
    end

    if ndims == 2
        @turbo for i in 1:dims[1], j in 2:dims[2]-2
            result[i, j] = pfa * (field[i, j+1] - field[i, j]) +
                           pfb * (field[i, j+2] - field[i, j-1])
        end
    elseif ndims == 3
        @turbo for i in 1:dims[1], j in 2:dims[2]-2, k in 1:dims[3]
            result[i, j, k] = pfa * (field[i, j+1, k] - field[i, j, k]) +
                              pfb * (field[i, j+2, k] - field[i, j-1, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end


function ddy_up_6th(field::AbstractArray{<:AbstractFloat}, dy::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf6.a / dy
    pfb = pf6.b / dy
    pfc = pf6.c / dy

    if dims[2] < 6
        return result  # Return zeros for insufficient data
    end
    if ndims == 2
        @turbo for i in 1:dims[1], j in 3:dims[2]-3
            result[i, j] = pfa * (field[i, j+1] - field[i, j]) +
                           pfb * (field[i, j+2] - field[i, j-1]) +
                           pfc * (field[i, j+3] - field[i, j-2])
        end
    elseif ndims == 3
        @turbo for i in 1:dims[1], j in 3:dims[2]-3, k in 1:dims[3]
            result[i, j, k] = pfa * (field[i, j+1, k] - field[i, j, k]) +
                              pfb * (field[i, j+2, k] - field[i, j-1, k]) +
                              pfc * (field[i, j+3, k] - field[i, j-2, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end

#----------------------------------------------------------------------------------------------------------------------

#---------------------- ddz_up --------------------------------------------------------------------------------------------

function ddz_up(field::AbstractArray{<:AbstractFloat}, dz::AbstractFloat, order::Int)

    #------------- for 1D arrays use ddx 
    if length(size(field)) == 1
        return ddx_up(field, dz, order)
    end 

    #------------- for 2D arrays use ddy 
    if length(size(field)) == 2
        return ddy_up(field, dz, order)
    end 

    if order == 2
        return ddz_up_2nd(field, dz)
    elseif order == 4
        return ddz_up_4th(field, dz)
    elseif order == 6
        return ddz_up_6th(field, dz)
    else
        throw(ArgumentError("Unsupported order: $order"))
    end
end

function ddz_up_2nd(field::AbstractArray{<:AbstractFloat}, dz::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf2.a / dz

    if dims[3] < 3
        return result  # Return zeros for insufficient data
    end

    if ndims == 3
        @turbo for i in 1:dims[1], j in 1:dims[2], k in 1:dims[3]-1
            result[i, j, k] = pfa * (field[i, j, k+1] - field[i, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end


function ddz_up_4th(field::AbstractArray{<:AbstractFloat}, dz::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf4.a / dz
    pfb = pf4.b / dz

    if dims[3] < 4
        return result  # Return zeros for insufficient data
    end

    if ndims == 3
        @turbo for i in 1:dims[1], j in 1:dims[2], k in 2:dims[3]-2
            result[i, j, k] = pfa * (field[i, j, k+1] - field[i, j, k]) +
                              pfb * (field[i, j, k+2] - field[i, j, k-1])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end


function ddz_up_6th(field::AbstractArray{<:AbstractFloat}, dz::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf6.a / dz
    pfb = pf6.b / dz
    pfc = pf6.c / dz

    if dims[3] < 6
        return result  # Return zeros for insufficient data
    end

    if ndims == 3
        @turbo for i in 1:dims[1], j in 1:dims[2], k in 3:dims[3]-3
            result[i, j, k] = pfa * (field[i, j, k+1] - field[i, j, k]) +
                              pfb * (field[i, j, k+2] - field[i, j, k-1]) +
                              pfc * (field[i, j, k+3] - field[i, j, k-2])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end

#----------------------------------------------------------------------------------------------------------------------

#---------------------- ddx_up --------------------------------------------------------------------------------------------
function ddx_dn(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
    if order == 2
        return ddx_dn_2nd(field, dx)
    elseif order == 4
        return ddx_dn_4th(field, dx)
    elseif order == 6
        return ddx_dn_6th(field, dx)
    else
        throw(ArgumentError("Unsupported order: $order"))
    end
end

function ddx_dn_2nd(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf2.a / dx
    if dims[1] < 3
        return result  # Return zeros for insufficient data
    end
    if ndims == 1
        @turbo for i in 2:dims[1]
            result[i] = pfa * (field[i] - field[i-1])
        end
    elseif ndims == 2
        @turbo for i in 2:dims[1], j in 1:dims[2]
            result[i, j] = pfa * (field[i, j] - field[i-1, j])
        end
        
    elseif ndims == 3
        @turbo for i in 2:dims[1], j in 1:dims[2], k in 1:dims[3]
            result[i, j, k] = pfa * (field[i, j, k] - field[i-1, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end


function ddx_dn_4th(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf4.a / dx
    pfb = pf4.b / dx
    
    if dims[1] < 4
        return result  # Return zeros for insufficient data
    end
    
    if ndims == 1
        @turbo for i in 3:dims[1]-1
            result[i] = pfa * (field[i] - field[i-1]) +
                        pfb * (field[i+1] - field[i-2])
        end
    elseif ndims == 2
        @turbo for i in 3:dims[1]-1, j in 1:dims[2]
            result[i, j] = pfa * (field[i, j] - field[i-1, j]) +
                           pfb * (field[i+1, j] - field[i-2, j])
        end
    elseif ndims == 3
        @turbo for i in 3:dims[1]-1, j in 1:dims[2], k in 1:dims[3]
            result[i, j, k] = pfa * (field[i, j, k] - field[i-1, j, k]) +
                              pfb * (field[i+1, j, k] - field[i-2, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end


function ddx_dn_6th(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf6.a / dx
    pfb = pf6.b / dx
    pfc = pf6.c / dx

    if dims[1] < 6
        return result  # Return zeros for insufficient data
    end

    if ndims == 1
        @turbo for i in 4:dims[1]-2
            result[i] = pfa * (field[i  ] - field[i-1]) +
                        pfb * (field[i+1] - field[i-2]) +
                        pfc * (field[i+2] - field[i-3])
        end
    elseif ndims == 2
        @turbo for i in 4:dims[1]-2, j in 1:dims[2]
            result[i, j] = pfa * (field[i  , j] - field[i-1, j]) +
                           pfb * (field[i+1, j] - field[i-2, j]) +
                           pfc * (field[i+2, j] - field[i-3, j])
        end
    elseif ndims == 3
        @turbo for i in 4:dims[1]-2, j in 1:dims[2], k in 1:dims[3]
            result[i, j, k] = pfa * (field[i  , j, k] - field[i-1, j, k]) +
                              pfb * (field[i+1, j, k] - field[i-2, j, k]) +
                              pfc * (field[i+2, j, k] - field[i-3, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end

#----------------------------------------------------------------------------------------------------------------------


#---------------------- ddy_up --------------------------------------------------------------------------------------------

function ddy_dn(field::AbstractArray{<:AbstractFloat}, dy::AbstractFloat, order::Int)


    #------------- for 1D arrays use ddx 
    if length(size(field)) == 1
        return ddx_dn(field, dy, order)
    end 


    if order == 2
        return ddy_dn_2nd(field, dy)
    elseif order == 4
        return ddy_dn_4th(field, dy)
    elseif order == 6
        return ddy_dn_6th(field, dy)
    else
        throw(ArgumentError("Unsupported order: $order"))
    end
end

function ddy_dn_2nd(field::AbstractArray{<:AbstractFloat}, dy::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf2.a / dy

    if dims[2] < 3
        return result  # Return zeros for insufficient data
    end

    if ndims == 2
        @turbo for i in 1:dims[1], j in 2:dims[2]
            result[i, j] = pfa * (field[i, j] - field[i, j-1])
        end
    elseif ndims == 3
        @turbo for i in 1:dims[1], j in 2:dims[2], k in 1:dims[3]
            result[i, j, k] = pfa * (field[i, j, k] - field[i, j-1, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end


function ddy_dn_4th(field::AbstractArray{<:AbstractFloat}, dy::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf4.a / dy
    pfb = pf4.b / dy

    if dims[2] < 4
        return result  # Return zeros for insufficient data
    end

    if ndims == 2
        @turbo for i in 1:dims[1], j in 3:dims[2]-1
            result[i, j] = pfa * (field[i, j  ] - field[i, j-1]) +
                           pfb * (field[i, j+1] - field[i, j-2])
        end
    elseif ndims == 3
        @turbo for i in 1:dims[1], j in 3:dims[2]-1, k in 1:dims[3]
            result[i, j, k] = pfa * (field[i, j  , k] - field[i, j-1, k]) +
                              pfb * (field[i, j+1, k] - field[i, j-2, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end


function ddy_dn_6th(field::AbstractArray{<:AbstractFloat}, dy::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf6.a / dy
    pfb = pf6.b / dy
    pfc = pf6.c / dy

    if dims[2] < 6
        return result  # Return zeros for insufficient data
    end
    if ndims == 2
        @turbo for i in 1:dims[1], j in 4:dims[2]-2
            result[i, j] = pfa * (field[i, j  ] - field[i, j-1]) +
                           pfb * (field[i, j+1] - field[i, j-2]) +
                           pfc * (field[i, j+2] - field[i, j-3])
        end
    elseif ndims == 3
        @turbo for i in 1:dims[1], j in 4:dims[2]-2, k in 1:dims[3]
            result[i, j, k] = pfa * (field[i, j  , k] - field[i, j-1, k]) +
                              pfb * (field[i, j+1, k] - field[i, j-2, k]) +
                              pfc * (field[i, j+2, k] - field[i, j-3, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end

#----------------------------------------------------------------------------------------------------------------------

#---------------------- ddz_up --------------------------------------------------------------------------------------------

function ddz_dn(field::AbstractArray{<:AbstractFloat}, dz::AbstractFloat, order::Int)

    #------------- for 1D arrays use ddx 
    if length(size(field)) == 1
        return ddx_dn(field, dz, order)
    end 

    #------------- for 2D arrays use ddy 
    if length(size(field)) == 2
        return ddy_dn(field, dz, order)
    end 

    if order == 2
        return ddz_dn_2nd(field, dz)
    elseif order == 4
        return ddz_dn_4th(field, dz)
    elseif order == 6
        return ddz_dn_6th(field, dz)
    else
        throw(ArgumentError("Unsupported order: $order"))
    end
end

function ddz_dn_2nd(field::AbstractArray{<:AbstractFloat}, dz::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf2.a / dz

    if dims[3] < 3
        return result  # Return zeros for insufficient data
    end

    if ndims == 3
        @turbo for i in 1:dims[1], j in 1:dims[2], k in 2:dims[3]
            result[i, j, k] = pfa * (field[i, j, k] - field[i, j, k-1])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end


function ddz_dn_4th(field::AbstractArray{<:AbstractFloat}, dz::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf4.a / dz
    pfb = pf4.b / dz

    if dims[3] < 4
        return result  # Return zeros for insufficient data
    end

    if ndims == 3
        @turbo for i in 1:dims[1], j in 1:dims[2], k in 3:dims[3]-1
            result[i, j, k] = pfa * (field[i, j, k  ] - field[i, j, k-1]) +
                              pfb * (field[i, j, k+1] - field[i, j, k-2])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end


function ddz_dn_6th(field::AbstractArray{<:AbstractFloat}, dz::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    pfa = pf6.a / dz
    pfb = pf6.b / dz
    pfc = pf6.c / dz

    if dims[3] < 6
        return result  # Return zeros for insufficient data
    end

    if ndims == 3
        @turbo for i in 1:dims[1], j in 1:dims[2], k in 4:dims[3]-2
            result[i, j, k] = pfa * (field[i, j, k  ] - field[i, j, k-1]) +
                              pfb * (field[i, j, k+1] - field[i, j, k-2]) +
                              pfc * (field[i, j, k+2] - field[i, j, k-3])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end

#----------------------------------------------------------------------------------------------------------------------

end
