module CPUVectorOps
#---------------- Module containing CPU implementations of vector operations ---------------------------------------------
# expects the following layout:
#                                      1st dimension: x-direction   
#                                      2nd dimension: y-direction
#                                      3rd dimension: z-direction
#       1D arrays which call ddy or ddz will fall back to ddx
#       2D arrays which call ddz will fall back to ddy 
#--------------------
#       The following functions are implemented:
#                 upwind derivative: ddx_up, ddy_up, ddz_up
#                 downwind derivative: ddx_dn, ddy_dn, ddz_dn
#--------------------
#       2nd, 4th and 6th order methods are implemented for each function.
#       These functions should not be exported, but rather called by the exported functions
#       Each function must take 'order' as an argument to specify the order of the derivative
#       order signature should be given as integer:
#                       order = 2, 4, 6
#--------------------
#       Multiple dispatch is used to enable two types of function calls:
#          1. f(input, kwargs) -> This will create new output array and return it
#          2. f!(input, output, kwargs) -> This will write the result to the output array with no return value
# ------------------------------------------------------------------------------------------------------------------------
#-------------------------- Exported function signatures ----------------------------------------------------------------
export ddx_up, ddy_up, ddz_up, ddx_dn, ddy_dn, ddz_dn,
       curl_up, curl_dn
#----------------------------------------------------------------------------------------------------------------------
#----------------- import prefacts for derivatives, interpolation, and Laplace operator -------------------------------
include("Prefactors.jl")
using LoopVectorization
#-------------------------------------------------------------------------------------------------------------------

function curl_up(field::AbstractArray{<:AbstractFloat, 4}, dx::AbstractFloat, dy::AbstractFloat, dz::AbstractFloat, order::Int)
    # Ensure the order is valid
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 

    dims = size(field)  # dims should be (3, Nx, Ny, Nz)
    
    if (dims[1] != 3) 
        throw(ArgumentError("Curl_up: Expects first component by contant 3 components"))
    end 

    result = zeros(eltype(field), dims)

    x_in = @view field[1,:,:,:]
    y_in = @view field[2,:,:,:]
    z_in = @view field[3,:,:,:]

    # Use views for better performance
    res_x = @view result[1,:,:,:]
    res_y = @view result[2,:,:,:]
    res_z = @view result[3,:,:,:]


    curl_x_up(y_in, z_in, res_x, dy, dz, order)
    curl_y_up(x_in, z_in, res_y, dx, dz, order)
    curl_z_up(x_in, y_in, res_z, dx, dy, order)

    return result
end

function curl_up(field::AbstractArray{<:AbstractFloat, 4}, result::AbstractArray{<:AbstractFloat, 4}, dx::AbstractFloat, dy::AbstractFloat, dz::AbstractFloat, order::Int)
    # Ensure the order is valid
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 

    dims = size(field)  # dims should be (3, Nx, Ny, Nz)
    
    if (dims[1] != 3) 
        throw(ArgumentError("Curl_up: Expects first component by contant 3 components"))
    end 

    x_in = @view field[1,:,:,:]
    y_in = @view field[2,:,:,:]
    z_in = @view field[3,:,:,:]

    # Use views for better performance
    res_x = @view result[1,:,:,:]
    res_y = @view result[2,:,:,:]
    res_z = @view result[3,:,:,:]


    curl_x_up(y_in, z_in, res_x, dy, dz, order)
    curl_y_up(x_in, z_in, res_y, dx, dz, order)
    curl_z_up(x_in, y_in, res_z, dx, dy, order)
end

function curl_x_up(y::AbstractArray{<:AbstractFloat},z::AbstractArray{<:AbstractFloat}, out::AbstractArray{<:AbstractFloat}, dy::AbstractFloat, dz::AbstractFloat,order::Int)
    out .= ddy_up(z,dy,order) .- ddz_up(y,dz,order)
end
function curl_y_up(x::AbstractArray{<:AbstractFloat}, z::AbstractArray{<:AbstractFloat}, out::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, dz::AbstractFloat, order::Int)
    out .= ddz_up(x, dz, order) .- ddx_up(z, dx, order)
end

function curl_z_up(x::AbstractArray{<:AbstractFloat}, y::AbstractArray{<:AbstractFloat}, out::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, dy::AbstractFloat, order::Int)
    out .= ddx_up(y, dx, order) .- ddy_up(x, dy, order)
end

function curl_dn(field::AbstractArray{<:AbstractFloat, 4}, dx::AbstractFloat, dy::AbstractFloat, dz::AbstractFloat, order::Int)
    # Ensure the order is valid
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 

    dims = size(field)  # dims should be (3, Nx, Ny, Nz)
    
    if (dims[1] != 3) 
        throw(ArgumentError("Curl_up: Expects first component by contant 3 components"))
    end 

    result = zeros(eltype(field), dims)

    x_in = @view field[1,:,:,:]
    y_in = @view field[2,:,:,:]
    z_in = @view field[3,:,:,:]

    # Use views for better performance
    res_x = @view result[1,:,:,:]
    res_y = @view result[2,:,:,:]
    res_z = @view result[3,:,:,:]


    curl_x_dn(y_in, z_in, res_x, dy, dz, order)
    curl_y_dn(x_in, z_in, res_y, dx, dz, order)
    curl_z_dn(x_in, y_in, res_z, dx, dy, order)

    return result
end

function curl_dn(field::AbstractArray{<:AbstractFloat, 4}, result::AbstractArray{<:AbstractFloat, 4}, dx::AbstractFloat, dy::AbstractFloat, dz::AbstractFloat, order::Int)
    # Ensure the order is valid
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 

    dims = size(field)  # dims should be (3, Nx, Ny, Nz)
    
    if (dims[1] != 3) 
        throw(ArgumentError("Curl_up: Expects first component by contant 3 components"))
    end 

    x_in = @view field[1,:,:,:]
    y_in = @view field[2,:,:,:]
    z_in = @view field[3,:,:,:]

    # Use views for better performance
    res_x = @view result[1,:,:,:]
    res_y = @view result[2,:,:,:]
    res_z = @view result[3,:,:,:]


    curl_x_dn(y_in, z_in, res_x, dy, dz, order)
    curl_y_dn(x_in, z_in, res_y, dx, dz, order)
    curl_z_dn(x_in, y_in, res_z, dx, dy, order)
end

function curl_x_dn(y::AbstractArray{<:AbstractFloat}, z::AbstractArray{<:AbstractFloat}, out::AbstractArray{<:AbstractFloat}, dy::AbstractFloat, dz::AbstractFloat, order::Int)
    out .= ddy_dn(z, dy, order) .- ddz_dn(y, dz, order)
end

function curl_y_dn(x::AbstractArray{<:AbstractFloat}, z::AbstractArray{<:AbstractFloat}, out::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, dz::AbstractFloat, order::Int)
    out .= ddz_dn(x, dz, order) .- ddx_dn(z, dx, order)
end

function curl_z_dn(x::AbstractArray{<:AbstractFloat}, y::AbstractArray{<:AbstractFloat}, out::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, dy::AbstractFloat, order::Int)
    out .= ddx_dn(y, dx, order) .- ddy_dn(x, dy, order)
end
#---------------------------------------------------------------------------------------------------------------------






#---------------------- ddx_up --------------------------------------------------------------------------------------------
function ddx_up(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
    dims = size(field)
    result = zeros(eltype(field), dims)

    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 

    if order == 2 && dims[1] > 1
        ddx_up_2nd(field, result, dx)
    elseif order == 4 && dims[1] > 3
        ddx_up_4th(field, result, dx)
    elseif order == 6 && dims[1] > 5
        ddx_up_6th(field, result, dx)
    end
    return result
end

function ddx_up(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)

    dims = size(field)
    dims_res = size(result)

    if dims != dims_res
        throw(ArgumentError("ddx_up: Input and output arrays must have the same dimensions"))
    end

    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 


    if order == 2 && dims[1] > 1
        ddx_up_2nd(field, result, dx)
    elseif order == 4 && dims[1] > 3
        ddx_up_4th(field, result, dx)
    elseif order == 6 && dims[1] > 5
        ddx_up_6th(field, result, dx)
    end
end

function ddx_up_2nd(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dx::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf2.a / dx
    if ndims == 1
        @turbo for i in 1:dims[1]-1
            result[i] = pfa * (field[i+1] - field[i])
        end
    elseif ndims == 2
        @turbo for j in 1:dims[2], i in 1:dims[1]-1
            result[i, j] = pfa * (field[i+1, j] - field[i, j])
        end
        
    elseif ndims == 3
        @turbo for k in 1:dims[3], j in 1:dims[2], i in 1:dims[1]-1
            result[i, j, k] = pfa * (field[i+1, j, k] - field[i, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end


function ddx_up_4th(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dx::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf4.a / dx
    pfb = pf4.b / dx
    
    if ndims == 1
        @turbo for i in 2:dims[1]-2
            result[i] = pfa * (field[i+1] - field[i]) +
                        pfb * (field[i+2] - field[i-1])
        end
    elseif ndims == 2
        @turbo for j in 1:dims[2], i in 2:dims[1]-2
            result[i, j] = pfa * (field[i+1, j] - field[i, j]) +
                           pfb * (field[i+2, j] - field[i-1, j])
        end
    elseif ndims == 3
        @turbo for k in 1:dims[3], j in 1:dims[2], i in 2:dims[1]-2
            result[i, j, k] = pfa * (field[i+1, j, k] - field[i, j, k]) +
                              pfb * (field[i+2, j, k] - field[i-1, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end


function ddx_up_6th(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dx::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf6.a / dx
    pfb = pf6.b / dx
    pfc = pf6.c / dx

    if ndims == 1
        @turbo for i in 3:dims[1]-3
            result[i] = pfa * (field[i+1] - field[i]) +
                        pfb * (field[i+2] - field[i-1]) +
                        pfc * (field[i+3] - field[i-2])
        end
    elseif ndims == 2
        @turbo for j in 1:dims[2], i in 3:dims[1]-3
            result[i, j] = pfa * (field[i+1, j] - field[i, j]) +
                           pfb * (field[i+2, j] - field[i-1, j]) +
                           pfc * (field[i+3, j] - field[i-2, j])
        end
    elseif ndims == 3
        @turbo for k in 1:dims[3], j in 1:dims[2], i in 3:dims[1]-3
            result[i, j, k] = pfa * (field[i+1, j, k] - field[i, j, k]) +
                              pfb * (field[i+2, j, k] - field[i-1, j, k]) +
                              pfc * (field[i+3, j, k] - field[i-2, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end

#----------------------------------------------------------------------------------------------------------------------


#---------------------- ddy_up --------------------------------------------------------------------------------------------

function ddy_up(field::AbstractArray{<:AbstractFloat}, dy::AbstractFloat, order::Int)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 
    #------------- for 1D arrays use ddx 
    if ndims == 1
        ddx_up(field, result, dy, order)
        return result
    end 


    if order == 2 && dims[2] > 1
        ddy_up_2nd(field, result, dy)
    elseif order == 4 && dims[2] > 3
        ddy_up_4th(field, result, dy)
    elseif order == 6 && dims[2] > 5
        ddy_up_6th(field, result, dy)
    end
    return result
end

function ddy_up(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dy::AbstractFloat, order::Int)
    dims = size(field)
    ndims = length(dims)
    dims_res = size(result)

    if dims != dims_res
        throw(ArgumentError("ddy_up: Input and output arrays must have the same dimensions"))
    end
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 
    #------------- for 1D arrays use ddx 
    if ndims == 1
        ddx_up(field, result, dy, order)
        return
    end 


    if order == 2 && dims[2] > 1
        ddy_up_2nd(field, result, dy)
    elseif order == 4 && dims[2] > 3
        ddy_up_4th(field, result, dy)
    elseif order == 6 && dims[2] > 5
        ddy_up_6th(field, result, dy)
    end
end

function ddy_up_2nd(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dy::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf2.a / dy

    if ndims == 2
        @turbo for j in 1:dims[2]-1, i in 1:dims[1]
            result[i, j] = pfa * (field[i, j+1] - field[i, j])
        end
    elseif ndims == 3
        @turbo for k in 1:dims[3], j in 1:dims[2]-1, i in 1:dims[1]
            result[i, j, k] = pfa * (field[i, j+1, k] - field[i, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end


function ddy_up_4th(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dy::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf4.a / dy
    pfb = pf4.b / dy

    if ndims == 2
        @turbo for j in 2:dims[2]-2, i in 1:dims[1]
            result[i, j] = pfa * (field[i, j+1] - field[i, j]) +
                           pfb * (field[i, j+2] - field[i, j-1])
        end
    elseif ndims == 3
        @turbo for  k in 1:dims[3], j in 2:dims[2]-2, i in 1:dims[1]
            result[i, j, k] = pfa * (field[i, j+1, k] - field[i, j, k]) +
                              pfb * (field[i, j+2, k] - field[i, j-1, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end


function ddy_up_6th(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dy::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf6.a / dy
    pfb = pf6.b / dy
    pfc = pf6.c / dy

    if ndims == 2
        @turbo for  j in 3:dims[2]-3, i in 1:dims[1]
            result[i, j] = pfa * (field[i, j+1] - field[i, j]) +
                           pfb * (field[i, j+2] - field[i, j-1]) +
                           pfc * (field[i, j+3] - field[i, j-2])
        end
    elseif ndims == 3
        @turbo for k in 1:dims[3], j in 3:dims[2]-3, i in 1:dims[1]
            result[i, j, k] = pfa * (field[i, j+1, k] - field[i, j, k]) +
                              pfb * (field[i, j+2, k] - field[i, j-1, k]) +
                              pfc * (field[i, j+3, k] - field[i, j-2, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end

#----------------------------------------------------------------------------------------------------------------------

#---------------------- ddz_up --------------------------------------------------------------------------------------------

function ddz_up(field::AbstractArray{<:AbstractFloat}, dz::AbstractFloat, order::Int)
    dims = size(field)
    result = zeros(eltype(field), dims)

    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 


    #------------- for 1D arrays use ddx 
    if length(size(field)) == 1
        ddx_up(field, result, dz, order)
        return result
    end 

    #------------- for 2D arrays use ddy 
    if length(size(field)) == 2
        ddy_up(field, result, dz, order)
        return result
    end 

    if order == 2 && dims[3] > 1
        ddz_up_2nd(field, result, dz)
    elseif order == 4 && dims[3] > 3
        ddz_up_4th(field, result, dz)
    elseif order == 6 && dims[3] > 5
        ddz_up_6th(field, result, dz)
    end
    return result
end

function ddz_up(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dz::AbstractFloat, order::Int)
    dims = size(field)
    ndims = length(dims)
    dims_res = size(result)
    if dims != dims_res
        throw(ArgumentError("ddz_up: Input and output arrays must have the same dimensions"))
    end
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 


    #------------- for 1D arrays use ddx 
    if ndims == 1
        ddx_up(field, result, dz, order)
        return 
    end 

    #------------- for 2D arrays use ddy 
    if ndims== 2
        ddy_up(field, result, dz, order)
        return 
    end 

    if order == 2 && dims[3] > 1
        ddz_up_2nd(field, result, dz)
    elseif order == 4 && dims[3] > 3
        ddz_up_4th(field, result, dz)
    elseif order == 6 && dims[3] > 5
        ddz_up_6th(field, result, dz)
    end
end

function ddz_up_2nd(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dz::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf2.a / dz

    if ndims == 3
        @turbo for k in 1:dims[3]-1, j in 1:dims[2], i in 1:dims[1]
            result[i, j, k] = pfa * (field[i, j, k+1] - field[i, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end


function ddz_up_4th(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dz::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf4.a / dz
    pfb = pf4.b / dz

    if ndims == 3
        @turbo for k in 2:dims[3]-2, j in 1:dims[2], i in 1:dims[1]
            result[i, j, k] = pfa * (field[i, j, k+1] - field[i, j, k]) +
                              pfb * (field[i, j, k+2] - field[i, j, k-1])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end


function ddz_up_6th(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dz::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf6.a / dz
    pfb = pf6.b / dz
    pfc = pf6.c / dz

    if ndims == 3
        @turbo for k in 3:dims[3]-3, j in 1:dims[2], i in 1:dims[1]
            result[i, j, k] = pfa * (field[i, j, k+1] - field[i, j, k]) +
                              pfb * (field[i, j, k+2] - field[i, j, k-1]) +
                              pfc * (field[i, j, k+3] - field[i, j, k-2])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end

#----------------------------------------------------------------------------------------------------------------------

#---------------------- ddx_up --------------------------------------------------------------------------------------------
function ddx_dn(field::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
    dims = size(field)
    result = zeros(eltype(field), dims)
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 
    if order == 2 && dims[1] > 1
        ddx_dn_2nd(field, result, dx)
    elseif order == 4 && dims[1] > 3
        ddx_dn_4th(field, result, dx)
    elseif order == 6 && dims[1] > 5
        ddx_dn_6th(field, result, dx)
    end
    return result
end

function ddx_dn(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dx::AbstractFloat, order::Int)
    dims = size(field)
    dims_res = size(result)
    if dims != dims_res
        throw(ArgumentError("ddx_dn: Input and output arrays must have the same dimensions"))
    end
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 
    if order == 2 && dims[1] > 1
        ddx_dn_2nd(field, result, dx)
    elseif order == 4 && dims[1] > 3
        ddx_dn_4th(field, result, dx)
    elseif order == 6 && dims[1] > 5
        ddx_dn_6th(field, result, dx)
    end
end

function ddx_dn_2nd(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dx::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf2.a / dx
    if ndims == 1
        @turbo for i in 2:dims[1]
            result[i] = pfa * (field[i] - field[i-1])
        end
    elseif ndims == 2
        @turbo for j in 1:dims[2], i in 2:dims[1]
            result[i, j] = pfa * (field[i, j] - field[i-1, j])
        end
        
    elseif ndims == 3
        @turbo for k in 1:dims[3], j in 1:dims[2], i in 2:dims[1]
            result[i, j, k] = pfa * (field[i, j, k] - field[i-1, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end


function ddx_dn_4th(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dx::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf4.a / dx
    pfb = pf4.b / dx
    
    if ndims == 1
        @turbo for i in 3:dims[1]-1
            result[i] = pfa * (field[i] - field[i-1]) +
                        pfb * (field[i+1] - field[i-2])
        end
    elseif ndims == 2
        @turbo for j in 1:dims[2], i in 3:dims[1]-1
            result[i, j] = pfa * (field[i, j] - field[i-1, j]) +
                           pfb * (field[i+1, j] - field[i-2, j])
        end
    elseif ndims == 3
        @turbo for k in 1:dims[3], j in 1:dims[2], i in 3:dims[1]-1
            result[i, j, k] = pfa * (field[i, j, k] - field[i-1, j, k]) +
                              pfb * (field[i+1, j, k] - field[i-2, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end


function ddx_dn_6th(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dx::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf6.a / dx
    pfb = pf6.b / dx
    pfc = pf6.c / dx

    if ndims == 1
        @turbo for i in 4:dims[1]-2
            result[i] = pfa * (field[i  ] - field[i-1]) +
                        pfb * (field[i+1] - field[i-2]) +
                        pfc * (field[i+2] - field[i-3])
        end
    elseif ndims == 2
        @turbo for  j in 1:dims[2], i in 4:dims[1]-2
            result[i, j] = pfa * (field[i  , j] - field[i-1, j]) +
                           pfb * (field[i+1, j] - field[i-2, j]) +
                           pfc * (field[i+2, j] - field[i-3, j])
        end
    elseif ndims == 3
        @turbo for  k in 1:dims[3], j in 1:dims[2], i in 4:dims[1]-2
            result[i, j, k] = pfa * (field[i  , j, k] - field[i-1, j, k]) +
                              pfb * (field[i+1, j, k] - field[i-2, j, k]) +
                              pfc * (field[i+2, j, k] - field[i-3, j, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end

#----------------------------------------------------------------------------------------------------------------------


#---------------------- ddy_up --------------------------------------------------------------------------------------------

function ddy_dn(field::AbstractArray{<:AbstractFloat}, dy::AbstractFloat, order::Int)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 

    #------------- for 1D arrays use ddx 
    if ndims== 1
        ddx_dn(field, result, dy, order)
        return result
    end 

    if order == 2 && dims[2] > 1
        ddy_dn_2nd(field, result, dy)
    elseif order == 4 && dims[2] > 3
        ddy_dn_4th(field, result, dy)
    elseif order == 6 && dims[2] > 5
        ddy_dn_6th(field, result, dy)
    end
    return result
end

function ddy_dn(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dy::AbstractFloat, order::Int)
    dims = size(field)
    ndims = length(dims)
    dims_res = size(result)
    if dims != dims_res
        throw(ArgumentError("ddy_dn: Input and output arrays must have the same dimensions"))
    end
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 

    #------------- for 1D arrays use ddx 
    if ndims== 1
        ddx_dn(field, result, dy, order)
        return
    end 

    if order == 2 && dims[2] > 1
        ddy_dn_2nd(field, result, dy)
    elseif order == 4 && dims[2] > 3
        ddy_dn_4th(field, result, dy)
    elseif order == 6 && dims[2] > 5
        ddy_dn_6th(field, result, dy)
    end
end

function ddy_dn_2nd(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dy::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf2.a / dy

    if ndims == 2
        @turbo for j in 2:dims[2], i in 1:dims[1]
            result[i, j] = pfa * (field[i, j] - field[i, j-1])
        end
    elseif ndims == 3
        @turbo for k in 1:dims[3], j in 2:dims[2], i in 1:dims[1]
            result[i, j, k] = pfa * (field[i, j, k] - field[i, j-1, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end


function ddy_dn_4th(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dy::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf4.a / dy
    pfb = pf4.b / dy


    if ndims == 2
        @turbo for j in 3:dims[2]-1, i in 1:dims[1]
            result[i, j] = pfa * (field[i, j  ] - field[i, j-1]) +
                           pfb * (field[i, j+1] - field[i, j-2])
        end
    elseif ndims == 3
        @turbo for  k in 1:dims[3], j in 3:dims[2]-1, i in 1:dims[1]
            result[i, j, k] = pfa * (field[i, j  , k] - field[i, j-1, k]) +
                              pfb * (field[i, j+1, k] - field[i, j-2, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end


function ddy_dn_6th(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dy::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf6.a / dy
    pfb = pf6.b / dy
    pfc = pf6.c / dy

    if ndims == 2
        @turbo for  j in 4:dims[2]-2,i in 1:dims[1]
            result[i, j] = pfa * (field[i, j  ] - field[i, j-1]) +
                           pfb * (field[i, j+1] - field[i, j-2]) +
                           pfc * (field[i, j+2] - field[i, j-3])
        end
    elseif ndims == 3
        @turbo for  k in 1:dims[3], j in 4:dims[2]-2, i in 1:dims[1]
            result[i, j, k] = pfa * (field[i, j  , k] - field[i, j-1, k]) +
                              pfb * (field[i, j+1, k] - field[i, j-2, k]) +
                              pfc * (field[i, j+2, k] - field[i, j-3, k])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end

#----------------------------------------------------------------------------------------------------------------------

#---------------------- ddz_up --------------------------------------------------------------------------------------------

function ddz_dn(field::AbstractArray{<:AbstractFloat}, dz::AbstractFloat, order::Int)
    dims = size(field)
    ndims = length(dims)
    result = zeros(eltype(field), dims)
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 

    #------------- for 1D arrays use ddx 
    if ndims== 1
        ddx_dn(field, result, dz, order)
        return result
    end 

    #------------- for 2D arrays use ddy 
    if ndims == 2
        ddy_dn(field, result, dz, order)
        return result
    end 

    if order == 2 && dims[3] > 1
        ddz_dn_2nd(field, result, dz)
    elseif order == 4 && dims[3] > 3
        ddz_dn_4th(field, result, dz)
    elseif order == 6 && dims[3] > 5
        ddz_dn_6th(field, result, dz)
    end
    return result
end


function ddz_dn(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dz::AbstractFloat, order::Int)
    dims = size(field)
    ndims = length(dims)
    dims_res = size(result)
    if dims != dims_res
        throw(ArgumentError("ddz_dn: Input and output arrays must have the same dimensions"))
    end
    if !(order == 2 || order == 4 || order == 6)
        throw(ArgumentError("Unsupported order: $order"))
    end 

    #------------- for 1D arrays use ddx 
    if ndims == 1
        ddx_dn(field, result, dz, order)
        return
    end 

    #------------- for 2D arrays use ddy 
    if ndims == 2
        ddy_dn(field, result, dz, order)
        return
    end 

    if order == 2 && dims[3] > 1
         ddz_dn_2nd(field, result, dz)
    elseif order == 4 && dims[3] > 3
        ddz_dn_4th(field, result, dz)
    elseif order == 6 && dims[3] > 5
        ddz_dn_6th(field, result, dz)
    end
end

function ddz_dn_2nd(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dz::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf2.a / dz


    if ndims == 3
        @turbo for k in 2:dims[3], j in 1:dims[2], i in 1:dims[1]
            result[i, j, k] = pfa * (field[i, j, k] - field[i, j, k-1])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end


function ddz_dn_4th(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dz::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf4.a / dz
    pfb = pf4.b / dz

 

    if ndims == 3
        @turbo for k in 3:dims[3]-1, j in 1:dims[2], i in 1:dims[1]
            result[i, j, k] = pfa * (field[i, j, k  ] - field[i, j, k-1]) +
                              pfb * (field[i, j, k+1] - field[i, j, k-2])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end


function ddz_dn_6th(field::AbstractArray{<:AbstractFloat}, result::AbstractArray{<:AbstractFloat}, dz::AbstractFloat)
    dims = size(field)
    ndims = length(dims)
    pfa = pf6.a / dz
    pfb = pf6.b / dz
    pfc = pf6.c / dz

    if ndims == 3
        @turbo for k in 4:dims[3]-2, j in 1:dims[2], i in 1:dims[1]
            result[i, j, k] = pfa * (field[i, j, k  ] - field[i, j, k-1]) +
                              pfb * (field[i, j, k+1] - field[i, j, k-2]) +
                              pfc * (field[i, j, k+2] - field[i, j, k-3])
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
end

#----------------------------------------------------------------------------------------------------------------------

end
