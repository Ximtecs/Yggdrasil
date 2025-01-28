module CPUVectorOps

include("Prefactors.jl")
using .Prefactors: pf2, pf4, pf6

export ddx_up  # Export only what's necessary

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
        for i in 2:dims[1]-1
            result[i] = pfa * (field[i+1] - field[i])
        end
    elseif ndims == 2
        for j in 2:dims[2]-1
            for i in 2:dims[1]-1
                result[i, j] = pfa * (field[i+1, j] - field[i-1, j])
            end
        end
    elseif ndims == 3
        for k in 2:dims[3]-1
            for j in 2:dims[2]-1
                for i in 2:dims[1]-1
                    result[i, j, k] = pfa * (field[i+1, j, k] - field[i-1, j, k])
                end
            end
        end
    else
        throw(ArgumentError("Unsupported array dimension: $ndims"))
    end
    return result
end

end
