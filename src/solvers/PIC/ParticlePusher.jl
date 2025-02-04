#export vay_pusher
"""
   vay_pusher(p, w, B_arr, E_arr, dt, q_m_ratio)

Performs the Vay particle pusher algorithm to update the momenta `p` of charged particles 
based on the provided electric (`E_arr`) and magnetic (`B_arr`) fields.

# Arguments
- `p::Array{Float64,2}`: (3, n) array of particle momenta.
- `w::Array{Float64,1}`: (n) array of particle weights
        w is optional, when not given, @simd will be used, otherwise ifstatement will check if w>0
- `B_arr::Array{Float64,2}`: (3, n) array of magnetic field values at particle positions.
- `E_arr::Array{Float64,2}`: (3, n) array of electric field values at particle positions.
- `dt::Float64`: Time step.
- `q_m_ratio::Float64`: Charge-to-mass ratio of the particles.

# Notes
- Assumes column-major storage (efficient for Julia).
- Uses `@inbounds` to disable bounds checking for better performance.
- Uses `@simd` to enable SIMD vectorization when possible.
"""
#function vay_pusher(p::Array{<:AbstractFloat,2}, 
#                     B_arr::Array{<:AbstractFloat,2}, E_arr::Array{<:AbstractFloat,2}, 
#                     dt::Float64, q_m_ratio::Float64, k_F::Float64, c :: Float64)

function vay_pusher(p::Array{<:AbstractFloat,2}, 
                    B_arr::Array{<:AbstractFloat,2}, 
                    E_arr::Array{<:AbstractFloat,2};
                    kwargs...)


    allowed_keys = [:dt, :q_m_ratio, :k_F, :c]#, :TRACE_STATE]
    for key in keys(kwargs)
        if key ∉ allowed_keys
            error("Invalid keyword argument: $key. Allowed keys: $(allowed_keys)")
        end
    end
    #---------------- default values ----------------#
    dt = get(kwargs, :dt, 1.0)
    q_m_ratio = get(kwargs, :q_m_ratio, 1.0)
    k_F = get(kwargs, :k_F, 1.0)
    c = get(kwargs, :c, 1.0)
    #------------------------------------------------#
    @trace "vay_pusher" begin


    n = size(p, 2)  # Number of particles

    qmt_h  = 0.5 * dt * q_m_ratio  # Precomputed factors
    qmt_hB = qmt_h  * k_F
    o_c2 = 1.0 / c^2.0

    @inbounds @simd for i in 1:n
        # Extract electric and magnetic fields
        E = @SVector [E_arr[1, i], E_arr[2, i], E_arr[3, i]]
        B = @SVector [B_arr[1, i], B_arr[2, i], B_arr[3, i]]

        # Proper velocity and velocity calculations
        u_i = @SVector [p[1, i], p[2, i], p[3, i]]
        o_gam = 1.0 / sqrt(1.0 + sum(u_i.^2 * o_c2))  # Assuming `pic_scaling%elm%o_c2 = 1.0`
        v_i = u_i * o_gam
        tau = qmt_hB * B

        # Equation (9) - Vay 2008
        u_dash = u_i + 2 * qmt_h * E + cross(v_i, tau)

        gam_dash_2 = 1 + sum(u_dash.^2) * o_c2
        tau_2 = sum(tau .* tau)
        sigma = gam_dash_2 - tau_2
        u_star = sum(u_dash .* tau)

        o_gam_i1 = sqrt(2.0) / sqrt(sigma + sqrt(sigma^2 + 4 * (tau_2 + u_star^2)))

        t = tau * o_gam_i1
        s = 1.0 / (1.0 + sum(t .* t))

        # Equation (12) - Vay 2008
        u_i1 = s * (u_dash + sum(u_dash .* t) * t + cross(u_dash, t))

        # Store updated momentum
        p[1, i] = u_i1[1]
        p[2, i] = u_i1[2]
        p[3, i] = u_i1[3]
    end
    end

    #if TRACE_STATE !== nothing
    #    trace_end(TRACE_STATE, "vay_pusher")
    #end

end

function vay_pusher(p::Array{<:AbstractFloat,2}, w::Array{<:AbstractFloat,1}, 
    B_arr::Array{<:AbstractFloat,2}, E_arr::Array{<:AbstractFloat,2}, 
    dt::Float64, q_m_ratio::Float64, k_F::Float64, c :: Float64)
    #trace_begin("vay_pusher")
    n = size(p, 2)  # Number of particles

    qmt_h  = 0.5 * dt * q_m_ratio  # Precomputed factors
    qmt_hB = qmt_h  * k_F
    o_c2 = 1.0 / c^2.0

    @inbounds for i in 1:n
        if w[i] < 0.0
            continue
        end
        # Extract electric and magnetic fields
        E = @SVector [E_arr[1, i], E_arr[2, i], E_arr[3, i]]
        B = @SVector [B_arr[1, i], B_arr[2, i], B_arr[3, i]]

        # Proper velocity and velocity calculations
        u_i = @SVector [p[1, i], p[2, i], p[3, i]]
        o_gam = 1.0 / sqrt(1.0 + sum(u_i.^2 * o_c2))  # Assuming `pic_scaling%elm%o_c2 = 1.0`
        v_i = u_i * o_gam
        tau = qmt_hB * B

        # Equation (9) - Vay 2008
        u_dash = u_i + 2 * qmt_h * E + cross(v_i, tau)

        gam_dash_2 = 1 + sum(u_dash.^2) * o_c2
        tau_2 = sum(tau .* tau)
        sigma = gam_dash_2 - tau_2
        u_star = sum(u_dash .* tau)

        o_gam_i1 = sqrt(2.0) / sqrt(sigma + sqrt(sigma^2 + 4 * (tau_2 + u_star^2)))

        t = tau * o_gam_i1
        s = 1.0 / (1.0 + sum(t .* t))

        # Equation (12) - Vay 2008
        u_i1 = s * (u_dash + sum(u_dash .* t) * t + cross(u_dash, t))

        # Store updated momentum
        p[1, i] = u_i1[1]
        p[2, i] = u_i1[2]
        p[3, i] = u_i1[3]
    end


    #trace_end("vay_pusher")
end

#end #
