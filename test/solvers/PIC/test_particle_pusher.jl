l                 = 720.0847332647693
n                 = 1.0e9
B                 = 61.60986870674753
T                 = 6.873175661794552e9

l_s           = 720.0847332647693
ρ_s     = 1.6893481169800001e-15
t_s             = 2.401944125174654e-8
temp_s      = 1.0

me_s          = 18.361526737600673
q_s           = 1.0
ϵ_s           = 1.0
μ_s           = 1.0


base = BaseUnits("HL"; ϵ_s=μ_s, μ_s=μ_s, me_s=me_s, q_s=q_s)

ds = 0.1
PPC = 128


scaling = ScalingHL(base, n, l, T, T, T, B, 
                    l_s, ρ_s, t_s, temp_s)
set_macro_particle_weights(scaling, ds, PPC)

c = scaling.c_code
k_F = scaling.k_F
q_m_ratio = scaling.e_code / scaling.m_e_code
dt = 0.25 * ds / scaling.c_code



@testset "vay const E for ParticlePusher" begin
    @trace "test" begin

    n_p = 10000

    w = ones(Float32 ,n_p)
    p = rand(Float32, 3, n_p) ./ 100 


    E = zeros(Float32, 3,n_p)
    B = zeros(Float32, 3,n_p)


    B0 = scaling.B_flux_code
    E0 = 1e-2 * B0
    E[1,:] .=  1.0 * E0
    E[2,:] .=  2.0 * E0
    E[3,:] .= -1.0 * E0

    n_steps = 100
    p_copies = []
    

    args = []


    args = (; dt=dt, q_m_ratio=q_m_ratio, k_F=k_F, c=c)#, TRACE_STATE=TRACE_STATE)  # Named tuple


    for i in 1:n_steps
        #vay_pusher(p, B, E, dt, q_m_ratio, k_F, c)  # In-place update
        vay_pusher(p, B, E; args...)  # In-place update

        if i % 10 == 0 
            push!(p_copies, copy(p))  # Store a snapshot every 10 steps
        end
    end
    
    # Compute differences in momentum every 10 steps
    diff = [p_copies[i] .- p_copies[i-1] for i in 2:length(p_copies)]
    
    # Theoretical momentum difference per 10 steps
    p_diff_theo = [1.0, 2.0, -1.0] .* E0 .* (10 * dt) .* q_m_ratio  # (3,)
    
    # Verify that every particle's momentum difference matches the theory
    for d in diff  # Loop over stored snapshots (each d is a (3, n_particles) array)
        @test all(isapprox.(d, p_diff_theo; rtol=1e-5))  # Check all (3, n_particles) entries
    end

    end
    print_timings(GLOBAL_TIMER)
 
        
end 