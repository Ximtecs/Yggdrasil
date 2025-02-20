struct ParticleContainer{T<:AbstractFloat}
    p :: Array{T, 3}  # velocities, should be (3, n)
    r :: Array{T, 3}  # positions, should be (3, n)     - between -0.5 and 0.5
    w :: Array{T, 2}  # weights, should be (n)
    q :: Vector{Int}  # integer position, common for all particles

    qm :: Float64  # charge to mass ratio


    #--------- bookkeeping variables ---------------------------------------
    nt         :: Int  # number of timeslots
    it         :: Int  # current timeslot
    new        :: Int  # new timeslot
    n_target   :: Int  # target number of particles
    n          :: Int  # current size of container

    
    m          :: Vector{Int}  # current number of particles
    first_dead :: Vector{Int}  # index of the first dead particle
    last_alive :: Vector{Int} # index of the last alive particle
    #-----------------------------------------------------------------------

    #--------------- Constructors --------------------------------------------
    function ParticleContainer(q::Vector{Int}, qm::Float64, n::Int, nt :: Int)
        if length(q) != 3
            throw(ArgumentError("q must be a 1D vector with 3 elements"))
        end

        n_target = n
        n        = n * 1.2 # 20% more than the target number

        p = zeros(Float64, 3, n, nt)
        r = zeros(Float64, 3, n, nt)
        w = -ones(Float64, n, nt)


        m         = zeros(Int, nt)
        first_dead = ones(Int, nt)
        last_alive = zeros(Int, nt)

        it = 1
        new_it = min(2, nt)


        new{Float64}(p, r, w, q, qm, 
                    nt, it, new_it,
                     n_target, n, m, first_dead, last_alive)
    end
    #--------------------------------------------------------------------------
end

function compress_particles!(pc::ParticleContainer, timeslot::Int)
    if timeslot < 1 || timeslot > pc.nt
        throw(ArgumentError("Invalid timeslot"))
    end

    while pc.first_dead[timeslot] <= pc.last_alive[timeslot]
        dead_idx = pc.first_dead[timeslot]
        alive_idx = pc.last_alive[timeslot]
        # Swap dead and alive particles
        pc.r[:, dead_idx, timeslot], pc.r[:, alive_idx, timeslot] = pc.r[:, alive_idx, timeslot], pc.r[:, dead_idx, timeslot]
        pc.p[:, dead_idx, timeslot], pc.p[:, alive_idx, timeslot] = pc.p[:, alive_idx, timeslot], pc.p[:, dead_idx, timeslot]
        pc.w[dead_idx, timeslot], pc.w[alive_idx, timeslot] = pc.w[alive_idx, timeslot], pc.w[dead_idx, timeslot]

        update_indices(pc, dead_idx, timeslot)
        update_indices(pc, alive_idx, timeslot)
    end

    # Check that compression has been done properly
    if pc.first_dead[timeslot] <= pc.last_alive[timeslot]
        throw(AssertionError("Compression failed: first_dead should be greater than last_alive"))
    end
end


function decrease_size!(pc::ParticleContainer, factor::Float64)
    new_n = Int(pc.n * factor)
    
    #----------------- first check if we can decrease the size -----------------
    for timeslot in 1:pc.nt
        if pc.m[timeslot] > 0.9 * new_n
            println("Timeslot $timeslot has more than 0.9 * new_n particles, skipping size decrease for this timeslot")
            return
        end
    end
    #---------------------------------------------------------------------------

    #-------------- chcek if compression is needed in order to decrease size ---
    for timeslot in 1:pc.nt
        if pc.last_alive[timeslot] > new_n
            compress_particles!(pc, timeslot)
        end
    end
    #-----------------------------------------------------------------------------

    new_p = zeros(eltype(pc.p), 3, new_n, pc.nt)
    new_r = zeros(eltype(pc.r), 3, new_n, pc.nt)
    new_w = -ones(eltype(pc.w), new_n, pc.nt)

    new_p[:, 1:new_n, :] .= pc.p[:, 1:new_n, :]
    new_r[:, 1:new_n, :] .= pc.r[:, 1:new_n, :]
    new_w[1:new_n, :] .= pc.w[1:new_n, :]

    pc.p = new_p
    pc.r = new_r
    pc.w = new_w
    pc.n = new_n
end


function increase_size!(pc::ParticleContainer, factor::Float64)
    new_n = Int(pc.n * factor)
    new_p = zeros(eltype(pc.p), 3, new_n, pc.nt)
    new_r = zeros(eltype(pc.r), 3, new_n, pc.nt)
    new_w = -ones(eltype(pc.w), new_n, pc.nt)

    new_p[:, 1:pc.n, :] .= pc.p
    new_r[:, 1:pc.n, :] .= pc.r
    new_w[1:pc.n, :] .= pc.w

    # Delete old arrays
    empty!(pc.p)
    empty!(pc.r)
    empty!(pc.w)

    pc.p = new_p
    pc.r = new_r
    pc.w = new_w
    pc.n = new_n
end




function add_particle!(pc::ParticleContainer, position::Vector{T}, velocity::Vector{T}, weight::T, timeslot::Int) where T<:AbstractFloat
    if length(position) != 3 || length(velocity) != 3
        throw(ArgumentError("position and velocity must be 3-element vectors"))
    end

    if pc.n >= pc.n_target
        throw(ArgumentError("Particle container is full"))
    end

    if timeslot < 1 || timeslot > pc.nt
        throw(ArgumentError("Invalid timeslot"))
    end

    i = pc.first_dead[timeslot]
    pc.r[:, i, timeslot] = position
    pc.p[:, i, timeslot] = velocity
    pc.w[i, timeslot] = weight

    update_indices(pc, i, timeslot)
    pc.m[timeslot] += 1
end

function remove_particle!(pc::ParticleContainer, index::Int, timeslot::Int)
    if index < 1 || index > pc.n
        throw(ArgumentError("Invalid particle index"))
    end

    if timeslot < 1 || timeslot > pc.nt
        throw(ArgumentError("Invalid timeslot"))
    end

    pc.w[index, timeslot] = -1.0  # Mark the particle as dead
    pc.p[index, timeslot] = 0.0
    pc.r[index, timeslot] = 0.0

    update_indices(pc, index, timeslot)
    pc.m[timeslot] -= 1
end




#---------------------- bookkeeping - keep track of last alive and first dead particles ---------------------
function update_indices!(pc::ParticleContainer, i::Int, nt :: Int)
    update_last_alive!(pc, i, nt)
    update_first_dead!(pc, i, nt)
end
function update_last_alive!(pc::ParticleContainer, i::Int, nt :: Int)
    if i > pc.last_alive[nt]
        if pc.w[i,nt] > 0.0
            pc.last_alive[nt] = i
        end 
    end

    if i == pc.last_alive[nt]
        if pc.w[i,nt] <= 0.0
            pc.last_alive[nt] = findlast(x -> x > 0.0, pc.w[:,nt]) #TODO - should this be optimized=
        end
    end
end
function update_first_dead!(pc::ParticleContainer, i::Int, nt :: Int)
    if i < pc.first_dead[nt]
        if pc.w[i,nt] <= 0.0
            pc.first_dead[nt] = i
        end
    end
    if i == pc.first_dead[nt]
        if pc.w[i,nt] > 0.0
            pc.first_dead[nt] = findfirst(x -> x <= 0.0, pc.w[:,nt]) #TODO - should this be optimized=
        end
    end
end
#-----------------------------------------------------------------------------------------------------------------