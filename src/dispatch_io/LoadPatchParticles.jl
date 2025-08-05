#---------- Functions to load particle data for specific patches -------------

# Internal helper: compute byte size of a particle patch on disk
function particle_patch_bytes(patch::Particles_NML)
    int_bytes = sizeof(Int32)
    float_bytes = sizeof(Float32)
    per_particle = 4*int_bytes + 8*float_bytes
    return per_particle * sum(patch.M)
end

# Internal helper: skip a particle patch in an open file
function skip_particles!(f::IO, patch::Particles_NML)
    seek(f, position(f) + particle_patch_bytes(patch))
end

#---------------- Load particles for a single patch ----------------
function load_patch_particles(Snapshot_meta::Snapshot_metadata, patch_ID::Int)
    n_patches = Snapshot_meta.N_PARTICLE_PATCHES
    DO_PARTICLES = Snapshot_meta.DO_PARTICLES
    if n_patches == 0 || !DO_PARTICLES
        error("Trying to include particles but no particles in snapshot")
    end

    particle_index = findfirst(p -> p.ID == patch_ID, Snapshot_meta.PARTICLES)
    if particle_index === nothing
        error("Patch ID $patch_ID not found in particle data")
    end
    particle_patch = Snapshot_meta.PARTICLES[particle_index]
    data_file = particle_patch.DATA_FILE
    N_SPECIES = particle_patch.N_SPECIES
    n_particles = particle_patch.M

    all_q = Vector{Array{Int32,2}}()
    all_r = Vector{Array{Float32,2}}()
    all_p = Vector{Array{Float32,2}}()
    all_w = Vector{Array{Float32,1}}()
    all_e = Vector{Array{Float32,1}}()
    all_nr = Vector{Array{Int32,1}}()
    all_pos = Vector{Array{Float32,2}}()
    all_ids = Vector{Array{Int32,1}}()

    for j in 1:N_SPECIES
        Mj = n_particles[j]
        push!(all_q,  zeros(Int32, 3, Mj))
        push!(all_r,  zeros(Float32, 3, Mj))
        push!(all_p,  zeros(Float32, 3, Mj))
        push!(all_w,  zeros(Float32, Mj))
        push!(all_e,  zeros(Float32, Mj))
        push!(all_nr, zeros(Int32, Mj))
        push!(all_pos, zeros(Float32, 3, Mj))
        push!(all_ids, zeros(Int32, Mj))
    end

    f = open(data_file, "r")
    for i in 1:particle_index-1
        prev_patch = Snapshot_meta.PARTICLES[i]
        if prev_patch.DATA_FILE == data_file
            skip_particles!(f, prev_patch)
        end
    end

    patch_index = findfirst(p -> p.ID == patch_ID, Snapshot_meta.PATCHES)
    patch = Snapshot_meta.PATCHES[patch_index]
    snapshot = Snapshot_meta.SNAPSHOT

    for j in 1:N_SPECIES
        q_data = all_q[j];  read!(f, q_data)
        r_data = all_r[j];  read!(f, r_data)
        p_data = all_p[j];  read!(f, p_data)
        w_data = all_w[j];  read!(f, w_data)
        e_data = all_e[j];  read!(f, e_data)
        nr_data = all_nr[j]; read!(f, nr_data)
        all_ids[j] .= Int32(patch_ID)
        global_pos = calc_global_particle_pos(q_data, r_data, patch, snapshot)
        all_pos[j] .= global_pos
    end
    close(f)

    for j in 1:N_SPECIES
        sorted_indices = sortperm(all_nr[j][:])
        all_q[j] = all_q[j][:, sorted_indices]
        all_r[j] = all_r[j][:, sorted_indices]
        all_p[j] = all_p[j][:, sorted_indices]
        all_w[j] = all_w[j][sorted_indices]
        all_e[j] = all_e[j][sorted_indices]
        all_nr[j] = all_nr[j][sorted_indices]
        all_ids[j] = all_ids[j][sorted_indices]
        all_pos[j] = all_pos[j][:, sorted_indices]
    end

    return all_q, all_r, all_p, all_w, all_e, all_nr, all_pos, all_ids
end
#--------------------------------------------------------------------------------

#---------------- Load particles for multiple patches ----------------
function load_patch_particles(Snapshot_meta::Snapshot_metadata, patch_IDs::Vector{Int}; as_list::Bool=false)
    if length(patch_IDs) == 1
        return load_patch_particles(Snapshot_meta, patch_IDs[1])
    end

    if as_list
        q_list = Vector{Any}()
        r_list = Vector{Any}()
        p_list = Vector{Any}()
        w_list = Vector{Any}()
        e_list = Vector{Any}()
        nr_list = Vector{Any}()
        pos_list = Vector{Any}()
        ids_list = Vector{Any}()
        for ID in patch_IDs
            q,r,p,w,e,nr,pos,ids = load_patch_particles(Snapshot_meta, ID)
            push!(q_list, q); push!(r_list, r); push!(p_list, p)
            push!(w_list, w); push!(e_list, e); push!(nr_list, nr)
            push!(pos_list, pos); push!(ids_list, ids)
        end
        return q_list, r_list, p_list, w_list, e_list, nr_list, pos_list, ids_list
    end

    patch_results = [load_patch_particles(Snapshot_meta, ID) for ID in patch_IDs]
    N_SPECIES = Snapshot_meta.PARTICLES[1].N_SPECIES
    total_particles = [0 for _ in 1:N_SPECIES]
    for res in patch_results
        q_patch = res[1]
        for j in 1:N_SPECIES
            total_particles[j] += size(q_patch[j],2)
        end
    end

    all_q = [zeros(Int32, 3, total_particles[j]) for j in 1:N_SPECIES]
    all_r = [zeros(Float32, 3, total_particles[j]) for j in 1:N_SPECIES]
    all_p = [zeros(Float32, 3, total_particles[j]) for j in 1:N_SPECIES]
    all_w = [zeros(Float32, total_particles[j]) for j in 1:N_SPECIES]
    all_e = [zeros(Float32, total_particles[j]) for j in 1:N_SPECIES]
    all_nr = [zeros(Int32, total_particles[j]) for j in 1:N_SPECIES]
    all_pos = [zeros(Float32, 3, total_particles[j]) for j in 1:N_SPECIES]
    all_ids = [zeros(Int32, total_particles[j]) for j in 1:N_SPECIES]

    particle_index = [1 for _ in 1:N_SPECIES]
    for res in patch_results
        q_patch, r_patch, p_patch, w_patch, e_patch, nr_patch, pos_patch, ids_patch = res
        for j in 1:N_SPECIES
            n = size(q_patch[j],2)
            inds = particle_index[j]:particle_index[j]+n-1
            all_q[j][:, inds] = q_patch[j]
            all_r[j][:, inds] = r_patch[j]
            all_p[j][:, inds] = p_patch[j]
            all_w[j][inds] = w_patch[j]
            all_e[j][inds] = e_patch[j]
            all_nr[j][inds] = nr_patch[j]
            all_pos[j][:, inds] = pos_patch[j]
            all_ids[j][inds] = ids_patch[j]
            particle_index[j] += n
        end
    end

    for j in 1:N_SPECIES
        sorted_indices = sortperm(all_nr[j][:])
        all_q[j] = all_q[j][:, sorted_indices]
        all_r[j] = all_r[j][:, sorted_indices]
        all_p[j] = all_p[j][:, sorted_indices]
        all_w[j] = all_w[j][sorted_indices]
        all_e[j] = all_e[j][sorted_indices]
        all_nr[j] = all_nr[j][sorted_indices]
        all_ids[j] = all_ids[j][sorted_indices]
        all_pos[j] = all_pos[j][:, sorted_indices]
    end

    return all_q, all_r, all_p, all_w, all_e, all_nr, all_pos, all_ids
end
#--------------------------------------------------------------------------------



function load_patch_particles(
    Snapshot_meta::Snapshot_metadata,
    patch::Patch_NML
)
    return load_patch_particles(Snapshot_meta, patch.ID)
end

function load_patch_particles(
    Snapshot_meta::Snapshot_metadata,
    patches::Vector{Patch_NML},
    as_list::Bool=false,
)
    IDs = [p.ID for p in patches]
    return load_patch_particles(Snapshot_meta, IDs; as_list=as_list)
end