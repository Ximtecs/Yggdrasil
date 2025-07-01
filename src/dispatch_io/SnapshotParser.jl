function int_to_fixed_length_string(number::Int, length::Int)
    number_str = string(number)
    return lpad(number_str, length, '0')
end

function parse_name_content_NML(content::String)
    matches = eachmatch(r"(?s)&.*?/", content)
    sections  = [match.match[2:end-3] for match in matches]

    params_list = []

    for (i, section) in enumerate(sections)
        lines = split(section, "\n")
        name = strip(lines[1])[1:end-4]
        content = join(lines[2:end], "\n")
        push!(params_list, (name, content))
    end

    return params_list
end
#--------------------------------------------------------------------------------------------------------

#----------- Creates a dictionary with parameters and values for a single section of an nml file ------
function parse_params(params::String)
    lines = split(params, "\n")
    param_dict = Dict{String, String}()
    current_key = ""
    for line in lines
        if contains(line, "=")
            key_value = split(line, "=")
            key = strip(key_value[1])
            value = replace(strip(key_value[2]), " " => "")
            param_dict[key] = value
            current_key = key
        elseif current_key != ""
            value = replace(strip(line), " " => "")
            param_dict[current_key] *= value
        end
    end
    return param_dict
end
#--------------------------------------------------------------------------------------------------------

function parse_SNAPSHOT_NML(params :: String)
    dict = parse_params(params)
    Snapshot_params = SNAPSHOT_NML(
        parse_value(dict["IOFORMAT"], Int),
        parse_value(dict["IOUT"], Int),
        parse_value(dict["TIME"], Float64),
        parse_value(dict["NTOTAL"], Int),
        parse_value(dict["BOX"], Vector{Float64}),
        parse_value(dict["LI"], Vector{Int}),
        parse_value(dict["UI"], Vector{Int}),
        parse_value(dict["NG"], Vector{Int}),
        parse_value(dict["GN"], Vector{Int}),
        parse_value(dict["N"], Vector{Int}),
        parse_value(dict["NV"], Int),
        parse_value(dict["MV"], Int),
        parse_value(dict["NT"], Int),
        parse_value(dict["GAMMA"], Float64),
        parse_value(dict["EOS_NAME"], String),
        parse_value(dict["OPACITY"], String),
        parse_value(dict["PERIODIC"], Vector{Bool}),
        parse_value(dict["GUARD_ZONES"], Bool),
        parse_value(dict["TIME_DERIVS"], Int),
        parse_value(dict["NO_MANS_LAND"], Bool),
        parse_value(dict["OMP_NTHREADS"], Int),
        parse_value(dict["MPI_SIZE"], Int),
        parse_value(dict["MESH_TYPE"], Int),
        parse_value(dict["MPI_DIMS"], Vector{Int}),
        parse_value(dict["REFINE_RATIO"], Int),
        parse_value(dict["ORIGIN"], Vector{Float64})
    )
    return Snapshot_params
end


function parse_IO_NML(params :: String)
    dict = parse_params(params)
    IO_params = IO_NML(
      parse_value(dict["FORMAT"], Int) ,
      parse_value(dict["NTOTAL"], Int) ,
      parse_value(dict["OUT_TIME"], Float64) ,
      parse_value(dict["GUARD_ZONES"], Bool) ,
      parse_value(dict["TIME_DERIVS"], Int) ,
      parse_value(dict["METHOD"], String) ,
      parse_value(dict["NML_VERSION"], Int) ,
      parse_value(dict["DO_GENERIC"], Bool) ,
      parse_value(dict["DO_PIC"], Bool) ,
      parse_value(dict["DO_PARTICLES"], Bool)
        )
    return IO_params
end

function parse_IDX_NML(params :: String)
    dict = parse_params(params)

    #---- we add 1 to get indices starting from 1 instead of 0 ----
    IDX_params = IDX_NML(
        parse_value(dict["D"], Int) + 1,
        parse_value(dict["E"], Int) + 1,
        parse_value(dict["ET"], Int) + 1,
        parse_value(dict["S"], Int) + 1,
        parse_value(dict["PX"], Int) + 1,
        parse_value(dict["PY"], Int) + 1,
        parse_value(dict["PZ"], Int) + 1,
        parse_value(dict["BX"], Int) + 1,
        parse_value(dict["BY"], Int) + 1,
        parse_value(dict["BZ"], Int) + 1,
        parse_value(dict["QR"], Int) + 1,
        parse_value(dict["TT"], Int) + 1,
        parse_value(dict["PHI"], Int) + 1,
        parse_value(dict["P1"], Int) + 1,
        parse_value(dict["P2"], Int) + 1,
        parse_value(dict["P3"], Int) + 1,
        parse_value(dict["B1"], Int) + 1,
        parse_value(dict["B2"], Int) + 1,
        parse_value(dict["B3"], Int) + 1,
        parse_value(dict["EX"], Int) + 1,
        parse_value(dict["EY"], Int) + 1,
        parse_value(dict["EZ"], Int) + 1,
        parse_value(dict["JX"], Int) + 1,
        parse_value(dict["JY"], Int) + 1,
        parse_value(dict["JZ"], Int) + 1,
        parse_value(dict["RPHI"], Int) + 1
    )
    return IDX_params
end

function parse_NBOR_NML(params :: String)
    dict = parse_params(params)

    NBOR_params = NBOR_NML(
        parse_value(dict["PARENT_ID"], Int),
        parse_value(dict["NBOR_IDS"], Vector{Int})
    )

    return NBOR_params
end

function parse_PATCH_NML(params :: String, NBOR_params :: NBOR_NML, data_pos :: Int, data_file :: String)
    dict = parse_params(params)

    PATCH_params = Patch_NML(
        parse_value(dict["ID"], Int),
        parse_value(dict["POSITION"], Vector{Float64}),
        parse_value(dict["SIZE"], Vector{Float64}),
        parse_value(dict["LEVEL"], Int),
        parse_value(dict["DTIME"], Float64),
        parse_value(dict["ISTEP"], Int),
        parse_value(dict["DS"], Vector{Float64}),
        parse_value(dict["NCELL"], Vector{Int}),
        parse_value(dict["N"], Vector{Int}),
        parse_value(dict["NW"], Int),
        parse_value(dict["NV"], Int),
        parse_value(dict["VELOCITY"], Vector{Float64}),
        parse_value(dict["QUALITY"], Float64),
        parse_value(dict["MESH_TYPE"], Int),
        parse_value(dict["KIND"], String),
        parse_value(dict["ETYPE"], String),
        parse_value(dict["RECORD"], Int),
        parse_value(dict["RANK"], Int),
        parse_value(dict["IPOS"], Vector{Int}),
        parse_value(dict["COST"], Float64),
        parse_value(dict["CENTRE_NAT"], Vector{Float64}),
        parse_value(dict["LLC_NAT"], Vector{Float64}),
        parse_value(dict["EROT1"], Vector{Float64}),
        parse_value(dict["EROT2"], Vector{Float64}),
        parse_value(dict["EROT3"], Vector{Float64}),
        NBOR_params,
        data_pos,
        data_file
    )

    return PATCH_params
end 


function parse_PARTICLES_NML(params :: String, data_pos :: Int, data_file :: String)
    dict = parse_params(params)

    Particles_params = Particles_NML(
        parse_value(dict["ID"], Int),
        parse_value(dict["N_SPECIES"], Int),
        parse_value(dict["IS_ELECTRON"], Vector{Bool}),
        parse_value(dict["MASS"], Vector{Float64}),
        parse_value(dict["CHARGE"], Vector{Float64}),
        parse_value(dict["M"], Vector{Int}),
        data_pos,
        data_file
    )
    return Particles_params
end

#--------------- load and store information in the snapshot.nml file ------------------------------
function parse_snapshot_nml(file_path::String)
    content = read(file_path, String)

    params_list = parse_name_content_NML(content)


    IO_content = [content for (name,content) in params_list if name == "IO"][1]
    IDX_content = [content for (name,content) in params_list if name == "IDX"][1]
    SNAPSHOT_content = [content for (name,content) in params_list if name == "SNAPSHOT"][1]


    IO_params = parse_IO_NML(IO_content)
    IDX_params = parse_IDX_NML(IDX_content)
    Snapshot_params = parse_SNAPSHOT_NML(SNAPSHOT_content)

    return IO_params, IDX_params, Snapshot_params
end 
#--------------------------------------------------------------------------------------------------------



#--------------- parse a single section with both NBOR_NML and Patch_NML information ------------------------------
function parse_patch(patch_param :: String, nbor_param :: String, data_pos :: Int, data_file :: String)
    NBOR_params = parse_NBOR_NML(nbor_param)
    Patch_params = parse_PATCH_NML(patch_param, NBOR_params, data_pos, data_file)

    return Patch_params
end 
#--------------------------------------------------------------------------------------------------------





#----------------- parse all pathces in a patches.nml file ------------------------------
function parse_patches_nml(file_path::String, data_file::String)
    content = read(file_path, String)
    params_list = parse_name_content_NML(content)
    params_list = filter(x -> x[1] != "IDX", params_list) # remove all IDX_NML sections
    #params_list = params_list[2:end] # remove the IDX_NML section


    n_sections = size(params_list)[1]
    n_patches = Int(n_sections/2)
    patches_params = []

    n_pic_patches = 0
    data_pos = 1
    for i in 1:2:n_sections
        patch_param = params_list[i][2]
        nbor_param = params_list[i+1][2]

        PATCH_params = parse_patch(patch_param, nbor_param, data_pos, data_file)

        if PATCH_params.KIND == "PIC"
            n_pic_patches += 1
        end
        push!(patches_params, PATCH_params)
        data_pos = data_pos + 1
    end 

    return patches_params, n_patches, n_pic_patches
end
#--------------------------------------------------------------------------------------------------------

function parse_particles_nml(file_path::String, data_file::String)
    content = read(file_path, String)
    params_list = parse_name_content_NML(content)

    #------ number of patches with particle information ------
    n_sections = size(params_list)[1]
    #---------------------------------------------------------

    particles_params = []

    data_pos = 1
    for i in 1:n_sections
        particle_param = params_list[i][2]

        PARTICLE_paramms = parse_PARTICLES_NML(particle_param, data_pos, data_file)

        push!(particles_params, PARTICLE_paramms)

        data_pos = data_pos + 1
    end

    return particles_params, n_sections
end



function find_system(content)
    lines = split(content, "\n")
    for line in lines
        if occursin("SYSTEM", line)
            return line
        end
    end
end 

function parse_unit_system(file)
    content = read(file, String)
    system_line = find_system(content)
    system = split(system_line, "=")[2]
    system = system[2:end-2]
    system = strip(system) 
    return system 
end 


#-------------- parse all meta information from a snapshot folder ------------------------------
function read_snapshot(data_folder :: String, snap :: Int)

    snap_folder = data_folder * int_to_fixed_length_string(snap, 5) * "/"


    params_file = data_folder * "params.nml"
    system = parse_unit_system(params_file)

    snapshot_nml_file = snap_folder * "snapshot.nml"
    IO_params, IDX_params, Snapshot_params = parse_snapshot_nml(snapshot_nml_file)


    patches_params = []
    particles_params = Vector{Particles_NML}()


    N_PARTICLES= Vector{Int}()

    DO_PIC = IO_params.DO_PIC
    DO_PARTICLES = IO_params.DO_PARTICLES
    n_patches = 0
    n_pic_patches = 0
    n_particle_patches = 0

    if DO_PARTICLES    
        particle_folder = data_folder * "particles/" * int_to_fixed_length_string(snap, 5) * "/"
    else
        particle_folder = ""
    end


    for MPI_rank in 0:Snapshot_params.MPI_SIZE-1
        #patches_nml_file = snap_folder * "rank_" * int_to_fixed_length_string(MPI_rank, 5) * "_patches.nml"
        patches_nml_file = snap_folder * "snapshot_" * int_to_fixed_length_string(MPI_rank, 5) * ".nml"
        data_file = snap_folder * "snapshot_" * int_to_fixed_length_string(MPI_rank, 5) * ".dat"



        patches_params_rank, n_patches_rank, n_pic_patches_rank = parse_patches_nml(patches_nml_file,data_file)
        patches_params = vcat(patches_params, patches_params_rank)
        n_patches = n_patches + n_patches_rank
        n_pic_patches = n_pic_patches + n_pic_patches_rank

        if DO_PARTICLES && isdir(particle_folder)  # Check if the folder exists
            particles_nml_file = particle_folder * "rank_" * int_to_fixed_length_string(MPI_rank, 5) * "_particles.nml"
            particles_data_file = particle_folder * "rank_" * int_to_fixed_length_string(MPI_rank, 5) * ".data"
            
            # Only proceed if the expected files exist
            if isfile(particles_nml_file) && isfile(particles_data_file)
                particles_params_rank, n_particle_patches_rank = parse_particles_nml(particles_nml_file, particles_data_file)
                n_particle_patches += n_particle_patches_rank
                particles_params = vcat(particles_params, particles_params_rank)
            end
        end
        

        # if DO_PARTICLES
        #     particles_nml_file = particle_folder * "rank_" * int_to_fixed_length_string(MPI_rank, 5) * "_particles.nml"
        #     particles_data_file = particle_folder * "rank_" * int_to_fixed_length_string(MPI_rank, 5) * ".data"
        #     particles_params_rank, n_particle_patches_rank = parse_particles_nml(particles_nml_file, particles_data_file)
        #     n_particle_patches += n_particle_patches_rank
        #     particles_params = vcat(particles_params, particles_params_rank)
        # end
    end 

    level_min = 1000
    level_max = 0
    for patch in patches_params
        if patch.LEVEL < level_min
            level_min = patch.LEVEL
        end
        if patch.LEVEL > level_max
            level_max = patch.LEVEL
        end
    end 


    NV_PIC = 27
    if n_pic_patches > 0
        for patch in patches_params
            if patch.KIND == "PIC"
                NV_PIC = patch.NV
                break
            end
        end
    end
    NV_MHD = 8
    if n_pic_patches > 0
        for patch in patches_params
            if patch.KIND != "PIC"
                NV_MHD = patch.NV
                break
            end
        end
    end

    if length(particles_params) > 0
        N_SPECIES = particles_params[1].N_SPECIES
        for i in 1:N_SPECIES
            push!(N_PARTICLES, 0)
        end

        for particle in particles_params
            M = particle.M
            for i in 1:N_SPECIES
                N_PARTICLES[i] += M[i]
            end
        end

    else
        push!(N_PARTICLES, 0)
    end





    Snapshot_meta = Snapshot_metadata(IO_params, Snapshot_params, IDX_params,
                                     n_patches, n_pic_patches, patches_params, snap_folder,
                                     DO_PIC, DO_PARTICLES, NV_PIC, NV_MHD, N_PARTICLES, n_particle_patches,
                                     particle_folder, particles_params, system, level_min, level_max)

    return Snapshot_meta
end 
#--------------------------------------------------------------------------------------------------------
