function load_multiple_snapshots(data_folder :: String, i_start :: Int, i_end :: Int, stride :: Int, load_pic :: Bool; use_level::Union{Nothing, Int} = nothing)

    
    #---------------- count number of snapshots ------------
    n_snaps = Int(floor( (i_end - i_start) / stride) + 1)
    #------------------------------------------------------
    #------  Get initial snapshot and memory size  --------
    initial_snapshot = read_snapshot(data_folder, i_start);
    level = isnothing(use_level) ? initial_snapshot.LEVELMIN : use_level
    mem_size = get_mem_size(initial_snapshot, level)
    #------------------------------------------------------

    #------ Allocate memory for all snapshots -------------
    all_data = zeros(Float32,mem_size..., n_snaps)
    all_t = zeros(Float32,n_snaps)
    #------------------------------------------------------


    #-------- Loop through all snapshots and load them into the allocated memory ---------
    Threads.@threads for snap in i_start:stride:i_end
        iter = Int(((snap - i_start) / stride) + 1 )
        Snapshot_meta = read_snapshot(data_folder, snap);
        all_data[:,:,:,:,iter] = load_snapshot(Snapshot_meta, load_pic; use_level=level)
        all_t[iter] = Snapshot_meta.SNAPSHOT.TIME
    end
    #---------------------------------------------------------------------------------------

    return all_data, all_t
end

function load_multiple_snapshots(data_folder :: String, i_start :: Int, i_end :: Int, stride :: Int, llc ::Vector{Int}, urc:: Vector{Int}, load_pic :: Bool; use_level::Union{Nothing, Int} = nothing)

    
    #---------------- count number of snapshots ------------
    n_snaps = Int(floor( (i_end - i_start) / stride) + 1)
    #------------------------------------------------------
    #------  Get initial snapshot and memory size  --------
    initial_snapshot = read_snapshot(data_folder, i_start);
    level = isnothing(use_level) ? initial_snapshot.LEVELMIN : use_level
    mem_size = get_mem_size(initial_snapshot, level)
    #------------------------------------------------------

    mem_size_red = urc .- llc .+ 1
    mem_size[1:3] = mem_size_red

    #------ Allocate memory for all snapshots -------------
    all_data = zeros(Float32,mem_size..., n_snaps)
    all_t = zeros(Float32,n_snaps)
    #------------------------------------------------------


    #-------- Loop through all snapshots and load them into the allocated memory ---------
    Threads.@threads for snap in i_start:stride:i_end
        iter = Int(((snap - i_start) / stride) + 1 )
        Snapshot_meta = read_snapshot(data_folder, snap);
        all_data[:,:,:,:,iter] = load_snapshot(Snapshot_meta, load_pic, llc, urc; use_level=level)
        all_t[iter] = Snapshot_meta.SNAPSHOT.TIME
    end
    #---------------------------------------------------------------------------------------

    return all_data, all_t
end

