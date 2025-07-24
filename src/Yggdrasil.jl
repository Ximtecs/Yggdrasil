 module Yggdrasil
    
    using Printf
    using TimerOutputs
    using LinearAlgebra
    using LoopVectorization
    using StaticArrays
    using Interpolations
    using PyPlot
    using PyCall


    include("io/Timer_mod.jl")
    include("dispatch_io/DispatchIO.jl")
    include("vector_ops/VectorOps.jl")
    include("scaling/Scaling.jl")
    include("solvers/PIC/ParticlePusher.jl")
    include("processing/Processing.jl")
    include("plotting/Plotting.jl")


    #---------- folder ---------------------------------------- 
        #------- macros -------------------
        #----------------------------------
        #-------- vars --------------------
        #----------------------------------
        #------- structs ------------------
        #----------------------------------
        #------- functions ----------------
        #----------------------------------
    #---------------------------------------------------------- 


    #---------- IO --------------------------------------------- 
        #------- macros -------------------
    export @trace
        #----------------------------------
        #-------- vars --------------------
    export GLOBAL_TIMER, ENABLE_PROFILING
        #----------------------------------
        #------- structs ------------------
        #----------------------------------
        #------- functions ----------------
    export print_timings
        #----------------------------------
    #---------------------------------------------------------- 

    #---------- DISPATCH IO ---------------------------------------- 
        #------- macros -------------------
        #----------------------------------
        #-------- vars --------------------
        #----------------------------------
        #------- structs ------------------
        export IDX_NML, IO_NML, SNAPSHOT_NML, NBOR_NML, Patch_NML,
              Particles_NML, Snapshot_metadata
        #----------------------------------
        #------- functions ----------------
        export load_patch, load_snapshot, load_snapshot_particles
        export load_multiple_snapshots_particles, load_multiple_snapshots
        export read_snapshot
        export get_idx_value, get_xyz, get_ds
        export get_area
        export print_info
        export drop_unit_dims
        export find_patch
        export find_patches_line
        export find_patches_square
        export find_patches_cube
        export get_dispatch_scaling
        #----------------------------------
    #---------------------------------------------------------- 

    #---------- VectorOps -------------------------------------
        #------- macros -------------------
        #----------------------------------
        #-------- vars --------------------
        #----------------------------------
        #------- structs ------------------
    export Backend, CPU, GPU
        #----------------------------------
        #------- functions ----------------
    export ddx_dn, ddx_up, ddy_dn, ddy_up, ddz_dn, ddz_up, curl_dn, curl_up
        #----------------------------------
    #---------------------------------------------------------- 


    #---------- Scaling ---------------------------------------- 
        #------- macros -------------------
        #----------------------------------
        #-------- vars --------------------
        #----------------------------------
        #------- structs ------------------
    export BaseUnits, ScalingHL, ScalingCGS
    export code_units, known_real_units, real_unit_input, ratios_and_scaling
        #----------------------------------
        #------- functions ----------------
    export find_real_units_CGS, find_real_units_HL # to go from code to real units
    export get_scaleHL, get_scaleCGS # to get code unit from real units
    export set_macro_particle_weights # for setting particle per cell and cell size to Scaling 
    export print_basic_info, print_all_HL, print_all_CGS # for printing scaling info
        #----------------------------------
    #---------------------------------------------------------- 

    #---------- LineAnalysis ---------------------------------------- 
        #------- macros -------------------
        #----------------------------------
        #-------- vars --------------------
        #----------------------------------
        #------- structs ------------------
        #----------------------------------
        #------- functions ----------------
        export get_line, generate_line_plane
        export interpolate_to_line, interpolate_to_plane
        export mom_to_vel!, vel_to_mom!
        export rotate_aligned_vectors! 
        #----------------------------------
    #---------------------------------------------------------- 

    #---------- Plotting ---------------------------------------- 
        #------- macros -------------------
        #----------------------------------
        #-------- vars --------------------
        export default_color_palette
        #----------------------------------
        #------- structs ------------------
        #----------------------------------
        #------- functions ----------------
        export plot_heatmap, plot_heatmaps
        export plot_line, plot_lines
        #----------------------------------
    #---------------------------------------------------------- 

    #---------- solvers/PIC -----------------------------------
        #------- macros -------------------
        #----------------------------------
        #-------- vars --------------------
        #----------------------------------
        #------- structs ------------------
        #----------------------------------
        #------- functions ----------------
    export vay_pusher
        #----------------------------------
    #---------------------------------------------------------- 

end
