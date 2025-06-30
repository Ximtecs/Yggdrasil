 module Yggdrasil
    
    using Printf
    using TimerOutputs
    using LinearAlgebra
    using LoopVectorization
    using StaticArrays


     include("io/Timer_mod.jl")
    include("vector_ops/VectorOps.jl")
    include("scaling/Scaling.jl")
    include("solvers/PIC/ParticlePusher.jl")


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
