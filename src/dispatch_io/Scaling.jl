function get_dispatch_scaling(scaling_nml :: SCALING_NML)
    system = scaling_nml.SYSTEM

    n = scaling_nml.N_P_REAL

    l = scaling_nml.L_REAL
    T = scaling_nml.TEMP_REAL
    Te = scaling_nml.ELECTRON_TEMP_REAL
    Ti = scaling_nml.ION_TEMP_REAL
    B = scaling_nml.B_REAL


    l_scale = scaling_nml.L
    d_scale = scaling_nml.D
    t_scale = scaling_nml.T
    temp_scale = 1.0 # assume no temperature scaling

    ds = scaling_nml.DS
    PPC = scaling_nml.PER_CELL


    q_e_scale = scaling_nml.Q_E_SCALE
    m_e_scale = scaling_nml.M_E_SCALE
    mu_0_scale = scaling_nml.MU_0_SCALE
    eps_0_scale = scaling_nml.EPS_0_SCALE

    if system == "HL"
        base_units = BaseUnits("HL")
        scale_base_units(base_units, eps_0_scale, mu_0_scale, m_e_scale, q_e_scale)
        scaling = ScalingHL(base_units, n, l, T, Te, Ti, B,
                         l_scale, d_scale, t_scale, temp_scale)
        set_macro_particle_weights(scaling, ds, PPC)
        return scaling

    elseif system == "CGS"
        base_units = BaseUnits("CGS")
        scale_base_units(base_units, eps_0_scale, mu_0_scale, m_e_scale, q_e_scale)
        scaling = ScalingCGS(base_units, n, l, T, Te, Ti, B,
                         l_scale, d_scale, t_scale, temp_scale)
        set_macro_particle_weights(scaling, ds, PPC)
        return scaling
    else
        error("Unknown system: $system")
    end 
end 


function get_dispatch_scaling(snapshot :: Snapshot_metadata)
    scaling_nml = snapshot.SCALING
    return get_dispatch_scaling(scaling_nml)
end 