function get_scaleHL(r_s, real_units, c_code, ds, per_cell, scale_type)
    base = BaseUnits("HL")
    #------ calculate electron mass scaing based on ratio of electron to proton mass -----
    electron_mass_scaling = (base.m_p / base.m_e) / r_s.electron_proton_ratio
    #-------------------------------------------------------------------------------------
    scale_base_units(base, r_s.eps_0_scaling, 1.0, electron_mass_scaling, r_s.q_charge_scaling)
    
    
    mu = ( (base.m_e + base.m_p)  / 2. ) /  base.m_u# assuming only protons and electrons


    B = real_units.B_field
    T = real_units.T
    T_elec = T
    T_ion = T
    n_p = real_units.number_density 
    rho = real_units.mass_density

    if (n_p> 0.0) && (rho> 0.0) 
        error("Both number density and mass density are positive. Only one should be positive")
    end 

    if (n_p > 0.0) 
        rho = n_p * 2 * mu * base.m_u
    elseif (rho > 0.0) 
        n_p = rho / (2 * mu * base.m_u)
    else
        error("Both number density and mass density are negative. One should be positive")
    end 


    rho_scale = rho 
    T_scale= 1.0

    #------- scale type 1 we set timescale to electron plasma freqeuncy, then set c ------
    if (scale_type == 1) 
        e_plasma = base.e * sqrt(n_p / base.m_e)
        time_scale         = 2pi / e_plasma
        l_scale            = time_scale * base.c  /   c_code
    elseif (scale_type == 2) 
        p_plasma = base.e * sqrt(n_p / base.m_p)
        p_skin_depth = base.c / p_plasma

        l_scale = p_skin_depth
        time_scale = l_scale / base.c * c_code
    elseif (scale_type == 3) 
        e_plasma = base.e * sqrt(n_p / base.m_e)
        e_skin_depth = base.c / e_plasma

        l_scale = e_skin_depth
        time_scale = l_scale / base.c * c_code
    else
        error("scale type not implemented")
    end

    L = l_scale


    scaling = ScalingHL(base,
                    n_p, L, T, T_elec, T_ion, B,
                    l_scale, rho_scale, time_scale, T_scale)

    set_macro_particle_weights(scaling, ds, per_cell)
    return scaling
end

function get_scaleCGS(r_s, real_units, c_code, ds, per_cell, scale_type)
    base = BaseUnits("CGS")
    #------ calculate electron mass scaing based on ratio of electron to proton mass -----
    electron_mass_scaling = (base.m_p / base.m_e) / r_s.electron_proton_ratio
    #-------------------------------------------------------------------------------------
    scale_base_units(base, r_s.eps_0_scaling, 1.0, electron_mass_scaling, r_s.q_charge_scaling)
    
    
    mu = ( (base.m_e + base.m_p)  / 2. ) /  base.m_u# assuming only protons and electrons


    B = real_units.B_field
    T = real_units.T
    T_elec = T   # assuming T_elec = T
    T_ion = T    # assuming T_ion = T
    n_p = real_units.number_density 
    rho = real_units.mass_density

    #------------ check that only one of number density or mass density is positive ---------
    if (n_p> 0.0) && (rho> 0.0) 
        error("Both number density and mass density are positive. Only one should be positive")
    end 
    #-------------------------------------------------------------------------------------
    #-------------- get number density from mass density or vice versa -------------------
    if (n_p > 0.0) 
        rho = n_p * 2 * mu * base.m_u
    elseif (rho > 0.0) 
        n_p = rho / (2 * mu * base.m_u)
    else
        error("Both number density and mass density are negative. One should be positive")
    end 
    #-------------------------------------------------------------------------------------

    #------------ always set rho_Scale so rho=1 in code units ------------------------------
    rho_scale = rho 
    #-------------------------------------------------------------------------------------
    #------------- don't scale temperature for now ----------------------------------------
    T_scale= 1.0 
    #-------------------------------------------------------------------------------------

    #------- scale type 1 we set timescale to electron plasma freqeuncy, then set c ------
    if (scale_type == 1) 
        e_plasma = base.e * sqrt(4pi * n_p / base.m_e)
        time_scale         = 2pi / e_plasma
        l_scale            = time_scale * base.c  /   c_code
    #------ scale type 2 we set length scale to proton skin depth, then set time scale get correct c_code------
    elseif (scale_type == 2) 
        p_plasma = base.e * sqrt(4pi * n_p / base.m_p)
        p_skin_depth = base.c / p_plasma

        l_scale = p_skin_depth
        time_scale = l_scale / base.c * c_code
    #------ scale type 3 we set length scale to electron skin depth, then set time scale get correct c_code------
    elseif (scale_type == 3) 
        e_plasma = base.e * sqrt(4pi * n_p / base.m_e)
        e_skin_depth = base.c / e_plasma

        l_scale = e_skin_depth
        time_scale = l_scale / base.c * c_code
    else
        error("scale type not implemented")
    end

    L = l_scale


    scaling = ScalingCGS(base,
                    n_p, L, T, T_elec, T_ion, B,
                    l_scale, rho_scale, time_scale, T_scale)

    set_macro_particle_weights(scaling, ds, per_cell)
    return scaling
end