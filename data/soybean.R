soybean <- list(
    direct_modules = list(
        "BioCro:format_time",
        stomata_water_stress = "BioCro:stomata_water_stress_linear",
        specific_leaf_area = "BioCro:sla_linear",
        "BioCro:parameter_calculator",
        "BioCro:soybean_development_rate_calculator",
        leaf_water_stress = "BioCro:leaf_water_stress_exponential",
        partitioning_coefficients = "BioCro:partitioning_coefficient_logistic",
        "BioCro:soil_evaporation",
        solar_coordinates = "BioCro:solar_position_michalsky",
        "BioCro:shortwave_atmospheric_scattering",
        "BioCro:incident_shortwave_from_ground_par",
        "BioCro:height_from_lai",
        "BioCro:canopy_gbw_thornley",
        "BioCro:stefan_boltzmann_longwave",
        "BioCro:ten_layer_canopy_properties",
        canopy_photosynthesis = "BioCro:ten_layer_c3_canopy",
        "BioCro:ten_layer_canopy_integrator",
        partitioning_growth_calculator = "BioCro:partitioning_growth_calculator",
        "BioCro:senescence_coefficient_logistic",
        "BioCro:carbon_assimilation_to_biomass",
        "BioCro:maintenance_respiration_calculator"
    ),
    differential_modules = list(
        senescence = "BioCro:senescence_logistic",
        "BioCro:maintenance_respiration",
        "BioCro:partitioning_growth",
        soil_profile = "BioCro:two_layer_soil_profile",
        "BioCro:development_index",
        thermal_time = "BioCro:thermal_time_linear"
    ),
    ode_solver = list(
        type = 'boost_rkck54',
        output_step_size = 1.0,
        adaptive_rel_error_tol = 1e-4,
        adaptive_abs_error_tol = 1e-4,
        adaptive_max_steps = 200
    ),
    initial_values = list(
        Leaf               = 0.06312,
        Stem               = 0.00789,
        Root               = 0.00789,
        Grain              = 0.00001,
        Shell              = 0.00001,
        LeafLitter         = 0,
        RootLitter         = 0,
        StemLitter         = 0,
        soil_water_content = 0.32,
        cws1               = 0.32,
        cws2               = 0.32,
        DVI                = -1,
        TTc                = 0,
        Rhizome            = 0.0000001,
        RhizomeLitter      = 0
    ),
    parameters = list(
        # soil parameters (clay loam)
        soil_air_entry              = -2.6,
        soil_b_coefficient          = 5.2,
        soil_bulk_density           = 1.35,
        soil_clay_content           = 0.34,
        soil_field_capacity         = 0.32,
        soil_sand_content           = 0.32,
        soil_saturated_conductivity = 6.4e-05,
        soil_saturation_capacity    = 0.52,
        soil_silt_content           = 0.34,
        soil_wilting_point          = 0.2,

        # sla_linear module
        iSp                         = 3.5,
        Sp_thermal_time_decay       = 0,

        # parameter_calculator module
        alpha1                      = 0,
        alphab1                     = 0,
        LeafN                       = 2,
        LeafN_0                     = 2,
        Vcmax_at_25                 = 110,

        # soybean_development_rate_calculator module
        maturity_group              = 3,
        Tbase_emr                   = 10,
        TTemr_threshold             = 60,
        Rmax_emrV0                  = 0.1990,
        Tmin_emrV0                  = 5.0,
        Topt_emrV0                  = 31.5,
        Tmax_emrV0                  = 45.0,
        Tmin_R0R1                   = 5.0,
        Topt_R0R1                   = 31.5,
        Tmax_R0R1                   = 45.0,
        Tmin_R1R7                   = 0.0,
        Topt_R1R7                   = 21.5,
        Tmax_R1R7                   = 38.7,
        sowing_fractional_doy       = 0,

        # partitioning_coefficient_logistic module
        alphaLeaf                   = 24.7116,
        alphaRhizome                = 0,
        alphaRoot                   = 36.9670,
        alphaShell                  = 10.8835,
        alphaStem                   = 24.5764,
        betaLeaf                    = -19.2275,
        betaRhizome                 = -Inf,
        betaRoot                    = -40.1915,
        betaShell                   = -7.9549,
        betaStem                    = -18.3517,
        kRhizome_emr                = 0,
        kRhizome_emr_DVI            = 0,

        # soil_evaporation module
        rsec                        = 0.2,
        soil_clod_size              = 0.04,
        soil_reflectance            = 0.2,
        soil_transmission           = 0.01,
        specific_heat_of_air        = 1010,

        # solar_position_michalsky module
        lat                         = 40,
        longitude                   = -88,

        # shortwave_atmospheric_scattering module
        atmospheric_pressure        = 101325,
        atmospheric_transmittance   = 0.6,
        atmospheric_scattering      = 0.3,

        # incident_shortwave_from_ground_par module
        par_energy_fraction         = 0.5,
        par_energy_content          = 0.219,

        # height_from_lai module
        heightf                     = 6,

        # canopy_gbw_thornley module
        min_gbw_canopy              = 0.005,

        # carbon_asismilation_to_biomass module
        dry_biomass_per_carbon = 30.026,

        # stefan_boltzmann_longwave module
        emissivity_sky              = 1,

        # ten_layer_canopy_properties module
        chil                        = 0.81,
        k_diffuse                   = 0.7,
        kpLN                        = 0,
        leaf_reflectance_nir        = 0.42,
        leaf_reflectance_par        = 0.10,
        leaf_transmittance_nir      = 0.42,
        leaf_transmittance_par      = 0.05,
        lnfun                       = 0,

        # ten_layer_c3_canopy module (temperature response)
        Gstar_c    = 19.02,
        Gstar_Ea   = 37.83e3,
        Jmax_c     = 17.57,
        Jmax_Ea    = 43.54e3,
        Kc_c       = 38.05,
        Kc_Ea      = 79.43e3,
        Ko_c       = 20.30,
        Ko_Ea      = 36.38e3,
        phi_PSII_0 = 0.352,
        phi_PSII_1 = 0.022,
        phi_PSII_2 = -3.4e-4,
        RL_c       = 18.72,
        RL_Ea      = 46.39e3,
        theta_0    = 0.76,
        theta_1    = 0.018,
        theta_2    = -3.7e-4,
        Tp_c       = 19.77399,
        Tp_Ha      = 62.99e3,
        Tp_Hd      = 182.14e3,
        Tp_S       = 0.588e3,
        Vcmax_c    = 26.35,
        Vcmax_Ea   = 65.33e3,

        # ten_layer_c3_canopy module
        b0                          = 0.008,
        b1                          = 10.6,
        beta_PSII                   = 0.5,
        Catm                        = 372.59,
        electrons_per_carboxylation = 4.5,
        electrons_per_oxygenation   = 5.25,
        Gs_min                      = 1e-3,
        Jmax_at_25                  = 195,
        Jmax_at_25_mature           = 195,
        leafwidth                   = 0.1,
        O2                          = 210,
        RL_at_25                    = 1.28,
        sf_jmax                     = 0.2,
        Tp_at_25                    = 13,
        windspeed_height            = 5,

        # ten_layer_canopy_integrator module
        growth_respiration_fraction = 0,

        # partitioning_growth_calculator module
        grc_grain                       = 0.0,
        grc_leaf                        = 0.0,
        grc_rhizome                     = 0.0,
        grc_root                        = 0.00253,
        grc_shell                       = 0.0,
        grc_stem                        = 0.02256,

        # maintenance_respiration_calculator module
        mrc_grain                       = 0.0,
        mrc_leaf                        = 0.00036626,
        mrc_rhizome                     = 0.0,
        mrc_root                        = 0.00001017,
        mrc_shell                       = 0.0,
        mrc_stem                        = 0.00036626,

        # partitioning_growth module
        retrans                     = 0.9,
        retrans_rhizome             = 1.0,

        # senescence_coefficient_logistic module
        rateSeneLeaf                = 0.011955,
        rateSeneStem                = 0.000229,
        rateSeneRoot                = 0,
        rateSeneRhizome             = 0,
        alphaSeneLeaf               = 49.7562,
        alphaSeneStem               = 25.2319,
        alphaSeneRoot               = 10,
        alphaSeneRhizome            = 10,
        betaSeneLeaf                = -30.2365,
        betaSeneStem                = -16.3467,
        betaSeneRoot                = -10,
        betaSeneRhizome             = -10,

        # thermal_time_senescence_logistic module
        remobilization_fraction     = 0.6,

        # two_layer_soil_profile module
        soil_depth1                 = 0.0,
        soil_depth2                 = 2.5,
        soil_depth3                 = 10.0,
        wsFun                       = 2,
        hydrDist                    = 0,
        rfl                         = 0.2,
        rsdf                        = 0.44,
        phi1                        = 0.01,
        phi2                        = 1.5,

        # thermal_time_linear module
        tbase                       = 10,

        # litter_cover module
        km_leaf_litter              = 2.0,

        # system parameters
        timestep                    = 1
    )
)
