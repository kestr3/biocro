test_that('cumulative carbon and water dynamics can be added to the soybean model', {
    soybean_res <- expect_silent(
        with(soybean, {run_biocro(
            c(initial_values, list(
                canopy_assimilation = 0,
                canopy_gross_assimilation = 0,
                canopy_non_photorespiratory_CO2_release = 0,
                canopy_photorespiration = 0,
                canopy_transpiration = 0,
                Grain_gr = 0,
                Grain_mr = 0,
                Leaf_gr = 0,
                Leaf_mr = 0,
                Leaf_WS_loss = 0,
                Rhizome_gr = 0,
                Rhizome_mr = 0,
                Root_gr = 0,
                Root_mr = 0,
                Shell_gr = 0,
                Shell_mr = 0,
                soil_evaporation = 0,
                Stem_gr = 0,
                Stem_mr = 0,
                total_precip = 0,
                whole_plant_growth_respiration = 0
            )),
            within(parameters, {
                growth_respiration_fraction = 0.01
            }),
            soybean_weather[['2002']],
            direct_modules,
            c(differential_modules, list(
                'BioCro:cumulative_carbon_dynamics',
                'BioCro:cumulative_water_dynamics'
            ))
        )})
    )

    # Check that all the assimilated carbon (gross assimilation) is balanced by
    # the sum costs of photorespiration, non-photorespiratory CO2 release by the
    # leaf, growth respiration, maintenance respiration, tissue growth, and
    # litter formation.
    soybean_res$total_carbon_use <- with(soybean_res, {
        (Grain - soybean$initial_values$Grain) +
        (Leaf - soybean$initial_values$Leaf) +
        (LeafLitter - soybean$initial_values$LeafLitter) +
        (Rhizome - soybean$initial_values$Rhizome) +
        (RhizomeLitter - soybean$initial_values$RhizomeLitter) +
        (Root - soybean$initial_values$Root) +
        (RootLitter - soybean$initial_values$RootLitter) +
        (Shell - soybean$initial_values$Shell) +
        (Stem - soybean$initial_values$Stem) +
        (StemLitter - soybean$initial_values$StemLitter) +
        canopy_non_photorespiratory_CO2_release +
        canopy_photorespiration +
        Grain_gr +
        Grain_mr +
        Leaf_gr +
        Leaf_mr +
        Leaf_WS_loss +
        Rhizome_gr +
        Rhizome_mr +
        Root_gr +
        Root_mr +
        Shell_gr +
        Shell_mr +
        Stem_gr +
        Stem_mr +
        whole_plant_growth_respiration
    })

    expect_equal(
        soybean_res$canopy_gross_assimilation,
        soybean_res$total_carbon_use
    )
})

test_that('cumulative carbon and water dynamics can be added to the willow model', {
    willow_res <- expect_silent(
        with(willow, {run_biocro(
            c(initial_values, list(
                canopy_assimilation = 0,
                canopy_gross_assimilation = 0,
                canopy_non_photorespiratory_CO2_release = 0,
                canopy_photorespiration = 0,
                canopy_transpiration = 0,
                Grain_gr = 0,
                Grain_mr = 0,
                Leaf_gr = 0,
                Leaf_mr = 0,
                Leaf_WS_loss = 0,
                Rhizome_gr = 0,
                Rhizome_mr = 0,
                Root_gr = 0,
                Root_mr = 0,
                Shell_gr = 0,
                Shell_mr = 0,
                soil_evaporation = 0,
                Stem_gr = 0,
                Stem_mr = 0,
                total_precip = 0,
                whole_plant_growth_respiration = 0
            )),
            c(parameters, list(
                Grain_mr_rate = 0,
                Leaf_mr_rate = 0,
                Rhizome_mr_rate = 0,
                Root_mr_rate = 0,
                Shell_mr_rate = 0,
                Stem_mr_rate = 0
            )),
            weather[['2002']],
            direct_modules,
            c(differential_modules, list(
                'BioCro:cumulative_carbon_dynamics',
                'BioCro:cumulative_water_dynamics'
            ))
        )})
    )
})

test_that('cumulative carbon and water dynamics can be added to the miscanthus model', {
    willow_res <- expect_silent(
        with(miscanthus_x_giganteus, {run_biocro(
            c(initial_values, list(
                canopy_assimilation = 0,
                canopy_gross_assimilation = 0,
                canopy_non_photorespiratory_CO2_release = 0,
                canopy_photorespiration = 0,
                canopy_transpiration = 0,
                Grain_gr = 0,
                Grain_mr = 0,
                Leaf_gr = 0,
                Leaf_mr = 0,
                Leaf_WS_loss = 0,
                Rhizome_gr = 0,
                Rhizome_mr = 0,
                Root_gr = 0,
                Root_mr = 0,
                Shell_gr = 0,
                Shell_mr = 0,
                soil_evaporation = 0,
                Stem_gr = 0,
                Stem_mr = 0,
                total_precip = 0,
                whole_plant_growth_respiration = 0
            )),
            c(parameters, list(
                Grain_mr_rate = 0,
                Leaf_mr_rate = 0,
                Rhizome_mr_rate = 0,
                Root_mr_rate = 0,
                Shell_mr_rate = 0,
                Stem_mr_rate = 0
            )),
            weather[['2002']],
            direct_modules,
            c(differential_modules, list(
                'BioCro:cumulative_carbon_dynamics',
                'BioCro:cumulative_water_dynamics'
            ))
        )})
    )
})

test_that('all partitioning calculators produce the same outputs', {
    expect_equal(
        module_info('BioCro:no_leaf_resp_neg_assim_partitioning_growth_calculator', verbose = FALSE)$outputs,
        module_info('BioCro:no_leaf_resp_partitioning_growth_calculator', verbose = FALSE)$outputs
    )

    expect_equal(
        module_info('BioCro:no_leaf_resp_partitioning_growth_calculator', verbose = FALSE)$outputs,
        module_info('BioCro:partitioning_growth_calculator', verbose = FALSE)$outputs
    )

    expect_equal(
        module_info('BioCro:no_leaf_resp_neg_assim_partitioning_growth_calculator', verbose = FALSE)$outputs,
        module_info('BioCro:partitioning_growth_calculator', verbose = FALSE)$outputs
    )
})
