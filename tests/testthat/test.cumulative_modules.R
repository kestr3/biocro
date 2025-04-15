test_that('cumulative carbon and water dynamics can be added to crop models', {
    res <- expect_silent(
        with(soybean, {run_biocro(
            c(initial_values, list(
                canopy_assimilation = 0,
                canopy_gross_assimilation = 0,
                canopy_non_photorespiratory_CO2_release = 0,
                canopy_photorespiration = 0,
                canopy_transpiration = 0,
                soil_evaporation = 0,
                total_precip = 0
            )),
            parameters,
            soybean_weather[['2002']],
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
