#ifndef C3PHOTO_H
#define C3PHOTO_H

#include "photosynthesis_outputs.h"   // for photosynthesis_outputs
#include "c3_temperature_response.h"  // for c3_temperature_response_parameters

photosynthesis_outputs c3photoC(
    c3_temperature_response_parameters const tr_param,
    double const absorbed_ppfd,
    double const Tleaf,
    double const Tambient,
    double const RH,
    double const Vcmax_at_25,
    double const Jmax_at_25,
    double const TPU_rate_max,
    double const RL_at_25,
    double const bb0,
    double const bb1,
    double const Gs_min,
    double Ca,
    double const AP,
    double const O2,
    double const StomWS,
    double const electrons_per_carboxylation,
    double const electrons_per_oxygenation,
    double const beta_PSII,
    double const gbw);

double solc(double LeafT);
double solo(double LeafT);

double check_assim(
    double Ci,
    c3_temperature_response_parameters const tr_param,
    double const absorbed_ppfd,                // micromol / m^2 / s
    double const Tleaf,                        // degrees C
    double const Tambient,                     // degrees C
    double const RH,                           // dimensionless
    double const Vcmax_at_25,                  // micromol / m^2 / s
    double const Jmax_at_25,                   // micromol / m^2 / s
    double const TPU_rate_max,                 // micromol / m^2 / s
    double const RL_at_25,                     // micromol / m^2 / s
    double const b0,                           // mol / m^2 / s
    double const b1,                           // dimensionless
    double const Gs_min,                       // mol / m^2 / s
    double const Ca,                           // micromol / mol
    double const AP,                           // Pa (TEMPORARILY UNUSED)
    double const O2,                           // millimol / mol (atmospheric oxygen mole fraction)
    double const StomWS,                       // dimensionless
    double const electrons_per_carboxylation,  // self-explanatory units
    double const electrons_per_oxygenation,    // self-explanatory units
    double const beta_PSII,                    // dimensionless (fraction of absorbed light that reaches photosystem II)
    double const gbw                           // mol / m^2 / s
);
#endif
