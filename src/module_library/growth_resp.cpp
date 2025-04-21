#include "growth_resp.h"
#include "temperature_response_functions.h"  // for Q10_temperature_response

/**
 *  @brief Calculates respiratory losses associated with a particular base rate
 *         of biomass available for growth.
 *
 *  @param [in] base_rate The base rate of carbon production that does not
 *         include respiratory losses. Any units are acceptable, such as
 *         mol / m^2 / s or Mg / ha / hour.
 *
 *  @param [in] grc Growth respiration coefficient (dimensionless)
 *
 *  @param [in] temp Temperature (degrees C)
 *
 *  @return A rate of respiratory losses having the same units as `base_rate`
 *
 *  ### Model overview
 *
 *  The idea here is that when `A_base` is the rate of carbon being allocated to
 *  a tissue for growth, a fraction `f_g` of this carbon is lost to respiration
 *  rather than being converted into new biomass. In other words, the rate of
 *  respiratory losses `R_g` is given by
 *
 *  > `R_g = f_g * A_base (1)
 *
 *  The remaining carbon available for growth (`A_growth`) is given by
 *
 *  > `A_growth = A_base * (1 - f_g) (2)
 *
 *  In this function, the temperature dependence of the respiration coefficient
 *  is modeled using a simple "Q10" method where `Q10 = 2` and the base
 *  temperature is 0 degrees C, i.e.,
 *
 *  > `f_g = f_g_0 * 2^(T / 10)` (3)
 *
 *  Ever since the days of WIMOVAC, there has been an additional constraint that
 *  `A_growth` is clamped to be greater than or equal to zero. From
 *  Equation (2), we can see that `A_growth` would be negative whenever `A_base`
 *  is negative, since the fraction `f_g` must lie on [0, 1] by definition. So,
 *  this clamping operation is equivalent to setting `f_g` = 1 when `A_base` is
 *  negative.
 *
 *  In the code below, `base_rate` represents `A_base`, `grc0` and `grc`
 *  represent `f_g_0` and `f_g`, and `temp` represents the temperature `T`.
 *
 *  ### Sources
 *
 *  In Stephen Humphries's thesis, he describes the respiration model in the
 *  following way: "The respiration model of McCree (1970) modified according to
 *  Penning de Vries (1972) and Thornley (1970) is used here to predict plant
 *  respiration... The respiration associated with each plant structure is
 *  modified according to the temperature of the structure using a Q10
 *  approximation with a value of 2 as described by Spain and Keen (1992)."
 *
 *  I (EBL) believe these references to be to the following documents:
 *
 *  - [de Vries, F. P. "Respiration and growth" in "Crop processes in controlled
 *    environments" 327–347 (Academic Press, 1972)]
 *    (https://library.wur.nl/WebQuery/wurpubs/fulltext/218533)
 *
 *  - [Thornley, J. H. M. "Respiration, Growth and Maintenance in Plants" Nature
 *    227, 304–305 (1970)](https://doi.org/10.1038/227304b0)
 *
 *  - McCree, K. "An equation for the rate of respiration of whiteclover plants
 *    grown under controlled conditions" in: "Prediction and measurement
 *    of photosynthetic productivity" Wageningen: Centre for Agricultural
 *    Publishing and Documentation, 221-229 (1970)
 *
 *  - Keen, R. E. & Spain, J. D. "Computer simulation in biology. A BASIC
 *    introduction" (1992)
 *
 *  I believe the McCree, Thornley, and de Vries papers describe the
 *  relationship between gross assimilation and respiration, while the Keen &
 *  Spain book describes the temperature dependence. Unfortunately, I can't find
 *  an online version of that book so I can't be sure.
 *
 *  For some general discussions about respiration, see the following two
 *  sources:
 *
 *  - [Amthor, J. S. "The McCree–de Wit–Penning de Vries–Thornley Respiration
 *    Paradigms: 30 Years Later" Ann Bot 86, 1–20 (2000)]
 *    (https://doi.org/10.1006/anbo.2000.1175)
 *
 *  - [Amthor, J. S. "The role of maintenance respiration in plant growth."
 *    Plant, Cell & Environment 7, 561–569 (1984)]
 *    (https://doi.org/10.1111/1365-3040.ep11591833)
 *
 *  The Humphries thesis is also unfortunately not available online:
 *
 *  Humphries, S. "Will mechanistically rich models provide us with new insights
 *  into the response of plant production to climate change?: development and
 *  experiments with WIMOVAC: (Windows Intuitive Model of Vegetation response
 *  to Atmosphere & Climate Change)" (University of Essex, 2002)
 *
 *  Originally, this function was described as applying costs due to maintenance
 *  respiration. Yufeng [YH] has pointed out that since the cost is proportional
 *  to the growth rate, this is actually growth respiration. See these papers:
 *  - Apsim: (https://apsimdev.apsim.info/ApsimX/Documents/AgPastureScience.pdf)
 *  - Thornley, J. H. M. "Growth, maintenance and respiration: a re-interpretation."
 *    Annals of Botany 41.6 (1977): 1191-1203.
 */
double growth_resp(double const base_rate, double const grc0, double const temp)
{
    // Get the growth respiration coefficient at the current temperature
    double const grc{base_rate < 0 ? 1.0
                                   : grc0 * Q10_temperature_response(temp, 0.0)};  // dimensionless

    // Check for error conditions
    if (grc < 0.0 || grc > 1.0) {
        throw std::range_error("Thrown in growth_resp: grc is outside [0, 1].");
    }

    return base_rate * grc;
}
