#ifndef SOLAR_ZENITH_ANGLE_H
#define SOLAR_ZENITH_ANGLE_H

#include "../framework/state_map.h"
#include "../framework/module.h"
//#include "..framework/constants.h" 

namespace standardBML
{
/**
 * @class solar_zenith_angle
 * 
 * @brief Calculates the solar zenith angle using a simple model. A major
 * shortcoming of this model is that solar noon always occurs at 12:00 PM.
 */
class solar_zenith_angle : public direct_module
{
   public:
    solar_zenith_angle(
        state_map const& input_quantities,
        state_map* output_quantities)
        : direct_module{},

          // Get references to input parameters
          lat(get_input(input_quantities, "lat")),
          fractional_doy(get_input(input_quantities, "fractional_doy")),

          // Get pointers to output parameters
          cosine_zenith_angle_op(get_op(output_quantities, "cosine_zenith_angle"))
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "solar_zenith_angle"; }


   private:
    // References to input parameters
    double const& lat;
    double const& fractional_doy;
    // Pointers to output parameters
    double* cosine_zenith_angle_op;
    // Main operation
    void do_operation() const;
};

string_vector solar_zenith_angle::get_inputs()
{
    return {
        "lat",      // degrees (North is positive)
        "fractional_doy",  // time expressed as a fractional day of year
    };
}

string_vector solar_zenith_angle::get_outputs()
{
    return {
        "cosine_zenith_angle"  // dimensionless
    };
}

    // copied from AuxBioCro.cpp from Yufeng's
double cos_zenith_angle(const double latitude, const int day_of_year,
                        const double hour_of_day)
{  
    constexpr double pi = 3.1415926;
    constexpr double radians_per_degree = pi/180;
    constexpr int solar_noon = 12;
    constexpr double radians_rotation_per_hour = 15 * radians_per_degree;
    constexpr double axial_tilt = 23.5 * radians_per_degree;

    const double phi = latitude * radians_per_degree;
    const int NDS = day_of_year + 10;

    const double omega = 360.0 * (NDS / 365.0) * radians_per_degree;

    const double delta = -axial_tilt * cos(omega);

    const double tau = (hour_of_day - solar_noon) * radians_rotation_per_hour;

    return sin(delta) * sin(phi) + cos(delta) * cos(phi) * cos(tau);
}

void solar_zenith_angle::do_operation() const
{
    // Unpack the doy and hour
    const double doy = floor(fractional_doy);
    const double hour = 24.0 * (fractional_doy - doy);

    // Update the output pointers
    update(cosine_zenith_angle_op, cos_zenith_angle(lat, doy, hour));
}

} // namespace standardBML
#endif
