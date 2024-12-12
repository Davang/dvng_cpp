/*!
 * \file    units.hpp
 * \version 1.0.0
 * \author  Davang
 * \date    2024/12/10
 * \copyright GNU GPL V3
 * \brief   header only SI units library
 */


#ifndef DVNG_UNITS_H
#define DVNG_UNITS_H

#include <cstdint>
#include <ratio>

namespace dvng::units
{

	/*!
	 *  \brief unit_t struct that holds a long double and 7 templates one for each quantity of the SI
	 *  \tparam T_T exponent of time.
	 *  \tparam L_T exponent of length.
	 *  \tparam M_T exponent of mass.
	 *  \tparam I_T exponent of electric current.
	 *  \tparam THETA_T exponent of thermodynamic temperature.
	 *  \tparam N_T exponent of amount of substance.
	 *  \tparam J_T exponent of luminous intensity.
	 */
	template< int32_t T_T, int32_t L_T, int32_t M_T, int32_t I_T, int32_t THETA_T, int32_t N_T, int32_t J_T >
	struct unit_t
	{
		static constexpr int32_t T = T_T;
		static constexpr int32_t L = L_T;
		static constexpr int32_t M = M_T;
		static constexpr int32_t I = I_T;
		static constexpr int32_t THTETA = THETA_T;
		static constexpr int32_t N = N_T;
		static constexpr int32_t J = J_T;

		long double m_value; //!< value holder

		inline constexpr unit_t() : m_value{ 0.0 } {}

		inline constexpr unit_t(const unit_t& t_unit) : m_value{ t_unit.m_value } {}

		inline constexpr unit_t(const unit_t&& t_unit) : m_value{ t_unit.m_value } {}

		inline constexpr explicit unit_t(const long double& t_value) : m_value{ t_value } {}

		inline constexpr unit_t operator=(const unit_t& Lhs)
		{
			if (this != &Lhs)
			{
				this->m_value = Lhs.m_value;
			}

			return *this;
		}

		inline constexpr  unit_t operator=(const long double& Lhs)
		{
			this->m_value = Lhs;
			return *this;
		}

		long double& operator()()
		{
			return m_value;
		}

		const long double& operator()() const
		{
			return m_value;
		}

		inline bool constexpr operator==(const unit_t& Lhs) const
		{
			return (this->m_value == Lhs.m_value);
		}

		inline bool constexpr operator!=(const unit_t& Lhs) const
		{
			return not (*this == Lhs);
		}

		inline bool constexpr operator< (const unit_t& Lhs) const
		{
			return (this->m_value < Lhs.m_value);
		}

		inline bool constexpr operator> (const unit_t& Lhs) const
		{
			return (Lhs < *this);
		}

		inline bool constexpr operator<= (const unit_t& Lhs) const
		{
			return not (*this > Lhs);
		}

		inline bool constexpr operator>= (const unit_t& Lhs) const
		{
			return not (*this < Lhs);
		}

		inline unit_t constexpr operator+(const unit_t& Lhs) const
		{
			return unit_t(this->m_value + Lhs.m_value);
		}

		inline unit_t constexpr  operator-(const unit_t& Lhs) const
		{
			return unit_t(this->m_value - Lhs.m_value);
		}

		template<  int32_t T_LHS, int32_t L_LHS, int32_t M_LHS, int32_t I_LHS, int32_t THETA_LHS, int32_t N_LHS, int32_t J_LHS >
		inline unit_t<  T_T + T_LHS, L_T + L_LHS, M_T + M_LHS, I_T + I_LHS, THETA_T + THETA_LHS, N_T + N_LHS, J_T + J_LHS>
			constexpr operator* (const unit_t< T_LHS, L_LHS, M_LHS, I_LHS, THETA_LHS, N_LHS, J_LHS>& Lhs) const
		{
			return
				unit_t< T_T + T_LHS, L_T + L_LHS, M_T + M_LHS, I_T + I_LHS, THETA_T + THETA_LHS, N_T + N_LHS, J_T + J_LHS>
				(this->m_value * Lhs.m_value);
		}


		template<  int32_t T_LHS, int32_t L_LHS, int32_t M_LHS, int32_t I_LHS, int32_t THETA_LHS, int32_t N_LHS, int32_t J_LHS >
		inline unit_t< T_T - T_LHS, L_T - L_LHS, M_T - M_LHS, I_T - I_LHS, THETA_T - THETA_LHS, N_T - N_LHS, J_T - J_LHS>
			constexpr  operator/(const unit_t<T_LHS, L_LHS, M_LHS, I_LHS, THETA_LHS, N_LHS, J_LHS>& Lhs) const
		{
			return
				unit_t< T_T - T_LHS, L_T - L_LHS, M_T - M_LHS, I_T - I_LHS, THETA_T - THETA_LHS, N_T - N_LHS, J_T - J_LHS>
				(this->m_value / Lhs.m_value);
		}

		inline unit_t constexpr operator+(const long double& Lhs) const
		{
			return unit_t(this->m_value + Lhs);
		}

		inline unit_t constexpr  operator-(const long double& Lhs) const
		{
			return unit_t(this->m_value - Lhs);
		}

		inline unit_t constexpr operator*(const long double& Lhs) const
		{
			return unit_t(this->m_value * Lhs);
		}

		inline unit_t constexpr  operator/(const long double& Lhs) const
		{
			return unit_t(this->m_value / Lhs);
		}

		inline unit_t constexpr operator-()
		{
			m_value = -m_value;
			return *this;
		}


		inline unit_t constexpr operator+=(const unit_t& Lhs)
		{
			this->m_value = this->m_value + Lhs.m_value;
			return *this;
		}

		inline unit_t constexpr  operator-=(const unit_t& Lhs)
		{
			this->m_value = this->m_value - Lhs.m_value;
			return *this;
		}

		inline unit_t constexpr operator+=(const long double& Lhs)
		{
			this->m_value = this->m_value + Lhs;
			return *this;
		}

		inline unit_t constexpr  operator-=(const long double& Lhs)
		{
			this->m_value = this->m_value - Lhs;
			return this;
		}

		inline unit_t constexpr operator*=(const long double& Lhs)
		{
			this->m_value = this->m_value * Lhs;
			return *this;
		}

		inline unit_t constexpr  operator/=(const long double& Lhs)
		{
			this->m_value = this->m_value / Lhs;
			return this;
		}
	};

	template< int32_t T_T, int32_t L_T, int32_t M_T, int32_t I_T, int32_t THETA_T, int32_t N_T, int32_t J_T >
	inline unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>
		constexpr operator+(const long double& Rhs, const unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>& Lhs)
	{
		return unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>(Rhs + Lhs.m_value);
	}

	template< int32_t T_T, int32_t L_T, int32_t M_T, int32_t I_T, int32_t THETA_T, int32_t N_T, int32_t J_T >
	inline unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T >
		constexpr operator-(const long double& Rhs, const unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>& Lhs)
	{
		return unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>(Rhs - Lhs.m_value);
	}
	template< int32_t T_T, int32_t L_T, int32_t M_T, int32_t I_T, int32_t THETA_T, int32_t N_T, int32_t J_T >
	inline unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>
		constexpr operator*(const long double& Rhs, const unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>& Lhs)
	{
		return unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>(Rhs * Lhs.m_value);
	}
	template< int32_t T_T, int32_t L_T, int32_t M_T, int32_t I_T, int32_t THETA_T, int32_t N_T, int32_t J_T >
	inline unit_t< -T_T, -L_T, -M_T, -I_T, -THETA_T, -N_T, -J_T>
		constexpr operator/(const long double& Rhs, const unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>& Lhs)
	{
		return unit_t< -T_T, -L_T, -M_T, -I_T, -THETA_T, -N_T, -J_T>(Rhs / Lhs.m_value);
	}

	template< int32_t T_T, int32_t L_T, int32_t M_T, int32_t I_T, int32_t THETA_T, int32_t N_T, int32_t J_T >
	inline bool constexpr operator==(const long double& Rhs, const unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>& Lhs)
	{
		return (Lhs == Rhs);
	}

	template< int32_t T_T, int32_t L_T, int32_t M_T, int32_t I_T, int32_t THETA_T, int32_t N_T, int32_t J_T >
	inline bool constexpr operator!=(const long double& Rhs, const unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>& Lhs)
	{
		return (Lhs != Rhs);
	}

	template< int32_t T_T, int32_t L_T, int32_t M_T, int32_t I_T, int32_t THETA_T, int32_t N_T, int32_t J_T >
	inline bool constexpr operator< (const long double& Rhs, const unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>& Lhs)
	{
		return (Lhs > Rhs);
	}

	template< int32_t T_T, int32_t L_T, int32_t M_T, int32_t I_T, int32_t THETA_T, int32_t N_T, int32_t J_T >
	inline bool constexpr operator> (const long double& Rhs, const unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>& Lhs)
	{
		return (Lhs < Rhs);
	}

	template< int32_t T_T, int32_t L_T, int32_t M_T, int32_t I_T, int32_t THETA_T, int32_t N_T, int32_t J_T >
	inline bool constexpr operator<= (const long double& Rhs, const unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>& Lhs)
	{
		return (Lhs >= Rhs);
	}

	template< int32_t T_T, int32_t L_T, int32_t M_T, int32_t I_T, int32_t THETA_T, int32_t N_T, int32_t J_T >
	inline bool constexpr operator>= (const long double& Rhs, const unit_t< T_T, L_T, M_T, I_T, THETA_T, N_T, J_T>& Lhs)
	{
		return (Lhs <= Rhs);
	}

	//helpers

	//base units

	using time_t = unit_t<1, 0, 0, 0, 0, 0, 0>;

	using length_t = unit_t<0, 1, 0, 0, 0, 0, 0>;

	using mass_t = unit_t<0, 0, 1, 0, 0, 0, 0>;

	using electric_current_t = unit_t<0, 0, 0, 1, 0, 0, 0>;

	using thermodinamic_temperature_t = unit_t<0, 0, 0, 0, 1, 0, 0>;

	using amount_of_substance_t = unit_t<0, 0, 0, 0, 0, 1, 0>;

	using luminous_intensity_t = unit_t<0, 0, 0, 0, 0, 0, 1>;

	//derived units
	using frequency_t = unit_t<-1, 0, 0, 0, 0, 0, 0>;

	using velocity_t = unit_t<-1, 1, 0, 0, 0, 0, 0>;

	using acceleration_t = unit_t<-2, 1, 0, 0, 0, 0, 0>;

	using area_t = unit_t<0, 2, 0, 0, 0, 0, 0>;

	using volume_t = unit_t<0, 3, 0, 0, 0, 0, 0>;

	using force_t = unit_t<-2, 1, 1, 0, 0, 0, 0>;

	using pressure_t = unit_t<-2, -1, 1, 0, 0, 0, 0>;

	using energy_t = unit_t<-2, 2, 1, 0, 0, 0, 0>;

	using power_t = unit_t<-3, 2, 1, 0, 0, 0, 0>;

	using electric_charge_t = unit_t<1, 0, 0, 1, 0, 0, 0>;

	using electric_potential_t = unit_t< -3, 2, 1, -1, 0, 0, 0>;

	using capacitance_t = unit_t< 4, -2, -1, 2, 0, 0, 0>;

	using resistance_t = unit_t< -3, 2, 1, -2, 0, 0, 0>;

	using electrical_conductance_t = unit_t< 3, -2, -1, 2, 0, 0, 0>;

	using magentic_flux_t = unit_t< -2, 2, 1, -1, 0, 0, 0>;

	using magentic_flux_density_t = unit_t< -2, 0, 1, -1, 0, 0, 0>;

	using inductance_t = unit_t< -2, 2, 1, -2, 0, 0, 0>;

	using density_t = unit_t< 0, -3, 1, 0, 0, 0, 0>;

	using specific_volume_t = unit_t< 0, 3, -1, 0, 0, 0, 0>;

	using current_density_t = unit_t< 0, -2, 0, 1, 0, 0, 0>;

	using magnetic_field_strength_t = unit_t< 0, -1, 0, 1, 0, 0, 0>;

	using concentration_t = unit_t< 0, -3, 0, 0, 0, 1, 0>;

	using angle_t = unit_t<0, 0, 0, 0, 0, 0, 0>;

	using solid_angle_t = unit_t<0, 0, 0, 0, 0, 0, 0>;

	using moment_force_t = unit_t<-2, 2, 1, 0, 0, 0, 0>;

	using angular_velocity_t = unit_t<-1, 0, 0, 0, 0, 0, 0>;

	using angular_acceleration_t = unit_t<-2, 0, 0, 0, 0, 0, 0>;

} // dvng::strs


/*!
 *  \brief dvng::strs string literals for units.
 */
namespace dvng::strs
{
	//si units 
	constexpr dvng::units::time_t operator""_s(const long double t_time) { return dvng::units::time_t(t_time); }
	constexpr dvng::units::length_t operator""_m(const long double t_length) { return dvng::units::length_t(t_length); }
	constexpr dvng::units::mass_t operator""_kg(const long double t_mass) { return dvng::units::mass_t(t_mass); }
	constexpr dvng::units::electric_current_t operator""_A(const long double t_electric_current) { return dvng::units::electric_current_t(t_electric_current); }
	constexpr dvng::units::thermodinamic_temperature_t operator""_K(const long double t_thermodinamic_temperature) { return dvng::units::thermodinamic_temperature_t(t_thermodinamic_temperature); }
	constexpr dvng::units::amount_of_substance_t operator""_mol(const long double t_amount_of_substance) { return dvng::units::amount_of_substance_t(t_amount_of_substance); }
	constexpr dvng::units::luminous_intensity_t operator""_cd(const long double t_luminous_intensity) { return dvng::units::luminous_intensity_t(t_luminous_intensity); }

	//derived si units
	constexpr dvng::units::frequency_t operator""_Hz(long double t_frequency) { return dvng::units::frequency_t(t_frequency); }
	constexpr dvng::units::velocity_t operator""_m_s(long double t_velocity) { return dvng::units::velocity_t(t_velocity); }
	constexpr dvng::units::acceleration_t operator""_m_s2(long double t_acceleration) { return dvng::units::acceleration_t(t_acceleration); }
	constexpr dvng::units::area_t operator""_m2(long double t_area) { return dvng::units::area_t(t_area); }
	constexpr dvng::units::volume_t operator""_m3(long double t_volume) { return dvng::units::volume_t(t_volume); }
	constexpr dvng::units::force_t operator""_N(long double t_force) { return dvng::units::force_t(t_force); }
	constexpr dvng::units::pressure_t operator""_Pa(long double t_pressure) { return dvng::units::pressure_t(t_pressure); }
	constexpr dvng::units::energy_t operator""_J(long double t_energy) { return dvng::units::energy_t(t_energy); }
	constexpr dvng::units::power_t operator""_W(long double t_power) { return dvng::units::power_t(t_power); }
	constexpr dvng::units::electric_charge_t operator""_Q(long double t_electric_charge) { return dvng::units::electric_charge_t(t_electric_charge); }
	constexpr dvng::units::electric_potential_t operator""_V(long double t_electric_potential) { return dvng::units::electric_potential_t(t_electric_potential); }
	constexpr dvng::units::capacitance_t operator""_C(long double t_capacitance) { return dvng::units::capacitance_t(t_capacitance); }
	constexpr dvng::units::resistance_t operator""_Ohm(long double t_resistance) { return dvng::units::resistance_t(t_resistance); }
	constexpr dvng::units::electrical_conductance_t operator""_S(long double t_electrical_conductance) { return dvng::units::electrical_conductance_t(t_electrical_conductance); }
	constexpr dvng::units::magentic_flux_t operator""_Wb(long double t_magentic_flux) { return dvng::units::magentic_flux_t(t_magentic_flux); }
	constexpr dvng::units::magentic_flux_density_t operator""_T(long double t_magentic_flux_density) { return dvng::units::magentic_flux_density_t(t_magentic_flux_density); }
	constexpr dvng::units::inductance_t operator""_H(long double t_inductance) { return dvng::units::inductance_t(t_inductance); }
	constexpr dvng::units::density_t operator""_kg_m3(long double t_density) { return dvng::units::density_t(t_density); }
	constexpr dvng::units::specific_volume_t operator""_m3_kg(long double t_specific_volume) { return dvng::units::specific_volume_t(t_specific_volume); }
	constexpr dvng::units::current_density_t operator""_A_m2(long double t_current_density) { return dvng::units::current_density_t(t_current_density); }
	constexpr dvng::units::magnetic_field_strength_t operator""_A_m(long double t_magnetic_field_strength) { return dvng::units::magnetic_field_strength_t(t_magnetic_field_strength); }
	constexpr dvng::units::concentration_t operator""_mol_m3(long double t_concentration) { return dvng::units::concentration_t(t_concentration); }
	constexpr dvng::units::angle_t operator""_rad(const long double t_angle) { return dvng::units::angle_t(t_angle); }
	constexpr dvng::units::solid_angle_t operator""_sr(const long double t_solid_angle) { return dvng::units::solid_angle_t(t_solid_angle); }
	constexpr dvng::units::moment_force_t operator""_Nm(const long double t_moment_force) { return dvng::units::moment_force_t(t_moment_force); }
	constexpr dvng::units::angular_velocity_t operator""_omega(const long double t_angular_velocity) { return dvng::units::angular_velocity_t(t_angular_velocity); }
	constexpr dvng::units::angular_acceleration_t operator""_alpha(const long double t_angular_acceleration) { return dvng::units::angular_acceleration_t(t_angular_acceleration); }

	//non si units literals
	constexpr dvng::units::time_t operator""_min(const long double t_time) { return dvng::units::time_t( 60.0 * t_time); }
	constexpr dvng::units::time_t operator""_h(const long double t_time) { return dvng::units::time_t(3'600.0 * t_time); }
	constexpr dvng::units::time_t operator""_day(const long double t_time) { return dvng::units::time_t(86'400.0 * t_time); }
	constexpr dvng::units::time_t operator""_week(const long double t_time) { return dvng::units::time_t(604'800.0 * t_time); }

	constexpr dvng::units::length_t operator""_in(const long double t_length) { return dvng::units::length_t(0.0254 * t_length); }
	constexpr dvng::units::length_t operator""_ft(const long double t_length) { return dvng::units::length_t(0.3048 * t_length); }
	constexpr dvng::units::length_t operator""_yd(const long double t_length) { return dvng::units::length_t(0.9144 * t_length); }
	constexpr dvng::units::length_t operator""_mi(const long double t_length) { return dvng::units::length_t(1'609.344 * t_length); }
	constexpr dvng::units::length_t operator""_nmi(const long double t_length) { return dvng::units::length_t(1852 * t_length); }

	//prefixes 
	constexpr long double operator_atto(const long double t_value) { return (t_value * std::atto::num); };
	constexpr long double operator_femto(const long double t_value) { return (t_value * std::femto::num); };
	constexpr long double operator_pico(const long double t_value) { return (t_value * std::pico::num); };
	constexpr long double operator_nano(const long double t_value) { return (t_value * std::nano::num); };
	constexpr long double operator_micro(const long double t_value) { return (t_value * std::micro::num); };
	constexpr long double operator_milli(const long double t_value) { return (t_value * std::milli::num); };
	constexpr long double operator_centi(const long double t_value) { return (t_value * std::centi::num); };
	constexpr long double operator_deci(const long double t_value) { return (t_value * std::deci::num); };
	constexpr long double operator_deca(const long double t_value) { return (t_value * std::deca::num); };
	constexpr long double operator_hecto(const long double t_value) { return (t_value * std::hecto::num); };
	constexpr long double operator_kilo(const long double t_value) { return (t_value * std::kilo::num); };
	constexpr long double operator_mega(const long double t_value) { return (t_value * std::mega::num); };
	constexpr long double operator_giga(const long double t_value) { return (t_value * std::giga::num); };
	constexpr long double operator_tera(const long double t_value) { return (t_value * std::tera::num); };
	constexpr long double operator_peta(const long double t_value) { return (t_value * std::peta::num); };
	constexpr long double operator_exa(const long double t_value) { return (t_value * std::exa::num); };
}// dvng::strs

#endif /* DVNG_UNITS_H */