module finite_nuclear_potential_mod
    use types, only: dp
    use dft_data, only: dft_data_t
    use constants, only: pi
    use utils, only: stop_error
    ! We will need to add USE statements for the modules containing the charge density subroutines
    ! and normalize_rho_log, nuc_pot_log_grid once they are modularized.
    ! For now, this might cause errors if they are not found or if implicit interfaces are assumed.
    implicit none
private
public :: finite_nuclear_potential

contains

subroutine finite_nuclear_potential(d)

    !Wrapper for the entire finite nuclear potential calculation
    !Uses charge distribution data to build a nuclear charge density distribution
    !Then calculates nuclear potential from that
    !Puts it onto the correct grid in the end

    !De Vries, H, De Jager, C W, and De Vries, C. Nuclear charge-density-distribution parameters 
    !from elastic electron scattering. United States: N. p., 1987. 
    !Web. doi:10.1016/0092-640X(87)90013-1.

    type(dft_data_t), intent(inout) :: d

    real(dp), dimension(size(d%R)) :: rho_log_grid
    real(dp), dimension(size(d%R)) :: normalized_rho_log_grid

    SELECT CASE(d%nuclear_model_type)
        CASE (1) ! harmonic_osc_nuc_rho
            call harmonic_oscillator_charge_density(d, rho_log_grid)
        CASE (2) ! mod_harmonic_osc_nuc_rho
            call mod_harmonic_oscillator_charge_density(d, rho_log_grid)
        CASE (3) ! fourier_bessel_nuc_rho
            call fourier_bessel_charge_density(d, rho_log_grid)
        CASE (4) ! sum_of_gaussian_nuc_rho
            call sum_of_gaussians_nuc_charge_density(d, rho_log_grid)
        CASE (5) ! two_param_fermi_nuc_rho
            call fermi2p_charge_density(d, rho_log_grid)
        CASE (6) ! three_param_fermi_nuc_rho
            call fermi3p_charge_density(d, rho_log_grid)
        CASE (7) ! three_param_gaussian_nuc_rho
            call gaussian3p_charge_density(d, rho_log_grid)
        CASE(8) ! constant_density
            d%nuclear_radius = d%rho_nuc_distro_parameters(1)
            call const_charge_density(d, rho_log_grid)
        CASE DEFAULT
            call stop_error("finite_nuclear_potential: Unknown nuclear_model_type")
    END SELECT

    call normalize_rho_log(d, rho_log_grid, normalized_rho_log_grid)
    call nuc_pot_log_grid(d, normalized_rho_log_grid)
    
end subroutine finite_nuclear_potential

end module finite_nuclear_potential_mod