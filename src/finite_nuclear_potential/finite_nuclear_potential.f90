subroutine finite_nuclear_potential(d)

    !CURRENTLY ONLY FERMI 2PF IS FUNCTIONAL WITH DFTATOM

    !Wrapper for the entire finite nuclear potential calculation
    !Uses charge distribution data to build a nuclear charge density distribution
    !Then calculates nuclear potential from that
    !Puts it onto the correct grid in the end

    !De Vries, H, De Jager, C W, and De Vries, C. Nuclear charge-density-distribution parameters 
    !from elastic electron scattering. United States: N. p., 1987. 
    !Web. doi:10.1016/0092-640X(87)90013-1.


    use types, only: dp
    use dft_data, only: dft_data_t
    use constants, only: pi

    implicit none

    type(dft_data_t), intent(inout) :: d

    real(dp), dimension(size(d%R)) :: rho_log_grid
    real(dp), dimension(size(d%R)) :: normalized_rho_log_grid
    real(dp) :: radius_rms
    real(dp) :: a
    real(dp) :: c_param  ! Renamed from c to avoid conflict with d%c (speed of light)
    real(dp) :: alpha
    real(dp) :: z_param  ! Renamed from z to avoid conflict with d%Z (atomic number)
    real(dp) :: w

    ! Ensure pointers in d are associated if they are to be used without allocation here
    ! For fb_coefficients and sog_parameters, they are allocated based on nuclear_model_type

    SELECT CASE(d%nuclear_model_type)
        CASE (1) ! harmonic_osc_nuc_rho
            ! Parameters radius_rms, a, alpha are now extracted inside the subroutine from d
            ! radius_rms = d%rho_nuc_distro_parameters(1)
            ! a = d%rho_nuc_distro_parameters(2)
            ! alpha = d%rho_nuc_distro_parameters(3)

            call harmonic_oscillator_charge_density(d, rho_log_grid)


        CASE (2) ! mod_harmonic_osc_nuc_rho
            ! Parameters radius_rms, a, alpha are now extracted inside the subroutine from d
            ! radius_rms = d%rho_nuc_distro_parameters(1)
            ! a = d%rho_nuc_distro_parameters(2)
            ! alpha = d%rho_nuc_distro_parameters(3)

            call mod_harmonic_oscillator_charge_density(d, rho_log_grid)

        CASE (3) ! fourier_bessel_nuc_rho
            ! fb_coefficients should be allocated and filled in drivers
            call fourier_bessel_charge_density(d, rho_log_grid)
            ! The last element of fb_coefficients is the cutoff radius
            if (.not. associated(d%fb_coefficients) .or. size(d%fb_coefficients) < 1) then
                 call stop_error("finite_nuclear_potential: d%fb_coefficients not properly set for fourier_bessel before nuclear_radius assignment.")
            end if
            d%nuclear_radius = d%fb_coefficients(size(d%fb_coefficients))


        CASE (4) ! sum_of_gaussian_nuc_rho
            ! sog_parameters should be allocated and filled in drivers
            call sum_of_gaussians_nuc_charge_density(d, rho_log_grid)


        CASE (5) ! two_param_fermi_nuc_rho
            ! Parameters are taken from d%rho_nuc_distro_parameters by the subroutine
            ! c_param = d%rho_nuc_distro_parameters(1) ! Was c in original, index 2
            ! z_param = d%rho_nuc_distro_parameters(2) ! Was z in original, index 3

            call fermi2p_charge_density(d, rho_log_grid)


        CASE (6) ! three_param_fermi_nuc_rho
            ! Parameters c_param, z_param, w are now extracted inside the subroutine from d
            ! c_param = d%rho_nuc_distro_parameters(1)
            ! z_param = d%rho_nuc_distro_parameters(2)
            ! w       = d%rho_nuc_distro_parameters(3)

            call fermi3p_charge_density(d, rho_log_grid)



        CASE (7) ! three_param_gaussian_nuc_rho
            ! Parameters c_param, z_param, w are now extracted inside the subroutine from d
            ! c_param = d%rho_nuc_distro_parameters(1)
            ! z_param = d%rho_nuc_distro_parameters(2)
            ! w       = d%rho_nuc_distro_parameters(3)

            call gaussian3p_charge_density(d, rho_log_grid)

        CASE(8) ! constant_density
            if (.not. associated(d%rho_nuc_distro_parameters) .or. size(d%rho_nuc_distro_parameters) < 1) then
                 call stop_error("finite_nuclear_potential: d%rho_nuc_distro_parameters not properly set for constant_density nuclear_radius.")
            end if
            d%nuclear_radius = d%rho_nuc_distro_parameters(1) ! Set d%nuclear_radius before the call
            call const_charge_density(d, rho_log_grid)

        CASE DEFAULT
            call stop_error("finite_nuclear_potential: Unknown nuclear_model_type")
    END SELECT

    call normalize_rho_log(d, rho_log_grid, normalized_rho_log_grid)

    ! nuc_pot is now d%V_coulomb
    call nuc_pot_log_grid(d, normalized_rho_log_grid)

    
end subroutine finite_nuclear_potential