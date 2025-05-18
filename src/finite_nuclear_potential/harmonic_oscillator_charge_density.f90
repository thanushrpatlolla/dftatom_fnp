function ho_rho(param_a, param_alpha, r_val)
    use types, only: dp
    implicit none
    real(dp), intent(in) :: param_a, param_alpha, r_val
    real(dp) :: ho_rho

    ho_rho = (1.0_dp + param_alpha * (r_val / param_a)**2) * dexp(-1.0_dp * (r_val / param_a)**2)

end function ho_rho

subroutine harmonic_oscillator_charge_density(d, rho_out)
    use types, only: dp
    use dft_data, only: dft_data_t
    use utils, only: stop_error
    implicit none

    type(dft_data_t), intent(inout) :: d
    real(dp), dimension(size(d%R)), intent(out) :: rho_out

    real(dp) :: param_a, param_alpha, param_radius_rms 
    integer :: i_r
    real(dp) :: cutoff_rho_val
    logical :: cutoff_reached
    real(dp) :: r0 = 0.0_dp
    real(dp), external :: ho_rho

    param_radius_rms = d%rho_nuc_distro_parameters(2) ! Convention from finite_nuclear_potential
    param_a          = d%rho_nuc_distro_parameters(3)
    param_alpha      = d%rho_nuc_distro_parameters(4)

    cutoff_reached = .false.

    cutoff_rho_val = d%min_nuc_rho_cutoff * ho_rho(param_a, param_alpha, r0)
    rho_out = 0.0_dp

    do i_r = 1, size(d%R)
        if (cutoff_reached) then
            rho_out(i_r) = 0.0_dp
        else
            rho_out(i_r) = ho_rho(param_a, param_alpha, d%R(i_r))
            if (rho_out(i_r) < cutoff_rho_val) then
                cutoff_reached = .true.
                d%nuclear_radius = d%R(i_r) 
            endif
        endif
    end do

end subroutine harmonic_oscillator_charge_density