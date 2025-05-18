function sog_rho_internal(d, r_val)
    use types, only: dp
    use constants, only: pi
    use dft_data, only: dft_data_t
    implicit none

    type(dft_data_t), intent(in) :: d
    real(dp), intent(in) :: r_val
    real(dp) :: sog_rho_internal

    integer :: i_sog_param
    real(dp) :: r_rms_local, rp_local, gamma_local, a_i, r_i_local, q_i_local
    integer :: num_sog_terms_local

    num_sog_terms_local = size(d%sog_parameters, dim=1) - 1

    r_rms_local = d%sog_parameters(num_sog_terms_local + 1, 1)
    rp_local    = d%sog_parameters(num_sog_terms_local + 1, 2)
    
    gamma_local = rp_local * ((2.0_dp/3.0_dp)**0.5_dp)

    sog_rho_internal = 0.0_dp
    do i_sog_param = 1, num_sog_terms_local
        r_i_local = d%sog_parameters(i_sog_param, 1)
        q_i_local = d%sog_parameters(i_sog_param, 2)
        a_i = (real(d%Z, dp) * q_i_local) / (2.0_dp * (pi**1.5_dp) * (gamma_local**3) * (1.0_dp + (2.0_dp * r_i_local**2) / (gamma_local**2)))
        sog_rho_internal = sog_rho_internal + a_i * (dexp(-1.0_dp * ((r_val - r_i_local) / gamma_local)**2) + dexp(-1.0_dp * ((r_val + r_i_local) / gamma_local)**2))
    end do

end function sog_rho_internal

subroutine sum_of_gaussians_nuc_charge_density(d, rho_out)
    use types, only: dp
    use dft_data, only: dft_data_t
    use utils, only: stop_error
    implicit none

    type(dft_data_t), intent(inout) :: d
    real(dp), dimension(size(d%R)), intent(out) :: rho_out

    integer :: i_r
    real(dp) :: cutoff_rho_val
    logical :: cutoff_reached
    real(dp) :: r0 = 0.0_dp
    real(dp), external :: sog_rho_internal

    cutoff_reached = .false.
    d%nuclear_radius = d%R(size(d%R))

    cutoff_rho_val = d%min_nuc_rho_cutoff * sog_rho_internal(d, r0)

    do i_r = 1, size(d%R)
        if (cutoff_reached) then
            rho_out(i_r) = 0.0_dp
        else
            rho_out(i_r) = sog_rho_internal(d, d%R(i_r))
            if (rho_out(i_r) < cutoff_rho_val) then
                cutoff_reached = .true.
                d%nuclear_radius = d%R(i_r)
            endif
        endif
    end do

end subroutine sum_of_gaussians_nuc_charge_density