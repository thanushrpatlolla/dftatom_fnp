function gaussian3p_rho_internal(param_c, param_z, param_w, r_val)
    use types, only: dp
    implicit none
    real(dp), intent(in) :: param_c, param_z, param_w, r_val
    real(dp) :: gaussian3p_rho_internal

    gaussian3p_rho_internal = (1.0_dp + param_w * (r_val**2) / param_c**2)&
                             / (1.0_dp + dexp(((r_val**2) - (param_c**2)) / (param_z**2)))

end function gaussian3p_rho_internal

subroutine gaussian3p_charge_density(d, rho_out)
    use types, only: dp
    use dft_data, only: dft_data_t
    use utils, only: stop_error
    implicit none

    type(dft_data_t), intent(inout) :: d
    real(dp), dimension(size(d%R)), intent(out) :: rho_out

    real(dp) :: param_c, param_z, param_w
    integer :: i_r
    real(dp) :: cutoff_rho_val
    logical :: cutoff_reached
    real(dp) :: r0 = 0.0_dp
    real(dp), external :: gaussian3p_rho_internal

    if (size(d%rho_nuc_distro_parameters) < 3) then
        call stop_error("gaussian3p_charge_density: Not enough parameters for 3P Gaussian model.")
    end if
    param_c = d%rho_nuc_distro_parameters(2)
    param_z = d%rho_nuc_distro_parameters(3)
    param_w = d%rho_nuc_distro_parameters(4)

    cutoff_reached = .false.

    cutoff_rho_val = d%min_nuc_rho_cutoff * gaussian3p_rho_internal(param_c, param_z, param_w, r0)
    rho_out = 0.0_dp

    do i_r = 1, size(d%R)
        if (cutoff_reached) then
            rho_out(i_r) = 0.0_dp
        else
            rho_out(i_r) = gaussian3p_rho_internal(param_c, param_z, param_w, d%R(i_r))
            if (rho_out(i_r) < cutoff_rho_val) then 
                cutoff_reached = .true.
                d%nuclear_radius = d%R(i_r)
            endif
        endif
    end do

end subroutine gaussian3p_charge_density