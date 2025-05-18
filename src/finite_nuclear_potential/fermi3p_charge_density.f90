function fermi3p_rho(param_c, param_z, param_w, r_val)
    use types, only: dp
    implicit none
    real(dp), intent(in) :: param_c, param_z, param_w, r_val
    real(dp) :: fermi3p_rho

    fermi3p_rho = (1.0_dp + param_w * (r_val**2) / param_c**2) / (1.0_dp + dexp((r_val - param_c) / param_z))

end function fermi3p_rho

subroutine fermi3p_charge_density(d, rho_out)
    use types, only: dp
    use dft_data, only: dft_data_t
    use utils, only: stop_error
    implicit none

    type(dft_data_t), intent(inout) :: d ! Changed to inout
    real(dp), dimension(size(d%R)), intent(out) :: rho_out
    ! real(dp), intent(out) :: nuclear_radius_out ! Removed

    real(dp) :: param_c, param_z, param_w
    ! real(dp) :: r_val_local ! Removed
    ! real(dp) :: current_rho_val ! Removed
    integer :: i_r
    real(dp) :: cutoff_rho_val
    logical :: cutoff_reached
    real(dp) :: r0 = 0.0_dp
    real(dp), external :: fermi3p_rho

    param_c = d%rho_nuc_distro_parameters(2)
    param_z = d%rho_nuc_distro_parameters(3)
    param_w = d%rho_nuc_distro_parameters(4)

    cutoff_reached = .false.
    d%nuclear_radius = d%R(size(d%R)) ! Initialize d%nuclear_radius

    cutoff_rho_val = d%min_nuc_rho_cutoff * fermi3p_rho(param_c, param_z, param_w, r0)

    rho_out = 0.0_dp

    do i_r = 1, size(d%R)
        if (cutoff_reached) then
            rho_out(i_r) = 0.0_dp
        else
            rho_out(i_r) = fermi3p_rho(param_c, param_z, param_w, d%R(i_r))
            if (rho_out(i_r) < cutoff_rho_val) then
                cutoff_reached = .true.
                d%nuclear_radius = d%R(i_r)
            endif
        endif
    end do

end subroutine fermi3p_charge_density
