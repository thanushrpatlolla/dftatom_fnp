function fermi2p_rho(fermi_2p_c, fermi_2p_z, r_val)
    use types, only: dp ! Use the project-wide precision
    implicit none
    real(dp), intent(in) :: fermi_2p_c
    real(dp), intent(in) :: fermi_2p_z
    real(dp), intent(in) :: r_val
    real(dp) :: fermi2p_rho

    fermi2p_rho = 1.0_dp / (1.0_dp + dexp((r_val - fermi_2p_c) / fermi_2p_z))

end function fermi2p_rho


subroutine fermi2p_charge_density(d, rho_out)
    use types, only: dp
    use dft_data, only: dft_data_t
    implicit none

    type(dft_data_t), intent(in) :: d
    real(dp), dimension(size(d%R)), intent(out) :: rho_out

    ! Parameters for Fermi 2P model are expected in d%rho_nuc_distro_parameters
    real(dp) :: fermi_2p_c, fermi_2p_z

    integer :: i_r
    logical :: cutoff_reached
    real(dp) :: r0 = 0.0_dp
    real(dp), external :: fermi2p_rho ! Function is defined above

    ! Extract model parameters from d
    ! Assuming they are stored as: c = params(1), z = params(2)
    if (size(d%rho_nuc_distro_parameters) < 2) then
        call stop_error("fermi2p_charge_density: Not enough parameters in d%rho_nuc_distro_parameters")
    end if
    fermi_2p_c = d%rho_nuc_distro_parameters(2)
    fermi_2p_z = d%rho_nuc_distro_parameters(3)

    cutoff_reached = .false.

    ! Calculate the cutoff threshold using min_nuc_rho_cutoff from d
    cutoff_rho_val = d%min_nuc_rho_cutoff * fermi2p_rho(fermi_2p_c, fermi_2p_z, r0)

    do i_r = 1, size(d%R)
        if (cutoff_reached) then
            rho_out(i_r) = 0.0_dp
        else
            rho_out(i_r) = fermi2p_rho(fermi_2p_c, fermi_2p_z, d%R(i_r))
            if (current_rho_val < cutoff_rho_val) then
                cutoff_reached = .true.
                d%nuclear_radius = d%R(i_r)
            endif
        endif
    end do

end subroutine fermi2p_charge_density