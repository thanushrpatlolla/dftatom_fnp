subroutine const_charge_density(d, rho_out)
    use types, only: dp
    use dft_data, only: dft_data_t
    implicit none

    type(dft_data_t), intent(in) :: d
    real(dp), dimension(size(d%R)), intent(out) :: rho_out

    integer :: i_r

    rho_out = 0.0_dp

    do i_r = 1, size(d%R)
        if (d%R(i_r) <= d%nuclear_radius) then
            rho_out(i_r) = 1.0_dp
        else
            rho_out(i_r) = 0.0_dp
        endif
    end do

end subroutine const_charge_density
