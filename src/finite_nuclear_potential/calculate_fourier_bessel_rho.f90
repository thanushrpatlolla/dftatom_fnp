subroutine fourier_bessel_charge_density(d, rho_out)
    use types, only: dp
    use dft_data, only: dft_data_t
    use constants, only: pi
    use utils, only: stop_error
    implicit none

    type(dft_data_t), intent(in) :: d
    real(dp), dimension(size(d%R)), intent(out) :: rho_out

    integer :: num_fb_terms
    real(dp) :: cutoff_radius_local
    integer :: i_r, v
    real(dp) :: j0_input_arg


    ! The last element of fb_coefficients is the cutoff radius.
    ! The preceding elements are the coefficients themselves.
    num_fb_terms = size(d%fb_coefficients) - 1
    cutoff_radius_local = d%fb_coefficients(num_fb_terms + 1)

    rho_out = 0.0_dp ! Initialize

    do i_r = 1, size(d%R)
        if (d%R(i_r) < cutoff_radius_local) then
            do v = 1, num_fb_terms
                j0_input_arg = dble(v) * pi * d%R(i_r) / cutoff_radius_local
                rho_out(i_r) = rho_out(i_r) + d%fb_coefficients(v) * sin(j0_input_arg) / j0_input_arg
            end do
        else
            rho_out(i_r) = 0.0_dp ! Outside cutoff radius
        endif
    end do


end subroutine fourier_bessel_charge_density




