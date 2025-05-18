subroutine normalize_rho_log(d, rho_in, normalized_rho_out)
    use types, only: dp
    use dft_data, only: dft_data_t
    use constants, only: pi
    use utils, only: stop_error ! For stop_error
    implicit none

    type(dft_data_t), intent(in) :: d
    real(dp), dimension(size(d%R)), intent(in) :: rho_in
    real(dp), dimension(size(d%R)), intent(out) :: normalized_rho_out

    real(dp) :: integral_val
    real(dp) :: alpha
    integer :: i
    integer :: grid_size

    grid_size = size(d%R)
    if (grid_size < 4) then
        call stop_error("normalize_rho_log: Grid size must be at least 4 for this integration scheme.")
    end if

    integral_val = 0.0_dp
    alpha = log(d%R(2) / d%R(1)) ! Assuming d%R is not empty and d%R(1) > 0

    ! Using a specific integration scheme. Ensure grid_size is sufficient.
    integral_val = integral_val + (1.0_dp/2.0_dp) * 4.0_dp * pi * alpha * &
        (rho_in(2) * d%R(2)**3 + rho_in(1) * d%R(1)**3)
    integral_val = integral_val + (1.0_dp/12.0_dp) * 4.0_dp * pi * alpha * &
        (5.0_dp * rho_in(3) * d%R(3)**3 + &
         8.0_dp * rho_in(2) * d%R(2)**3 - rho_in(1) * d%R(1)**3)

    do i = 4, grid_size
        integral_val = integral_val + (1.0_dp/24.0_dp) * 4.0_dp * pi * alpha * &
            (9.0_dp  * rho_in(i)   * d%R(i)**3 + &
             19.0_dp * rho_in(i-1) * d%R(i-1)**3 - &
             5.0_dp  * rho_in(i-2) * d%R(i-2)**3 + &
             rho_in(i-3) * d%R(i-3)**3)
    end do

    do i = 1, grid_size
        normalized_rho_out(i) = (rho_in(i) * d%Z) / integral_val
    end do

end subroutine normalize_rho_log