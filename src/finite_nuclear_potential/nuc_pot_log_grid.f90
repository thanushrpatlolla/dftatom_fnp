subroutine nuc_pot_log_grid(d, normalized_rho_in)

    use types, only: dp
    use dft_data, only: dft_data_t
    use constants, only: pi
    use utils, only: stop_error ! For stop_error

    implicit none

    type(dft_data_t), intent(inout) :: d ! Reads d%R, d%nuclear_radius, d%Z; Writes d%V_coulomb
    real(dp), dimension(size(d%R)), intent(in) :: normalized_rho_in

    integer :: grid_size
    real(dp), dimension(size(d%R)) :: int1, int2
    integer :: i
    real(dp) :: alpha

    grid_size = size(d%R)

    if (grid_size < 4) then
        call stop_error("nuc_pot_log_grid: Grid size must be at least 4 for this integration scheme.")
    end if
    if (associated(d%V_coulomb) .and. size(d%V_coulomb) /= grid_size) then
        call stop_error("nuc_pot_log_grid: d%V_coulomb is associated but has incorrect size.")
    else if (.not. associated(d%V_coulomb)) then
        call stop_error("nuc_pot_log_grid: d%V_coulomb is not associated.")
    end if

    alpha = log(d%R(2) / d%R(1))

    int1 = 0.0_dp
    int2 = 0.0_dp

    ! Integration for int1 (integral of rho * r^2)
    int1(1) = 0.0_dp ! By definition, integral up to the first point
    int1(2) = int1(1) + (1.0_dp/2.0_dp) * 4.0_dp * pi * alpha * &
        (normalized_rho_in(2) * d%R(2)**3 + normalized_rho_in(1) * d%R(1)**3)
    int1(3) = int1(2) + (1.0_dp/12.0_dp) * 4.0_dp * pi * alpha * &
        (5.0_dp * normalized_rho_in(3) * d%R(3)**3 + &
         8.0_dp * normalized_rho_in(2) * d%R(2)**3 - normalized_rho_in(1) * d%R(1)**3)

    ! Integration for int2 (integral of rho * r)
    int2(1) = 0.0_dp ! By definition
    int2(2) = int2(1) + (1.0_dp/2.0_dp) * 4.0_dp * pi * alpha * &
        (normalized_rho_in(2) * d%R(2)**2 + normalized_rho_in(1) * d%R(1)**2)
    int2(3) = int2(2) + (1.0_dp/12.0_dp) * 4.0_dp * pi * alpha * &
        (5.0_dp * normalized_rho_in(3) * d%R(3)**2 + &
         8.0_dp * normalized_rho_in(2) * d%R(2)**2 - normalized_rho_in(1) * d%R(1)**2)

    do i = 4, grid_size
        int1(i) = int1(i-1) + (1.0_dp/24.0_dp) * 4.0_dp * pi * alpha * &
            (9.0_dp  * normalized_rho_in(i)   * d%R(i)**3 + &
             19.0_dp * normalized_rho_in(i-1) * d%R(i-1)**3 - &
             5.0_dp  * normalized_rho_in(i-2) * d%R(i-2)**3 + &
             normalized_rho_in(i-3) * d%R(i-3)**3)

        int2(i) = int2(i-1) + (1.0_dp/24.0_dp) * 4.0_dp * pi * alpha * &
            (9.0_dp  * normalized_rho_in(i)   * d%R(i)**2 + &
             19.0_dp * normalized_rho_in(i-1) * d%R(i-1)**2 - &
             5.0_dp  * normalized_rho_in(i-2) * d%R(i-2)**2 + &
             normalized_rho_in(i-3) * d%R(i-3)**2)
    end do

    do i = 1, grid_size
        if (d%R(i) <= d%nuclear_radius) then
            ! V(r) = - (1/r) * integral_0_r (rho(r') 4*pi*r'^2) dr' - integral_r_R (rho(r') 4*pi*r') dr'
            ! Here int1(i) is integral_0_r (rho(r') r'^2 dr' * 4pi), so first term is -int1(i)/r
            ! And int2(grid_size) is integral_0_R_max (rho(r') r' dr' * 4pi)
            ! And int2(i) is integral_0_r (rho(r') r' dr' * 4pi)
            ! So second term is -(int2(grid_size) - int2(i))
            d%V_coulomb(i) = -int1(i) / d%R(i) - (int2(grid_size) - int2(i))
        else
            ! Outside the nuclear radius, potential is -Z/r
            d%V_coulomb(i) = -d%Z / d%R(i)
        endif
    end do

end subroutine nuc_pot_log_grid