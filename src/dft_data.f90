module dft_data

! Contains the 'dft_data_t' type used in the DFT routines.
! This data type stores mesh, potential, atomic configuration, orbitals
! and other parameters of the DFT problem.

use types, only: dp
implicit none
private
public dft_data_t

type dft_data_t
    real(dp), dimension(:), pointer :: R, Rp, V_coulomb, V_h, V_xc, &
        e_xc, V_tot, rho
    real(dp), dimension(:, :), pointer :: orbitals
    real(dp) :: reigen_eps, alpha, c
    integer :: Z, scf_iter, reigen_max_iter
    ! If .true., we are solving the Dirac equation
    logical :: dirac
    ! Triples of (n, l, f), where "n", "l" are quantum numbers and "f" is the
    ! occupation number:
    integer, dimension(:), pointer :: no, lo, so
    real(dp), dimension(:), pointer :: fo
    real(dp), pointer :: ks_energies(:), Emax_init(:), Emin_init(:)
    ! Total energy (and its parts):
    real(dp) :: Ekin, Ecoul, Exc, Eenuc, Etot
    logical :: perturb ! use perturbation acceleration in the eigenproblem


    ! Parameters for finite nuclear potential model
    logical :: fnp !if true use a finite_nuclear_potential, V_coulomb is modified
    integer :: nuclear_model_type ! Type of nuclear model
    real(dp) :: nuclear_radius ! radius of the nucleus when using finite nuclear potential
    real(dp), dimension(:), pointer :: rho_nuc_distro_parameters !(size based on model)
    real(dp), dimension(:), pointer :: fb_coefficients           !(size based on model)
    real(dp), dimension(:,:), pointer :: sog_parameters          !(size based on model)
    ! Default cutoff for nuclear charge density
    real(dp) :: min_nuc_rho_cutoff = 1.0e-8_dp
end type

end module
