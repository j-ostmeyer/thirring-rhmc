module comms5
  use params
  implicit none

contains
  pure subroutine update_halo_5(size5, Array)
    !
    integer, intent(in) :: size5
    complex(dp), intent(inout) :: Array(kthird, 0:ksizex_l + 1, 0:ksizey_l + 1, 0:ksizet_l + 1, size5)
    !
    Array(:, 0, :, :, :) = Array(:, ksizex_l, :, :, :)
    Array(:, ksizex_l + 1, :, :, :) = Array(:, 1, :, :, :)
    Array(:, :, 0, :, :) = Array(:, :, ksizey_l, :, :)
    Array(:, :, ksizey_l + 1, :, :) = Array(:, :, 1, :, :)
    Array(:, :, :, 0, :) = Array(:, :, :, ksizet_l, :)
    Array(:, :, :, ksizet_l + 1, :) = Array(:, :, :, 1, :)
    !
    return
    !
  end subroutine update_halo_5
  module comms5
