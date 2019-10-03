module reductions
  implicit none

contains

  function reduce_real_dp5d(input_arr) result(r)
    use params
    implicit none
    real(sp) :: input_arr(:, :, :, :, :)
    real(dp) :: r

    integer :: input_shape(5)
    integer :: idxs1
    integer :: idxs2
    integer :: idxs3
    integer :: idxs4
    integer :: idxs5

    input_shape = shape(input_arr)
    r = 0
    do idxs1 = 1, input_shape(1)
      do idxs2 = 1, input_shape(2)
        do idxs3 = 1, input_shape(3)
          do idxs4 = 1, input_shape(4)
            do idxs5 = 1, input_shape(5)
              r = r + input_arr(idxs1, idxs2, idxs3, idxs4, idxs5)
            enddo
          enddo
        enddo
      enddo
    enddo
  end function reduce_real_dp5d
end module reductions
