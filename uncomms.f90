module comms
  use params
  implicit none
 
contains
!***********************************************************************
!   Update boundary terms
!***********************************************************************
  pure subroutine complete_halo_update_4(size4, Array)
!     
    integer, intent(in) :: size4
    complex(dp), intent(inout) :: Array(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, size4)
!
    Array(0,:,:,:) = Array(ksizex_l,:,:,:)
    Array(ksizex_l+1,:,:,:) = Array(1,:,:,:)
    Array(:,0,:,:) = Array(:,ksizey_l,:,:)
    Array(:,ksizey_l+1,:,:) = Array(:,1,:,:)
    Array(:,:,0,:) = Array(:,:,ksizet_l,:)
    Array(:,:,ksizet_l+1,:) = Array(:,:,1,:)
!
    return
!      
  end subroutine complete_halo_update_4
!***********************************************************************
  pure subroutine complete_halo_update_4_real(size4, Array)
!     
    integer, intent(in) :: size4
    real, intent(inout) :: Array(0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, size4)
!
    Array(0,:,:,:) = Array(ksizex_l,:,:,:)
    Array(ksizex_l+1,:,:,:) = Array(1,:,:,:)
    Array(:,0,:,:) = Array(:,ksizey_l,:,:)
    Array(:,ksizey_l+1,:,:) = Array(:,1,:,:)
    Array(:,:,0,:) = Array(:,:,ksizet_l,:)
    Array(:,:,ksizet_l+1,:) = Array(:,:,1,:)
!      
    return
!      
  end subroutine complete_halo_update_4_real
!***********************************************************************
  pure subroutine complete_halo_update_5(size5, Array)
!     
    integer, intent(in) :: size5
    complex(dp), intent(inout) :: Array(kthird, 0:ksizex_l+1, 0:ksizey_l+1, 0:ksizet_l+1, size5)
!
    Array(:,0,:,:,:) = Array(:,ksizex_l,:,:,:)
    Array(:,ksizex_l+1,:,:,:) = Array(:,1,:,:,:)
    Array(:,:,0,:,:) = Array(:,:,ksizey_l,:,:)
    Array(:,:,ksizey_l+1,:,:) = Array(:,:,1,:,:)
    Array(:,:,:,0,:) = Array(:,:,:,ksizet_l,:)
    Array(:,:,:,ksizet_l+1,:) = Array(:,:,:,1,:)
!      
    return
!      
  end subroutine complete_halo_update_5
!***********************************************************************
  pure subroutine complete_halo_update_6(size5, size6, Array)
!     
    integer, intent(in) :: size5, size6
    complex(dp), intent(inout) :: Array(kthird, 0:ksizex_l+1, 0:ksizey_l+1, &
         &                              0:ksizet_l+1, size5, size6)
!
    Array(:,0,:,:,:,:) = Array(:,ksizex_l,:,:,:,:)
    Array(:,ksizex_l+1,:,:,:,:) = Array(:,1,:,:,:,:)
    Array(:,:,0,:,:,:) = Array(:,:,ksizey_l,:,:,:)
    Array(:,:,ksizey_l+1,:,:,:) = Array(:,:,1,:,:,:)
    Array(:,:,:,0,:,:) = Array(:,:,:,ksizet_l,:,:)
    Array(:,:,:,ksizet_l+1,:,:) = Array(:,:,:,1,:,:)
!      
    return
!      
  end subroutine complete_halo_update_6
end module comms
