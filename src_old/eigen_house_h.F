      module eigen_house_h_mod
      implicit NONE

         real(8), pointer        ::  u0_z(:), v0_z(:)
         real(8), pointer        ::  u1_z(:), v1_z(:)
         real(8), pointer        ::  tsave(:, :)
         complex(8), pointer     ::  u0_z_h(:), v0_z_h(:)
         complex(8), pointer     ::  u1_z_h(:), v1_z_h(:)
         complex(8), pointer     ::  tsave_h(:, :)
         integer                 ::  offset1, offset2
         integer                 ::  offset3, offset4

      end module eigen_house_h_mod

