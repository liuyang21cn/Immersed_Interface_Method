c-----------------------------------------------------------------------
c
      subroutine singular_call

      call jc_firstsecond
      call jc_pressure0
      call mac_distribute
      call correction_interpolate
      call field_interpolate
      call correction_strain
      call uv_strain
      call dudv_surface
      call correction_velocity

      return
      end


c-----------------------------------------------------------------------
