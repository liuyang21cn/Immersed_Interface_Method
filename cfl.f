c-----------------------------------------------------------------------
c
      subroutine cfl
      include 'parameter.inc'
      include 'field.inc'
      real*8 umax,vmax,dtc,dtv

      if(icfl.eq.0) then
        dt=dt0
        t=t+dt
        return
      endif

      umax=0.0d0
      do j=0,ny+1
        do i=0,nx
          umax=max(umax,abs(u(i,j)))
        enddo
      enddo

      vmax=0.0d0
      do j=0,ny
        do i=0,nx+1
          vmax=max(vmax,abs(v(i,j)))
        enddo
      enddo

      if(umax.eq.0.0.and.vmax.eq.0.0) then
        dtc=100.0d0
      else
        dtc=cflc/(umax*dx1+vmax*dy1)
      endif
      dtv=cflv*Re/(dx2+dy2)

      dt=min(dtc,dtv,dtcfl)
      dt1=1.0d0/dt
      dt2=dt1/dt

      t=t+dt

      return
      end


c-----------------------------------------------------------------------
