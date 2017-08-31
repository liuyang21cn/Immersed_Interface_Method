c-----------------------------------------------------------------------
c
      subroutine initial
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'old.inc'
      integer ni,nj
      real*8 tc

      t=0.0d0

      ! setup flapping movement here 
      if(move.eq.1) then
        t0=t
        l=1
        tc=pi/0.4d0
    
        theta(l)=0.75d0*pi+0.25d0*pi*(dsin(0.8d0*t0)*
     .                      (1.0d0-dexp(-t0/tc)))
        xsc(l)=1.25d0*(dcos(0.8d0*t0)+1.0d0)*dcos(pi/3.0d0)
        ysc(l)=1.25d0*(dcos(0.8d0*t0)+1.0d0)*dsin(pi/3.0d0)

        thetat(l)=0.25d0*pi*0.8d0*dcos(0.8d0*t0)*(1.0d0-dexp(-t0/tc))+
     .            0.25d0*pi*(dsin(0.8d0*t0))*(dexp(-t0/tc)/tc)
        xsct(l)=1.25d0*(-0.8d0*dsin(0.8d0*t0))*dcos(pi/3.0d0)
        ysct(l)=1.25d0*(-0.8d0*dsin(0.8d0*t0))*dsin(pi/3.0d0)

        thetatt(l)=-0.25d0*pi*0.64d0*dsin(0.8d0*t0)*
     .                              (1.0d0-dexp(-t0/tc))+
     .              0.5d0*pi*(0.8d0*dcos(0.8d0*t0))*(dexp(-t0/tc)/tc)-
     .              0.25d0*pi*(dsin(0.8d0*t0))*(dexp(-t0/tc)/tc/tc)
        xsctt(l)=1.25d0*(-0.64d0*dcos(0.8d0*t0))*dcos(pi/3.0d0)
        ysctt(l)=1.25d0*(-0.64d0*dcos(0.8d0*t0))*dsin(pi/3.0d0)
      endif

      call surface_parametrization
 
      do l=1,ms
        do m=0,ns
          us(m,l)=xsct(l)-thetat(l)*(ys(m,l)-ysc(l))
          vs(m,l)=ysct(l)+thetat(l)*(xs(m,l)-xsc(l))
        enddo
      enddo

      do j=0,ny+1
        do i=0,nx
            u(i,j)=1.0d0
        enddo
      enddo

      do j=0,ny
        v(0,j)=0.0d0
        do i=1,nx+1
            v(i,j)=0.0d0
        enddo
      enddo

      do j=0,ny+1
        do i=0,nx+1
          p(i,j)=0.0d0
        enddo
      enddo

      do j=0,ny+1
        do i=0,nx+1
          d(i,j)=0.0d0
          dn(i,j)=0.0d0
        enddo
      enddo

      if(lw_pbc.eq.0.and.le_pbc.eq.0.and.
     .   ls_pbc.eq.0.and.ln_pbc.eq.0) then
        ni=nx-1
        nj=ny-1
        call vrffti(ni,wsavei)
        call vrffti(nj,wsavej)
      endif

      if(lw_pbc.eq.1.and.le_pbc.eq.1.and.
     .   ls_pbc.eq.1.and.ln_pbc.eq.1) then
        ni=nx-2
        nj=ny-2
        call vsinti(ni,wsavei)
        call vsinti(nj,wsavej)
      endif

      if(lw_pbc.eq.2.and.le_pbc.eq.2.and.
     .   ls_pbc.eq.2.and.ln_pbc.eq.2) then
        ni=nx
        nj=ny
        call vcosti(ni,wsavei)
        call vcosti(nj,wsavej)
      endif

      if(lw_pbc.eq.1.and.le_pbc.eq.2.and.
     .   ls_pbc.eq.2.and.ln_pbc.eq.2) then
        ni=nx-1
        nj=ny
        call vsinqi(ni,wsavei)
        call vcosti(nj,wsavej)
      endif      

      if(lw_pbc.eq.1.and.le_pbc.eq.2.and.
     .   ls_pbc.eq.0.and.ln_pbc.eq.0) then
        ni=nx-1
        nj=ny-1
        call vsinqi(ni,wsavei)
        call vrffti(nj,wsavej)
      endif

      call velocity_initial
      call velocity_reset

      return
      end


c-----------------------------------------------------------------------
