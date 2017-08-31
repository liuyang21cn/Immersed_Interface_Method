c-----------------------------------------------------------------------  
c
      subroutine surface_move(krk,fac)
      include 'parameter.inc'
      include 'surface.inc'
      integer krk
      real*8 fac,tc

      if(krk.eq.1) t0=t0+0.5d0*dt
      if(krk.eq.3) t0=t0+0.5d0*dt

      l=1
      tc=pi/0.4d0

      theta(l)=0.75*pi+0.25d0*pi*(dsin(0.8d0*t0)*
     .                      (1.0d0-dexp(-t0/tc)))
      xsc(l)=1.25d0*(dcos(0.8d0*t0)+1.0d0)*dcos(pi/3.0d0)
      ysc(l)=1.25d0*(dcos(0.8d0*t0)+1.0d0)*dsin(pi/3.0d0)

      thetat(l)=0.25d0*pi*0.8d0*dcos(0.8d0*t0)*(1.0d0-dexp(-t0/tc))+
     .          0.25d0*pi*(dsin(0.8d0*t0))*(dexp(-t0/tc)/tc)
      xsct(l)=1.25d0*(-0.8d0*dsin(0.8d0*t0))*dcos(pi/3.0d0)
      ysct(l)=1.25d0*(-0.8d0*dsin(0.8d0*t0))*dsin(pi/3.0d0)

      thetatt(l)=-0.25d0*pi*0.64d0*dsin(0.8d0*t0)*
     .                             (1.0d0-dexp(-t0/tc))+
     .          0.5d0*pi*(0.8d0*dcos(0.8d0*t0))*(dexp(-t0/tc)/tc)-
     .          0.25d0*pi*(dsin(0.8d0*t0))*(dexp(-t0/tc)/tc/tc)
      xsctt(l)=1.25d0*(-0.64d0*dcos(0.8d0*t0))*dcos(pi/3.0d0)
      ysctt(l)=1.25d0*(-0.64d0*dcos(0.8d0*t0))*dsin(pi/3.0d0)

      DO l=1,ms

      do m=0,ns-1
        xs(m,l)=xsc(l)+xs0(m,l)*dcos(theta(l))-ys0(m,l)*dsin(theta(l))
        ys(m,l)=ysc(l)+xs0(m,l)*dsin(theta(l))+ys0(m,l)*dcos(theta(l))
        us(m,l)=xsct(l)-thetat(l)*(ys(m,l)-ysc(l))
        vs(m,l)=ysct(l)+thetat(l)*(xs(m,l)-xsc(l))
      enddo
      xs(ns,l)=xs(0,l)
      ys(ns,l)=ys(0,l)
      us(ns,l)=us(0,l)
      vs(ns,l)=vs(0,l)

      ENDDO

      return
      end


c-----------------------------------------------------------------------
