c-----------------------------------------------------------------------
c
      subroutine correction_velocity
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ie,ie1,ic,ic1,je,je1,jc,jc1
      real*8 sx,sy,hp,hm,signx,signy
      real*8 ui,vi,djc,ddjc


      DO l=1,ms

      do i=0,ncxe(l)
        signy=falfaxe(3,i,l)
        sy=falfaxe(1,i,l)
        ie=ixe(i,l)
        je=jexe(i,l)
        je1=je+1
        jc=jcxe(i,l)
        jc1=jc+1
        if(jc.eq.je) then
          hp=yc(jc1)-sy
          hm=sy-ye(je)
          ui=falfaxe(13,i,l)
          vi=falfaxe(14,i,l)

          djc=ui*vjcxe(2,i,l)+vi*ujcxe(2,i,l)
          ddjc=ui*vjcxe(4,i,l)+vi*ujcxe(4,i,l)+2.0d0*
     .         (falfaxe(4,i,l)*falfaxe(6,i,l)-
     .          falfaxe(5,i,l)*falfaxe(7,i,l))

          ucdy(i,l)=-dy1*(djc*(ye(je)-sy)+
     .              0.5d0*ddjc*(ye(je)-sy)**2.0d0)
          ucdyy(1,i,l)=-dy2*(ujcxe(2,i,l)*(yc(jc1)-sy)+
     .                 0.5d0*ujcxe(4,i,l)*(yc(jc1)-sy)**2.0d0)
          ucdyy(2,i,l)=dy2*(ujcxe(2,i,l)*(yc(jc)-sy)+
     .                0.5d0*ujcxe(4,i,l)*(yc(jc)-sy)**2.0d0)
        else
          hp=ye(je1)-sy
          hm=sy-yc(jc)
          ui=falfaxe(13,i,l)
          vi=falfaxe(14,i,l)

          djc=ui*vjcxe(2,i,l)+vi*ujcxe(2,i,l)
          ddjc=ui*vjcxe(4,i,l)+vi*ujcxe(4,i,l)+2.0d0*
     .         (falfaxe(4,i,l)*falfaxe(6,i,l)-
     .          falfaxe(5,i,l)*falfaxe(7,i,l))

          ucdy(i,l)=-dy1*(djc*(ye(je1)-sy)+
     .             0.5d0*ddjc*(ye(je1)-sy)**2.0d0)
          ucdyy(1,i,l)=-dy2*(ujcxe(2,i,l)*(yc(jc1)-sy)+
     .                 0.5d0*ujcxe(4,i,l)*(yc(jc1)-sy)**2.0d0)
          ucdyy(2,i,l)=dy2*(ujcxe(2,i,l)*(yc(jc)-sy)+
     .                0.5d0*ujcxe(4,i,l)*(yc(jc)-sy)**2.0d0)
        endif
      enddo

      do i=0,ncxc(l)
        signy=falfaxc(3,i,l)
        sy=falfaxc(1,i,l)
        ic=ixc(i,l)
        je=jexc(i,l)
        je1=je+1
        jc=jcxc(i,l)
        jc1=jc+1
        if(jc.eq.je) then
          hp=yc(jc1)-sy
          hm=sy-ye(je)
          ui=falfaxc(13,i,l)
          vi=falfaxc(14,i,l)

          djc=2.0d0*vi*vjcxc(2,i,l)
          ddjc=2.0d0*vi*vjcxc(4,i,l)+2.0d0*
     .         (falfaxc(10,i,l)*falfaxc(10,i,l)-
     .          falfaxc(11,i,l)*falfaxc(11,i,l))

          vcvdy(i,l)=-dy1*(djc*(yc(jc1)-sy)+
     .              0.5d0*ddjc*(yc(jc1)-sy)**2.0d0)
          vcdyy(1,i,l)=-dy2*(vjcxc(2,i,l)*(ye(je1)-sy)+
     .                 0.5d0*vjcxc(4,i,l)*(ye(je1)-sy)**2.0d0)
          vcdyy(2,i,l)=dy2*(vjcxc(2,i,l)*(ye(je)-sy)+
     .                0.5d0*vjcxc(4,i,l)*(ye(je)-sy)**2.0d0)
        else
          hp=ye(je1)-sy
          hm=sy-yc(jc)
          ui=falfaxc(13,i,l)
          vi=falfaxc(14,i,l)

          djc=2.0d0*vi*vjcxc(2,i,l)
          ddjc=2.0d0*vi*vjcxc(4,i,l)+2.0d0*
     .         (falfaxc(10,i,l)*falfaxc(10,i,l)-
     .          falfaxc(11,i,l)*falfaxc(11,i,l))

          vcvdy(i,l)=-dy1*(djc*(yc(jc)-sy)+
     .              0.5d0*ddjc*(yc(jc)-sy)**2.0d0)
          vcdyy(1,i,l)=-dy2*(vjcxc(2,i,l)*(ye(je1)-sy)+
     .                 0.5d0*vjcxc(4,i,l)*(ye(je1)-sy)**2.0d0)
          vcdyy(2,i,l)=dy2*(vjcxc(2,i,l)*(ye(je)-sy)+
     .                0.5d0*vjcxc(4,i,l)*(ye(je)-sy)**2.0d0)
        endif
      enddo

      do j=0,ncye(l)
        signx=falfaye(2,j,l)
        sx=falfaye(1,j,l)
        je=jye(j,l)
        ie=ieye(j,l)
        ie1=ie+1
        ic=icye(j,l)
        ic1=ic+1
        if(ic.eq.ie) then
          hp=xc(ic1)-sx
          hm=sx-xe(ie)
          ui=falfaye(13,j,l)
          vi=falfaye(14,j,l)

          djc=ui*vjcye(1,j,l)+vi*ujcye(1,j,l)
          ddjc=ui*vjcye(3,j,l)+vi*ujcye(3,j,l)+2.0d0*
     .         (falfaye(4,j,l)*falfaye(6,j,l)-
     .          falfaye(5,j,l)*falfaye(7,j,l))

          vcdx(j,l)=-dx1*(djc*(xe(ie)-sx)+
     .             0.5d0*ddjc*(xe(ie)-sx)**2.0d0)
          vcdxx(1,j,l)=-dx2*(vjcye(1,j,l)*(xc(ic1)-sx)+
     .                 0.5d0*vjcye(3,j,l)*(xc(ic1)-sx)**2.0d0)
          vcdxx(2,j,l)=dx2*(vjcye(1,j,l)*(xc(ic)-sx)+
     .                0.5d0*vjcye(3,j,l)*(xc(ic)-sx)**2.0d0)
        else
          hp=xe(ie1)-sx
          hm=sx-xc(ic)
          ui=falfaye(13,j,l)
          vi=falfaye(14,j,l)

          djc=ui*vjcye(1,j,l)+vi*ujcye(1,j,l)
          ddjc=ui*vjcye(3,j,l)+vi*ujcye(3,j,l)+2.0d0*
     .         (falfaye(4,j,l)*falfaye(6,j,l)-
     .          falfaye(5,j,l)*falfaye(7,j,l))

          vcdx(j,l)=-dx1*(djc*(xe(ie1)-sx)+
     .             0.5d0*ddjc*(xe(ie1)-sx)**2.0d0)
          vcdxx(1,j,l)=-dx2*(vjcye(1,j,l)*(xc(ic1)-sx)+
     .                 0.5d0*vjcye(3,j,l)*(xc(ic1)-sx)**2.0d0)
          vcdxx(2,j,l)=dx2*(vjcye(1,j,l)*(xc(ic)-sx)+
     .                0.5d0*vjcye(3,j,l)*(xc(ic)-sx)**2.0d0)
        endif
      enddo

      do j=0,ncyc(l)
        signx=falfayc(2,j,l)
        sx=falfayc(1,j,l)
        jc=jyc(j,l)
        ie=ieyc(j,l)
        ie1=ie+1
        ic=icyc(j,l)
        ic1=ic+1
        if(ic.eq.ie) then
          hp=xc(ic1)-sx
          hm=sx-xe(ie)
          ui=falfayc(13,j,l)
          vi=falfayc(14,j,l)

          djc=2.0d0*ui*ujcyc(1,j,l)
          ddjc=2.0d0*ui*ujcyc(3,j,l)+2.0d0*
     .         (falfayc(4,j,l)*falfayc(4,j,l)-
     .          falfayc(5,j,l)*falfayc(5,j,l))

          ucudx(j,l)=-dx1*(djc*(xc(ic1)-sx)+
     .              0.5d0*ddjc*(xc(ic1)-sx)**2.0d0)
          ucdxx(1,j,l)=-dx2*(ujcyc(1,j,l)*(xe(ie1)-sx)+
     .                 0.5d0*ujcyc(3,j,l)*(xe(ie1)-sx)**2.0d0)
          ucdxx(2,j,l)=dx2*(ujcyc(1,j,l)*(xe(ie)-sx)+
     .                0.5d0*ujcyc(3,j,l)*(xe(ie)-sx)**2.0d0)
        else
          hp=xe(ie1)-sx
          hm=sx-xc(ic)
          ui=falfayc(13,j,l)
          vi=falfayc(14,j,l)

          djc=2.0d0*ui*ujcyc(1,j,l)
          ddjc=2.0d0*ui*ujcyc(3,j,l)+2.0d0*
     .         (falfayc(4,j,l)*falfayc(4,j,l)-
     .          falfayc(5,j,l)*falfayc(5,j,l))

          ucudx(j,l)=-dx1*(djc*(xc(ic)-sx)+
     .              0.5d0*ddjc*(xc(ic)-sx)**2.0d0)
          ucdxx(1,j,l)=-dx2*(ujcyc(1,j,l)*(xe(ie1)-sx)+
     .                 0.5d0*ujcyc(3,j,l)*(xe(ie1)-sx)**2.0d0)
          ucdxx(2,j,l)=dx2*(ujcyc(1,j,l)*(xe(ie)-sx)+
     .                0.5d0*ujcyc(3,j,l)*(xe(ie)-sx)**2.0d0)
        endif
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
