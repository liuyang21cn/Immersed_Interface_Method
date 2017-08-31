c-----------------------------------------------------------------------
c
      subroutine correction_strain
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ie,ie1,ic,je,je1,jc
      real*8 sx,sy

      DO l=1,ms

      do i=0,ncxc(l)
        sy=falfaxc(1,i,l)
        jc=jcxc(i,l)
        je=jexc(i,l)
        je1=je+1
        if(jc.eq.je) then
          pcudy(i,l)=-dy1*(ujcxc(2,i,l)*(ye(je)-sy)+
     .               0.5d0*ujcxc(4,i,l)*(ye(je)-sy)**2.0d0)
          pcvdy(i,l)=-dy1*(vjcxc(2,i,l)*(ye(je)-sy)+
     .               0.5d0*vjcxc(4,i,l)*(ye(je)-sy)**2.0d0)
        else
          pcudy(i,l)=-dy1*(ujcxc(2,i,l)*(ye(je1)-sy)+
     .               0.5d0*ujcxc(4,i,l)*(ye(je1)-sy)**2.0d0)
          pcvdy(i,l)=-dy1*(vjcxc(2,i,l)*(ye(je1)-sy)+
     .               0.5d0*vjcxc(4,i,l)*(ye(je1)-sy)**2.0d0)
        endif
      enddo

      do j=0,ncyc(l)
        sx=falfayc(1,j,l)
        ic=icyc(j,l)
        ie=ieyc(j,l)
        ie1=ie+1
        if(ic.eq.ie) then
          pcudx(j,l)=-dx1*(ujcyc(1,j,l)*(xe(ie)-sx)+
     .               0.5d0*ujcyc(3,j,l)*(xe(ie)-sx)**2.0d0)
          pcvdx(j,l)=-dx1*(vjcyc(1,j,l)*(xe(ie)-sx)+
     .               0.5d0*vjcyc(3,j,l)*(xe(ie)-sx)**2.0d0)
        else
          pcudx(j,l)=-dx1*(ujcyc(1,j,l)*(xe(ie1)-sx)+
     .               0.5d0*ujcyc(3,j,l)*(xe(ie1)-sx)**2.0d0)
          pcvdx(j,l)=-dx1*(vjcyc(1,j,l)*(xe(ie1)-sx)+
     .               0.5d0*vjcyc(3,j,l)*(xe(ie1)-sx)**2.0d0)
        endif
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
