c-----------------------------------------------------------------------
c
      subroutine correction_pressure
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ie,ic,ic1,je,jc,jc1
      real*8 sx,sy

      DO l=1,ms

      do i=0,ncxc(l)
        sy=falfaxc(1,i,l)
        je=jexc(i,l)
        jc=jcxc(i,l)
        jc1=jc+1
        pcdyy(1,i,l)=-dy2*(pjcxc(6,i,l)+
     .                     pjcxc(2,i,l)*(yc(jc1)-sy)+
     .               0.5d0*pjcxc(4,i,l)*(yc(jc1)-sy)**2.0d0)
        pcdyy(2,i,l)=dy2*(pjcxc(6,i,l)+
     .                    pjcxc(2,i,l)*(yc(jc)-sy)+
     .              0.5d0*pjcxc(4,i,l)*(yc(jc)-sy)**2.0d0)
        if(jc.eq.je) then
          vcpdy(i,l)=-dy1*(pjcxc(6,i,l)+
     .                     pjcxc(2,i,l)*(yc(jc1)-sy)+
     .               0.5d0*pjcxc(4,i,l)*(yc(jc1)-sy)**2.0d0)
        else 
          vcpdy(i,l)=-dy1*(pjcxc(6,i,l)+
     .                     pjcxc(2,i,l)*(yc(jc)-sy)+
     .               0.5d0*pjcxc(4,i,l)*(yc(jc)-sy)**2.0d0)
        endif
      enddo

      do j=0,ncyc(l)
        sx=falfayc(1,j,l)
        ie=ieyc(j,l)
        ic=icyc(j,l)
        ic1=ic+1
        pcdxx(1,j,l)=-dx2*(pjcyc(5,j,l)+
     .                     pjcyc(1,j,l)*(xc(ic1)-sx)+
     .               0.5d0*pjcyc(3,j,l)*(xc(ic1)-sx)**2.0d0)
        pcdxx(2,j,l)=dx2*(pjcyc(5,j,l)+
     .                    pjcyc(1,j,l)*(xc(ic)-sx)+
     .              0.5d0*pjcyc(3,j,l)*(xc(ic)-sx)**2.0d0)
        if(ic.eq.ie) then
          ucpdx(j,l)=-dx1*(pjcyc(5,j,l)+
     .                     pjcyc(1,j,l)*(xc(ic1)-sx)+
     .               0.5d0*pjcyc(3,j,l)*(xc(ic1)-sx)**2.0d0)
        else
          ucpdx(j,l)=-dx1*(pjcyc(5,j,l)+
     .                     pjcyc(1,j,l)*(xc(ic)-sx)+
     .               0.5d0*pjcyc(3,j,l)*(xc(ic)-sx)**2.0d0)
        endif
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
