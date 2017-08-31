c-----------------------------------------------------------------------
c
      subroutine correction_interpolate
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer is,js
      real*8 sx,sy

      DO l=1,ms

      do j=0,ncye(l) 
        sx=falfaye(1,j,l)       
        if(ieye(j,l).eq.icye(j,l)) then
          ivxeye(j,l)=ieye(j,l)
          is=icye(j,l)+1
          vcixeye(j,l)=-0.5d0*(vjcye(1,j,l)*(xc(is)-sx)+
     .                   0.5d0*vjcye(3,j,l)*(xc(is)-sx)**2.0d0)
        else
          ivxeye(j,l)=ieye(j,l)+1
          is=icye(j,l)
          vcixeye(j,l)=0.5d0*(vjcye(1,j,l)*(xc(is)-sx)+
     .                  0.5d0*vjcye(3,j,l)*(xc(is)-sx)**2.0d0)
        endif
      enddo

      do j=0,ncyc(l)
        sx=falfayc(1,j,l)
        if(ieyc(j,l).eq.icyc(j,l)) then
          iuxcyc(j,l)=icyc(j,l)+1
          is=ieyc(j,l)
          ucixcyc(j,l)=0.5d0*(ujcyc(1,j,l)*(xe(is)-sx)+
     .                  0.5d0*ujcyc(3,j,l)*(xe(is)-sx)**2.0d0)
        else
          iuxcyc(j,l)=icyc(j,l)
          is=ieyc(j,l)+1
          ucixcyc(j,l)=-0.5d0*(ujcyc(1,j,l)*(xe(is)-sx)+
     .                   0.5d0*ujcyc(3,j,l)*(xe(is)-sx)**2.0d0)
        endif
      enddo

      do i=0,ncxc(l)
        sy=falfaxc(1,i,l)
        if(jexc(i,l).eq.jcxc(i,l)) then
          jvycxc(i,l)=jcxc(i,l)+1
          js=jexc(i,l)
          vciycxc(i,l)=0.5d0*(vjcxc(2,i,l)*(ye(js)-sy)+
     .                  0.5d0*vjcxc(4,i,l)*(ye(js)-sy)**2.0d0)
          juyexc(i,l)=jexc(i,l)
          js=jcxc(i,l)+1
          uciyexc(i,l)=-0.5d0*(ujcxc(2,i,l)*(yc(js)-sy)+
     .                   0.5d0*ujcxc(4,i,l)*(yc(js)-sy)**2.0d0)
        else
          jvycxc(i,l)=jcxc(i,l)
          js=jexc(i,l)+1
          vciycxc(i,l)=-0.5d0*(vjcxc(2,i,l)*(ye(js)-sy)+
     .                   0.5d0*vjcxc(4,i,l)*(ye(js)-sy)**2.0d0)
          juyexc(i,l)=jexc(i,l)+1
          js=jcxc(i,l)
          uciyexc(i,l)=0.5d0*(ujcxc(2,i,l)*(yc(js)-sy)+
     .                 0.5d0*ujcxc(4,i,l)*(yc(js)-sy)**2.0d0)
        endif
      enddo

      do i=0,ncxe(l)
        sy=falfaxe(1,i,l)
        if(jexe(i,l).eq.jcxe(i,l)) then
          juyexe(i,l)=jexe(i,l)
          js=jcxe(i,l)+1
          uciyexe(i,l)=-0.5d0*(ujcxe(2,i,l)*(yc(js)-sy)+
     .                   0.5d0*ujcxe(4,i,l)*(yc(js)-sy)**2.0d0)
          jvycxe(i,l)=jcxe(i,l)+1
          js=jexe(i,l)
          vciycxe(i,l)=0.5d0*(vjcxe(2,i,l)*(ye(js)-sy)+
     .                  0.5d0*vjcxe(4,i,l)*(ye(js)-sy)**2.0d0)
        else
          juyexe(i,l)=jexe(i,l)+1
          js=jcxe(i,l)
          uciyexe(i,l)=0.5d0*(ujcxe(2,i,l)*(yc(js)-sy)+
     .                  0.5d0*ujcxe(4,i,l)*(yc(js)-sy)**2.0d0)
          jvycxe(i,l)=jcxe(i,l)
          js=jexe(i,l)+1
          vciycxe(i,l)=-0.5d0*(vjcxe(2,i,l)*(ye(js)-sy)+
     .                   0.5d0*vjcxe(4,i,l)*(ye(js)-sy)**2.0d0)
        endif
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
