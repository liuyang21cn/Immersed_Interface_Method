c-----------------------------------------------------------------------
c
      subroutine initial_reset
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      real*8 xx,yy

      DO l=1,ms

      do j=1,ny
        do i=1,nx-1
          xx=xe(i)-xsc(l)
          yy=yc(j)-ysc(l)
          if(IOxeyc(i,j).eq.l) then
            u(i,j)=xsct(l)-thetat(l)*yy
          endif
        enddo
      enddo

      do j=1,ny-1
        do i=1,nx
          xx=xc(i)-xsc(l)
          yy=ye(j)-ysc(l)
          if(IOxcye(i,j).eq.l) then
            v(i,j)=ysct(l)+thetat(l)*xx
          endif
        enddo
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------

