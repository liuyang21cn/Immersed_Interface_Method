c-----------------------------------------------------------------------
c
      subroutine divergence_reset
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'rhsp.inc'
      integer ic,ic1,jc,jc1
      real*8 xx,yy,ang,signx,signy

      DO l=1,ms

      do j=1,ny
        do i=1,nx
          xx=xc(i)-xsc(l)
          yy=yc(j)-ysc(l)
          if(IOxcyc(i,j).eq.l) then
            d(i,j)=0.0d0
          endif
        enddo
      enddo

      do i=0,ncxc(l)
        signy=falfaxc(3,i,l)
        ic=ixc(i,l)
        jc=jcxc(i,l)
        jc1=jc+1
        if(signy.ge.0.0d0) then
          d(ic,jc1)=2.0d0*d(ic,jc1+1)-d(ic,jc1+2)
        else
          d(ic,jc)=2.0d0*d(ic,jc-1)-d(ic,jc-2)
        endif
      enddo

      do j=0,ncyc(l)
        signx=falfayc(2,j,l)
        jc=jyc(j,l)
        ic=icyc(j,l)
        ic1=ic+1
        if(signx.ge.0.0d0) then
          d(ic1,jc)=2.0d0*d(ic1+1,jc)-d(ic1+2,jc)
        else
          d(ic,jc)=2.0d0*d(ic-1,jc)-d(ic-2,jc)
        endif
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------

