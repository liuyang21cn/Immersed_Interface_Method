c-----------------------------------------------------------------------
c
      subroutine animation
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      do j=1,ny
        do i=1,nx
          o(i,j)=(vec(i,j)-vec(i-1,j))*dx1-
     .           (uce(i,j)-uce(i,j-1))*dy1
        enddo
      enddo
      if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxc(l)
            o(ixc(i,l),jexc(i,l)+1)=o(ixc(i,l),jexc(i,l)+1)+pcudy(i,l)
          enddo
          do j=0,ncyc(l)
            o(ieyc(j,l)+1,jyc(j,l))=o(ieyc(j,l)+1,jyc(j,l))+pcvdx(j,l)
          enddo
        enddo
      endif

      do j=1,ny
        write(26,400)(p(i,j),i=1,nx)
        write(36,400)(ucc(i,j),i=1,nx)
        write(46,400)(vcc(i,j),i=1,nx)
        write(56,400)(o(i,j),i=1,nx)
      enddo
      do m=0,ns
        write(76,400)(xs(m,l),ys(m,l),l=1,ms)
      enddo
400   format(1x,1000e16.6)

      return
      end


c-----------------------------------------------------------------------
