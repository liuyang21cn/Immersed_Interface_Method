c-----------------------------------------------------------------------
c
      subroutine mesh
      include 'parameter.inc'
      include 'field.inc'
 
      dx=xl/dble(nx)
      hdx=0.5d0*dx
      dx1=1.0d0/dx
      dx2=dx1*dx1
      dy=yl/dble(ny)
      hdy=0.5d0*dy
      dy1=1.0d0/dy
      dy2=dy1*dy1

      do i=-2,2*nx+2
        x(i)=x0+hdx*dble(i)
      enddo
      do j=-2,2*ny+2
        y(j)=y0+hdy*dble(j)
      enddo

      xe(-1)=x(-2)
      do i=0,nx+1
        xe(i)=x(2*i)
        xc(i)=xe(i)-hdx
      enddo
      ye(-1)=y(-2)
      do j=0,ny+1
        ye(j)=y(2*j)
        yc(j)=ye(j)-hdy
      enddo

      return
      end


c-----------------------------------------------------------------------
