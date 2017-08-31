c-----------------------------------------------------------------------
c
      subroutine surface_property
      include 'parameter.inc'
      include 'surface.inc'
      real*8 gacobi,coord(2,0:ns)

      DO l=1,ms

        do m=0,ns
          coord(1,m)=xs(m,l)
          coord(2,m)=ys(m,l)
        enddo

        do m=0,ns-1
          taox(m,l)=coord(1,m+1)-coord(1,m)
          taoy(m,l)=coord(2,m+1)-coord(2,m)
          gacobi=dsqrt(taox(m,l)*taox(m,l)+taoy(m,l)*taoy(m,l))
          taox(m,l)=taox(m,l)/gacobi
          taoy(m,l)=taoy(m,l)/gacobi
         enddo
         taox(ns,l)=taox(0,l)
         taoy(ns,l)=taoy(0,l)

        ENDDO

      return
      end


c-----------------------------------------------------------------------
