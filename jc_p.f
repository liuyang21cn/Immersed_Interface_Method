c-----------------------------------------------------------------------
c
      subroutine jc_p
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer many,nany,ic,jc,id,jd
      parameter(many=3,nany=3)
      real*8 ds,xx,yy,xn,yn,gacobi,foo,pn,pplus,pminus
      real*8 pp(nany),xa(many),ya(many),xb(many),yb(many)
      real*8 r(1,ns),w(1,ns),uu(nany),vv(nany),so
      real*8 dudnjcn, dvdnjcn,su,sv

      ds=1.01d0*sqrt(dx*dx+dy*dy)
      DO l=1,ms
        do m=0,ns-1
           xn=taoy(m,l)
           yn=-taox(m,l)
           gacobi=dsqrt(xn*xn+yn*yn)
           xn=xn/gacobi
           yn=yn/gacobi

           su=us(m,l)
           sv=vs(m,l)
        do n=1,nany
           xx=xs(m,l)+dble(n)*ds*xn
           yy=ys(m,l)+dble(n)*ds*yn
           i=int((xx-x0)/hdx)
           j=int((yy-y0)/hdy)
           if(mod(i,2).eq.0) then
           ic=i/2
        else
           ic=(i+1)/2
           endif
        if(mod(j,2).eq.0) then
        jc=j/2
        else
        jc=(j+1)/2
        endif

       id=int(sign(1.0d0,xn))
       jd=int(sign(1.0d0,yn))
       if(id.lt.0) then
        ic=ic+1
       endif
       if(jd.lt.0) then
       jc=jc+1
       endif

       ! interpolate pp(n)
       do i=1,many
          do j=1,many
          xa(j)=yc(jc+jd*(j-1))
          ya(j)=p(ic+id*(i-1),jc+jd*(j-1))
          enddo
       xb(i)=xc(ic+id*(i-1))
       call interpolate(xa,ya,many,yy,yb(i),foo)
       enddo
      call interpolate(xb,yb,many,xx,pp(n),foo)
      enddo  !enddo n=1,nany

      pplus=3.0d0*pp(1)-3.0d0*pp(2)+pp(3)
      pminus=-xsctt(l)*xs(m,l)-ysctt(l)*ys(m,l)+0.5d0*thetat(l)*
     .        thetat(l)*((xs(m,l)-xsc(l))*(xs(m,l)-xsc(l))+
     .        (ys(m,l)-ysc(l))*(ys(m,l)-ysc(l)))+0.0d0
       pjc(5,m,l)=pplus-pminus
       pjc(6,m,l)=pplus-pminus
       enddo

      do m=0,ns-1
      if(m.eq.0) then
      xn= taoy(ns-1,l)
      yn=-taox(ns-1,l)
       else
      xn= taoy(m-1,l)
       yn=-taox(m-1,l)
        endif
      gacobi=dsqrt(xn*xn+yn*yn)
       xn=xn/gacobi
      yn=yn/gacobi

      su=us(m,l)
      sv=vs(m,l)
      do n=1,nany
       xx=xs(m,l)+dble(n)*ds*xn
         yy=ys(m,l)+dble(n)*ds*yn
         i=int((xx-x0)/hdx)
       j=int((yy-y0)/hdy)
        if(mod(i,2).eq.0) then
       ic=i/2
        else
        ic=(i+1)/2
       endif
       if(mod(j,2).eq.0) then
          jc=j/2
        else
           jc=(j+1)/2
          endif

         id=int(sign(1.0d0,xn))
           jd=int(sign(1.0d0,yn))
          if(id.lt.0) then
         ic=ic+1
          endif
           if(jd.lt.0) then
           jc=jc+1
          endif

! interpolate pp(n)
         do i=1,many
             do j=1,many
          xa(j)=yc(jc+jd*(j-1))
          ya(j)=p(ic+id*(i-1),jc+jd*(j-1))
          enddo
           xb(i)=xc(ic+id*(i-1))
            call interpolate(xa,ya,many,yy,yb(i),foo)
           enddo
          call interpolate(xb,yb,many,xx,pp(n),foo)
          enddo !enddo n=1,nany

          pplus=3.0d0*pp(1)-3.0d0*pp(2)+pp(3)
          pminus=-xsctt(l)*xs(m,l)-ysctt(l)*ys(m,l)+0.5d0*thetat(l)*
     .            thetat(l)*((xs(m,l)-xsc(l))*(xs(m,l)-xsc(l))+
     .            (ys(m,l)-ysc(l))*(ys(m,l)-ysc(l)))+0.0d0
          pjc(5,m,l)=0.5d0*(pplus-pminus+pjc(5,m,l))
            pjc(6,m,l)=pjc(5,m,l)
          enddo
           pjc(5,ns,l)=pjc(5,0,l)
           pjc(6,ns,l)=pjc(6,0,l)
          ENDDO

      return
      end


c-----------------------------------------------------------------------
