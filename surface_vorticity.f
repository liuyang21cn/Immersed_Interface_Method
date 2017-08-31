c-----------------------------------------------------------------------
c     
      subroutine surface_vorticity
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer many,nany,ic,jc,id,jd
      parameter(many=3,nany=3)
      real*8 ds,xx,yy,xn,yn,gacobi,foo
      real*8 ww(nany),xa(many),ya(many),xb(many),yb(many)

      ds=1.01d0*sqrt(dx*dx+dy*dy)

      DO l=1,ms
         do m=0,ns-1
            xn= taoy(m,l)
            yn=-taox(m,l)
            gacobi=dsqrt(xn*xn+yn*yn)
            xn=xn/gacobi
            yn=yn/gacobi
            
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
                     ya(j)=o(ic+id*(i-1),jc+jd*(j-1))
                  enddo
                  xb(i)=xc(ic+id*(i-1))
                  call interpolate(xa,ya,many,yy,yb(i),foo)
               enddo
               call interpolate(xb,yb,many,xx,ww(n),foo)
            enddo               !enddo n=1,nany
            
            swo(m,l)=3.0d0*ww(1)-3.0d0*ww(2)+ww(3)

         enddo
         swo(ns,l)=swo(0,l)
         swo(ns,l)=swo(0,l)

         
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
                     ya(j)=o(ic+id*(i-1),jc+jd*(j-1))
                  enddo
                  xb(i)=xc(ic+id*(i-1))
                  call interpolate(xa,ya,many,yy,yb(i),foo)
               enddo
               call interpolate(xb,yb,many,xx,ww(n),foo)
            enddo               !enddo n=1,nany
            
            swo(m,l)=0.5d0*(swo(m,l)+3.0d0*ww(1)-3.0d0*ww(2)+ww(3))
         enddo
         swo(ns,l)=swo(0,l)
         swo(ns,l)=swo(0,l)

         ENDDO

      return
      end

c-----------------------------------------------------------------------
