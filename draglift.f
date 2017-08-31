c-----------------------------------------------------------------------
c     
      subroutine draglift
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer many,nany,ic,jc,id,jd
      parameter(many=3,nany=3)
      real*8 ds,xx,yy,xn,yn,gacobi,foo,pplus,pminus
      real*8 pp(nany),xa(many),ya(many),xb(many),yb(many)
      real*8 uu(nany),vv(nany)
      real*8 dudnjcn,dvdnjcn,su,sv
      real*8 dist,fxL(0:ns,ms),fyL(0:ns,ms)

      ds=1.01d0*sqrt(dx*dx+dy*dy)

      DO l=1,ms
         do m=0,ns-1
            xn= taoy(m,l)
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
            enddo               !enddo n=1,nany
            
            dudnjcn= -thetat(l)*yn
            dvdnjcn=  thetat(l)*xn
            pplus=3.0d0*pp(1)-3.0d0*pp(2)+pp(3)
            fxL(m,l)=-pplus*xn+1.0d0/re*
     .           (ujc(1,m,l)*xn+ujc(2,m,l)*yn+dudnjcn)
            fyL(m,l)=-pplus*yn+1.0d0/re*
     .           (vjc(1,m,l)*xn+vjc(2,m,l)*yn+dvdnjcn)

         enddo
         fxL(ns,l)=fxL(0,l)
         fyL(ns,l)=fyL(0,l)

         
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
            enddo               !enddo n=1,nany
            
            dudnjcn= -thetat(l)*yn
            dvdnjcn=  thetat(l)*xn
            pplus=3.0d0*pp(1)-3.0d0*pp(2)+pp(3)
            fxL(m,l)=0.5d0*(fxL(m,l)+(-pplus*xn+1.0d0/re*
     .           (ujc(1,m,l)*xn+ujc(2,m,l)*yn+dudnjcn)))
            fyL(m,l)=0.5d0*(fyL(m,l)+(-pplus*yn+1.0d0/re*
     .           (vjc(1,m,l)*xn+vjc(2,m,l)*yn+dvdnjcn)))
         enddo
         fxL(ns,l)=fxL(0,l)
         fyL(ns,l)=fyL(0,l)

         cd(l)=0.0d0
         cl(l)=0.0d0
         do m=0,ns-1
         dist=dsqrt((xs(m+1,l)-xs(m,l))*(xs(m+1,l)-xs(m,l))+
     .        (ys(m+1,l)-ys(m,l))*(ys(m+1,l)-ys(m,l)))
         cd(l)=cd(l)+(fxL(m,l)+fxL(m+1,l))*dist
         cl(l)=cl(l)+(fyL(m,l)+fyL(m+1,l))*dist
      enddo

      ENDDO

      return
      end

c-----------------------------------------------------------------------
