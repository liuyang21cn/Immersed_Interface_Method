c-----------------------------------------------------------------------
      subroutine jc_firstsecond
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      real*8 gacobi2,gacobi,r1,r2,r3,tx,ty,xn,yn,ft,fn
      real*8 duxjc,duyjc,dduxjc,dduyjc,dduxyjc,ddpxjc0,ddpxjc1
      real*8 dvxjc,dvyjc,ddvxjc,ddvyjc,ddvxyjc,ddpyjc0,ddpyjc1
      real*8 dpxjc,dpyjc     
      real*8 dudnjc,dudnjcn,dudnjcp,dvdnjc,dvdnjcn,dvdnjcp
      real*8 ddudnp,ddvdnp,qx,qy 
      real*8 deltaujc1(0:ns,ms),deltavjc1(0:ns,ms)
      real*8 deltaujc2(0:ns,ms),deltavjc2(0:ns,ms)
      real*8 ujc2(6,0:ns,ms),vjc2(6,0:ns,ms),pjc2(8,0:ns,ms)
      ! dudxdt: derivative of dudx over tao 
      real*8 dudxdt,dvdxdt,dpdxdt,dudydt,dvdydt,dpdydt

      ! variables used in interpolation of dudnjc ddudnjc
      real*8 ds,su,sv,xx,yy,uu(3),vv(3),dist,foo
      integer many,nany,ie,ic,je,jc,id,jd,iu,ju,iv,jv,io,jo
      parameter(many=3,nany=3)
      real*8 xa(many),ya(many),xb(many),yb(many)

      ! used in interpolation of
      ds=1.01d0*sqrt(dx*dx+dy*dy)
      DO l=1,ms
         do m=0,ns-1
            ! to test cylinder, set curvture 1/R
            tx=taox(m,l)
            ty=taoy(m,l)

            gacobi=dsqrt(tx*tx+ty*ty)
            gacobi2=tx*tx+ty*ty

            ! normal direction
            !rotate direction correctly!            
            xn=ty
            yn=-tx
            xn=xn/gacobi
            yn=yn/gacobi

            ! velocity of lagrangian point
            su=us(m,l)
            sv=vs(m,l)
            ! 1.Compute jump condition of dudx, dudy, dpdx, dpdy
            ! 1.1 jump condtion of dudn
            ! 1.1.1 dudn jump conditions on negative side
            ! dudnjcn= thetat X n
            
            dudnjcn = -thetat(l)*yn
            dvdnjcn =  thetat(l)*xn

            ! 1.1.2 dudn jump conditions on positive side
            ! ============ put this part into a subroutine ============
            ! a. find 3 points on normal direction s1, s2, s3

            do n=1,nany

               ! coordinates of interpolation point
               xx=xs(m,l)+dble(n)*ds*xn
               yy=ys(m,l)+dble(n)*ds*yn
               ! indices of interpolation point
               i=int((xx-x0)/hdx)
               j=int((yy-y0)/hdy)
               if(mod(i,2).eq.0) then
                  ie=i/2
                  ic=ie
               else
                  ic=(i+1)/2
                  ie=ic-1
               endif
               if(mod(j,2).eq.0) then
                  je=j/2
                  jc=je
               else
                  jc=(j+1)/2
                  je=jc-1
               endif
               iu=ie
               ju=jc
               iv=ic
               jv=je

               id=int(sign(1.d0,xn))
               jd=int(sign(1.d0,yn))
               if(id.lt.0.0d0) then
                  iu=iu+1
                  iv=iv+1
               endif
               if(jd.lt.0.0d0) then
                  ju=ju+1
                  jv=jv+1
               endif

               ! interpolation of uu(n)
               do j=1,many
                  do i=1,many
                     xa(i)=xe(iu+id*(i-1))
                     ya(i)=u(iu+id*(i-1),ju+jd*(j-1))
                  enddo
                  xb(j)=yc(ju+jd*(j-1))
                  call interpolate(xa,ya,many,xx,yb(j),foo)
               enddo
               call interpolate(xb,yb,many,yy,uu(n),foo)
               
               ! interpolation of vv(n)
               do i=1,many
                  do j=1,many
                     xa(j)=ye(jv+jd*(j-1))
                     ya(j)=v(iv+id*(i-1),jv+jd*(j-1))
                  enddo
                  xb(i)=xc(iv+id*(i-1))
                  call interpolate(xa,ya,many,yy,yb(i),foo)
               enddo
               call interpolate(xb,yb,many,xx,vv(n),foo)               
            ! end do n=1,nany
            enddo
            
            ! b. dudn jump conditions on positive side
            dudnjcp=(4.0d0*uu(1)-uu(2)-3.0d0*su)/ds/2.0d0
            dvdnjcp=(4.0d0*vv(1)-vv(2)-3.0d0*sv)/ds/2.0d0

            ! c. dudn jump conditions  
            dudnjc = dudnjcp - dudnjcn
            dvdnjc = dvdnjcp - dvdnjcn

            ! 1.2 delta u jump conditions
            ! 1.2.1 curvature
            ddudnp=(-5.0d0*uu(1)+4.0d0*uu(2)-1.0d0*uu(3)+2.0d0*su)/ds/ds
            ddvdnp=(-5.0d0*vv(1)+4.0d0*vv(2)-1.0d0*vv(3)+2.0d0*sv)/ds/ds
                                
            ! 1.2.3 ddudn jump conditions
            deltaujc1(m,l) = ddudnp + curv(m,l)*dudnjc
            deltavjc1(m,l) = ddvdnp + curv(m,l)*dvdnjc

            ! 1.3 dudx dvdx dudy dvdy
            duxjc = -ty/(tx*yn-ty*xn)*dudnjc
            dvxjc = -ty/(tx*yn-ty*xn)*dvdnjc
            duyjc =  tx/(tx*yn-ty*xn)*dudnjc
            dvyjc =  tx/(tx*yn-ty*xn)*dvdnjc 
            
            ! 1.4 dpdx dpdy jump conditions
            ! 1.4.1 body force qx, qy
            qx =  thetatt(l)*(ys(m,l)-ysc(l))
            qy = -thetatt(l)*(xs(m,l)-xsc(l))
            
            ! 1.4.2 dpdx, dpdy jump conditions
            dpxjc = 1.0d0/re*deltaujc1(m,l) + qx
            dpyjc = 1.0d0/re*deltavjc1(m,l) + qy

            ujc(1,m,l)=duxjc
            vjc(1,m,l)=dvxjc
            pjc(1,m,l)=dpxjc

            ujc(2,m,l)=duyjc
            vjc(2,m,l)=dvyjc
            pjc(2,m,l)=dpyjc

         enddo
         deltaujc1(ns,l)=deltaujc1(0,l)
         deltavjc1(ns,l)=deltavjc1(0,l)
         ujc(1,ns,l)=ujc(1,0,l)
         vjc(1,ns,l)=vjc(1,0,l)
         pjc(1,ns,l)=pjc(1,0,l)
         ujc(2,ns,l)=ujc(2,0,l)
         vjc(2,ns,l)=vjc(2,0,l)
         pjc(2,ns,l)=pjc(2,0,l)

         do m=0,ns-1
            tx=taox(m,l)
            ty=taoy(m,l)
            gacobi2=tx*tx+ty*ty
            ! 2.  second order jump conditions
            ! 2.1 ddudx ddudy ddvdx ddvdy jump conditions 
            ! 2.1.1 dudxdt, dvdxdt, dudydt, dxdydt
            ! dudxdt = (dudxB-dudxA)/|AB|
            !===== indices of ujc, vjc, pjc =====

            dist=dsqrt((ys(m+1,l)-ys(m,l))*(ys(m+1,l)-ys(m,l))+
     .                 (xs(m+1,l)-xs(m,l))*(xs(m+1,l)-xs(m,l)))
            dudxdt = (ujc(1,m+1,l) - ujc(1,m,l))/dist
            dvdxdt = (vjc(1,m+1,l) - vjc(1,m,l))/dist
            dpdxdt = (pjc(1,m+1,l) - pjc(1,m,l))/dist

            dudydt = (ujc(2,m+1,l) - ujc(2,m,l))/dist
            dvdydt = (vjc(2,m+1,l) - vjc(2,m,l))/dist
            dpdydt = (pjc(2,m+1,l) - pjc(2,m,l))/dist

            ! 2.1.2 dduxjc dduyjc dduxyjc  
            dduxjc=(tx*dudxdt-ty*dudydt+ty*ty*deltaujc1(m,l))/gacobi2
            ddvxjc=(tx*dvdxdt-ty*dvdydt+ty*ty*deltavjc1(m,l))/gacobi2
            
            dduyjc=(ty*dudydt-tx*dudxdt+tx*tx*deltaujc1(m,l))/gacobi2
            ddvyjc=(ty*dvdydt-tx*dvdxdt+tx*tx*deltavjc1(m,l))/gacobi2
            
            dduxyjc=(ty*dudxdt+tx*dudydt-tx*ty*deltaujc1(m,l))/gacobi2
            ddvxyjc=(ty*dvdxdt+tx*dvdydt-tx*ty*deltavjc1(m,l))/gacobi2
            
            ! 2.2 ddpdx ddpdy jump conditions
            r3 = 0.0d0
            ddpxjc0  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
            ddpyjc0  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

            r3 = 1.0d0
            ddpxjc1  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
            ddpyjc1  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

            ujc(3,m,l)=dduxjc
            vjc(3,m,l)=ddvxjc
            pjc(3,m,l)=ddpxjc1
            ujc(5,m,l)=dduxyjc
            vjc(5,m,l)=ddvxyjc
            pjc(7,m,l)=ddpxjc0

            ujc(4,m,l)=dduyjc
            vjc(4,m,l)=ddvyjc
            pjc(4,m,l)=ddpyjc1
            ujc(6,m,l)=dduxyjc
            vjc(6,m,l)=ddvxyjc
            pjc(8,m,l)=ddpyjc0
         enddo
         ujc(3,ns,l)=ujc(3,0,l)
         vjc(3,ns,l)=vjc(3,0,l)
         pjc(3,ns,l)=pjc(3,0,l)
         ujc(5,ns,l)=ujc(5,0,l)
         vjc(5,ns,l)=vjc(5,0,l)
         pjc(7,ns,l)=pjc(7,0,l)

         ujc(4,ns,l)=ujc(4,0,l)
         vjc(4,ns,l)=vjc(4,0,l)
         pjc(4,ns,l)=pjc(4,0,l)
         ujc(6,ns,l)=ujc(6,0,l)
         vjc(6,ns,l)=vjc(6,0,l)
         pjc(8,ns,l)=pjc(8,0,l)
      ENDDO

      ! reverse direction
      DO l=1,ms
         do m=ns,1,-1
            ! to test cylinder, set curvture 1/R
            tx=-taox(m-1,l)
            ty=-taoy(m-1,l)

            gacobi=dsqrt(tx*tx+ty*ty)
            gacobi2=tx*tx+ty*ty

            ! normal direction
            !rotate direction correctly!            
            xn=-ty
            yn=tx
            xn=xn/gacobi
            yn=yn/gacobi

            ! velocity of lagrangian point
            su=us(m,l)
            sv=vs(m,l)
       
            ! 1.Compute jump condition of dudx, dudy, dpdx, dpdy
            ! 1.1 jump condtion of dudn
            ! 1.1.1 dudn jump conditions on negative side
            ! dudnjcn= thetat X n
            
            dudnjcn = -thetat(l)*yn
            dvdnjcn =  thetat(l)*xn

            ! 1.1.2 dudn jump conditions on positive side
            ! ============ put this part into a subroutine ============
            ! a. find 3 points on normal direction s1, s2, s3

            do n=1,nany

               ! coordinates of interpolation point
               xx=xs(m,l)+dble(n)*ds*xn
               yy=ys(m,l)+dble(n)*ds*yn
               ! indices of interpolation point
               i=int((xx-x0)/hdx)
               j=int((yy-y0)/hdy)
               if(mod(i,2).eq.0) then
                  ie=i/2
                  ic=ie
               else
                  ic=(i+1)/2
                  ie=ic-1
               endif
               if(mod(j,2).eq.0) then
                  je=j/2
                  jc=je
               else
                  jc=(j+1)/2
                  je=jc-1
               endif
               iu=ie
               ju=jc
               iv=ic
               jv=je

               id=int(sign(1.d0,xn))
               jd=int(sign(1.d0,yn))
               if(id.lt.0.0d0) then
                  iu=iu+1
                  iv=iv+1
               endif
               if(jd.lt.0.0d0) then
                  ju=ju+1
                  jv=jv+1
               endif

               ! interpolation of uu(n)
               do j=1,many
                  do i=1,many
                     xa(i)=xe(iu+id*(i-1))
                     ya(i)=u(iu+id*(i-1),ju+jd*(j-1))
                  enddo
                  xb(j)=yc(ju+jd*(j-1))
                  call interpolate(xa,ya,many,xx,yb(j),foo)
               enddo
               call interpolate(xb,yb,many,yy,uu(n),foo)
               
               ! interpolation of vv(n)
               do i=1,many
                  do j=1,many
                     xa(j)=ye(jv+jd*(j-1))
                     ya(j)=v(iv+id*(i-1),jv+jd*(j-1))
                  enddo
                  xb(i)=xc(iv+id*(i-1))
                  call interpolate(xa,ya,many,yy,yb(i),foo)
               enddo
               call interpolate(xb,yb,many,xx,vv(n),foo)               
            ! end do n=1,nany
            enddo
            
            ! b. dudn jump conditions on positive side
            dudnjcp=(4.0d0*uu(1)-uu(2)-3.0d0*su)/ds/2.0d0
            dvdnjcp=(4.0d0*vv(1)-vv(2)-3.0d0*sv)/ds/2.0d0

            ! c. dudn jump conditions  
            dudnjc = dudnjcp - dudnjcn
            dvdnjc = dvdnjcp - dvdnjcn

            ! 1.2 delta u jump conditions
            ! 1.2.1 curvature
            ddudnp=(-5.0d0*uu(1)+4.0d0*uu(2)-1.0d0*uu(3)+2.0d0*su)/ds/ds
            ddvdnp=(-5.0d0*vv(1)+4.0d0*vv(2)-1.0d0*vv(3)+2.0d0*sv)/ds/ds
                                
            ! 1.2.3 ddudn jump conditions
            deltaujc2(m,l) = ddudnp + curv(m,l)*dudnjc
            deltavjc2(m,l) = ddvdnp + curv(m,l)*dvdnjc

            ! 1.3 dudx dvdx dudy dvdy
            duxjc = -ty/(tx*yn-ty*xn)*dudnjc
            dvxjc = -ty/(tx*yn-ty*xn)*dvdnjc
            duyjc =  tx/(tx*yn-ty*xn)*dudnjc
            dvyjc =  tx/(tx*yn-ty*xn)*dvdnjc 
            
            ! 1.4 dpdx dpdy jump conditions
            ! 1.4.1 body force qx, qy
            qx =  thetatt(l)*(ys(m,l)-ysc(l))
            qy = -thetatt(l)*(xs(m,l)-xsc(l))
            
            ! 1.4.2 dpdx, dpdy jump conditions
            dpxjc = 1.0d0/re*deltaujc2(m,l) + qx
            dpyjc = 1.0d0/re*deltavjc2(m,l) + qy

            ujc2(1,m,l)=duxjc
            vjc2(1,m,l)=dvxjc
            pjc2(1,m,l)=dpxjc

            ujc2(2,m,l)=duyjc
            vjc2(2,m,l)=dvyjc
            pjc2(2,m,l)=dpyjc

         enddo
         deltaujc2(0,l)=deltaujc2(ns,l)
         deltavjc2(0,l)=deltavjc2(ns,l)
         ujc2(1,0,l)=ujc2(1,ns,l)
         vjc2(1,0,l)=vjc2(1,ns,l)
         pjc2(1,0,l)=pjc2(1,ns,l)
         ujc2(2,0,l)=ujc2(2,ns,l)
         vjc2(2,0,l)=vjc2(2,ns,l)
         pjc2(2,0,l)=pjc2(2,ns,l)

         do m=ns,1,-1
            tx=-taox(m-1,l)
            ty=-taoy(m-1,l)
            gacobi2=tx*tx+ty*ty
            ! 2.  second order jump conditions
            ! 2.1 ddudx ddudy ddvdx ddvdy jump conditions 
            ! 2.1.1 dudxdt, dvdxdt, dudydt, dxdydt
            ! dudxdt = (dudxB-dudxA)/|AB|
            !===== indices of ujc, vjc, pjc =====
            dist=dsqrt((ys(m-1,l)-ys(m,l))*(ys(m-1,l)-ys(m,l))+
     .                 (xs(m-1,l)-xs(m,l))*(xs(m-1,l)-xs(m,l)))
            dudxdt = (ujc2(1,m-1,l) - ujc2(1,m,l))/dist
            dvdxdt = (vjc2(1,m-1,l) - vjc2(1,m,l))/dist
            dpdxdt = (pjc2(1,m-1,l) - pjc2(1,m,l))/dist

            dudydt = (ujc2(2,m-1,l) - ujc2(2,m,l))/dist
            dvdydt = (vjc2(2,m-1,l) - vjc2(2,m,l))/dist
            dpdydt = (pjc2(2,m-1,l) - pjc2(2,m,l))/dist
            ! 2.1.2 dduxjc dduyjc dduxyjc  
            dduxjc=(tx*dudxdt-ty*dudydt+ty*ty*deltaujc2(m,l))/gacobi2
            ddvxjc=(tx*dvdxdt-ty*dvdydt+ty*ty*deltavjc2(m,l))/gacobi2
            
            dduyjc=(ty*dudydt-tx*dudxdt+tx*tx*deltaujc2(m,l))/gacobi2
            ddvyjc=(ty*dvdydt-tx*dvdxdt+tx*tx*deltavjc2(m,l))/gacobi2
            
            dduxyjc=(ty*dudxdt+tx*dudydt-tx*ty*deltaujc2(m,l))/gacobi2
            ddvxyjc=(ty*dvdxdt+tx*dvdydt-tx*ty*deltavjc2(m,l))/gacobi2
            
            ! 2.2 ddpdx ddpdy jump conditions
            r3 = 0.0d0
            ddpxjc0  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
            ddpyjc0  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

            r3 = 1.0d0
            ddpxjc1  = (tx*dpdxdt - ty*dpdydt + ty*ty*r3)/gacobi2
            ddpyjc1  = (ty*dpdydt - tx*dpdxdt + tx*tx*r3)/gacobi2

            ujc2(3,m,l)=dduxjc
            vjc2(3,m,l)=ddvxjc
            pjc2(3,m,l)=ddpxjc1
            ujc2(5,m,l)=dduxyjc
            vjc2(5,m,l)=ddvxyjc
            pjc2(7,m,l)=ddpxjc0

            ujc2(4,m,l)=dduyjc
            vjc2(4,m,l)=ddvyjc
            pjc2(4,m,l)=ddpyjc1
            ujc2(6,m,l)=dduxyjc
            vjc2(6,m,l)=ddvxyjc
            pjc2(8,m,l)=ddpyjc0
         enddo
         ujc2(3,0,l)=ujc2(3,ns,l)
         vjc2(3,0,l)=vjc2(3,ns,l)
         pjc2(3,0,l)=pjc2(3,ns,l)
         ujc2(5,0,l)=ujc2(5,ns,l)
         vjc2(5,0,l)=vjc2(5,ns,l)
         pjc2(7,0,l)=pjc2(7,ns,l)

         ujc2(4,0,l)=ujc2(4,ns,l)
         vjc2(4,0,l)=vjc2(4,ns,l)
         pjc2(4,0,l)=pjc2(4,ns,l)
         ujc2(6,0,l)=ujc2(6,ns,l)
         vjc2(6,0,l)=vjc2(6,ns,l)
         pjc2(8,0,l)=pjc2(8,ns,l)
      ENDDO
      DO l=1,ms
        do m=0,ns
          ujc(1,m,l)=(ujc(1,m,l)+ujc2(1,m,l))/2.0d0
          ujc(2,m,l)=(ujc(2,m,l)+ujc2(2,m,l))/2.0d0
          ujc(3,m,l)=(ujc(3,m,l)+ujc2(3,m,l))/2.0d0
          ujc(4,m,l)=(ujc(4,m,l)+ujc2(4,m,l))/2.0d0
          ujc(5,m,l)=(ujc(5,m,l)+ujc2(5,m,l))/2.0d0
          ujc(6,m,l)=(ujc(6,m,l)+ujc2(6,m,l))/2.0d0
         
          vjc(1,m,l)=(vjc(1,m,l)+vjc2(1,m,l))/2.0d0
          vjc(2,m,l)=(vjc(2,m,l)+vjc2(2,m,l))/2.0d0
          vjc(3,m,l)=(vjc(3,m,l)+vjc2(3,m,l))/2.0d0
          vjc(4,m,l)=(vjc(4,m,l)+vjc2(4,m,l))/2.0d0
          vjc(5,m,l)=(vjc(5,m,l)+vjc2(5,m,l))/2.0d0
          vjc(6,m,l)=(vjc(6,m,l)+vjc2(6,m,l))/2.0d0
         
          pjc(1,m,l)=(pjc(1,m,l)+pjc2(1,m,l))/2.0d0
          pjc(2,m,l)=(pjc(2,m,l)+pjc2(2,m,l))/2.0d0
          pjc(3,m,l)=(pjc(3,m,l)+pjc2(3,m,l))/2.0d0
          pjc(4,m,l)=(pjc(4,m,l)+pjc2(4,m,l))/2.0d0
          pjc(5,m,l)=(pjc(5,m,l)+pjc2(5,m,l))/2.0d0
          pjc(6,m,l)=(pjc(6,m,l)+pjc2(6,m,l))/2.0d0
          pjc(7,m,l)=(pjc(7,m,l)+pjc2(7,m,l))/2.0d0
          pjc(8,m,l)=(pjc(8,m,l)+pjc2(8,m,l))/2.0d0
        enddo
      ENDDO

      return
      end


c-----------------------------------------------------------------------
