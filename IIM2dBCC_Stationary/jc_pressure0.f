c-----------------------------------------------------------------------
c
      subroutine jc_pressure0
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      ! coordinates of middle point
      real*8 Mx(0:ns,ms), My(0:ns,ms)
      ! variables for princilple jump conditions
      real*8 tx,ty,tx2,ty2,xn,yn,gacobi
      real*8 dudnjcp,dudnjcn,dudnjc,tx1,ty1
      real*8 dvdnjcp,dvdnjcn,dvdnjc
      real*8 ddudnp,ddvdnp,deltaujc,deltavjc,qx,qy
      real*8 temp1, temp2,temp,radius,ssl
      real*8 dpxjcM(0:ns,ms),dpyjcM(0:ns,ms)
      
      ! variables for interpolation
      real*8 ds,ds1,ds2,su,sv,xx,yy,uu(3),vv(3),dist,foo
      integer many,nany,ie,ic,je,jc,id,jd,iu,ju,iv,jv,io,jo
      parameter (many=3,nany=3)
      real*8 xa(many),ya(many),xb(many),yb(many)

      ! matrix variables
      real*8 a1(0:ns-1),b1(0:ns-1),c1(0:ns-1)
      real*8 d1(0:ns-1),x1(0:ns-1)


      ds=1.01d0*sqrt(dx*dx+dy*dy)
      DO l=1,ms
         do m=0,ns-1

         ! normal direction
         xn=taoy(m,l)
         yn=-taox(m,l)
         gacobi=dsqrt(xn*xn+yn*yn)         
         xn=xn/gacobi
         yn=yn/gacobi
       
         ! middle point coordinates
         Mx(m,l)=(xs(m+1,l)+xs(m,l))/2.0d0
         My(m,l)=(ys(m+1,l)+ys(m,l))/2.0d0

         ! velocity at middle point
         su = (us(m+1,l)+us(m,l))/2.0d0
         sv = (vs(m+1,l)+vs(m,l))/2.0d0

         ! 1.Compute jump condition of dux,duy,dvx,dvy at middle point
         ! 1.1 jump condition of dudn
         ! 1.1.1 dudn jump conditions on negative side
         dudnjcn = -thetat(l)*yn
         dvdnjcn =  thetat(l)*xn

         ! 1.1.2 dudn jump conditions on positive side
         ! a. interpolating 3 points on normal direction

         do n=1,nany
            
            ! coordinates of interpolation point
            xx=Mx(m,l)+dble(n)*ds*xn
            yy=My(m,l)+dble(n)*ds*yn
            ! indices of interpolation point
            i=int((xx-x0)/hdx)
            j=int((yy-y0)/hdx)
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
        ! cylinder curv is 1/radius
        
        ! 1.2.2 ddudn jump conditions on positive side
        ddudnp=(-5.0d0*uu(1)+4.0d0*uu(2)-1.0d0*uu(3)+2.0d0*su)/ds/ds
        ddvdnp=(-5.0d0*vv(1)+4.0d0*vv(2)-1.0d0*vv(3)+2.0d0*sv)/ds/ds

        ! 1.2.3 delta u jump conditions
        deltaujc = ddudnp+curv(m,l)*dudnjc
        deltavjc = ddvdnp+curv(m,l)*dvdnjc

        ! 1.4 dpdx dpdy jump conditions
        ! 1.4.1 body force qx qy
        qx =  thetatt(l)*(My(m,l)-ysc(l))
        qy = -thetatt(l)*(Mx(m,l)-xsc(l))

        ! 1.4.2 dpdxM, dpdyM jump conditions
        dpxjcM(m,l) = 1.0d0/re*deltaujc + qx
        dpyjcM(m,l) = 1.0d0/re*deltavjc + qy 

        enddo ! end vertex loop
        
        dpxjcM(ns,l)=dpxjcM(0,l)
        dpyjcM(ns,l)=dpyjcM(0,l)

        ! 2. Assemble pressure matrix
        ! 2.1 Assemble rhs
        do m=0,ns-1

        tx = taox(m,l)
        ty = taoy(m,l)
        gacobi=dsqrt(tx*tx+ty*ty)
        tx=tx/gacobi
        ty=ty/gacobi

        ! a. compute distance of two panels
        ds1 = dsqrt((xs(m+1,l)-xs(m,l))*(xs(m+1,l)-xs(m,l))+
     .              (ys(m+1,l)-ys(m,l))*(ys(m+1,l)-ys(m,l)))

        if(m.eq.0) then
          ds2=dsqrt((xs(ns,l)-xs(ns-1,l))*(xs(ns,l)-xs(ns-1,l))+
     .              (ys(ns,l)-ys(ns-1,l))*(ys(ns,l)-ys(ns-1,l)))
        else
          ds2=dsqrt((xs(m,l)-xs(m-1,l))*(xs(m,l)-xs(m-1,l))+
     .              (ys(m,l)-ys(m-1,l))*(ys(m,l)-ys(m-1,l)))
        endif

        ! b. compute rhs by simpson rule
        temp1=(pjc(1,m,l)+4.0d0*dpxjcM(m,l)+pjc(1,m+1,l))*tx+
     .        (pjc(2,m,l)+4.0d0*dpyjcM(m,l)+pjc(2,m+1,l))*ty
        temp1=temp1*ds1/6.0d0

        if(m.eq.0) then
           tx2=-taox(ns-1,l)
           ty2=-taoy(ns-1,l)
           gacobi=dsqrt(tx2*tx2+ty2*ty2)
           tx2=tx2/gacobi
           ty2=ty2/gacobi
           temp2=(pjc(1,0,l)+4.0d0*dpxjcM(ns-1,l)+pjc(1,ns-1,l))*tx2+
     .           (pjc(2,0,l)+4.0d0*dpyjcM(ns-1,l)+pjc(2,ns-1,l))*ty2
           temp2=temp2*ds2/6.0d0
        else
           tx2=-taox(m-1,l)
           ty2=-taoy(m-1,l)
           gacobi=dsqrt(tx2*tx2+ty2*ty2)
           tx2=tx2/gacobi
           ty2=ty2/gacobi
           temp2=(pjc(1,m,l)+4.0d0*dpxjcM(m-1,l)+pjc(1,m-1,l))*tx2+
     .           (pjc(2,m,l)+4.0d0*dpyjcM(m-1,l)+pjc(2,m-1,l))*ty2
           temp2=temp2*ds2/6.0d0
        endif

        d1(m)=temp1+temp2
        ! end do m=0,ns-1
        enddo 

        ! 2.2 solve this matrix problem
        ! set up matrix A
        a1(0:ns-1)=1.0d0
        b1(0:ns-1)=-2.0d0
        c1(0:ns-1)=1.0d0

        ! set last x is 0
        x1(ns-1)=0
        ! change matrix to upper diagonal
        c1(0)=c1(0)/b1(0)
        d1(0)=d1(0)/b1(0)
        do m=1,ns-3
           temp=b1(m)-a1(m)*c1(m-1)
           c1(m)=c1(m)/temp
           d1(m)=(d1(m)-a1(m)*d1(m-1))/temp
        enddo
           d1(ns-2)=(d1(ns-2)-a1(ns-2)*d1(ns-3))/
     .              (b1(ns-2)-a1(ns-2)*c1(ns-3))

        ! back substitute
        x1(ns-2)=d1(ns-2)
        do m=ns-3,0,-1
           x1(m)=d1(m)-c1(m)*x1(m+1)
        enddo

        do m=0,ns-1
        pjc(5,m,l)=x1(m)
        pjc(6,m,l)=x1(m)
        enddo
        pjc(5,ns,l)=pjc(5,0,l)
        pjc(6,ns,l)=pjc(6,0,l)


      ENDDO ! End objects loop
      return
      end


c-----------------------------------------------------------------------
