c-----------------------------------------------------------------------
c
      subroutine euler_link
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer iu,id,ju,jd
      real*8 xsm0,xsm1,ysm0,ysm1,tangent,signx,signy
      real*8 f(4)
     
 
      IOxcyc=0
      IOxcye=0
      IOxeyc=0
 
      DO l=1,ms

      ncxe(l)=0
      ncxc(l)=0
      ncye(l)=0
      ncyc(l)=0

      do m=0,ns-1

        xsm0=xs(m,l)
        xsm1=xs(m+1,l)
        ysm0=ys(m,l)
        ysm1=ys(m+1,l)
        
        if(xsm0.gt.xsm1) then
          iu=int((xsm0-x0)/hdx)
          id=int((xsm1-x0)/hdx)+1
          if(iu.ge.id) then
            do i=iu,id,-1
              signx=sign(1.0d0,ysm1-ysm0)
              signy=sign(1.0d0,xsm0-xsm1)
              if(m.eq.0) then
                tangent=signy+sign(1.0d0,xs(ns-1,l)-xs(ns,l))
              else
                tangent=signy+sign(1.0d0,xs(m-1,l)-xsm0)
              endif
              if(x(i).ne.xsm0.or.tangent.ne.0.0) then
                f(2)=ysm0+(ysm1-ysm0)*(x(i)-xsm0)/(xsm1-xsm0)
                f(3)=us(m,l)+(us(m+1,l)-us(m,l))*(x(i)-xsm0)/(xsm1-xsm0)
                f(4)=vs(m,l)+(vs(m+1,l)-vs(m,l))*(x(i)-xsm0)/(xsm1-xsm0)
                j=int((f(2)-y0)/hdy)
                if(f(2).eq.y(j)) j=j-max(0,int(signy))
                if(mod(i,2).eq.0) then
                  ixe(ncxe(l),l)=i/2
                  if(mod(j,2).eq.0) then
                    jexe(ncxe(l),l)=j/2
                    jcxe(ncxe(l),l)=j/2
                  else
                    jexe(ncxe(l),l)=(j-1)/2
                    jcxe(ncxe(l),l)=(j+1)/2
                  endif
                  IOxeyc(ixe(ncxe(l),l),0:jcxe(ncxe(l),l))=l-
     .            IOxeyc(ixe(ncxe(l),l),0:jcxe(ncxe(l),l))
                  falfaxe(0,ncxe(l),l)=m     
                  falfaxe(1,ncxe(l),l)=f(2)
                  falfaxe(2,ncxe(l),l)=signx
                  falfaxe(3,ncxe(l),l)=signy
                  falfaxe(13,ncxe(l),l)=f(3)
                  falfaxe(14,ncxe(l),l)=f(4)
                  falfaxe(15,ncxe(l),l)=x(i)
                  ncxe(l)=ncxe(l)+1
                else
                  ixc(ncxc(l),l)=(i+1)/2
                  if(mod(j,2).eq.0) then
                    jexc(ncxc(l),l)=j/2
                    jcxc(ncxc(l),l)=j/2
                  else
                    jexc(ncxc(l),l)=(j-1)/2
                    jcxc(ncxc(l),l)=(j+1)/2
                  endif
                  IOxcye(ixc(ncxc(l),l),0:jexc(ncxc(l),l))=l-
     .            IOxcye(ixc(ncxc(l),l),0:jexc(ncxc(l),l))
                  IOxcyc(ixc(ncxc(l),l),0:jcxc(ncxc(l),l))=l-
     .            IOxcyc(ixc(ncxc(l),l),0:jcxc(ncxc(l),l))
                  falfaxc(0,ncxc(l),l)=m
                  falfaxc(1,ncxc(l),l)=f(2)
                  falfaxc(2,ncxc(l),l)=signx
                  falfaxc(3,ncxc(l),l)=signy
                  falfaxc(13,ncxc(l),l)=f(3)
                  falfaxc(14,ncxc(l),l)=f(4)
                  falfaxc(15,ncxc(l),l)=x(i)
                  ncxc(l)=ncxc(l)+1
                endif
              endif
            enddo
          endif
        elseif(xsm1.gt.xsm0) then
          id=int((xsm0-x0)/hdx)+1
          iu=int((xsm1-x0)/hdx)
          if(iu.ge.id) then
            do i=id,iu
              signx=sign(1.0d0,ysm1-ysm0)
              signy=sign(1.0d0,xsm0-xsm1)
              if(m.eq.ns-1) then
                tangent=signy+sign(1.0d0,xs(0,l)-xs(1,l))
              else
                tangent=signy+sign(1.0d0,xsm1-xs(m+2,l))
              endif
              if(x(i).ne.xsm1.or.tangent.ne.0.0) then
                f(2)=ysm0+(ysm1-ysm0)*(x(i)-xsm0)/(xsm1-xsm0)
                f(3)=us(m,l)+(us(m+1,l)-us(m,l))*(x(i)-xsm0)/(xsm1-xsm0)
                f(4)=vs(m,l)+(vs(m+1,l)-vs(m,l))*(x(i)-xsm0)/(xsm1-xsm0)
                j=int((f(2)-y0)/hdy)
                if(f(2).eq.y(j)) j=j-max(0,int(signy))
                if(mod(i,2).eq.0) then
                  ixe(ncxe(l),l)=i/2
                  if(mod(j,2).eq.0) then
                    jexe(ncxe(l),l)=j/2
                    jcxe(ncxe(l),l)=j/2
                  else
                    jexe(ncxe(l),l)=(j-1)/2
                    jcxe(ncxe(l),l)=(j+1)/2
                  endif
                  IOxeyc(ixe(ncxe(l),l),0:jcxe(ncxe(l),l))=l-
     .            IOxeyc(ixe(ncxe(l),l),0:jcxe(ncxe(l),l))
                  falfaxe(0,ncxe(l),l)=m
                  falfaxe(1,ncxe(l),l)=f(2)
                  falfaxe(2,ncxe(l),l)=signx
                  falfaxe(3,ncxe(l),l)=signy
                  falfaxe(13,ncxe(l),l)=f(3)
                  falfaxe(14,ncxe(l),l)=f(4)
                  falfaxe(15,ncxe(l),l)=x(i)
                  ncxe(l)=ncxe(l)+1
                else
                  ixc(ncxc(l),l)=(i+1)/2
                  if(mod(j,2).eq.0) then
                    jexc(ncxc(l),l)=j/2
                    jcxc(ncxc(l),l)=j/2
                  else
                    jexc(ncxc(l),l)=(j-1)/2
                    jcxc(ncxc(l),l)=(j+1)/2
                  endif
                  IOxcye(ixc(ncxc(l),l),0:jexc(ncxc(l),l))=l-
     .            IOxcye(ixc(ncxc(l),l),0:jexc(ncxc(l),l))
                  IOxcyc(ixc(ncxc(l),l),0:jcxc(ncxc(l),l))=l-
     .            IOxcyc(ixe(ncxc(l),l),0:jcxc(ncxc(l),l))
                  falfaxc(0,ncxc(l),l)=m
                  falfaxc(1,ncxc(l),l)=f(2)
                  falfaxc(2,ncxc(l),l)=signx
                  falfaxc(3,ncxc(l),l)=signy
                  falfaxc(13,ncxc(l),l)=f(3)
                  falfaxc(14,ncxc(l),l)=f(4)
                  falfaxc(15,ncxc(l),l)=x(i)
                  ncxc(l)=ncxc(l)+1
                endif
              endif
            enddo
          endif
        endif

        if(ysm0.gt.ysm1) then
          ju=int((ysm0-y0)/hdy)
          jd=int((ysm1-y0)/hdy)+1
          if(ju.ge.jd) then
            do j=ju,jd,-1
              signx=sign(1.0d0,ysm1-ysm0)
              signy=sign(1.0d0,xsm0-xsm1)
              if(m.eq.0) then
                tangent=signx+sign(1.0d0,ys(ns,l)-ys(ns-1,l))
              else
                tangent=signx+sign(1.0d0,ysm0-ys(m-1,l))
              endif
              if(y(j).ne.ysm0.or.tangent.ne.0.0) then
                f(1)=xsm0+(xsm1-xsm0)*(y(j)-ysm0)/(ysm1-ysm0)
                f(3)=us(m,l)+(us(m+1,l)-us(m,l))*(y(j)-ysm0)/(ysm1-ysm0)
                f(4)=vs(m,l)+(vs(m+1,l)-vs(m,l))*(y(j)-ysm0)/(ysm1-ysm0)
                i=int((f(1)-x0)/hdx)
                if(f(1).eq.x(i)) i=i-max(0,int(signx))
                if(mod(j,2).eq.0) then
                  jye(ncye(l),l)=j/2
                  if(mod(i,2).eq.0) then
                    ieye(ncye(l),l)=i/2
                    icye(ncye(l),l)=i/2
                  else
                    ieye(ncye(l),l)=(i-1)/2
                    icye(ncye(l),l)=(i+1)/2
                  endif
                  falfaye(0,ncye(l),l)=m
                  falfaye(1,ncye(l),l)=f(1)
                  falfaye(2,ncye(l),l)=signx
                  falfaye(3,ncye(l),l)=signy
                  falfaye(13,ncye(l),l)=f(3)
                  falfaye(14,ncye(l),l)=f(4)
                  falfaye(15,ncye(l),l)=y(j)
                  ncye(l)=ncye(l)+1
                else
                  jyc(ncyc(l),l)=(j+1)/2
                  if(mod(i,2).eq.0) then
                    ieyc(ncyc(l),l)=i/2
                    icyc(ncyc(l),l)=i/2
                  else
                    ieyc(ncyc(l),l)=(i-1)/2
                    icyc(ncyc(l),l)=(i+1)/2
                  endif
                  falfayc(0,ncyc(l),l)=m
                  falfayc(1,ncyc(l),l)=f(1)
                  falfayc(2,ncyc(l),l)=signx
                  falfayc(3,ncyc(l),l)=signy
                  falfayc(13,ncyc(l),l)=f(3)
                  falfayc(14,ncyc(l),l)=f(4)
                  falfayc(15,ncyc(l),l)=y(j)
                  ncyc(l)=ncyc(l)+1
                endif
              endif
            enddo
          endif
        elseif(ysm1.gt.ysm0) then
          jd=int((ysm0-y0)/hdy)+1
          ju=int((ysm1-y0)/hdy)
          if(ju.ge.jd) then
            do j=jd,ju
              signx=sign(1.0d0,ysm1-ysm0)
              signy=sign(1.0d0,xsm0-xsm1)
              if(m.eq.ns-1) then
                tangent=signx+sign(1.0d0,ys(1,l)-ys(0,l))
              else
                tangent=signx+sign(1.0d0,ys(m+2,l)-ysm1)
              endif
              if(y(j).ne.ysm1.or.tangent.ne.0.0) then
                f(1)=xsm0+(xsm1-xsm0)*(y(j)-ysm0)/(ysm1-ysm0)
                f(3)=us(m,l)+(us(m+1,l)-us(m,l))*(y(j)-ysm0)/(ysm1-ysm0)
                f(4)=vs(m,l)+(vs(m+1,l)-vs(m,l))*(y(j)-ysm0)/(ysm1-ysm0)
                i=int((f(1)-x0)/hdx)
                if(f(1).eq.x(i)) i=i-max(0,int(signx))
                if(mod(j,2).eq.0) then
                  jye(ncye(l),l)=j/2
                  if(mod(i,2).eq.0) then
                    ieye(ncye(l),l)=i/2
                    icye(ncye(l),l)=i/2
                  else
                    ieye(ncye(l),l)=(i-1)/2
                    icye(ncye(l),l)=(i+1)/2
                  endif
                  falfaye(0,ncye(l),l)=m
                  falfaye(1,ncye(l),l)=f(1)
                  falfaye(2,ncye(l),l)=signx
                  falfaye(3,ncye(l),l)=signy
                  falfaye(13,ncye(l),l)=f(3)
                  falfaye(14,ncye(l),l)=f(4)
                  falfaye(15,ncye(l),l)=y(j)
                  ncye(l)=ncye(l)+1
                else
                  jyc(ncyc(l),l)=(j+1)/2
                  if(mod(i,2).eq.0) then
                    ieyc(ncyc(l),l)=i/2
                    icyc(ncyc(l),l)=i/2
                  else
                    ieyc(ncyc(l),l)=(i-1)/2
                    icyc(ncyc(l),l)=(i+1)/2
                  endif
                  falfayc(0,ncyc(l),l)=m
                  falfayc(1,ncyc(l),l)=f(1)
                  falfayc(2,ncyc(l),l)=signx
                  falfayc(3,ncyc(l),l)=signy
                  falfayc(13,ncyc(l),l)=f(3)
                  falfayc(14,ncyc(l),l)=f(4)
                  falfayc(15,ncyc(l),l)=y(j)
                  ncyc(l)=ncyc(l)+1
                endif
              endif
            enddo
          endif
        endif

      enddo

      ncxe(l)=ncxe(l)-1
      ncxc(l)=ncxc(l)-1
      ncye(l)=ncye(l)-1
      ncyc(l)=ncyc(l)-1

      ENDDO

      return
      end


c-----------------------------------------------------------------------
