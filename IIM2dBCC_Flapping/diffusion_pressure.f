c-----------------------------------------------------------------------
c
      subroutine diffusion_pressure(gama,krk,fac)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'rhsp.inc'
      include 'old.inc'
      integer ie,ic,ic1,je,jc,jc1,krk
      real*8 signx,signy,sx,sy,gama,fac
      real*8 f(6),xx,yy,xsm0,xsm1,ysm0,ysm1
      real*8 pw(ny),pe(ny),ps(nx),pn(nx)

      call jc_firstsecond
      call jc_p

c mac_distribute
      
      DO l=1,ms

      do i=0,ncxc(l) 
        m=falfaxc(0,i,l)
        xx=falfaxc(15,i,l)
        xsm0=xs(m,l)
        ysm0=ys(m,l)
        xsm1=xs(m+1,l)
        ysm1=ys(m+1,l)

        signx=falfaxc(2,i,l)
        signy=falfaxc(3,i,l)

        do n=1,6
          f(n)=pjc(n,m,l)+(pjc(n,m+1,l)-pjc(n,m,l))*
     .         (xx-xsm0)/(xsm1-xsm0)
        enddo
        do n=1,6
          pjcxc(n,i,l)=f(n)
        enddo

        pjcxc(1,i,l)=signx*pjcxc(1,i,l)
        pjcxc(3,i,l)=signx*pjcxc(3,i,l)
        pjcxc(5,i,l)=signx*pjcxc(5,i,l)

        pjcxc(2,i,l)=signy*pjcxc(2,i,l)
        pjcxc(4,i,l)=signy*pjcxc(4,i,l)
        pjcxc(6,i,l)=signy*pjcxc(6,i,l)

      enddo

      do j=0,ncyc(l)
        m=falfayc(0,i,l)
        yy=falfayc(15,i,l)
        xsm0=xs(m,l)
        ysm0=ys(m,l)
        xsm1=xs(m+1,l)
        ysm1=ys(m+1,l)
        signx=falfayc(2,j,l)
        signy=falfayc(3,j,l)


        do n=1,6
          f(n)=pjc(n,m,l)+(pjc(n,m+1,l)-pjc(n,m,l))*
     .         (yy-ysm0)/(ysm1-ysm0)
        enddo
        do n=1,6
          pjcyc(n,j,l)=f(n)
        enddo

        pjcyc(1,j,l)=signx*pjcyc(1,j,l)
        pjcyc(3,j,l)=signx*pjcyc(3,j,l)
        pjcyc(5,j,l)=signx*pjcyc(5,j,l)

        pjcyc(2,j,l)=signy*pjcyc(2,j,l)
        pjcyc(4,j,l)=signy*pjcyc(4,j,l)
        pjcyc(6,j,l)=signy*pjcyc(6,j,l)

      enddo

      ENDDO

c correction_difference

      DO l=1,ms

      do i=0,ncxc(l)
        sy=falfaxc(1,i,l)
        je=jexc(i,l)
        jc=jcxc(i,l)
        jc1=jc+1
        pcdyy(1,i,l)=-dy2*(pjcxc(6,i,l)+
     .                     pjcxc(2,i,l)*(yc(jc1)-sy)+
     .               0.5d0*pjcxc(4,i,l)*(yc(jc1)-sy)**2.0d0)
        pcdyy(2,i,l)=dy2*(pjcxc(6,i,l)+
     .                    pjcxc(2,i,l)*(yc(jc)-sy)+
     .               0.5d0*pjcxc(4,i,l)*(yc(jc)-sy)**2.0d0)
        if(jc.eq.je) then
          vcpdy(i,l)=-dy1*(pjcxc(6,i,l)+
     .                     pjcxc(2,i,l)*(yc(jc1)-sy)+
     .               0.5d0*pjcxc(4,i,l)*(yc(jc1)-sy)**2.0d0)
        else 
          vcpdy(i,l)=-dy1*(pjcxc(6,i,l)+
     .                     pjcxc(2,i,l)*(yc(jc)-sy)+
     .               0.5d0*pjcxc(4,i,l)*(yc(jc)-sy)**2.0d0)
        endif
      enddo

      do j=0,ncyc(l)
        sx=falfayc(1,j,l)
        ie=ieyc(j,l)
        ic=icyc(j,l)
        ic1=ic+1
        pcdxx(1,j,l)=-dx2*(pjcyc(5,j,l)+
     .                     pjcyc(1,j,l)*(xc(ic1)-sx)+
     .               0.5d0*pjcyc(3,j,l)*(xc(ic1)-sx)**2.0d0)
        pcdxx(2,j,l)=dx2*(pjcyc(5,j,l)+
     .                    pjcyc(1,j,l)*(xc(ic)-sx)+
     .               0.5d0*pjcyc(3,j,l)*(xc(ic)-sx)**2.0d0)
        if(ic.eq.ie) then
          ucpdx(j,l)=-dx1*(pjcyc(5,j,l)+
     .                     pjcyc(1,j,l)*(xc(ic1)-sx)+
     .               0.5d0*pjcyc(3,j,l)*(xc(ic1)-sx)**2.0d0)
        else
          ucpdx(j,l)=-dx1*(pjcyc(5,j,l)+
     .                     pjcyc(1,j,l)*(xc(ic)-sx)+
     .               0.5d0*pjcyc(3,j,l)*(xc(ic)-sx)**2.0d0)
        endif
      enddo

      ENDDO

c  pressure

      do j=1,ny
        do i=1,nx
          d(i,j)=dx1*(u(i,j)-u(i-1,j))+dy1*(v(i,j)-v(i,j-1))
        enddo
      enddo
      
      if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxc(l)
            d(ixc(i,l),jexc(i,l)+1)=d(ixc(i,l),jexc(i,l)+1)+pcvdy(i,l)
          enddo
          do j=0,ncyc(l)
            d(ieyc(j,l)+1,jyc(j,l))=d(ieyc(j,l)+1,jyc(j,l))+pcudx(j,l)
          enddo
         enddo
       endif
       
      call divergence_reset
      
      do j=1,ny
        do i=1,nx
        prhs(i,j)=
     .             dn(i,j)/(fac*dt)
     .            +(dx2*Re1)*(d(i+1,j)-2.0d0*d(i,j)+d(i-1,j))
     .            +(dy2*Re1)*(d(i,j+1)-2.0d0*d(i,j)+d(i,j-1))
     .            -2.0d0*(d(i,j)*d(i,j)+
     .             (ucc(i,j)/dx)*0.5d0*(d(i+1,j)-d(i-1,j))+
     .             (vcc(i,j)/dy)*0.5d0*(d(i,j+1)-d(i,j-1)))
     .            -2.0d0*(dvx(i,j)*duy(i,j)-dux(i,j)*dvy(i,j))
        enddo
       enddo
     
       if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxc(l)
            prhs(ixc(i,l),jcxc(i,l))=prhs(ixc(i,l),jcxc(i,l))-
     .                               pcdyy(1,i,l)
            prhs(ixc(i,l),jcxc(i,l)+1)=prhs(ixc(i,l),jcxc(i,l)+1)-
     .                                 pcdyy(2,i,l)
          enddo
          do j=0,ncyc(l)
            prhs(icyc(j,l),jyc(j,l))=prhs(icyc(j,l),jyc(j,l))-
     .                               pcdxx(1,j,l)
            prhs(icyc(j,l)+1,jyc(j,l))=prhs(icyc(j,l)+1,jyc(j,l))-
     .                                 pcdxx(2,j,l)
          enddo
        enddo
      endif
                   
      do j=1,ny
        do i=1,nx
          o(i,j)=p(i,j)
          prhs(i,j)=prhs(i,j)-gama*p(i,j)
        enddo
      enddo

      call pbc_fft(pw,pe,ps,pn,krk)
      call helmholtz_fft(pw,pe,ps,pn,gama,
     .                   lw_pbc,le_pbc,ls_pbc,ln_pbc,1)


      return
      end


c-----------------------------------------------------------------------
