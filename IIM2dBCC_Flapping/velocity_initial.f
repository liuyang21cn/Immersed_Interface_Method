c-----------------------------------------------------------------------
c
      subroutine velocity_initial
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'rhsp.inc'
      include 'old.inc'
      real*8 ds,xx,yy,grad
      real*8 signx,signy
      real*8 f(6),sx,sy
      real*8 pw(ny),pe(ny),ps(nx),pn(nx)
      real*8 xsm0,xsm1,ysm0,ysm1
      integer ic,ic1,ie,jc,jc1,je

      call surface_property
      call euler_link

      ds=1.01d0*dsqrt(dx*dx+dy*dy)

c jc_firstsecond
      call jc_firstsecond
c jc_pressure0
      call jc_pressure0 
c mac_distribute
      
      DO l=1,ms

      do i=0,ncxc(l) 
        m=falfaxc(0,i,l)
        xx=falfaxc(15,i,l)
        xsm0=xs(m,l)
        xsm1=xs(m+1,l)
        ysm0=ys(m,l)
        ysm1=ys(m+1,l)
        
        do n=1,6
        pjcxc(n,i,l)=pjc(n,m,l)+(pjc(n,m+1,l)-pjc(n,m,l))*
     .            (xx-xsm0)/(xsm1-xsm0)
        enddo

        pjcxc(1,i,l)=signx*pjcxc(1,i,l)
        pjcxc(3,i,l)=signx*pjcxc(3,i,l)
        pjcxc(5,i,l)=signx*pjcxc(5,i,l)

        pjcxc(2,i,l)=signy*pjcxc(2,i,l)
        pjcxc(4,i,l)=signy*pjcxc(4,i,l)
        pjcxc(6,i,l)=signy*pjcxc(6,i,l)

      enddo

      do j=0,ncyc(l)
        m=falfayc(0,j,l)
        yy=falfayc(15,j,l)
        xsm0=xs(m,l)
        xsm1=xs(m+1,l)
        ysm0=ys(m,l)
        ysm1=ys(m+1,l)
        do n=1,6
        pjcyc(n,j,l)=pjc(n,m,l)+(pjc(n,m+1,l)-pjc(n,m,l))*
     .            (yy-ysm0)/(ysm1-ysm0)
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
     .              0.5d0*pjcxc(4,i,l)*(yc(jc)-sy)**2.0d0)
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
     .              0.5d0*pjcyc(3,j,l)*(xc(ic)-sx)**2.0d0)
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
          prhs(i,j)=0.0d0
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
        pw(j)=0.0d0
        pe(j)=0.0d0
      enddo
      do i=1,nx
        ps(i)=0.0d0
        pn(i)=0.0d0
      enddo

      call poisson_fft(pw,pe,ps,pn,2,2,2,2,1)

c correction:

      do j=1,ny
        do i=1,nx-1
          grad=-dx1*(p(i+1,j)-p(i,j))
          u(i,j)=u(i,j)+grad
        enddo
      enddo

      do j=1,ny-1
        do i=1,nx
          grad=-dy1*(p(i,j+1)-p(i,j))
          v(i,j)=v(i,j)+grad
        enddo
      enddo

      if(isingular.eq.1) then
        do l=1,ms
          do j=0,ncyc(l)
            u(icyc(j,l),jyc(j,l))=u(icyc(j,l),jyc(j,l))
     .                            -ucpdx(j,l)
          enddo
          do i=0,ncxc(l)
            v(ixc(i,l),jcxc(i,l))=v(ixc(i,l),jcxc(i,l))
     .                            -vcpdy(i,l)
          enddo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
