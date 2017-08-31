c-----------------------------------------------------------------------
c
      subroutine pbc_fft(pw,pe,ps,pn,krk)
      include 'parameter.inc'
      include 'field.inc'
      integer krk
      real*8 pw(ny),pe(ny),ps(nx),pn(nx)
      real*8 o1,o2,r1,r2,a,b,foo,rr

      o1=1.0d0
      o2=-1.0d0
      r1=0.5d0
      r2=2.0d0
      a=(o2*r2*r2-o1*r1*r1)/(r2*r2-r1*r1)
      b=(o1-o2)*r1*r1*r2*r2/(r2*r2-r1*r1)

c  periodic
      do j=1,ny
        pw(j)=0.0d0
        pe(j)=0.0d0
      enddo
      do i=1,nx
        ps(i)=0.0d0
        pn(i)=0.0d0
      enddo

c  dirichlet 
      if(lw_pbc.eq.1) then
        do j=1,ny
          pw(j)=-(0.125d0*(yc(j)*yc(j)-xc(1)*xc(1))+
     .             0.25d0*(yc(j)*xc(1)))/
     .           (yc(j)*yc(j)+xc(1)*xc(1))**2.0d0
          pw(j)=0.0d0
        enddo
      endif

      if(le_pbc.eq.1) then
        do j=1,ny
          pe(j)=0.0d0
        enddo
      endif

      if(ls_pbc.eq.1) then
        do i=1,nx
          ps(i)=0.0d0
        enddo
      endif

      if(ln_pbc.eq.1) then
        do i=1,nx
          pn(i)=0.0d0
        enddo
      endif

c  newmann 
      if(lw_pbc.eq.2) then
        do j=1,ny
          pw(j)=(dx2/Re)*(ucc(2,j)-ucc(1,j)
     .         +(dx/dy)*(vee(0,j)-vee(0,j-1)))
     .         +(dy2/Re)*(ucc(1,j+1)-2.0d0*ucc(1,j)+ucc(1,j-1))
     .         -dx1*(u(1,j)*u(1,j)-u(0,j)*u(0,j))
     .         -dy1*(uce(1,j)*v(1,j)-uce(1,j-1)*v(1,j-1))
c           rr=xc(1)*xc(1)+yc(j)*yc(j)
c           pw(j)=(a*a+b*b/(rr*rr)+2.0d0*a*b/rr)*xc(1)
        enddo
      endif

      if(le_pbc.eq.2) then
        do j=1,ny
          pe(j)=(dx2/Re)*(ucc(nx-1,j)-ucc(nx,j)
     .         -(dx/dy)*(vee(nx,j)-vee(nx,j-1)))
     .         +(dy2/Re)*(ucc(nx,j+1)-2.0d0*ucc(nx,j)+ucc(nx,j-1))
     .         -dx1*(u(nx,j)*u(nx,j)-u(nx-1,j)*u(nx-1,j))
     .         -dy1*(uce(nx,j)*v(nx,j)-uce(nx,j-1)*v(nx,j-1))
c           rr=xc(nx)*xc(nx)+yc(j)*yc(j)
c           pe(j)=(a*a+b*b/(rr*rr)+2.0d0*a*b/rr)*xc(nx)
        enddo
      endif

      if(ls_pbc.eq.2) then
        do i=1,nx
          ps(i)=(dy2/Re)*(vcc(i,2)-vcc(i,1)
     .         +(dy/dx)*(uee(i,0)-uee(i-1,0)))
     .         +(dx2/Re)*(vcc(i+1,1)-2.0d0*vcc(i,1)+vcc(i-1,1))
     .         -dy1*(v(i,1)*v(i,1)-v(i,0)*v(i,0))
     .         -dx1*(vec(i,1)*u(i,1)-vec(i-1,1)*u(i-1,1))
c           rr=xc(i)*xc(i)+yc(1)*yc(1)
c           ps(i)=(a*a+b*b/(rr*rr)+2.0d0*a*b/rr)*yc(1)
        enddo
      endif

      if(ln_pbc.eq.2) then
        do i=1,nx
          pn(i)=(dy2/Re)*(vcc(i,ny-1)-vcc(i,ny)
     .         -(dy/dx)*(uee(i,ny)-uee(i-1,ny)))
     .         +(dx2/Re)*(vcc(i+1,ny)-2.0d0*vcc(i,ny)+vcc(i-1,ny))
     .         -dy1*(v(i,ny)*v(i,ny)-v(i,ny-1)*v(i,ny-1))
     .         -dx1*(vec(i,ny)*u(i,ny)-vec(i-1,ny)*u(i-1,ny))
c           rr=xc(i)*xc(i)+yc(ny)*yc(ny)
c           pn(i)=(a*a+b*b/(rr*rr)+2.0d0*a*b/rr)*yc(ny)
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
