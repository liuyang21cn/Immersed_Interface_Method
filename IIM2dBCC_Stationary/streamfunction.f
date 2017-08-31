c-----------------------------------------------------------------------
c
      subroutine streamfunction
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer le,lw,ls,ln
      real*8 phe(ny),phw(ny),phs(nx),phn(nx)

      do j=1,ny
        do i=1,nx
          o(i,j)=(vec(i,j)-vec(i-1,j))/dx-
     .           (uce(i,j)-uce(i,j-1))/dy
        enddo
      enddo
      if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxc(l)
            o(ixc(i,l),jexc(i,l)+1)=o(ixc(i,l),jexc(i,l)+1)-pcudy(i,l)
          enddo
          do j=0,ncyc(l)
            o(ieyc(j,l)+1,jyc(j,l))=o(ieyc(j,l)+1,jyc(j,l))+pcvdx(j,l)
          enddo
        enddo
      endif

      lw=2
      le=2
      do j=1,ny
        if(le.eq.2) phe(j)=vcc(1,j)
        if(lw.eq.2) phw(j)=vcc(nx,j)
        if(le.eq.1) phe(j)=0.0d0
        if(lw.eq.1) phw(j)=yc(j)*
     .                     (1.0d0-0.25d0/(xc(1)*xc(1)+yc(j)*yc(j)))
      enddo
      
      ls=2
      ln=2
      do i=1,nx
        if(ls.eq.2) phs(i)=-ucc(i,1)
        if(ln.eq.2) phn(i)=-ucc(i,ny)
        if(ls.eq.1) phs(i)=0.0d0
        if(ln.eq.1) phn(i)=0.0d0
      enddo

      call poisson_fft(phw,phe,phs,phn,lw,le,ls,ln,0)

      return
      end


c-----------------------------------------------------------------------
