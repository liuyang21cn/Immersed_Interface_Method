c-----------------------------------------------------------------------
c
      subroutine vbc(krk,fac)
      include 'parameter.inc'
      include 'field.inc'
      integer krk
      real*8 fac

c  periodic
      if(lw_vbc.eq.0.and.le_vbc.eq.0) then
        do j=0,ny
          v(0,j)=v(nx-1,j)
          v(nx+1,j)=v(2,j)
        enddo
      endif

      if(ls_vbc.eq.0.and.ln_vbc.eq.0) then
        do i=0,nx+1
          v(i,0)=v(i,ny-1)
          v(i,ny)=v(i,1)
        enddo
      endif

c  dirichlet
      if(lw_vbc.eq.1) then
        do j=0,ny
          v(0,j)=2.0d0*0.0d0-v(2,j)
        enddo
      endif

      if(le_vbc.eq.1) then
        do j=0,ny
          v(nx+1,j)=2.0d0*0.0d0-v(nx-1,j)
        enddo
      endif

      if(ls_vbc.eq.1) then
        do i=0,nx+1
          v(i,0)=2.0d0*0.0d0-v(i,1)
        enddo
      endif

      if(ln_vbc.eq.1) then
        do i=0,nx+1
          v(i,ny)=2.0d0*0.0d0-v(i,ny-1)
        enddo
      endif

c  newmann
      if(lw_vbc.eq.2) then
        do j=0,ny
          v(0,j)=v(2,j)
        enddo
      endif

      if(le_vbc.eq.2) then
        do j=0,ny
          v(nx+1,j)=v(nx-1,j)
        enddo
      endif

      if(ls_vbc.eq.2) then
        do i=0,nx+1
          v(i,0)=v(i,1)-dy*0.0d0
        enddo
      endif

      if(ln_vbc.eq.2) then
        do i=0,nx+1
          v(i,ny)=v(i,ny-1)+dy*0.0d0
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
