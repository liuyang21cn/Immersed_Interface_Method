c-----------------------------------------------------------------------
c
      subroutine helmholtz_fft(pw,pe,ps,pn,gama,lw,le,ls,ln,iflag)
      include 'parameter.inc'
      integer ni,nj,lw,le,ls,ln,iflag
      real*8 gama
      real*8 pw(ny),pe(ny),ps(nx),pn(nx)

      if(lw.eq.0.and.le.eq.0.and.ls.eq.0.and.ln.eq.0) then
        ni=nx-1
        nj=ny-1
        call hbc0000(pw,pe,ps,pn,gama,ni,nj,iflag)
      endif

      if(lw.eq.1.and.le.eq.1.and.ls.eq.1.and.ln.eq.1) then
        ni=nx-2
        nj=ny-2
        call hbc1111(pw,pe,ps,pn,gama,ni,nj,iflag)
      endif

      if(lw.eq.2.and.le.eq.2.and.ls.eq.2.and.ln.eq.2) then
        ni=nx
        nj=ny
        call hbc2222(pw,pe,ps,pn,gama,ni,nj,iflag)
      endif

      if(lw.eq.1.and.le.eq.2.and.ls.eq.2.and.ln.eq.2) then
        ni=nx-1
        nj=ny
        call hbc1222(pw,pe,ps,pn,gama,ni,nj,iflag)
      endif      

      if(lw.eq.1.and.le.eq.2.and.ls.eq.0.and.ln.eq.0) then
        ni=nx-1
        nj=ny-1
        call hbc1200(pw,pe,ps,pn,gama,ni,nj,iflag)
      endif

      return
      end


c-----------------------------------------------------------------------
c
      subroutine hbc1200(pw,pe,ps,pn,gama,ni,nj,iflag)
      include 'parameter.inc'
      include 'field.inc'
      include 'rhsp.inc'
      include 'fft1200.inc'
      integer ni,nj,jd,iflag
      real*8 beta,pai,fac,gama
      real*8 pw(ny),pe(ny),ps(nx),pn(nx)
      parameter(pai=
     .3.1415926535897932384626433832795028841971693993751058209749446D0)

      beta=dx/dy

      if(iflag.eq.1) then
        do i=1,ni
          do j=1,nj
            ri(j,i)=dx*dx*prhs(i+1,j)
          enddo
        enddo
      else
        do i=1,ni
          do j=1,nj
            ri(j,i)=dx*dx*o(i+1,j)
          enddo
        enddo
      endif

      do j=1,nj
        ri(j,1)=ri(j,1)-pw(j)
        ri(j,ni)=ri(j,ni)-2.0d0*dx*pe(j)
      enddo

      call vsinqf(nj,ni,ri,wi,nj,wsavei)

      do j=1,nj
        do i=1,ni
          rj(i,j)=ri(j,i)
        enddo
      enddo

      call vrfftf(ni,nj,rj,wj,ni,wsavej)

      jd=int(nj/2)
      if((nj-jd*2).eq.1) jd=jd+1 

      do i=1,ni
        rj(i,1)=0.5d0*rj(i,j)/(dcos(dble(2*i-1)*pai/dble(2*ni))-1.0d0)
        if(nj.eq.2*jd) then
           rj(i,nj)=0.5d0*rj(i,nj)
     .      /(dcos(dble(2*i-1)*pai/dble(2*ni))-1.0d0-2.0d0*beta*beta)
        endif
        do j=1,jd-1
          fac=0.5d0/(dcos(dble(2*i-1)*pai/dble(2*ni))+beta*beta*
     .               dcos(dble(2*j)*pai/dble(nj))-(1.0d0+beta*beta))
          rj(i,2*j)=fac*rj(i,2*j)
          rj(i,2*j+1)=fac*rj(i,2*j+1)
        enddo
      enddo

      call vrfftb(ni,nj,rj,wj,ni,wsavej)

      do j=1,nj
        do i=1,ni
          ri(j,i)=rj(i,j)
        enddo
      enddo

      call vsinqb(nj,ni,ri,wi,nj,wsavei)

      if(iflag.eq.1) then
        do j=1,ny-1
          p(1,j)=pw(j)
          do i=2,nx
            p(i,j)=ri(j,i-1)
          enddo
        enddo
        do i=1,nx
          p(i,ny)=p(i,1)
        enddo
      else
        do j=1,ny-1
          ph(1,j)=pw(j)
          do i=2,nx
            ph(i,j)=ri(j,i-1)
          enddo
        enddo
        do i=1,nx
          ph(i,ny)=ph(i,1)
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
c
      subroutine hbc1222(pw,pe,ps,pn,gama,ni,nj,iflag)
      include 'parameter.inc'
      include 'field.inc'
      include 'rhsp.inc'
      include 'fft1222.inc'
      integer ni,nj,iflag
      real*8 beta,pai,gama
      real*8 pw(ny),pe(ny),ps(nx),pn(nx)
      parameter(pai=
     .3.1415926535897932384626433832795028841971693993751058209749446D0)

      beta=dx/dy

      if(iflag.eq.1) then
        do i=1,ni
          do j=1,nj
            ri(j,i)=dx*dx*prhs(i+1,j)
          enddo
        enddo
       else
        do i=1,ni
          do j=1,nj
            ri(j,i)=dx*dx*o(i+1,j)
          enddo
        enddo
       endif

      do j=1,nj
        ri(j,1)=ri(j,1)-pw(j)
        ri(j,ni)=ri(j,ni)-2.0d0*dx*pe(j)
      enddo
      do i=1,ni
        ri(1,i)=ri(1,i)+beta*beta*2.0d0*dy*ps(i+1)
        ri(nj,i)=ri(nj,i)-beta*beta*2.0d0*dy*pn(i+1)
      enddo

      call vsinqf(nj,ni,ri,wi,nj,wsavei)

      do j=1,nj
        do i=1,ni
          rj(i,j)=ri(j,i)
        enddo
      enddo

      call vcost(ni,nj,rj,wj,ni,wsavej)

      do j=1,nj
        do i=1,ni
          rj(i,j)=0.5d0*rj(i,j)/
     .            (dcos(dble(2*i-1)*pai/dble(2*ni))+beta*beta*
     .             dcos(dble(j-1)*pai/dble(nj-1))-
     .             (1.0d0+beta*beta+0.5d0*gama*dx*dx))
        enddo
      enddo

      call vcost(ni,nj,rj,wj,ni,wsavej)

      do j=1,nj
        do i=1,ni
          ri(j,i)=rj(i,j)
        enddo
      enddo

      call vsinqb(nj,ni,ri,wi,nj,wsavei)

      if(iflag.eq.1) then
        do j=1,ny
          p(1,j)=pw(j)
          do i=2,nx
            p(i,j)=ri(j,i-1)
          enddo
        enddo
      else
        do j=1,ny
          ph(1,j)=pw(j)
          do i=2,nx
            ph(i,j)=ri(j,i-1)
          enddo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
c
      subroutine hbc2222(pw,pe,ps,pn,gama,ni,nj,iflag)
      include 'parameter.inc'
      include 'field.inc'
      include 'rhsp.inc'
      include 'fft2222.inc'
      integer ni,nj,iflag
      real*8 beta,pai,gama
      real*8 pw(ny),pe(ny),ps(nx),pn(nx)
      parameter(pai=
     .3.1415926535897932384626433832795028841971693993751058209749446D0)

      beta=dx/dy

      if(iflag.eq.1) then
        do i=1,ni
          do j=1,nj
            ri(j,i)=dx*dx*prhs(i,j)
          enddo
        enddo
      else
        do i=1,ni
          do j=1,nj
            ri(j,i)=dx*dx*o(i,j)
          enddo
        enddo
      endif

      do j=1,nj
        ri(j,1)=ri(j,1)+2.0d0*dx*pw(j)
        ri(j,ni)=ri(j,ni)-2.0d0*dx*pe(j)
      enddo

      do i=1,ni
        ri(1,i)=ri(1,i)+beta*beta*2.0d0*dy*ps(i)
        ri(nj,i)=ri(nj,i)-beta*beta*2.0d0*dy*pn(i)
      enddo

      call vcost(nj,ni,ri,wi,nj,wsavei)

      do j=1,nj
        do i=1,ni
          rj(i,j)=ri(j,i)
        enddo
      enddo

      call vcost(ni,nj,rj,wj,ni,wsavej)

      do j=1,nj
        do i=1,ni
          if(i.eq.1.and.j.eq.1) then
            rj(i,j)=0.0d0
          else
            rj(i,j)=0.5d0*rj(i,j)/
     .              (dcos(dble(i-1)*pai/dble(ni-1))+beta*beta*
     .               dcos(dble(j-1)*pai/dble(nj-1))-
     .               (1.0d0+beta*beta+0.5d0*gama*dx*dx))
          endif
        enddo
      enddo

      call vcost(ni,nj,rj,wj,ni,wsavej)

      do j=1,nj
        do i=1,ni
          ri(j,i)=rj(i,j)
        enddo
      enddo

      call vcost(nj,ni,ri,wi,nj,wsavei)

      if(iflag.eq.1) then
        do j=1,ny
          do i=1,nx
            p(i,j)=ri(j,i)
          enddo
        enddo
      else
        do j=1,ny
          do i=1,nx
            ph(i,j)=ri(j,i)
          enddo
        enddo
      endif
    
      return
      end


c-----------------------------------------------------------------------
c
      subroutine hbc1111(pw,pe,ps,pn,gama,ni,nj,iflag)
      include 'parameter.inc'
      include 'field.inc'
      include 'rhsp.inc'
      include 'fft1111.inc'
      integer ni,nj,iflag
      real*8 beta,pai,gama
      real*8 pw(ny),pe(ny),ps(nx),pn(nx)
      parameter(pai=
     .3.1415926535897932384626433832795028841971693993751058209749446D0)

      beta=dx/dy

      if(iflag.eq.1) then
        do i=1,ni
          do j=1,nj
            ri(j,i)=dx*dx*prhs(i+1,j+1)
          enddo
        enddo
      else
        do i=1,ni
          do j=1,nj
            ri(j,i)=dx*dx*o(i+1,j+1)
          enddo
        enddo
      endif

      do j=1,nj
        ri(j,1)=ri(j,1)-pw(j+1)
        ri(j,ni)=ri(j,ni)-pe(j+1)
      enddo
      do i=1,ni
        ri(1,i)=ri(1,i)-beta*beta*ps(i+1)
        ri(nj,i)=ri(nj,i)-beta*beta*pn(i+1)
      enddo

      call vsint(nj,ni,ri,wi,nj,wsavei)

      do j=1,nj
        do i=1,ni
          rj(i,j)=ri(j,i)
        enddo
      enddo

      call vsint(ni,nj,rj,wj,ni,wsavej)

      do j=1,nj
        do i=1,ni
          rj(i,j)=0.5d0*rj(i,j)/
     .            (dcos(dble(i)*pai/dble(ni+1))+beta*beta*
     .             dcos(dble(j)*pai/dble(nj+1))-
     .             (1.0d0+beta*beta+0.5d0*gama*dx*dx))
        enddo
      enddo

      call vsint(ni,nj,rj,wj,ni,wsavej)

      do j=1,nj
        do i=1,ni
          ri(j,i)=rj(i,j)
        enddo
      enddo

      call vsint(nj,ni,ri,wi,nj,wsavei)

      if(iflag.eq.1) then
        do j=1,ny
          p(1,j)=pw(j)
          p(nx,j)=pe(j)
        enddo
        do i=1,nx
          p(i,1)=ps(i)
          p(i,ny)=pn(i)
        enddo
        do j=2,ny-1
          do i=2,nx-1
            p(i,j)=ri(j-1,i-1)
          enddo
        enddo
      else
        do j=1,ny
          ph(1,j)=pw(j)
          ph(nx,j)=pe(j)
        enddo
        do i=1,nx
          ph(i,1)=ps(i)
          ph(i,ny)=pn(i)
        enddo
        do j=2,ny-1
          do i=2,nx-1
            ph(i,j)=ri(j-1,i-1)
          enddo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
c
      subroutine hbc0000(pw,pe,ps,pn,gama,ni,nj,iflag)
      include 'parameter.inc'
      include 'field.inc'
      include 'rhsp.inc'
      include 'fft0000.inc'
      integer ni,nj,id,jd,iflag
      real*8 beta,pai,fac,gama
      real*8 pw(ny),pe(ny),ps(nx),pn(nx)
      parameter(pai=
     .3.1415926535897932384626433832795028841971693993751058209749446D0)

      beta=dx/dy

      if(iflag.eq.1) then
        do j=1,nj
          do i=1,ni
            ri(j,i)=dx*dx*prhs(i,j)
          enddo
        enddo
      else
        do j=1,nj
          do i=1,ni
            ri(j,i)=dx*dx*o(i,j)
          enddo
        enddo
      endif

      call vrfftf(nj,ni,ri,wi,nj,wsavei)

      do j=1,nj
        do i=1,ni
          rj(i,j)=ri(j,i)
        enddo
      enddo

      call vrfftf(ni,nj,rj,wj,ni,wsavej)

      id=int(ni/2)
      jd=int(nj/2)
      if((ni-id*2).eq.1) id=id+1
      if((nj-jd*2).eq.1) jd=jd+1 

      rj(1,1)=0.0d0

      do j=1,jd-1
        fac=0.5d0/(beta*beta*(dcos(dble(2*j)*pai/dble(nj))-1.0d0))
        rj(1,2*j)=fac*rj(1,2*j)
        rj(1,2*j+1)=fac*rj(1,2*j+1)
        fac=0.5d0/
     .      (beta*beta*(dcos(dble(2*j)*pai/dble(nj))-1.0d0)-2.0d0)
        rj(ni,2*j)=fac*rj(ni,2*j)
        rj(ni,2*j+1)=fac*rj(ni,2*j+1)
      enddo
      if(nj.eq.2*jd) rj(1,nj)=-0.25d0*rj(1,nj)/(beta*beta)

      do i=1,id-1
        fac=0.5/(dcos(dble(2*i)*pai/dble(ni))-1.0d0)
        rj(2*i,1)=fac*rj(2*i,1)
        rj(2*i+1,1)=fac*rj(2*i+1,1)
        fac=0.5d0/
     .      (dcos(dble(2*i)*pai/dble(ni))-1.0d0-2.0d0*beta*beta)
        rj(2*i,nj)=fac*rj(2*i,nj)
        rj(2*i+1,nj)=fac*rj(2*i+1,nj)
      enddo
      if(ni.eq.2*id) rj(ni,1)=-0.25d0*rj(ni,1)

      if(ni.eq.2*id.and.nj.eq.2*jd) then
        rj(ni,nj)=-0.25d0*rj(ni,nj)/(1.0d0+beta*beta)
      endif

      do j=1,jd-1
        do i=1,id-1
            fac=0.5d0/(dcos(dble(2*i)*pai/dble(ni))+beta*beta*
     .                 dcos(dble(2*j)*pai/dble(nj))-(1.0d0+beta*beta))
            rj(2*i,2*j)=fac*rj(2*i,2*j)
            rj(2*i,2*j+1)=fac*rj(2*i,2*j+1)
            rj(2*i+1,2*j)=fac*rj(2*i+1,2*j)
            rj(2*i+1,2*j+1)=fac*rj(2*i+1,2*j+1)
        enddo
      enddo

      call vrfftb(ni,nj,rj,wj,ni,wsavej)

      do j=1,nj
        do i=1,ni
          ri(j,i)=rj(i,j)
        enddo
      enddo

      call vrfftb(nj,ni,ri,wi,nj,wsavei)

      if(iflag.eq.1) then
        do j=1,ny-1
          do i=1,nx-1
            p(i,j)=ri(j,i)
          enddo
        enddo
        do j=1,ny-1
          p(nx,j)=p(1,j)
        enddo
        do i=1,nx
          p(i,ny)=p(i,1)
        enddo
      else
        do j=1,ny-1
          do i=1,nx-1
            ph(i,j)=ri(j,i)
          enddo
        enddo
        do j=1,ny-1
          ph(nx,j)=ph(1,j)
        enddo
        do i=1,nx
          ph(i,ny)=ph(i,1)
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
