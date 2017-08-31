c-----------------------------------------------------------------------
c
      subroutine data_plot
      include 'parameter.inc'
      include 'field.inc'
      include 'surface.inc'
      include 'rhsp.inc'

      integer ih,jh,iref,igap,iend,jgap,jend,ii,jj,mark1,mark2,iexact
      real*8 foo,unorm,vnorm,pnorm1,pnorm2,cp1,cp2,xcr,ycr,d0,pmax,pmin
      real*8 au(0:nx,0:ny+1),av(0:nx+1,0:ny),ap(0:nx+1,0:ny+1)
      real*8 o1,o2,r1,r2,a,b,rr

      o1=1.0d0
      o2=-1.0d0
      r1=0.5d0
      r2=2.0d0
      a=(o2*r2*r2-o1*r1*r1)/(r2*r2-r1*r1)
      b=(o1-o2)*r1*r1*r2*r2/(r2*r2-r1*r1)


      call field_interpolate
      call streamfunction
      call surface_vorticity

      open(unit=29,file='DAT/xe.dat',status='unknown')
      open(unit=39,file='DAT/ye.dat',status='unknown')
      do i=0,nx
        write(29,200)xe(i)
      enddo
      do j=0,ny
        write(39,200)ye(j)
      enddo
      close(29)
      close(39)
      open(unit=49,file='DAT/xc.dat',status='unknown')
      open(unit=59,file='DAT/yc.dat',status='unknown')
      do i=1,nx
        write(49,200)xc(i)
      enddo
      do j=1,ny
        write(59,200)yc(j)
      enddo
      close(49)
      close(59)
200   format(1x,10e16.6)

      ih=(0.0d0-x0)/dx
      jh=(0.0d0-y0)/dy
      open(unit=19,file='DAT/gy.dat',status='unknown')
      open(unit=29,file='DAT/gx.dat',status='unknown')
      do j=1,ny
        write(19,200)yc(j),ucc(ih,j),vcc(ih,j),p(ih,j),o(ih,j)
      enddo
      do i=1,nx
        write(29,200)xc(i),ucc(i,jh),vcc(i,jh),p(i,jh),o(i,jh)
      enddo
      close(19)
      close(29)

      open(unit=28,file='DAT/rhsp.dat',status='unknown')
      open(unit=38,file='DAT/p.dat',status='unknown')
      open(unit=48,file='DAT/d.dat',status='unknown')
      open(unit=58,file='DAT/uc.dat',status='unknown')
      open(unit=68,file='DAT/vc.dat',status='unknown')
      open(unit=78,file='DAT/wo.dat',status='unknown')
      open(unit=88,file='DAT/ph.dat',status='unknown')
      open(unit=98,file='DAT/swo.dat',status='unknown')
      do j=1,ny
        write(28,400)(prhs(i,j),i=1,nx)
        write(38,400)(p(i,j),i=1,nx)
        write(48,400)(d(i,j),i=1,nx)
        write(58,400)(ucc(i,j),i=1,nx)
        write(68,400)(vcc(i,j),i=1,nx)
        write(78,400)(o(i,j),i=1,nx)
        write(88,400)(-ph(i,j),i=1,nx)
      enddo
      do l=1,ms
        write(98,400)(swo(m,l),m=0,ns)
      enddo
      close(28)
      close(38)
      close(48)
      close(58)
      close(68)
      close(78)
      close(88)
      close(98)
400   format(1x,2000e16.6)

      return
      end


c-----------------------------------------------------------------------
