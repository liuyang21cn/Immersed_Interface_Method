c  parameter.inc
      implicit none

      integer i,j,k,l,m,n,nstart,nend,nstep,nx,ny,ms,ns

      integer icfl,iread,iwrite,iplot,itout,ianimation
      integer isingular,move,IterSum,iIter
      integer lw_ubc,le_ubc,ls_ubc,ln_ubc
      integer lw_vbc,le_vbc,ls_vbc,ln_vbc
      integer lw_pbc,le_pbc,ls_pbc,ln_pbc
 
      real*8 pi,Re,Re1,cflc,cflv,dt,dtcfl,dt0,dt1,dt2,t,tout,t0
      real*8 x0,y0,xl,yl,dx,dy,dx1,dx2,dy1,dy2,hdx,hdy

      parameter(pi=
     .3.1415926535897932384626433832795028841971693993751058209749446d0)
      
      ! if iIter=0, it's fft solver, if =1, it's helmholtz solver
      parameter(Re=40.0d0,Re1=1.0d0/Re)
      parameter(x0=-10.0d0,y0=-7.5d0)
      parameter(xl=40.0d0,yl=15.0d0)
      parameter(nx=1200+1,ny=450+1,ms=1,ns=256)


      parameter(isingular=0,move=0,iIter=0)
      parameter(lw_ubc=1,le_ubc=2,ls_ubc=1,ln_ubc=1)
      parameter(lw_vbc=1,le_vbc=2,ls_vbc=1,ln_vbc=1)
      parameter(lw_pbc=2,le_pbc=2,ls_pbc=2,ln_pbc=2)

      common/run1/nstep,nstart,nend,
     .           iread,iwrite,itout,iplot,ianimation,icfl
      common/run2/cflc,cflv,dt,dtcfl,dt0,dt1,dt2,t,tout,t0,
     .           dx,dy,dx1,dy1,dx2,dy2,hdx,hdy
