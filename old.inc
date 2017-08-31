c  old.inc
      real*8 un(0:nx,0:ny+1),vn(0:nx+1,0:ny),dn(nx,ny)
      real*8 xsn(0:ns,ms),ysn(0:ns,ms),usrk(4,0:ns,ms),vsrk(4,0:ns,ms)
      real*8 urk(0:nx,0:ny+1,4),vrk(0:nx+1,0:ny,4)

      common/old/un,vn,dn,xsn,ysn,usrk,vsrk,urk,vrk
