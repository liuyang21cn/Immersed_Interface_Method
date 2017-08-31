c-----------------------------------------------------------------------
c
      subroutine data_plot
      include 'parameter.inc'
      include 'field.inc'
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

      open(unit=38,file='DAT/p.dat',status='unknown')
      open(unit=48,file='DAT/d.dat',status='unknown')
      open(unit=58,file='DAT/uc.dat',status='unknown')
      open(unit=68,file='DAT/vc.dat',status='unknown')
      open(unit=78,file='DAT/wo.dat',status='unknown')
      open(unit=88,file='DAT/ph.dat',status='unknown')
      do j=1,ny
        write(38,400)(p(i,j),i=1,nx)
        write(48,400)(d(i,j),i=1,nx)
        write(58,400)(ucc(i,j),i=1,nx)
        write(68,400)(vcc(i,j),i=1,nx)
        write(78,400)(o(i,j),i=1,nx)
        write(88,400)(-ph(i,j),i=1,nx)
      enddo
      close(38)
      close(48)
      close(58)
      close(68)
      close(78)
      close(88)
400   format(1x,2000e16.6)

      iref=0
      igap=2**2
      iend=igap*(nx-1)+1
      jgap=2**2
      jend=jgap*(ny-1)+1
      if(iref.eq.1.or.iref.eq.-1) then
        open(unit=17,file='DAT/au.dat',status='unknown')
        open(unit=27,file='DAT/av.dat',status='unknown')
        open(unit=37,file='DAT/ap.dat',status='unknown')
      endif
      if(iref.eq.1) then
        do j=1,ny
          write(17,500)(ucc(i,j),i=1,nx)
          write(27,500)(vcc(i,j),i=1,nx)
          write(37,500)(p(i,j),i=1,nx)
        enddo
      endif
      if(iref.eq.-1) then
        do j=1,jend
           do i=1,iend
             if(mod(j-1,jgap).eq.0.and.mod(i-1,igap).eq.0) then
               jj=(j-1)/jgap+1
               ii=(i-1)/igap+1
               read(17,500)au(ii,jj)
               read(27,500)av(ii,jj)
               read(37,500)ap(ii,jj)
             else
               read(17,500)foo
               read(27,500)foo
               read(37,500)foo
             endif
           enddo
         enddo
      endif
      if(iref.eq.1.or.iref.eq.-1) then
        close(17)
        close(27)
        close(37)
      endif

      iexact=0
      mark1=0
      mark2=0
      d0=0.0d0
      if(iexact.eq.1) then
        do j=1,ny
          do i=1,nx
            xcr=xc(i)-d0*dsin(t)
            ycr=yc(j)-d0*dsin(t)
            if((xcr*xcr+ycr*ycr).gt.r1*r1*(1.0d0+0.00d0)) then
              if(mark1.eq.0) then
                cp1=-(0.5d0*a*a*(xcr*xcr+ycr*ycr)-
     .                0.5d0*b*b/(xcr*xcr+ycr*ycr)+
     .                a*b*dlog(xcr*xcr+ycr*ycr))-
     .                d0*dsin(t)*(xc(i)+yc(j))+p(i,j)
                mark1=1
              endif
              au(i,j)=(a+
     .                b/(xcr*xcr+ycr*ycr))*(-ycr)+
     .                d0*dcos(t)
              av(i,j)=(a+
     .                b/(xcr*xcr+ycr*ycr))*(xcr)+
     .                d0*dcos(t)
              ap(i,j)=(0.5d0*a*a*(xcr*xcr+ycr*ycr)-
     .                 0.5d0*b*b/(xcr*xcr+ycr*ycr)+
     .                 a*b*dlog(xcr*xcr+ycr*ycr))+
     .                d0*dsin(t)*(xc(i)+yc(j))+cp1
c            elseif((xcr*xcr+ycr*ycr).lt.r1*r1*(1.0d0-0.01d0)) then
c              if(mark2.eq.0) then
c                cp2=-0.5d0*(xcr*xcr+ycr*ycr)-
c     .               d0*dsin(t)*(xc(i)+yc(j))+p(i,j)
c                mark2=1
c              endif
c              au(i,j)=-ycr+d0*dcos(t)
c              av(i,j)=xcr+d0*dcos(t)
c              ap(i,j)=0.5d0*(xcr*xcr+ycr*ycr)+
c     .                d0*dsin(t)*(xc(i)+yc(j))+cp2
            else
              au(i,j)=ucc(i,j)
              av(i,j)=vcc(i,j)
              ap(i,j)=p(i,j)
            endif
          enddo
        enddo
      endif

      unorm=0.0d0
      vnorm=0.0d0
      pnorm1=0.0d0
      pnorm2=0.0d0
      pmax=0.0d0
      pmin=0.0d0
      do j=1,ny
        do i=1,nx
          if(iexact.eq.1) then
            ucc(i,j)=ucc(i,j)-au(i,j)
            vcc(i,j)=vcc(i,j)-av(i,j)
            p(i,j)=p(i,j)-ap(i,j)
          endif
          unorm=max(unorm,abs(ucc(i,j)))
          vnorm=max(vnorm,abs(vcc(i,j)))
          pnorm1=max(pnorm1,abs(p(i,j)))
          pmax=max(pmax,p(i,j))
          pmin=min(pmin,p(i,j))
        enddo
      enddo
      pnorm2=0.5d0*(pmax-pmin)

      write(*,*)
      write(*,*)'unorm  = ',unorm
      write(*,*)'vnorm  = ',vnorm
      write(*,*)'pnorm1 = ',pnorm1
      write(*,*)'pnorm2 = ',pnorm2

      iref=-1
      if(iref.eq.-1) then
        open(unit=19,file='DAT/eu.dat',status='unknown')
        open(unit=29,file='DAT/ev.dat',status='unknown')
        open(unit=39,file='DAT/ep.dat',status='unknown')
        do j=1,ny
          write(19,400)(ucc(i,j),i=1,nx)
          write(29,400)(vcc(i,j),i=1,nx)  
          write(39,400)(p(i,j),i=1,nx)
        enddo
        close(19)
        close(29)
        close(39)
      endif
500   format(1x,e36.18)

      return
      end


c-----------------------------------------------------------------------
