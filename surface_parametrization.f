c-----------------------------------------------------------------------
c
      subroutine surface_parametrization
      include 'parameter.inc'
      include 'surface.inc'
      character*1 fnext
      character*32 fname

      open(unit=13,file='DAT/shape0.dat',status='unknown')
     
c 1: cylinder; 2: square; 3: rectangle; 4: rounded plate

      DO l=1,ms
        
        write(unit=fnext,fmt='(i1)') jectob(l)
        fname='DAT/Vertices/object'//fnext//'.run'
        open(unit=31,file=fname,status='old')
        rewind(31) 
        read(31,*)
        do m=0,ns
          read(31,*)xs0(m,l),ys0(m,l),curv(m,l)
        enddo
        close(31)
    
        do m=0,ns
          xs(m,l)=xsc(l)+
     +             xs0(m,l)*dcos(theta(l))-ys0(m,l)*dsin(theta(l))
          ys(m,l)=ysc(l)+
     +             xs0(m,l)*dsin(theta(l))+ys0(m,l)*dcos(theta(l))
        enddo
        xs(ns,l)=xs(0,l)
        ys(ns,l)=ys(0,l)

        do m=0,ns
          write(13,100)xs(m,l),ys(m,l)
        enddo

      ENDDO

      close(13)
100   format(1x,30e16.7)

      return
      end


c-----------------------------------------------------------------------
