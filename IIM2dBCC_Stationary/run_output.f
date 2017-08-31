c-----------------------------------------------------------------------
c
      subroutine run_output(tstart,tend,total,day1,now1,day2,now2)
      include 'parameter.inc'
      include 'surface.inc'
      character*8 day1,day2
      integer*4 nday1(3), now1(3),nday2(3),now2(3)
      real*4 total
      real*8 tstart,tend

      close(16)
      close(26)
      close(36)
      close(46)
      close(56)
      close(66)
      close(76)
      close(86)
!      close(96)
!      close(97)
!      close(98)

      if(iIter.eq.1) then
        close(87) !IterSum
      endif

      open(unit=8,file='DAT/output.run',status='unknown')
      rewind 8
      write(8,*)
      write(8,*)'Re             =       ',Re
      write(8,*)'x0             =       ',x0
      write(8,*)'xl             =       ',xl
      write(8,*)'y0             =       ',y0
      write(8,*)'yl             =       ',yl
      write(8,*)
      write(8,*)'nx             =       ',nx
      write(8,*)'ny             =       ',ny
      write(8,*)'ns             =       ',ns
      write(8,*)'ms             =       ',ms
      write(8,*)'dx             =       ',dx
      write(8,*)'dy             =       ',dy
      write(8,*)
      write(8,*)'icfl           =       ',icfl
      write(8,*)'cflc           =       ',cflc
      write(8,*)'cflv           =       ',cflv
      write(8,*)'nstart         =       ',nstart
      write(8,*)'nend           =       ',nend
      write(8,*)'nstep          =       ',nstep
      write(8,*)'tstart         =       ',tstart
      write(8,*)'tend           =       ',tend
      write(8,*)'dtcfl          =       ',dtcfl
      write(8,*)'dt0            =       ',dt0
      write(8,*)
      write(8,*)'isingular      =       ',isingular
      write(8,*)'move           =       ',move
      write(8,*)
      write(8,*)'total time     =       ',total/3600.0d0,' hours'
      write(8,*)'beginning time :'
c     write(8, 1000) nday1(2), nday1(1), nday1(3), now1
      write(8,*)'Date: ',day1
      write(8,2000) now1
      write(8,*)'ending time    :'
c     write(8, 1000) nday2(2), nday2(1), nday2(3), now2
      write(8,*)'Date: ',day2
      write(8,2000) now2
 1000 format (' Date ', i2.2, '/', i2.2, '/', i4.4, '; Time ',
     .         i2.2, ':', i2.2, ':', i2.2 )
 2000 format (' Time: ',i2.2, ':', i2.2, ':', i2.2 )
      close(8)

      return
      end


c-----------------------------------------------------------------------

