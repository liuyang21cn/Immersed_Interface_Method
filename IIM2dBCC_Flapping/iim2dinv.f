c----------------------------------------------------------------------
c
      program main
      include 'parameter.inc'
      include 'field.inc'
      include 'surface.inc'
      character*8 day1,day2
      integer*4 nday1(3), now1(3),nday2(3),now2(3)
      real*4 total,etime,cost(2)
      real*8 tn,tmp,tstart,tend

c     call idate(nday1)
      call date_and_time(day1)
      call itime(now1)
      call run_input

      call mesh
      call initial

      nstart=0
      if(iread.eq.1) then
        call data_read
      endif
      nstart=nstart+1
      nend=nstart+nstep-1
      tstart=t
      tmp=t

      do n=nstart,nend
        call cfl
        write(*,*)
        write(*,*)'! n = ',n,' t = ',t
        call rk
        if(itout.eq.1.and.mod(n,1).eq.0) then
          call time_output
        endif
        if(ianimation.eq.1) then
          tn=t-dt
          if(tmp.ge.tn.and.tmp.lt.t)then
            write(*,*)'  !! animation'
            tmp=tmp+tout
            call animation
          endif
        endif
      enddo

      tend=t

      if(iwrite.eq.1) then
        call data_write
      endif

      if(iplot.eq.1) then
        call surface_plot
        call data_plot
      endif

c     call idate(nday2)
      call date_and_time(day2)
      call itime(now2)
      total=etime(cost)
      call run_output(tstart,tend,total,day1,now1,day2,now2)

      write(*,*)
      write(*,*)'-- beginning time --'
c     write(*, 1000) nday1(2), nday1(1), nday1(3), now1
      write(*,*)'Date: ',day1
      write(*,2000) now1
      write(*,*)'-- ending time --'
c     write(*, 1000) nday2(2), nday2(1), nday2(3), now2
      write(*,*)'Date: ',day2
      write(*,2000) now2
      write(*,*)'user time   = ',cost(1)/3600.0d0,' hours'
      write(*,*)'system time = ',cost(2)/3600.0d0,' hours'
      write(*,*)'total time  = ',total/3600.0d0,' hours'
 1000 format (' Date ', i2.2, '/', i2.2, '/', i4.4, '; Time ',
     .         i2.2, ':', i2.2, ':', i2.2 )
 2000 format (' Time: ',i2.2, ':', i2.2, ':', i2.2 )

      stop
      end


c-----------------------------------------------------------------------
