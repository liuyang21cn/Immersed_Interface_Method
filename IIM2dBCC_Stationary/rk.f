c-----------------------------------------------------------------------
c
      subroutine rk
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'old.inc'
      integer iflag,krk,iter
      real*8 fac,err,tol,gama

c  krk=1:
      krk=1
      fac=0.5d0

      call surface_property
      call euler_link
      if(isingular.eq.1) then
        call singular_call
      else
        call field_interpolate
        call uv_strain
      endif
   
      call old_save
    
       if(isingular.eq.1) then   
         call jc_pressure
         call correction_pressure
       endif
       call pressure(krk,fac)

      if(iIter.eq.1) then 
!       do j=1,ny
!         do i=1,nx
!           p(i,j)=0.0d0
!         enddo
!       enddo

      IterSum=0           
      iter=0
      gama=1.0d-2
      tol=dx*dx
      err=1.0d6
9     continue
      if(err.gt.tol) then
        call diffusion_pressure(gama,krk,fac)
        err=0.0d0
        do j=1,ny
          do i=1,nx
            err=max(err,abs(p(i,j)-o(i,j)))
           enddo
        enddo
        iter=iter+1
        write(*,*)'iter=',iter,'err=',err
        goto 9
       endif
       IterSum=IterSum+iter
       endif

      call rhs(krk)
1     call pgrad(krk)

      if(move.eq.1) then
        call surface_move(krk,fac)
      endif

      do j=1,ny
        do i=1,nx-1
          u(i,j)=un(i,j)+fac*dt*urk(i,j,1)
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          v(i,j)=vn(i,j)+fac*dt*vrk(i,j,1)
        enddo
      enddo 
      call velocity_reset
      call ubc(krk,fac)
      call vbc(krk,fac)

c  krk=2:
      krk=2
      fac=0.5d0

      call surface_property
      call euler_link
      if(isingular.eq.1) then
        call singular_call
      else
        call field_interpolate
        call uv_strain
      endif

      if(isingular.eq.1) then
        call jc_pressure
        call correction_pressure
      endif
      call pressure(krk,fac)

      if(iIter.eq.1) then
!      do j=1,ny
!        do i=1,nx
!          p(i,j)=0.0d0
!        enddo
!      enddo
 
      iter=0
      gama=1.0d-2
      tol=dx*dx
      err=1.0d6
8     continue
      if(err.gt.tol) then
!        call jc_p
        call diffusion_pressure(gama,krk,fac)
        err=0.0d0
        do j=1,ny
          do i=1,nx
            err=max(err,abs(p(i,j)-o(i,j)))
           enddo
        enddo
        iter=iter+1
        write(*,*)'iter=',iter,'err=',err
        goto 8
       endif

      IterSum=IterSum+iter
      endif

      call rhs(krk)
2     call pgrad(krk)

      if(move.eq.1) then
        call surface_move(krk,fac)
      endif

      do j=1,ny
        do i=1,nx-1
          u(i,j)=un(i,j)+fac*dt*urk(i,j,2)
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          v(i,j)=vn(i,j)+fac*dt*vrk(i,j,2)
        enddo
      enddo
      call velocity_reset
      call ubc(krk,fac)
      call vbc(krk,fac)

c  krk=3:
      krk=3
      fac=1.0d0

      call surface_property
      call euler_link
      if(isingular.eq.1) then
        call singular_call
      else
        call field_interpolate
        call uv_strain
      endif

      if(isingular.eq.1) then
        call jc_pressure
        call correction_pressure
      endif
      call pressure(krk,fac)

!       do j=1,ny
!        do i=1,nx
!          p(i,j)=0.0d0
!        enddo
!      enddo
 
      if(iIter.eq.1) then
      iter=0
      gama=1.0d-2
      tol=dx*dx
      err=1.0d6
7     continue
      if(err.gt.tol) then
!        call jc_p
        call diffusion_pressure(gama,krk,fac)
        err=0.0d0
        do j=1,ny
          do i=1,nx
            err=max(err,abs(p(i,j)-o(i,j)))
           enddo
        enddo
        iter=iter+1
        write(*,*)'iter=',iter,'err=',err
        goto 7
      endif

      IterSum=IterSum+iter
      endif


      call rhs(krk)
3     call pgrad(krk)

      if(move.eq.1) then
        call surface_move(krk,fac)
      endif

      do j=1,ny
        do i=1,nx-1
          u(i,j)=un(i,j)+fac*dt*urk(i,j,3)
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          v(i,j)=vn(i,j)+fac*dt*vrk(i,j,3)
        enddo
      enddo
      call velocity_reset
      call ubc(krk,fac)
      call vbc(krk,fac)

c  krk=4:
      krk=4
      fac=1.0d0

      call surface_property
      call euler_link
      if(isingular.eq.1) then
        call singular_call
      else
        call field_interpolate
        call uv_strain
      endif

      if(isingular.eq.1) then
        call jc_pressure
        call correction_pressure
      endif
      call pressure(krk,fac)

!      do j=1,ny
!        do i=1,nx
!          p(i,j)=0.0d0
!        enddo
!      enddo
 
      if(iIter.eq.1) then
      iter=0
      gama=1.0d-2
      tol=dx*dx
      err=1.0d6
6     continue
      if(err.gt.tol) then
!        call jc_p
        call diffusion_pressure(gama,krk,fac)
        err=0.0d0
        do j=1,ny
          do i=1,nx
            err=max(err,abs(p(i,j)-o(i,j)))
           enddo
        enddo
        iter=iter+1
        write(*,*)'iter=',iter,'err=',err
        goto 6
      endif
      IterSum=IterSum+iter
      write(*,*)'IterSum=',IterSum
      endif

      call rhs(krk)
4     call pgrad(krk)

      if(move.eq.1) then
        call surface_move(krk,fac)
      endif

      do j=1,ny
        do i=1,nx-1
          u(i,j)=un(i,j)+(dt/6.0d0)*(urk(i,j,1)+
     .                   2.0d0*(urk(i,j,2)+urk(i,j,3))+urk(i,j,4))
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          v(i,j)=vn(i,j)+(dt/6.0d0)*(vrk(i,j,1)+
     .                   2.0d0*(vrk(i,j,2)+vrk(i,j,3))+vrk(i,j,4))
        enddo
      enddo
      call velocity_reset
      call ubc(krk,fac)
      call vbc(krk,fac)

      ! compute drag & lift
      return
      end


c-----------------------------------------------------------------------
