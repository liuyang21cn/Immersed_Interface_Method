c-----------------------------------------------------------------------
c
      subroutine rhs(krk)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'old.inc'
      integer krk
      real*8 conv,grad,visc

      do j=1,ny
        do i=1,nx-1
          conv=-dx1*(ucc(i+1,j)*ucc(i+1,j)-ucc(i,j)*ucc(i,j))
     .         -dy1*(uee(i,j)*vee(i,j)-uee(i,j-1)*vee(i,j-1))
          grad=0.0d0
          visc=(dx2*Re1)*(u(i+1,j)-2.0d0*u(i,j)+u(i-1,j))+
     .         (dy2*Re1)*(u(i,j+1)-2.0d0*u(i,j)+u(i,j-1))
          urk(i,j,krk)=conv+grad+visc
        enddo
      enddo

      if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxe(l)
            urk(ixe(i,l),jexe(i,l)+1,krk)=urk(ixe(i,l),jexe(i,l)+1,krk)-
     .                                    ucdy(i,l)
            urk(ixe(i,l),jcxe(i,l),krk)=urk(ixe(i,l),jcxe(i,l),krk)+
     .                                  ucdyy(1,i,l)*Re1
            urk(ixe(i,l),jcxe(i,l)+1,krk)=urk(ixe(i,l),jcxe(i,l)+1,krk)+
     .                                    ucdyy(2,i,l)*Re1
          enddo
          do j=0,ncyc(l)
            urk(icyc(j,l),jyc(j,l),krk)=urk(icyc(j,l),jyc(j,l),krk)-
     .                                  ucudx(j,l)-ucpdx(j,l)
            urk(ieyc(j,l),jyc(j,l),krk)=urk(ieyc(j,l),jyc(j,l),krk)+
     .                                  ucdxx(1,j,l)*Re1
            urk(ieyc(j,l)+1,jyc(j,l),krk)=urk(ieyc(j,l)+1,jyc(j,l),krk)+
     .                                    ucdxx(2,j,l)*Re1
          enddo
        enddo
      endif
 
      do j=1,ny-1
        do i=1,nx
          conv=-dy1*(vcc(i,j+1)*vcc(i,j+1)-vcc(i,j)*vcc(i,j))
     .         -dx1*(uee(i,j)*vee(i,j)-uee(i-1,j)*vee(i-1,j))
          grad=0.0d0
          visc=(dx2*Re1)*(v(i+1,j)-2.0d0*v(i,j)+v(i-1,j))+
     .         (dy2*Re1)*(v(i,j+1)-2.0d0*v(i,j)+v(i,j-1))
          vrk(i,j,krk)=conv+grad+visc
        enddo
      enddo

      if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxc(l)
            vrk(ixc(i,l),jcxc(i,l),krk)=vrk(ixc(i,l),jcxc(i,l),krk)-
     .                                  vcvdy(i,l)-vcpdy(i,l)
            vrk(ixc(i,l),jexc(i,l),krk)=vrk(ixc(i,l),jexc(i,l),krk)+
     .                                  vcdyy(1,i,l)*Re1
            vrk(ixc(i,l),jexc(i,l)+1,krk)=vrk(ixc(i,l),jexc(i,l)+1,krk)+
     .                                    vcdyy(2,i,l)*Re1
          enddo
          do j=0,ncye(l)
            vrk(ieye(j,l)+1,jye(j,l),krk)=vrk(ieye(j,l)+1,jye(j,l),krk)-
     .                                    vcdx(j,l)
            vrk(icye(j,l),jye(j,l),krk)=vrk(icye(j,l),jye(j,l),krk)+
     .                                  vcdxx(1,j,l)*Re1
            vrk(icye(j,l)+1,jye(j,l),krk)=vrk(icye(j,l)+1,jye(j,l),krk)+
     .                                    vcdxx(2,j,l)*Re1
          enddo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
