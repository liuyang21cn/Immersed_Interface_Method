c-----------------------------------------------------------------------
c
      subroutine field_interpolate
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      do j=0,ny
        do i=0,nx
          uee(i,j)=0.5d0*(u(i,j)+u(i,j+1))
          vee(i,j)=0.5d0*(v(i,j)+v(i+1,j))
        enddo
      enddo
      if(isingular.eq.1) then
        do l=1,ms
          do j=0,ncye(l)
            vee(ivxeye(j,l),jye(j,l))=vee(ivxeye(j,l),jye(j,l))+
     .                                  vcixeye(j,l)
          enddo
        enddo
      endif

      do j=0,ny+1
        do i=1,nx
          ucc(i,j)=0.5d0*(u(i-1,j)+u(i,j))
        enddo
      enddo
      if(isingular.eq.1) then
        do l=1,ms
          do j=0,ncyc(l)
            ucc(iuxcyc(j,l),jyc(j,l))=ucc(iuxcyc(j,l),jyc(j,l))+
     .                                  ucixcyc(j,l)
          enddo
        enddo
      endif

      do j=0,ny
        do i=1,nx
          uce(i,j)=0.5d0*(ucc(i,j)+ucc(i,j+1))
        enddo
      enddo

      do j=1,ny
        do i=0,nx+1
          vcc(i,j)=0.5d0*(v(i,j-1)+v(i,j))
        enddo
      enddo

      do j=1,ny
        do i=0,nx
          vec(i,j)=0.5d0*(vee(i,j-1)+vee(i,j))
        enddo
      enddo

      if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxe(l)
            uee(ixe(i,l),juyexe(i,l))=uee(ixe(i,l),juyexe(i,l))+
     .                                  uciyexe(i,l)
            vec(ixe(i,l),jvycxe(i,l))=vec(ixe(i,l),jvycxe(i,l))+
     .                                  vciycxe(i,l)
          enddo
          do i=0,ncxc(l)
            uce(ixc(i,l),juyexc(i,l))=uce(ixc(i,l),juyexc(i,l))+
     .                                  uciyexc(i,l)
            vcc(ixc(i,l),jvycxc(i,l))=vcc(ixc(i,l),jvycxc(i,l))+
     .                                  vciycxc(i,l)
          enddo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
