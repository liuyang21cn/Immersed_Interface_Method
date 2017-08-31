c-----------------------------------------------------------------------
c
      subroutine uv_strain
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'rhsp.inc'

      do j=1,ny
        do i=1,nx
          dux(i,j)=dx1*(u(i,j)-u(i-1,j))
          duy(i,j)=dy1*(uce(i,j)-uce(i,j-1))
          dvx(i,j)=dx1*(vec(i,j)-vec(i-1,j))
          dvy(i,j)=dy1*(v(i,j)-v(i,j-1))
        enddo
      enddo

      if(isingular.eq.1) then
        do l=1,ms
          do i=0,ncxc(l)
            duy(ixc(i,l),jexc(i,l)+1)=duy(ixc(i,l),jexc(i,l)+1)+
     .                                pcudy(i,l)
            dvy(ixc(i,l),jexc(i,l)+1)=dvy(ixc(i,l),jexc(i,l)+1)+
     .                                pcvdy(i,l)
          enddo
          do j=0,ncyc(l)
            dux(ieyc(j,l)+1,jyc(j,l))=dux(ieyc(j,l)+1,jyc(j,l))+
     .                                pcudx(j,l)
            dvx(ieyc(j,l)+1,jyc(j,l))=dvx(ieyc(j,l)+1,jyc(j,l))+
     .                                pcvdx(j,l)
          enddo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
