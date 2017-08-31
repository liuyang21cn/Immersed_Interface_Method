c-----------------------------------------------------------------------
c
      subroutine old_save
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'rhsp.inc'
      include 'old.inc'

      do j=0,ny+1
        do i=0,nx
          un(i,j)=u(i,j)
        enddo
      enddo

      do j=0,ny
        do i=0,nx+1
          vn(i,j)=v(i,j)
        enddo
      enddo

      do j=1,ny
        do i=1,nx
          dn(i,j)=dux(i,j)+dvy(i,j)
        enddo
      enddo

      DO l=1,ms      

      do m=0,ns
        xsn(m,l)=xs(m,l)
        ysn(m,l)=ys(m,l)
      enddo

      ncxen(l)=ncxe(l)
      ncxcn(l)=ncxc(l)
      ncyen(l)=ncye(l)
      ncycn(l)=ncyc(l)

      do i=0,ncxe(l)
        ixen(i,l)=ixe(i,l)
        jexen(i,l)=jexe(i,l)
        jcxen(i,l)=jcxe(i,l)
        falfaxen(1,i,l)=falfaxe(1,i,l)
        falfaxen(2,i,l)=falfaxe(13,i,l)
        falfaxen(3,i,l)=falfaxe(14,i,l)
        do n=1,2
          ujcxen(n,i,l)=ujcxe(n,i,l)
          vjcxen(n,i,l)=vjcxe(n,i,l)
        enddo
      enddo

      do i=0,ncxc(l)
        ixcn(i,l)=ixc(i,l)
        jcxcn(i,l)=jcxc(i,l)
        jexcn(i,l)=jexc(i,l)
        falfaxcn(1,i,l)=falfaxc(1,i,l)
        falfaxcn(2,i,l)=falfaxc(13,i,l)
        falfaxcn(3,i,l)=falfaxc(14,i,l)
        do n=1,2
          ujcxcn(n,i,l)=ujcxc(n,i,l)
          vjcxcn(n,i,l)=vjcxc(n,i,l)
        enddo
      enddo

      do j=0,ncye(l)
        jyen(j,l)=jye(j,l)
        ieyen(j,l)=ieye(j,l)
        icyen(j,l)=icye(j,l)
        falfayen(1,j,l)=falfaye(1,j,l)
        falfayen(2,j,l)=falfaye(13,j,l)
        falfayen(3,j,l)=falfaye(14,j,l)
        do n=1,2
          ujcyen(n,j,l)=ujcye(n,j,l)
          vjcyen(n,j,l)=vjcye(n,j,l)
        enddo
      enddo

      do j=0,ncyc(l)
        jycn(j,l)=jyc(j,l)
        icycn(j,l)=icyc(j,l)
        ieycn(j,l)=ieyc(j,l)
        falfaycn(1,j,l)=falfayc(1,j,l)
        falfaycn(2,j,l)=falfayc(13,j,l)
        falfaycn(3,j,l)=falfayc(14,j,l)
        do n=1,2
          ujcycn(n,j,l)=ujcyc(n,j,l)
          vjcycn(n,j,l)=vjcyc(n,j,l)
        enddo
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------

