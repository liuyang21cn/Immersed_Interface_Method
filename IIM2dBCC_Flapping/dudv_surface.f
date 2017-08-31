c-----------------------------------------------------------------------
c
      subroutine dudv_surface
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'rhsp.inc'
      integer is,js
      real*8 uu,vu,ud,vd,ur,vr,ul,vl
      real*8 signx,signy,sx,sy,hp,hm
      real*8 duxjc,duyjc,dduxjc,dduyjc,dduyxjc,dduxyjc
      real*8 dvxjc,dvyjc,ddvxjc,ddvyjc,ddvyxjc,ddvxyjc
      real*8 duu,dud,dvu,dvd
      real*8 duyp,duym,dvyp,dvym,duxp,duxm,dvxp,dvxm

      DO l=1,ms

      do i=0,ncxe(l)

        duyjc=ujcxe(2,i,l)
        dvyjc=vjcxe(2,i,l)
        dduyjc=ujcxe(4,i,l)
        ddvyjc=vjcxe(4,i,l)
        sy=falfaxe(1,i,l)        

        if(jexe(i,l).eq.jcxe(i,l)) then
          is=ixe(i,l)
          js=jexe(i,l)
          uu=uee(is,js+1)
          vu=vee(is,js+1)
          ud=uee(is,js)
          vd=vee(is,js)
          duu=(uu-ud)*dy1-dy1*(duyjc*(ye(js)-sy)+
     .                  0.5d0*dduyjc*(ye(js)-sy)**2.0d0)
          dvu=(vu-vd)*dy1-dy1*(dvyjc*(ye(js)-sy)+
     .                  0.5d0*ddvyjc*(ye(js)-sy)**2.0d0)

          js=jcxe(i,l)+1
          uu=u(is,js)
          vu=vec(is,js)
          ud=u(is,js-1)
          vd=vec(is,js-1)
          dud=(uu-ud)*dy1-dy1*(duyjc*(yc(js)-sy)+
     .                  0.5d0*dduyjc*(yc(js)-sy)**2.0d0)
          dvd=(vu-vd)*dy1-dy1*(dvyjc*(yc(js)-sy)+
     .                  0.5d0*ddvyjc*(yc(js)-sy)**2.0d0)

          hp=yc(js)-sy
          hm=hdy-hp

          duyp=(hp*dud+hm*duu)/hdy+
     .          hp*duyjc/hdy-hp*hm*dduyjc/hdy
          duym=(hp*dud+hm*duu)/hdy-
     .          hm*duyjc/hdy-hp*hm*dduyjc/hdy
          dvyp=(hp*dvd+hm*dvu)/hdy+
     .          hp*dvyjc/hdy-hp*hm*ddvyjc/hdy
          dvym=(hp*dvd+hm*dvu)/hdy-
     .          hm*dvyjc/hdy-hp*hm*ddvyjc/hdy
        else
          is=ixe(i,l)
          js=jcxe(i,l)
          uu=u(is,js+1)
          vu=vec(is,js+1)
          ud=u(is,js)
          vd=vec(is,js)
          duu=(uu-ud)*dy1-dy1*(duyjc*(yc(js)-sy)+
     .                  0.5d0*dduyjc*(yc(js)-sy)**2.0d0)
          dvu=(vu-vd)*dy1-dy1*(dvyjc*(yc(js)-sy)+
     .                  0.5d0*ddvyjc*(yc(js)-sy)**2.0d0)

          js=jexe(i,l)+1
          uu=uee(is,js)
          vu=vee(is,js)
          ud=uee(is,js-1)
          vd=vee(is,js-1)
          dud=(uu-ud)*dy1-dy1*(duyjc*(ye(js)-sy)+
     .                  0.5d0*dduyjc*(ye(js)-sy)**2.0d0)
          dvd=(vu-vd)*dy1-dy1*(dvyjc*(ye(js)-sy)+
     .                  0.5d0*ddvyjc*(ye(js)-sy)**2.0d0)

          hp=ye(js)-sy
          hm=hdy-hp

          duyp=(hp*dud+hm*duu)/hdy+
     .          hp*duyjc/hdy-hp*hm*dduyjc/hdy
          duym=(hp*dud+hm*duu)/hdy-
     .          hm*duyjc/hdy-hp*hm*dduyjc/hdy
          dvyp=(hp*dvd+hm*dvu)/hdy+
     .          hp*dvyjc/hdy-hp*hm*ddvyjc/hdy
          dvym=(hp*dvd+hm*dvu)/hdy-
     .          hm*dvyjc/hdy-hp*hm*ddvyjc/hdy
        endif

        falfaxe(4,i,l)=duyp
        falfaxe(5,i,l)=duym
        falfaxe(6,i,l)=dvyp
        falfaxe(7,i,l)=dvym
      enddo

      do j=0,ncye(l)

        duxjc=ujcye(1,j,l)
        dvxjc=vjcye(1,j,l)
        dduxjc=ujcye(3,j,l)
        ddvxjc=vjcye(3,j,l)
        sx=falfaye(1,j,l)

        if(ieye(j,l).eq.icye(j,l)) then
          is=ieye(j,l)
          js=jye(j,l)
          ur=uee(is+1,js)
          vr=vee(is+1,js)
          ul=uee(is,js)
          vl=vee(is,js)
          duu=(ur-ul)*dx1-dx1*(duxjc*(xe(is)-sx)+
     .                  0.5d0*dduxjc*(xe(is)-sx)**2.0d0)
          dvu=(vr-vl)*dx1-dx1*(dvxjc*(xe(is)-sx)+
     .                  0.5d0*ddvxjc*(xe(is)-sx)**2.0d0)

          is=icye(j,l)+1
          ur=uce(is,js)
          vr=v(is,js)
          ul=uce(is-1,js)
          vl=v(is-1,js)
          dud=(ur-ul)*dx1-dx1*(duxjc*(xc(is)-sx)+
     .                  0.5d0*dduxjc*(xc(is)-sx)**2.0d0)
          dvd=(vr-vl)*dx1-dx1*(dvxjc*(xc(is)-sx)+
     .                  0.5d0*ddvxjc*(xc(is)-sx)**2.0d0)

          hp=xc(is)-sx
          hm=hdx-hp

          duxp=(hp*dud+hm*duu)/hdx+
     .          hp*duxjc/hdx-hp*hm*dduxjc/hdx
          duxm=(hp*dud+hm*duu)/hdx-
     .          hm*duxjc/hdx-hp*hm*dduxjc/hdx
          dvxp=(hp*dvd+hm*dvu)/hdx+
     .          hp*dvxjc/hdx-hp*hm*ddvxjc/hdx
          dvxm=(hp*dvd+hm*dvu)/hdx-
     .          hm*dvxjc/hdx-hp*hm*ddvxjc/hdx
        else
          is=icye(j,l)
          js=jye(j,l)
          ur=uce(is+1,js)
          vr=v(is+1,js)
          ul=uce(is,js)
          vl=v(is,js)
          duu=(ur-ul)*dx1-dx1*(duxjc*(xc(is)-sx)+
     .                  0.5d0*dduxjc*(xc(is)-sx)**2.0d0)
          dvu=(vr-vl)*dx1-dx1*(dvxjc*(xc(is)-sx)+
     .                  0.5d0*ddvxjc*(xc(is)-sx)**2.0d0)

          is=ieye(j,l)+1
          ur=uee(is,js)
          vr=vee(is,js)
          ul=uee(is-1,js)
          vl=vee(is-1,js)
          dud=(ur-ul)*dx1-dx1*(duxjc*(xe(is)-sx)+
     .                  0.5d0*dduxjc*(xe(is)-sx)**2.0d0)
          dvd=(vr-vl)*dx1-dx1*(dvxjc*(xe(is)-sx)+
     .                  0.5d0*ddvxjc*(xe(is)-sx)**2.0d0)

          hp=xe(is)-sx
          hm=hdx-hp

          duxp=(hp*dud+hm*duu)/hdx+
     .          hp*duxjc/hdx-hp*hm*dduxjc/hdx
          duxm=(hp*dud+hm*duu)/hdx-
     .          hm*duxjc/hdx-hp*hm*dduxjc/hdx
          dvxp=(hp*dvd+hm*dvu)/hdx+
     .          hp*dvxjc/hdx-hp*hm*ddvxjc/hdx
          dvxm=(hp*dvd+hm*dvu)/hdx-
     .          hm*dvxjc/hdx-hp*hm*ddvxjc/hdx
        endif

        falfaye(4,j,l)=duxp
        falfaye(5,j,l)=duxm
        falfaye(6,j,l)=dvxp
        falfaye(7,j,l)=dvxm
      enddo

      do i=0,ncxc(l)

        signx=falfaxc(2,i,l)
        signy=falfaxc(3,i,l)
        sy=falfaxc(1,i,l)

        duxjc=signy*signx*ujcxc(1,i,l)
        dvxjc=signy*signx*vjcxc(1,i,l)
        dduxjc=signy*signx*ujcxc(3,i,l)
        ddvxjc=signy*signx*vjcxc(3,i,l)
        dduyxjc=signy*signx*ujcxc(5,i,l)
        ddvyxjc=signy*signx*vjcxc(5,i,l)

        duyjc=ujcxc(2,i,l)
        dvyjc=vjcxc(2,i,l)
        dduyjc=ujcxc(4,i,l)
        ddvyjc=vjcxc(4,i,l)
        dduxyjc=ujcxc(6,i,l)
        ddvxyjc=vjcxc(6,i,l)

        if(jexc(i,l).eq.jcxc(i,l)) then
          is=ixc(i,l)
          js=jexc(i,l)
          uu=uce(is,js+1)
          vu=v(is,js+1)
          ud=uce(is,js)
          vd=v(is,js)
          duu=(uu-ud)*dy1-dy1*(duyjc*(ye(js)-sy)+
     .                  0.5d0*dduyjc*(ye(js)-sy)**2.0d0)
          dvu=(vu-vd)*dy1-dy1*(dvyjc*(ye(js)-sy)+
     .                  0.5d0*ddvyjc*(ye(js)-sy)**2.0d0)

          js=jcxc(i,l)+1
          uu=ucc(is,js)
          vu=vcc(is,js)
          ud=ucc(is,js-1)
          vd=vcc(is,js-1)
          dud=(uu-ud)*dy1-dy1*(duyjc*(yc(js)-sy)+
     .                  0.5d0*dduyjc*(yc(js)-sy)**2.0d0)
          dvd=(vu-vd)*dy1-dy1*(dvyjc*(yc(js)-sy)+
     .                  0.5d0*ddvyjc*(yc(js)-sy)**2.0d0)

          hp=yc(js)-sy
          hm=hdy-hp

          duyp=(hp*dud+hm*duu)/hdy+
     .          hp*duyjc/hdy-hp*hm*dduyjc/hdy
          duym=(hp*dud+hm*duu)/hdy-
     .          hm*duyjc/hdy-hp*hm*dduyjc/hdy
          dvyp=(hp*dvd+hm*dvu)/hdy+
     .          hp*dvyjc/hdy-hp*hm*ddvyjc/hdy
          dvym=(hp*dvd+hm*dvu)/hdy-
     .          hm*dvyjc/hdy-hp*hm*ddvyjc/hdy          
        else
          is=ixc(i,l)
          js=jcxc(i,l)
          uu=ucc(is,js+1)
          vu=vcc(is,js+1)
          ud=ucc(is,js)
          vd=vcc(is,js)
          duu=(uu-ud)*dy1-dy1*(duyjc*(yc(js)-sy)+
     .                  0.5d0*dduyjc*(yc(js)-sy)**2.0d0)
          dvu=(vu-vd)*dy1-dy1*(dvyjc*(yc(js)-sy)+
     .                  0.5d0*ddvyjc*(yc(js)-sy)**2.0d0)

          js=jexc(i,l)+1
          uu=uce(is,js)
          vu=v(is,js)
          ud=uce(is,js-1)
          vd=v(is,js-1)
          dud=(uu-ud)*dy1-dy1*(duyjc*(ye(js)-sy)+
     .                  0.5d0*dduyjc*(ye(js)-sy)**2.0d0)
          dvd=(vu-vd)*dy1-dy1*(dvyjc*(ye(js)-sy)+
     .                  0.5d0*ddvyjc*(ye(js)-sy)**2.0d0)

          hp=ye(js)-sy
          hm=hdy-hp

          duyp=(hp*dud+hm*duu)/hdy+
     .          hp*duyjc/hdy-hp*hm*dduyjc/hdy
          duym=(hp*dud+hm*duu)/hdy-
     .          hm*duyjc/hdy-hp*hm*dduyjc/hdy
          dvyp=(hp*dvd+hm*dvu)/hdy+
     .          hp*dvyjc/hdy-hp*hm*ddvyjc/hdy
          dvym=(hp*dvd+hm*dvu)/hdy-
     .          hm*dvyjc/hdy-hp*hm*ddvyjc/hdy
        endif

        is=ixc(i,l)
        js=jcxc(i,l)+1
        hp=yc(js)-sy
        hm=dy-hp
        duu=dux(is,js)
        dvu=dvx(is,js)
        dud=dux(is,js-1)
        dvd=dvx(is,js-1)
        duxp=(hp*dud+hm*duu)*dy1+
     .        hp*duxjc*dy1-hp*hm*dduxyjc*dy1
        duxm=(hp*dud+hm*duu)*dy1-
     .        hm*duxjc*dy1-hp*hm*dduxyjc*dy1
        dvxp=(hp*dvd+hm*dvu)*dy1+
     .        hp*dvxjc*dy1-hp*hm*ddvxyjc*dy1
        dvxm=(hp*dvd+hm*dvu)*dy1-
     .        hm*dvxjc*dy1-hp*hm*ddvxyjc*dy1

        falfaxc(4,i,l)=duxp
        falfaxc(5,i,l)=duxm
        falfaxc(6,i,l)=dvxp
        falfaxc(7,i,l)=dvxm
        falfaxc(8,i,l)=duyp
        falfaxc(9,i,l)=duym
        falfaxc(10,i,l)=dvyp
        falfaxc(11,i,l)=dvym
      enddo

      do j=0,ncyc(l)

        signx=falfayc(2,j,l)
        signy=falfayc(3,j,l)
        sx=falfayc(1,j,l)

        duxjc=ujcyc(1,j,l)
        dvxjc=vjcyc(1,j,l)
        dduxjc=ujcyc(3,j,l)
        ddvxjc=vjcyc(3,j,l)
        dduyxjc=ujcyc(5,j,l)
        ddvyxjc=vjcyc(5,j,l)

        duyjc=signx*signy*ujcyc(2,j,l)
        dvyjc=signx*signy*vjcyc(2,j,l)
        dduyjc=signx*signy*ujcyc(4,j,l)
        ddvyjc=signx*signy*vjcyc(4,j,l)
        dduxyjc=signx*signy*ujcyc(6,j,l)
        ddvxyjc=signx*signy*vjcyc(6,j,l)

        if(ieyc(j,l).eq.icyc(j,l)) then
          is=ieyc(j,l)
          js=jyc(j,l)
          ur=u(is+1,js)
          vr=vec(is+1,js)
          ul=u(is,js)
          vl=vec(is,js)
          duu=(ur-ul)*dx1-dx1*(duxjc*(xe(is)-sx)+
     .                  0.5d0*dduxjc*(xe(is)-sx)**2.0d0)
          dvu=(vr-vl)*dx1-dx1*(dvxjc*(xe(is)-sx)+
     .                  0.5d0*ddvxjc*(xe(is)-sx)**2.0d0)

          is=icyc(j,l)+1
          ur=ucc(is,js)
          vr=vcc(is,js)
          ul=ucc(is-1,js)
          vl=vcc(is-1,js)
          dud=(ur-ul)*dx1-dx1*(duxjc*(xc(is)-sx)+
     .                  0.5d0*dduxjc*(xc(is)-sx)**2.0d0)
          dvd=(vr-vl)*dx1-dx1*(dvxjc*(xc(is)-sx)+
     .                  0.5d0*ddvxjc*(xc(is)-sx)**2.0d0)

          hp=xc(is)-sx
          hm=hdx-hp

          duxp=(hp*dud+hm*duu)/hdx+
     .          hp*duxjc/hdx-hp*hm*dduxjc/hdx
          duxm=(hp*dud+hm*duu)/hdx-
     .          hm*duxjc/hdx-hp*hm*dduxjc/hdx
          dvxp=(hp*dvd+hm*dvu)/hdx+
     .          hp*dvxjc/hdx-hp*hm*ddvxjc/hdx
          dvxm=(hp*dvd+hm*dvu)/hdx-
     .          hm*dvxjc/hdx-hp*hm*ddvxjc/hdx
        else
          is=icyc(j,l)
          js=jyc(j,l)
          ur=ucc(is+1,js)
          vr=vcc(is+1,js)
          ul=ucc(is,js)
          vl=vcc(is,js)
          duu=(ur-ul)*dx1-dx1*(duxjc*(xc(is)-sx)+
     .                  0.5d0*dduxjc*(xc(is)-sx)**2.0d0)
          dvu=(vr-vl)*dx1-dx1*(dvxjc*(xc(is)-sx)+
     .                  0.5d0*ddvxjc*(xc(is)-sx)**2.0d0)

          is=ieyc(j,l)+1
          ur=u(is,js)
          vr=vec(is,js)
          ul=u(is-1,js)
          vl=vec(is-1,js)
          dud=(ur-ul)*dx1-dx1*(duxjc*(xe(is)-sx)+
     .                  0.5d0*dduxjc*(xe(is)-sx)**2.0d0)
          dvd=(vr-vl)*dx1-dx1*(dvxjc*(xe(is)-sx)+
     .                  0.5d0*ddvxjc*(xe(is)-sx)**2.0d0)

          hp=xe(is)-sx
          hm=hdx-hp

          duxp=(hp*dud+hm*duu)/hdx+
     .          hp*duxjc/hdx-hp*hm*dduxjc/hdx
          duxm=(hp*dud+hm*duu)/hdx-
     .          hm*duxjc/hdx-hp*hm*dduxjc/hdx
          dvxp=(hp*dvd+hm*dvu)/hdx+
     .          hp*dvxjc/hdx-hp*hm*ddvxjc/hdx
          dvxm=(hp*dvd+hm*dvu)/hdx-
     .          hm*dvxjc/hdx-hp*hm*ddvxjc/hdx
        endif

        js=jyc(j,l)
        is=icyc(j,l)+1
        hp=xc(is)-sx
        hm=dx-hp
        duu=duy(is,js)
        dvu=dvy(is,js)
        dud=duy(is-1,js)
        dvd=dvy(is-1,js)
        duyp=(hp*dud+hm*duu)*dx1+
     .        hp*duyjc*dx1-hp*hm*dduyxjc*dx1
        duym=(hp*dud+hm*duu)*dx1-
     .        hm*duyjc*dx1-hp*hm*dduyxjc*dx1
        dvyp=(hp*dvd+hm*dvu)*dx1+
     .        hp*dvyjc*dx1-hp*hm*ddvyxjc*dx1
        dvym=(hp*dvd+hm*dvu)*dx1-
     .        hm*dvyjc*dx1-hp*hm*ddvyxjc*dx1

        falfayc(4,j,l)=duxp
        falfayc(5,j,l)=duxm
        falfayc(6,j,l)=dvxp
        falfayc(7,j,l)=dvxm
        falfayc(8,j,l)=duyp
        falfayc(9,j,l)=duym
        falfayc(10,j,l)=dvyp
        falfayc(11,j,l)=dvym
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
