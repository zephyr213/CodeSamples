      subroutine forc3lj(id,it)
c
c   ***  subroutine for doing force calculation with linked cells
c
      include 'hellv6.inc'

	dimension strbx0(3,3),strbxre0(3,3),strbxsr0(3,3)
	data nxold /0/, nyold /0/, nzold /0/

	epot=0.0d0



	do i=1,3
	do j=1,3
	 strbx(i,j)=0.0d0
	 strbxre(i,j)=0.0d0
	 strbxsr(i,j)=0.0d0
	enddo
	enddo
c

      do 10 n=1,natomz
c
c      if (mod(n,2).eq.0) sx(n)=sx(n)+0.001

       if (sx(n).ge. 0.5d0) then
            sx(n)=sx(n)-1.d0*ipbcx
            rx0(n)=rx0(n)-1.d0*ipbcx*h11
       else if (sx(n).lt.-0.5d0) then
            sx(n)=sx(n)+1.d0*ipbcx
            rx0(n)=rx0(n)+1.d0*ipbcx*h11
       end if
       if (sy(n).ge. 0.5d0) then
            sy(n)=sy(n)-1.d0*ipbcy
            ry0(n)=ry0(n)-1.d0*ipbcy*h22
       else if (sy(n).lt.-0.5d0) then
            sy(n)=sy(n)+1.d0*ipbcy
            ry0(n)=ry0(n)+1.d0*ipbcy*h22
       end if
       if (sz(n).ge. 0.5d0) then
           sz(n)=sz(n)-1.d0*ipbcz
           if(it.eq.1) rz0(n)=rz0(n)-1.d0*ipbcz*h33
       else if (sz(n).lt.-0.5d0) then
            sz(n)=sz(n)+1.d0*ipbcz
            if (it.eq.1) rz0(n)=rz0(n)+1.d0*ipbcz*h33
       end if

c      ntop(n)=0
c      nplus(n)=0
c      nminus(n)=0
       stress(n,1)=0.0d0
       stress(n,2)=0.0d0
       stress(n,3)=0.0d0
       stress(n,4)=0.0d0
       stress(n,5)=0.0d0
       stress(n,6)=0.0d0
       str1_sr(1,n)=0d0
       str2_sr(1,n)=0d0
       str3_sr(1,n)=0d0
       str4_sr(1,n)=0d0
       str5_sr(1,n)=0d0
       str6_sr(1,n)=0d0
       str1_re(1,n)=0d0
       str2_re(1,n)=0d0
       str3_re(1,n)=0d0
       str4_re(1,n)=0d0
       str5_re(1,n)=0d0
       str6_re(1,n)=0d0
c
10    continue
c
c
c call cell setup if system has changed shape enough to require
c greater or fewer number of cells
c
	
      size=rsr(1)
      if (rsrmax.gt.size) size=rsrmax
c
c call cell setup if system has changed shape enough to require
c greater or fewer number of cells
c
      rrx=dsqrt(h11*h11+h21*h21+h31*h31)
      rry=dsqrt(h12*h12+h22*h22+h32*h32)
      rrz=dsqrt(h13*h13+h23*h23+h33*h33)
      nnz=nproc*int(rrz/size/dfloat(nproc))

      nnx=int(rrx/size)
      nny=int(rry/size)
      if (nnx.lt.3) nnx=3
      if (nny.lt.3) nny=3
      cell_size_X=rrx/nnx
      cell_size_Y=rry/nny

      if(nnz.lt.2) goto 1313
      if(rrx.le.(2.0*size)) goto 1313
      if(rry.le.(2.0*size)) goto 1313
      cell_size_Z=rrz/nnz
      if(cell_size_Z.lt.size) goto 1313
      goto 777

1313  continue     !  Report an error and exit...
      if(mynod.eq.0) then
        write(6,*)'Box size ERROR!'
        write(6,*)'Change number of nodes or system size'
        write(6,*)'nnx=',nnx,' nny=',nny,' nnz=',nnz,' ncell=',ncell
        write(6,*)'cell_size_X=',cell_size_X,' Y=',cell_size_Y,
     1' Z=',cell_size_Z,' Inter.dist.=',size
      endif
        call MPI_FINALIZE(ierr)
        stop

777   continue    ! safe to continue...


      nflag=abs(nnx-nxold)+abs(nny-nyold)+abs(nnz-nzold)

	t1=secnds(0.0)
      if (nflag.gt.0)  call cell_setup_ion
	t2=secnds(0.0)

c	write(6,*) 'mynod,cellsetup=',mynod,t2-t1

      nxold=nnx
      nyold=nny
      nzold=nnz
c
c

c even though this is stillinger-weber, same finder routine can
c be d

	t1=secnds(0.0)
        call finder_tersoff(it)
	t2=secnds(0.0)
c	write(6,*) 'mynod,finder time=',mynod,t2-t1


	do i=1,3
	do j=1,3
	 strbx(i,j)=0.0d0
	 strbxre(i,j)=0.0d0
	 strbxsr(i,j)=0.0d0
	enddo
	enddo
c

      do  n=1,natomz
	      epot1(n)=0.0d0
	      epot2(n)=0.0d0
       ere(1,n)=0d0
       esr(1,n)=0d0
       fxs(n)=0.0d0
       fys(n)=0.0d0
       fzs(n)=0.0d0
       fxre(1,n)=0d0
       fyre(1,n)=0d0
       fzre(1,n)=0d0
      enddo

c

	t1=secnds(0.0)
	call lj_use(it,negtot)
	t2=secnds(0.0)

c	write(6,*) 'mynod,stillweb=',mynod,t2-t1

	etott1=0.0d0
	etott2=0.0d0
	etott3=0.0d0
	etott4=0.0d0
	do n=1,natomz
		vijx=h(1,1)*x1(n)+h(1,2)*y1(n)+h(1,3)*z1(n)
		vijy=h(2,1)*x1(n)+h(2,2)*y1(n)+h(2,3)*z1(n)
		vijz=h(3,1)*x1(n)+h(3,2)*y1(n)+h(3,3)*z1(n)

		vijx=vijx/dt
		vijy=vijy/dt
		vijz=vijz/dt

	      if(rz0(n).lt.0.0d0) then
		qtotx1=qtotx1+vijx*epot1(n)*e00
		qtoty1=qtoty1+vijy*epot1(n)*e00
		qtotz1=qtotz1+vijz*epot1(n)*e00

		qtotx3=qtotx3+vijx*epot2(n)*e00
		qtoty3=qtoty3+vijy*epot2(n)*e00
		qtotz3=qtotz3+vijz*epot2(n)*e00

		etott1=etott1+epot1(n)
		etott3=etott3+epot2(n)
	      endif

	      if(rz0(n).ge.0.0d0) then
		qtotx2=qtotx2+vijx*epot1(n)*e00
		qtoty2=qtoty2+vijy*epot1(n)*e00
		qtotz2=qtotz2+vijz*epot1(n)*e00

		qtotx4=qtotx4+vijx*epot2(n)*e00
		qtoty4=qtoty4+vijy*epot2(n)*e00
		qtotz4=qtotz4+vijz*epot2(n)*e00

		etott2=etott2+epot1(n)
		etott4=etott4+epot2(n)
	      endif

		if(n.le.natoms) then
			amasm=amass(imass(n))
			ekint=vijx**2+vijy**2+vijz**2
			ekint=(amasm/2.0d0)*ekint
			if(rz0(n).lt.0.0d0) then
			 qtotx1=qtotx1+ekint*vijx*e00
			 qtoty1=qtoty1+ekint*vijy*e00
			 qtotz1=qtotz1+ekint*vijz*e00

			 qtotx3=qtotx3+ekint*vijx*e00
			 qtoty3=qtoty3+ekint*vijy*e00
			 qtotz3=qtotz3+ekint*vijz*e00

			 etott1=etott1+ekint
			 etott3=etott3+ekint
			endif

			if(rz0(n).ge.0.0d0) then
			 qtotx2=qtotx2+ekint*vijx*e00
			 qtoty2=qtoty2+ekint*vijy*e00
			 qtotz2=qtotz2+ekint*vijz*e00

			 qtotx4=qtotx4+ekint*vijx*e00
			 qtoty4=qtoty4+ekint*vijy*e00
			 qtotz4=qtotz4+ekint*vijz*e00

			 etott2=etott2+ekint
			 etott4=etott4+ekint
			endif

		endif
		enddo


c
	t1=secnds(0.0)
	call MPI_REDUCE(qtotx1,qtotx10,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(qtoty1,qtoty10,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(qtotz1,qtotz10,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)

	call MPI_REDUCE(qtotx2,qtotx20,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(qtoty2,qtoty20,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(qtotz2,qtotz20,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)

	call MPI_REDUCE(qtotx3,qtotx30,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(qtoty3,qtoty30,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(qtotz3,qtotz30,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)

	call MPI_REDUCE(qtotx4,qtotx40,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(qtoty4,qtoty40,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(qtotz4,qtotz40,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)

	call MPI_REDUCE(etott1,etott10,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(etott2,etott20,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(etott3,etott30,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(etott4,etott40,1,MPI_REAL8,MPI_SUM,
     1    0,MPI_COMM_WORLD,ierr)





	call MPI_BCAST(qtotx10,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(qtoty10,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(qtotz10,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(qtotx20,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(qtoty20,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(qtotz20,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(qtotx30,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(qtoty30,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(qtotz30,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(qtotx40,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(qtoty40,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(qtotz40,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

	call MPI_BCAST(etott10,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(etott20,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(etott30,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(etott40,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

	qtotx1=qtotx10
	qtoty1=qtoty10
	qtotz1=qtotz10

	qtotx2=qtotx20
	qtoty2=qtoty20
	qtotz2=qtotz20

	qtotx3=qtotx30
	qtoty3=qtoty30
	qtotz3=qtotz30

	qtotx4=qtotx40
	qtoty4=qtoty40
	qtotz4=qtotz40

	etott1=etott10
	etott2=etott20
	etott3=etott30
	etott4=etott40
	
	
	t2=secnds(0.0)

c	write(6,*) 'mynod,reduces=',mynod,t2-t1

	  do 202 n=1,natomz
      strbxre(1,1)=strbxre(1,1)+(stress(n,1))/dh
      strbxre(2,2)=strbxre(2,2)+(stress(n,2))/dh
      strbxre(3,3)=strbxre(3,3)+(stress(n,3))/dh
      strbxre(1,2)=strbxre(1,2)+(stress(n,4))/dh
      strbxre(1,3)=strbxre(1,3)+(stress(n,5))/dh
      strbxre(2,3)=strbxre(2,3)+(stress(n,6))/dh
      fx=fxre(1,n)
      fy=fyre(1,n)
      fz=fzre(1,n)
      fck=1d0/amass(imass(n))
      fxs(n)=fx*fck
      fys(n)=fy*fck
      fzs(n)=fz*fck

      epot=epot+ere(1,n)


      forn=dabs(fxs(n))+dabs(fys(n))+dabs(fzs(n))
      if (jpr.ge.2.and.n.le.ndet) then
      write(ifile,203) mynod,n,esr(1,n),ere(1,n),
     1  ere(1,n)+esr(1,n),imass(n)
c      write(ifile,305) mynod,n,fxt(1,n),fyt(1,n),fzt(1,n)
      endif
202   continue


c even though this is stillinger weber, same fshift routine can be
c used.

	t1=secnds(0.0)
	call fshift_tersoff(ichk)
	t2=secnds(0.0)

c	write(6,*) 'mynod,fshift time=',mynod,t2-t1
	

	strbx(1,1)=strbxre(1,1)
	strbx(2,2)=strbxre(2,2)
	strbx(3,3)=strbxre(3,3)
	strbx(1,2)=strbxre(1,2)
	strbx(1,3)=strbxre(1,3)
	strbx(2,3)=strbxre(2,3)
	strbx(2,1)=strbxre(1,2)
	strbx(3,1)=strbxre(1,3)
	strbx(3,2)=strbxre(2,3)

	call MPI_REDUCE(strbx,strbx0,9,MPI_REAL8,MPI_SUM,
     1   0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(strbx0,9,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

	do k=1,3
	do l=1,3
	 strbx(k,l)=strbx0(k,l)
	enddo
	enddo

	
  203 format(1x,i1,2x,i4,3(3x,e11.4),4(2x,i3))
  305 format(1x,i2,2x,i4,3(3x,e12.5))

c
      if (jpr.eq.5) then
      if (mynod.eq.0) then
	write(ifile,*) 'mynod=',mynod
      write(ifile,*) 'total potential energy=',epot
      write(ifile,*) 'total real energy=',ereal
      write(ifile,*) 'total short energy=',eshort
      write(ifile,*) 'total number of atoms=',nsum
      eee=epot/nsum
c     write(ifile,*) 'average energy per atom=',eee
      endif
      endif

      call MPI_REDUCE(epot,epot0,1,MPI_REAL8,MPI_SUM,
     1      0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(epot0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

	epot=epot0
	ecoh=ecoh+epot	
	ecoh=epot
	ereal0=0
	eshort0=0
		erealc=erealc+ereal0
	eshortc=eshortc+eshort0

      if (jpr.ge.1) then
      write(ifile,599) strbx
 599  format(1x,'stress-box',3(/,1x,3(e15.8,2x)))
      write(ifile,598) strbxre
 598  format(1x,'stress-box-real',3(/,1x,3(e15.8,2x)))
      write(ifile,597) strbxsr
 597  format(1x,'stress-box-short',3(/,1x,3(e15.8,2x)))
      endif
c

      return
      end
c
c-----------------------------------------------------------------------
c

      subroutine lj_use(it,negtot)
      include 'hellv6.inc'
      dimension pa(2)
      dimension epsilon(2),sigma(2)
      double precision dlj1,dlj6,ljrc,ljvc,ljvc1,ljvc2,ljvc6,ljdvr

	negtot=0
c
      pi=dacos(-1.0d0)
      twopi=2.0d0*dacos(-1.0d0)


	qtotx1=0.0d0
	qtoty1=0.0d0
	qtotz1=0.0d0

	qtotx2=0.0d0
	qtoty2=0.0d0
	qtotz2=0.0d0

	qtotx3=0.0d0
	qtoty3=0.0d0
	qtotz3=0.0d0

	qtotx4=0.0d0
	qtoty4=0.0d0
	qtotz4=0.0d0

	pa(1)=4.0d0
	ps(1)=1.49d0*3.616d0/2.314d0    ! rcut

	sigma(1)=2.314d0
	epsilon(1)=0.167d0


	pa(1)=pa(1)*epsilon(1)
	pl1(1)=pl1(1)*epsilon(1)

	rss=sigma(1)*(ps(1))/e00
	rss=rss*rss

	ljrc=ps(1)
	ljvc1=1.0d0/ljrc
	ljvc2=ljvc1*ljvc1
	ljvc6=ljvc2*ljvc2*ljvc2
	ljvc=pa(1)*(ljvc2-1.0d0)*ljvc2
        ljdvr=pa(1)*(1.0d0/sigma(1))*(-12.0d0*ljvc6+6.0d0)*ljvc1*ljvc6


	call matinv(h,hi,dh)


       do 5000 icell=ncstart1,ncend-1


c	nntotal=nntotal+nnabor
	nnabor=0
	nnpairs=0


        iz=icell/ncell_in_layer
        iy=mod(icell,ncell_in_layer)/nnx
        ix=mod(mod(icell,ncell_in_layer),nnx)
        ind=mod(iz+iZ_shift,nnz)*ncell_in_layer
     1  + mod(iy+nny,nny)*nnx + mod(ix+nnx,nnx)
        nr=ihead(ind)
c
c
c
1000   if (nr.gt.0) then

	if(nr.gt.natoms) write(6,*) 'mynod,nr,natoms=',mynod,nr,natoms
		call flush(6)

       imn=imass(nr)
c 
       sxn=sx(nr)
       syn=sy(nr)
       szn=sz(nr)



       neg=0

c first inner loop
c loop over all neighboring cells to
c obtain neighbors of ion nr.

       nabor=0; iyn=iy+nny; ixn=ix+nnx
       do 4000 kk=-1,1
        kzn=mod(iz+kk+iZ_shift,nnz)*ncell_in_layer
       do 4000 jj=-1,1
        jyn= kzn + mod(iyn+jj,nny)*nnx
       do 4000 ii=-1,1
        jcell= jyn + mod(ixn+ii,nnx)
        l=ihead(jcell)
c        nabor=nabor+1

3000     if (l.ne.0) then


       if (l.eq.nr) go to 310

	if(l.gt.natomz) write(6,*) 'l,mynod,natomz',l,mynod,natomz
		call flush(6)

c
       iml=imass(l)
       inty=intype(imn,iml)
c
       sx0=-sxn+sx(l)
       sy0=-syn+sy(l)
       sz0=-szn+sz(l)
       sx0=sx0-dfloat(nint(sx0))
       sy0=sy0-dfloat(nint(sy0))
       sz0=sz0-dfloat(nint(sz0))
c
       xsr=h11*sx0+h12*sy0+h13*sz0
       ysr=h21*sx0+h22*sy0+h23*sz0
       zsr=h31*sx0+h32*sy0+h33*sz0

       xsr2=xsr*xsr
       ysr2=ysr*ysr
       zsr2=zsr*zsr

       rsq=xsr2+ysr2+zsr2
       if(rsq.ge.rss) goto 310


	neg=neg+1
	r1=dsqrt(rsq)*e00
 	rx(neg)=xsr*e00
	ry(neg)=ysr*e00
	rz(neg)=zsr*e00
	rxt(neg)=rx(neg)/r1
	ryt(neg)=ry(neg)/r1
	rzt(neg)=rz(neg)/r1
	r22(neg)=rsq*e00**2
	rr(neg)=r1
	list2(neg)=l

	if(r1/sigma(1).ge.ps(1)) write(6,*) 'nr,mynod,r1/sigma=',nr,mynod,r1/sigma(1)
		call flush(6)


310	l=llist(l)
	goto 3000
	endif
4000   continue
c
c


	negtot=negtot+neg

	vix=h(1,1)*x1(nr)+h(1,2)*y1(nr)+h(1,3)*z1(nr)
	viy=h(2,1)*x1(nr)+h(2,2)*y1(nr)+h(2,3)*z1(nr)
	viz=h(3,1)*x1(nr)+h(3,2)*y1(nr)+h(3,3)*z1(nr)

c	if (mynod.eq.0) write(6,*)neg,e00
	do 305 jx=1,neg
	 j=list2(jx)
	 iml=imass(j)

	vjx=h(1,1)*x1(j)+h(1,2)*y1(j)+h(1,3)*z1(j)
	vjy=h(2,1)*x1(j)+h(2,2)*y1(j)+h(2,3)*z1(j)
	vjz=h(3,1)*x1(j)+h(3,2)*y1(j)+h(3,3)*z1(j)


c = energy, forces.

	 dlj1=(rr(jx)/sigma(1))**(-1.0d0)
	 dlj6=dlj1*dlj1*dlj1*dlj1*dlj1*dlj1
	 elm1a=dlj6*(dlj6-1.0d0)*pa(1)-ljvc-ljdvr*(rr(jx)/sigma(1)-ljrc)
	 vija=elm1a

	 ere(1,nr)=ere(1,nr)+vija/4.0d0
	 ere(1,j)=ere(1,j)+vija/4.0d0

	 epot1(nr)=epot1(nr)+vija/4.0d0
	 epot1(j)=epot1(j)+vija/4.0d0
	 epot2(nr)=epot2(nr)+vija/4.0d0
	 epot2(j)=epot2(j)+vija/4.0d0

	 elmd=-12.0d0*dlj6*dlj6*dlj1+6.0d0*dlj6*dlj1
	 elmd=elmd*(1.0d0/sigma(1))*pa(1)-ljdvr

	 fix=-elmd*rxt(jx)
	 fiy=-elmd*ryt(jx)
	 fiz=-elmd*rzt(jx)
	 fix=fix*e00
	 fiy=fiy*e00
	 fiz=fiz*e00
	 fxre(1,nr)=fxre(1,nr)-fix/2.0d0
	 fyre(1,nr)=fyre(1,nr)-fiy/2.0d0
	 fzre(1,nr)=fzre(1,nr)-fiz/2.0d0
	 
	 fxre(1,j)=fxre(1,j)+fix/2.0d0
	 fyre(1,j)=fyre(1,j)+fiy/2.0d0
	 fzre(1,j)=fzre(1,j)+fiz/2.0d0

c heat current

	fac=(fix*vix+fiy*viy+fiz*viz)/dt

	facx=rx(jx)*fac/2.0d0
	facy=ry(jx)*fac/2.0d0
	facz=rz(jx)*fac/2.0d0

	if(sz(nr).lt.0.0d0) then
	 qtotx1=qtotx1+facx
	 qtoty1=qtoty1+facy
	 qtotz1=qtotz1+facz

	 qtotx3=qtotx3+facx
	 qtoty3=qtoty3+facy
	 qtotz3=qtotz3+facz
	endif

	if(sz(nr).ge.0.0d0) then
	 qtotx2=qtotx2+facx
	 qtoty2=qtoty2+facy
	 qtotz2=qtotz2+facz

	 qtotx4=qtotx4+facx
	 qtoty4=qtoty4+facy
	 qtotz4=qtotz4+facz
	endif

c stress

	 str1_re(1,jx)=rxt(jx)*rx(jx)*elmd/2.0d0
	 str2_re(1,jx)=ryt(jx)*ry(jx)*elmd/2.0d0
	 str3_re(1,jx)=rzt(jx)*rz(jx)*elmd/2.0d0
	 str4_re(1,jx)=rxt(jx)*ry(jx)*elmd/2.0d0
	 str5_re(1,jx)=rxt(jx)*rz(jx)*elmd/2.0d0
	 str6_re(1,jx)=ryt(jx)*rz(jx)*elmd/2.0d0


c compute stress tensor

	stress(jx,1)=stress(jx,1)-str1_re(1,jx)
	stress(jx,2)=stress(jx,2)-str2_re(1,jx)
	stress(jx,3)=stress(jx,3)-str3_re(1,jx)
	stress(jx,4)=stress(jx,4)-str4_re(1,jx)
	stress(jx,5)=stress(jx,5)-str5_re(1,jx)
	stress(jx,6)=stress(jx,6)-str6_re(1,jx)

305 	continue
	nr=llist(nr)
	goto 1000

	endif
5000    continue  ! end loop on cells for ion i.

c
c
c put the forces into s-space
c
      do 294 n=1,natomz
      fxx=fxre(1,n)*hi11+fyre(1,n)*hi12+fzre(1,n)*hi13
      fyy=fxre(1,n)*hi21+fyre(1,n)*hi22+fzre(1,n)*hi23
      fzz=fxre(1,n)*hi31+fyre(1,n)*hi32+fzre(1,n)*hi33
      fxre(1,n)=fxx
      fyre(1,n)=fyy
      fzre(1,n)=fzz
294   continue


c
      return
      end
c
c
