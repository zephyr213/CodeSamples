c
c version which mixes polarizations
c
c lattice for C diamond
c
c create a t=0 phonon wave packet
c
c ila=1 longitudinal packet
c
c ita=1 transverse packet

	implicit double precision (a-h,o-z)
	parameter(nx=2,ny=2,nz=6000)
	parameter(nrmax=192000)
	parameter(nions=8)
	parameter(ila=1,ita=0)
	dimension r(nions,3),rx(nrmax),ry(nrmax),rz(nrmax)
	dimension rx0(nrmax),ry0(nrmax),rz0(nrmax),uzi(nrmax)
	dimension x1(nrmax),y1(nrmax),z1(nrmax),ch(nrmax),uz(nrmax)
	dimension zr(2*nz+1,6,6),zi(2*nz+1,6,6),
     1          w(2*nz+1,6),gk(2*nz+1),weight(2*nz+1)
	dimension wkx(2*nz+1),wky(2*nz+1),wkz(2*nz+1)
	dimension unowrealx(nrmax)
	dimension unowrealy(nrmax)
	dimension unowrealz(nrmax)
	dimension unowimagx(nrmax)
	dimension unowimagy(nrmax)
	dimension unowimagz(nrmax)
	dimension unowimag(nrmax)
	integer imass(nions),iyttrium(nrmax),ivac(nrmax),igrain(nrmax)
	integer ilatt(12),im(nrmax),nplane(nrmax),im1(nions),imn(nrmax)
	complex cmp1,cmp2,ctheta,stheta
	complex unow,anow(2*nz+1,6),expt,vecnow,tfac,vecnowx,vecnowy
	complex unowx,unowy,xfac,yfac,unowx2,unowy2,unowz,vecnowz
	
	data(r(1,i),i=1,3)/ 0.0, 0.0, 0.0/
	data(r(2,i),i=1,3)/ 0.5, 0.5, 0.0/
	data(r(3,i),i=1,3)/ 0.5, 0.0, 0.5/
	data(r(4,i),i=1,3)/ 0.0, 0.5, 0.5/
	data(r(5,i),i=1,3)/ 0.25, 0.25, 0.25/
	data(r(6,i),i=1,3)/ 0.75, 0.75, 0.25/
	data(r(7,i),i=1,3)/ 0.75, 0.25, 0.75/
	data(r(8,i),i=1,3)/ 0.25, 0.75, 0.75/

	data(imass(i),i=1,8) /1,1,1,1,1,1,1,1/
	data(im1(i),i=1,8) /1,1,3,3,2,2,4,4/

	open(unit=8,file='inp.dat')
	open(unit=12,file='displace')
	open(unit=13,file='propagate')
	open(unit=15,file='structure')
	open(unit=66,file='freq.out')
	open(unit=67,file='eigen.out')
	open(unit=22,file='kpoints.dat')
	read(8,*)
	read(8,*) iz0,izwidth,envel0
	read(8,*)
	read(8,*) dt
	read(8,*) 
	read(8,*) angle


	dt0=2.7640186d-4

c dt is number of steps assuming dt0=2.7640186d-4
c make dt time in ps

	dt=dt*dt0

	

	write(12,*) nx*ny*nz*8

	cmp1=(1.0d0,0.0d0)
	cmp2=(0.0d0,1.0d0)
c center of wavepacket in real space
	r0z=-250.0d0

c read in kpoints.dat file

	read(22,*) nktot
	do nk=1,nktot
	 read(22,*) ix,iy,iz,wx,wy,wz
	 wkx(nk)=wx
	 wky(nk)=wy
	 wkz(nk)=wz
	enddo

c There should be phonon data for nz k-points.
c freq.out contains the frequencies
c eigen.out contains the eigenfunctions which describe displacements

	open(unit=66,file='freq.out')
	open(unit=67,file='eigen.out')

c fort.66 contains frequency data as

c      n(kpoint)       nm(modes 1 to 6)    omega(modes 1 to 6)   freq(THz)

c 
c fort.67 contains real and imag. parts in columns, six components for each
c mode, six modes per k-point

c gk is the k vector not including the 1/a

c first the frequencies

	twopi=8.0d0*datan(1.0d0)
	pi=0.5d0*twopi
	angle=angle*(pi/180.0d0)

	fac1=dcos(angle)
	fac2=dsin(angle)

	do iz=1,nktot
c	 gk(iz)=twopi*(iz-1)
c	 gk(iz)=gk(iz)/nz
	 gk(iz)=twopi*wkz(iz)
	 do nm=1,6
	  read(66,*) n1,n2,w(iz,nm),freq
 	  write(6,*) n1,n2,w(iz,nm),freq
	 enddo
	enddo

c now fort.67 for eigenfunctions

	do iz=1,nktot
	do nm=1,6
	 w1=0.0d0
	 w2=0.0d0
	 do mm=1,6
	  read(67,*) zr(iz,nm,mm),zi(iz,nm,mm)
	  if(mm.le.3) w1=w1+zr(iz,nm,mm)**2+zi(iz,nm,mm)**2
	  if(mm.gt.3) w2=w2+zr(iz,nm,mm)**2+zi(iz,nm,mm)**2
	 enddo
c	 write(6,*) iz,nm,w1,w2
	enddo
	enddo


c now decide the weighting of each k-point around iz0, width is
c given by izwidth

c if LA is desired, use third eigenfunction.
c if TA is desired, use both first and second eigenfunctions
c first place atoms at diamond lattice positions

	ntot=0
	do iz=nz,1,-1
        do n=nions,1,-1
	do ix=1,nx
	do iy=1,ny
	 dx=float(ix-nx/2)
	 dy=float(iy-ny/2)
	 dz=float(iz-nz/2)

	 	ntot=ntot+1

	  im(ntot)=imass(n)

	  rx(ntot)=r(n,1)+dx
	  ry(ntot)=r(n,2)+dy
	  rz(ntot)=r(n,3)+dz-0.875d0
	
	 rx0(ntot)=rx(ntot)
	 ry0(ntot)=ry(ntot)
	 rz0(ntot)=rz(ntot)

	  x1(ntot)=0.0d0
	  y1(ntot)=0.0d0
	  z1(ntot)=0.0d0

	  imn(ntot)=(iz-1)*4+im1(n)

	  if(n.gt.4) igrain(ntot)=2
	  if(n.le.4) igrain(ntot)=1

	  if(igrain(ntot).eq.1) rznow=rz0(ntot)
	  if(igrain(ntot).eq.2) rznow=rz0(ntot)-0.25d0

	  charge=0.0d0
	  ch(ntot)=charge
	  nplane(ntot)=iz
	  m=ntot

	  write(6,*) ntot,rx(m),ry(m),rz(m),ch(m),im(m),
     1   	   nplane(m),igrain(m)

	enddo
	enddo
	enddo
	enddo
c now find eigenvector components to create uz displacements

	do iz=1,nktot
	do nlambda=1,6
	 anow(iz,nlambda)=0.0d0
	enddo
	enddo

	do n=1,ntot
	    gkn=gk(iz0)
	    nt=igrain(n)
	    ntx=(nt-1)*3+1
	    nty=(nt-1)*3+2
	    ntz=(nt-1)*3+3
	    unowx=0.0d0
	    unowy=0.0d0
	    unowz=cmp1*zr(iz0,3,ntz)+cmp2*zi(iz0,3,ntz)
	    if(nt.eq.1) rznow=rz0(n)
	    if(nt.eq.2) rznow=rz0(n)-0.25d0
	    expt=dcos(gkn*(rznow-r0z))*cmp1+dsin(gkn*(rznow-r0z))*cmp2
	    width=real(izwidth)
	    envel=envel0*dexp(-0.10d0*(rznow-r0z)**2/width**2)*cmp1 
	    unowz=unowz*expt*envel
	    write(15,*) imn(n),real(unowx),real(unowy),real(unowz)
	do iz=1,nktot
	  nlambda=3
	    gkn=gk(iz)
	  expt=dcos(gkn*rznow)*cmp1-dsin(gkn*rznow)*cmp2
	  vecnowx=zr(iz,nlambda,ntx)*cmp1-zi(iz,nlambda,ntx)*cmp2
	  vecnowy=zr(iz,nlambda,nty)*cmp1-zi(iz,nlambda,nty)*cmp2
	  vecnowz=zr(iz,nlambda,ntz)*cmp1-zi(iz,nlambda,ntz)*cmp2
	  anow(iz,nlambda)=anow(iz,nlambda)+
     1      expt*(unowx*vecnowx+unowy*vecnowy+unowz*vecnowz)
	enddo
	enddo
	
	do iz=1,nktot
	amag1=0.0d0
	do nlambda=3,3
	  anreal=real(anow(iz,nlambda))
	  animag=aimag(anow(iz,nlambda))
	  amag=anreal**2+animag**2
	  amag1=amag1+amag
	enddo
	write(14,*) iz,dsqrt(amag1)
	enddo


	do n=1,ntot
	  unowx=0.0d0
	  unowy=0.0d0
	  unowz=0.0d0
	  nt=igrain(n)
	  if(nt.eq.1) rznow=rz0(n)
	  if(nt.eq.2) rznow=rz0(n)-0.25d0
	  ntx=(nt-1)*3+1
	  nty=(nt-1)*3+2
	  ntz=(nt-1)*3+3
 	  do iz=1,nktot
	  gkn=gk(iz)
	  do nlambda=3,3
	   expt=dcos(gkn*rznow)*cmp1+dsin(gkn*rznow)*cmp2
	   vecnowx=zr(iz,nlambda,ntx)*cmp1+zi(iz,nlambda,ntx)*cmp2
	   vecnowy=zr(iz,nlambda,nty)*cmp1+zi(iz,nlambda,nty)*cmp2
	   vecnowz=zr(iz,nlambda,ntz)*cmp1+zi(iz,nlambda,ntz)*cmp2
	   unowx=unowx+anow(iz,nlambda)*expt*vecnowx/(ntot/2)
	   unowy=unowy+anow(iz,nlambda)*expt*vecnowy/(ntot/2)
	   unowz=unowz+anow(iz,nlambda)*expt*vecnowz/(ntot/2)
	  enddo
	  enddo
	  unowrealx(n)=real(unowx)
	  unowrealy(n)=real(unowy)
	  unowrealz(n)=real(unowz)
	  unowimagx(n)=aimag(unowx)
	  unowimagy(n)=aimag(unowy)
	  unowimagz(n)=aimag(unowz)
	  write(12,*) imn(n),unowrealx(n),unowrealy(n),unowrealz(n)
	  write(15,*) rz0(n),unowrealx(n),unowrealy(n),unowrealz(n)
	enddo

c now obtain the velocities 


	do n=1,ntot
	  unowx=0.0d0
	  unowy=0.0d0
	  unowz=0.0d0
	  nt=igrain(n)
	  if(nt.eq.1) rznow=rz0(n)
	  if(nt.eq.2) rznow=rz0(n)-0.25d0
	  ntx=(nt-1)*3+1
	  nty=(nt-1)*3+2
	  ntz=(nt-1)*3+3
 	  do iz=1,nktot
	  gkn=gk(iz)
	  do nlambda=3,3
	   expt=-w(iz,nlambda)*dcos(gkn*rznow)*cmp2+
     1          w(iz,nlambda)*dsin(gkn*rznow)*cmp1
	   vecnowx=zr(iz,nlambda,ntx)*cmp1+zi(iz,nlambda,ntx)*cmp2
	   vecnowy=zr(iz,nlambda,nty)*cmp1+zi(iz,nlambda,nty)*cmp2
	   vecnowz=zr(iz,nlambda,ntz)*cmp1+zi(iz,nlambda,ntz)*cmp2
	   unowx=unowx+anow(iz,nlambda)*expt*vecnowx/(ntot/2)
	   unowy=unowy+anow(iz,nlambda)*expt*vecnowy/(ntot/2)
	   unowz=unowz+anow(iz,nlambda)*expt*vecnowz/(ntot/2)
	  enddo
	  enddo
	  vnowrealx=real(unowx)
	  vnowrealy=real(unowy)
	  vnowrealz=real(unowz)
	  write(12,*) imn(n),vnowrealx,vnowrealy,vnowrealz
	enddo

c now obtain derivatives of the velocities 


	do n=1,ntot
	  unowx=0.0d0
	  unowy=0.0d0
	  unowz=0.0d0
	  nt=igrain(n)
	  if(nt.eq.1) rznow=rz0(n)
	  if(nt.eq.2) rznow=rz0(n)-0.25d0
	  ntx=(nt-1)*3+1
	  nty=(nt-1)*3+2
	  ntz=(nt-1)*3+3
 	  do iz=1,nktot
	  gkn=gk(iz)
	  do nlambda=3,3
	   expt=-w(iz,nlambda)**2*dcos(gkn*rznow)*cmp1-
     1          w(iz,nlambda)**2*dsin(gkn*rznow)*cmp2
	   vecnowx=zr(iz,nlambda,ntx)*cmp1+zi(iz,nlambda,ntx)*cmp2
	   vecnowy=zr(iz,nlambda,nty)*cmp1+zi(iz,nlambda,nty)*cmp2
	   vecnowz=zr(iz,nlambda,ntz)*cmp1+zi(iz,nlambda,ntz)*cmp2
	   unowx=unowx+0.5d0*anow(iz,nlambda)*expt*vecnowx/(ntot/2)
	   unowy=unowy+0.5d0*anow(iz,nlambda)*expt*vecnowy/(ntot/2)
	   unowz=unowz+0.5d0*anow(iz,nlambda)*expt*vecnowz/(ntot/2)
	  enddo
	  enddo
	  vnowrealx=real(unowx)
	  vnowrealy=real(unowy)
	  vnowrealz=real(unowz)
	  write(12,*) imn(n),vnowrealx,vnowrealy,vnowrealz
	enddo

c now obtain second derivatives of the velocities 


	do n=1,ntot
	  unowx=0.0d0
	  unowy=0.0d0
	  unowz=0.0d0
	  nt=igrain(n)
	  if(nt.eq.1) rznow=rz0(n)
	  if(nt.eq.2) rznow=rz0(n)-0.25d0
	  ntx=(nt-1)*3+1
	  nty=(nt-1)*3+2
	  ntz=(nt-1)*3+3
 	  do iz=1,nktot
	  gkn=gk(iz)
	  do nlambda=3,3
	   expt=w(iz,nlambda)**3*dcos(gkn*rznow)*cmp2-
     1          w(iz,nlambda)**3*dsin(gkn*rznow)*cmp1
	   vecnowx=zr(iz,nlambda,ntx)*cmp1+zi(iz,nlambda,ntx)*cmp2
	   vecnowy=zr(iz,nlambda,nty)*cmp1+zi(iz,nlambda,nty)*cmp2
	   vecnowz=zr(iz,nlambda,ntz)*cmp1+zi(iz,nlambda,ntz)*cmp2
	   unowx=unowx+(1.0d0/6.0d0)*anow(iz,nlambda)*expt*vecnowx/(ntot/2)
	   unowy=unowy+(1.0d0/6.0d0)*anow(iz,nlambda)*expt*vecnowy/(ntot/2)
	   unowz=unowz+(1.0d0/6.0d0)*anow(iz,nlambda)*expt*vecnowz/(ntot/2)
	  enddo
	  enddo
	  vnowrealx=real(unowx)
	  vnowrealy=real(unowy)
	  vnowrealz=real(unowz)
	  write(12,*) imn(n),vnowrealx,vnowrealy,vnowrealz
	enddo

c now obtain  third derivative of the velocities 


	do n=1,ntot
	  unowx=0.0d0
	  unowy=0.0d0
	  unowz=0.0d0
	  nt=igrain(n)
	  if(nt.eq.1) rznow=rz0(n)
	  if(nt.eq.2) rznow=rz0(n)-0.25d0
	  ntx=(nt-1)*3+1
	  nty=(nt-1)*3+2
	  ntz=(nt-1)*3+3
 	  do iz=1,nktot
	  gkn=gk(iz)
	  do nlambda=3,3
	   expt=w(iz,nlambda)**4*dcos(gkn*rznow)*cmp1+
     1          w(iz,nlambda)**4*dsin(gkn*rznow)*cmp2
	   vecnowx=zr(iz,nlambda,ntx)*cmp1+zi(iz,nlambda,ntx)*cmp2
	   vecnowy=zr(iz,nlambda,nty)*cmp1+zi(iz,nlambda,nty)*cmp2
	   vecnowz=zr(iz,nlambda,ntz)*cmp1+zi(iz,nlambda,ntz)*cmp2
	   unowx=unowx+(1.0d0/24.0d0)*anow(iz,nlambda)*expt*vecnowx/(ntot/2)
	   unowy=unowy+(1.0d0/24.0d0)*anow(iz,nlambda)*expt*vecnowy/(ntot/2)
	   unowz=unowz+(1.0d0/24.0d0)*anow(iz,nlambda)*expt*vecnowz/(ntot/2)
	  enddo
	  enddo
	  vnowrealx=real(unowx)
	  vnowrealy=real(unowy)
	  vnowrealz=real(unowz)
	  write(12,*) imn(n),vnowrealx,vnowrealy,vnowrealz
	enddo


c now obtain fourth deriviative of the velocities 


	do n=1,ntot
	  unowx=0.0d0
	  unowy=0.0d0
	  unowz=0.0d0
	  nt=igrain(n)
	  if(nt.eq.1) rznow=rz0(n)
	  if(nt.eq.2) rznow=rz0(n)-0.25d0
	  ntx=(nt-1)*3+1
	  nty=(nt-1)*3+2
	  ntz=(nt-1)*3+3
 	  do iz=1,nktot
	  gkn=gk(iz)
	  do nlambda=3,3
	   expt=-w(iz,nlambda)**5*dcos(gkn*rznow)*cmp2+
     1          w(iz,nlambda)**5*dsin(gkn*rznow)*cmp1
	   vecnowx=zr(iz,nlambda,ntx)*cmp1+zi(iz,nlambda,ntx)*cmp2
	   vecnowy=zr(iz,nlambda,nty)*cmp1+zi(iz,nlambda,nty)*cmp2
	   vecnowz=zr(iz,nlambda,ntz)*cmp1+zi(iz,nlambda,ntz)*cmp2
	   unowx=unowx+(1.0d0/120.0d0)*anow(iz,nlambda)*expt*vecnowx/(ntot/2)
	   unowy=unowy+(1.0d0/120.0d0)*anow(iz,nlambda)*expt*vecnowy/(ntot/2)
	   unowz=unowz+(1.0d0/120.0d0)*anow(iz,nlambda)*expt*vecnowz/(ntot/2)
	  enddo
	  enddo
	  vnowrealx=real(unowx)
	  vnowrealy=real(unowy)
	  vnowrealz=real(unowz)
	  write(12,*) imn(n),vnowrealx,vnowrealy,vnowrealz
	enddo

c now evolve wave packet in time


	do it=1,1
	time=dt*it
	do n=1,ntot
	  unowx=0.0d0
	  unowy=0.0d0
	  unowz=0.0d0
	  nt=igrain(n)
	  if(nt.eq.1) rznow=rz0(n)
	  if(nt.eq.2) rznow=rz0(n)-0.25d0
	  ntx=(nt-1)*3+1
	  nty=(nt-1)*3+2
	  ntz=(nt-1)*3+3
 	  do iz=1,nktot
	  gkn=gk(iz)
	  do nlambda=3,3
	   expt=dcos(gkn*rznow)*cmp1+dsin(gkn*rznow)*cmp2
	   vecnowx=zr(iz,nlambda,ntx)*cmp1+zi(iz,nlambda,ntx)*cmp2
	   vecnowy=zr(iz,nlambda,nty)*cmp1+zi(iz,nlambda,nty)*cmp2
	   vecnowz=zr(iz,nlambda,ntz)*cmp1+zi(iz,nlambda,ntz)*cmp2
	   tfac=cmp1*dcos(-w(iz,nlambda)*time)+cmp2*dsin(-w(iz,nlambda)*time)
	   unowx=unowx+tfac*anow(iz,nlambda)*expt*vecnowx/(ntot/2)
	   unowy=unowy+tfac*anow(iz,nlambda)*expt*vecnowy/(ntot/2)
	   unowz=unowz+tfac*anow(iz,nlambda)*expt*vecnowz/(ntot/2)
	  enddo
	  enddo
	  unowrealx(n)=real(unowx)
	  unowrealy(n)=real(unowy)
	  unowrealz(n)=real(unowz)
	  write(13,*) rz0(n),unowrealx(n),unowrealy(n),unowrealz(n)
	enddo
	enddo

c Compute from unowrealx,unowrealy,unowrealz all of the normal-mode components.  This
c is partly a check on how truly orthogonal the phonon wave functions are.

	do n=1,ntot
	nt=igrain(n)
	if(nt.eq.1) rznow=rz0(n)
	if(nt.eq.2) rznow=rz0(n)-0.25d0
	ntx=(nt-1)*3+1
	nty=(nt-1)*3+2
	ntz=(nt-1)*3+3
	do iz=1,nktot
	do nlambda=3,3
	    gkn=gk(iz)
	    unowx=unowrealx(n)*cmp1+unowimagx(n)*cmp2
	    unowy=unowrealy(n)*cmp1+unowimagy(n)*cmp2
	    unowz=unowrealz(n)*cmp1+unowimagz(n)*cmp2
	  expt=dcos(gkn*rznow)*cmp1-dsin(gkn*rznow)*cmp2
	  vecnowx=zr(iz,nlambda,ntx)*cmp1-zi(iz,nlambda,ntx)*cmp2
	  vecnowy=zr(iz,nlambda,nty)*cmp1-zi(iz,nlambda,nty)*cmp2
	  vecnowz=zr(iz,nlambda,ntz)*cmp1-zi(iz,nlambda,ntz)*cmp2
	  anow(iz,nlambda)=anow(iz,nlambda)+
     1      expt*(unowx*vecnowx+unowy*vecnowy+unowz*vecnowz)
	 enddo
	enddo
	enddo

	do nlambda=1,6
	do iz=1,nktot
		anreal=real(anow(iz,nlambda))
		animag=aimag(anow(iz,nlambda))
		amag=anreal**2+animag**2
		write(50+nlambda,*) iz,dsqrt(amag)
	enddo
	enddo

	stop
 	end
	  

