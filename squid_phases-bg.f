c**********************************************************************
c*	Program for Chimera States in SQUIDS                          *
C**********************************************************************


c	Use gnuplot> set palette defined (1 'blue',2 'green',3 'red',4 'black') for chimera imageplot
c	or set palette model CMY rgbformulae 7,5,15 for r
c	set palette model CMY rgbformulae 7,5,16 for SI and DM
c	Chimeras are seen after a long time, say time=3500 for L=256, dt=0.02 and others parameters as in paper
	

	implicit none
	integer i,j,k,L,t,seed,iaa,iab,ip,ii,nbin,nsbin
	double precision time,dt,a,sum1,strength,dmsum,threshold
	double precision meanz,meanz2
	parameter (L=64,dt=0.02d0,time=5000.d0)
	parameter (nbin=16,nsbin=L/nbin,threshold=0.005d0)
c	Remember to check L in the f subroutine
	double precision phi(L),theta(L),lambda(L,L),lambdainv(L,L),f(L)
	double precision z(L),s(nbin),sd(nbin),avsd(nbin)
	double precision phiold(L),thetaold(L),f1(L),f2(L),f3(L),f4(L)
	double precision beta,gam,phiac,omega,lambda0,phir,avz,avz2
	double complex r

	open (1,file='squid_phase-bg')
	

c	Initializing Parameter values 
	beta=0.08d0
	gam=0.002d0
	lambda0=-0.05d0
	phiac=0.01d0
	phir=0.9d0
	omega=2.d0*dacos(-1.d0)/5.9d0

c	Specifying Lambda matrix
	do i=1,L
	do j=1,L
	if (i .eq. j) then
	lambda(i,j)=1.d0
	else
	lambda(i,j)=lambda0/((abs(i-j))**3.d0)
	end if
	end do
	end do

c	Inverting the Lambda Matrix
	call inverse(lambda,lambdainv,L)


c	Phase Diagram Loops
	do iaa=0,600,15
	do iab=0,600,15

c	write (*,*) iaa,iab

	gam=iaa*1.d-5
	beta=iab/1000.d0

c	Specifying Lambda matrix
	do i=1,L
	do j=1,L
	if (i .eq. j) then
	lambda(i,j)=1.d0
	else
	lambda(i,j)=lambda0/((abs(i-j))**3.d0)
	end if
	end do
	end do
	
c	Inverting the Lambda Matrix
	call inverse(lambda,lambdainv,L)

c	Initilizing Average Order Parameter

	avz=0.d0
	avz2=0.d0

c	Initializing average sd value for SI and DM
	do i=1,nbin
	avsd(i)=0.d0
	end do

c	Specifying Initial Conditions
	seed=12345
	do i=1,L
	seed=seed*12345
	phi(i)=phir*(-1.d0+2.d0*rand(seed))/2.d0
	theta(i)=0.d0
	phiold(i)=phi(i)
	thetaold(i)=theta(i)
	end do

	

c	Time Loop
	do t=1,int(time/dt)

c	if (mod(t,int(400.d0/dt)) .eq. 0) then
c	write (*,*) int(t*dt)
c	end if

c	RK Order 4 evolution, theta=d phi/d t

	call rkfn(beta,gam,lambdainv,phiac,omega,phiold,thetaold,f,
     &(t-1.d0)*dt)
c	i is Squid Index Loop
	do i=1,L
	f1(i)=f(i)
	end do

	call rkfn(beta,gam,lambdainv,phiac,omega,phiold+thetaold*dt
     &/2.d0,thetaold+f1*dt/2.d0,f,(t-1.d0)*dt+dt/2.d0)
	do i=1,L
	f2(i)=f(i)
	end do

	call rkfn(beta,gam,lambdainv,phiac,omega,phiold+thetaold*dt
     &/2.d0+f1*dt**2.d0/4.d0,thetaold+f2*dt/2.d0,f,(t-1.d0)*dt+dt/2.d0)
	do i=1,L
	f3(i)=f(i)
	end do

	call rkfn(beta,gam,lambdainv,phiac,omega,phiold+thetaold*dt+
     &f2*dt/2.d0,thetaold+f3*dt,f,(t-1.d0)*dt+dt)
	do i=1,L
	f4(i)=f(i)
	end do

	do i=1,L
	theta(i)=thetaold(i)+(f1(i)+2.d0*f2(i)+2.d0*f3(i)+f4(i))*dt/6.d0
       phi(i)=phiold(i)+thetaold(i)*dt+(f1(i)+f2(i)+f3(i))*dt**2.d0/6.d0
	end do

c	Remembering the current phi and phidot values for use in the next step

	do i=1,L
	phiold(i)=phi(i)
	thetaold(i)=theta(i)
	end do

c	Printing for the space-time plot, discarding initial transients
c	if ((t*dt .ge. 0.d0) .and. (mod(t,int(1.d0/dt)) .eq. 0)) then
c	do i=1,L
c	write (1,*) i,t*dt,phi(i)
c	end do
c	end if



c	Calculation of Global Order Parameter
c	z is complex order parameter, avz and avz2 are temporal averages taken in last 10000 integration steps

	if (t .ge. (int(time/dt)-9999)) then

	r=dcmplx(0.d0,0.d0)
	do i=1,L
	r=r+dcmplx(dcos(2.d0*dacos(-1.d0)*phi(i))/(L*1.d0),dsin(2.d0*
     &dacos(-1.d0)*phi(i))/(L*1.d0))
	end do

	avz=avz+dsqrt(real(r*conjg(r)))*1.d-4
	avz2=avz2+real(r*conjg(r))*1.d-4

C  Calculation of z (To have ordering)
          do i=1,L
          ip=i+1
          if (i.eq.L) then
	  ip=1
	  end if
          z(i)=phi(i)-phi(ip)
          enddo

C  Calculation for binned/local standard deviation
           ii=0
	   meanz=0.d0

	   do i=1,L
	   meanz=meanz+z(i)
	   end do

	   meanz=meanz/(L*1.d0)

           do j=1,nbin
           meanz2=0.d0
           do i=nsbin*ii+1,nsbin*ii+nsbin
           meanz2=meanz2+(z(i)-meanz)**2.d0
           enddo

           meanz2=meanz2/nsbin
           
           sd(j)=dsqrt(meanz2)
           ii=ii+1
           enddo

	   do i=1,nbin
	   avsd(i)=avsd(i)+sd(i)*1.d-4
	   end do

	end if

c	End of Time Loop
	end do

c	Calculation of SI and DM

	   do i=1,nbin
           if(avsd(i).lt.threshold)then
           s(i)=1.d0
           else
           s(i)=0.d0
           endif
           enddo

C   Strength of Incoherence
           sum1=0.d0
           do i=1,nbin
           sum1=sum1+s(i)
           enddo
           strength=1.d0-sum1/nbin

C    Discontinuity measure (chimera/multichimera)
           dmsum=0.d0
           do i=1,nbin
           j=i+1
           if(i.eq.nbin) then
	   j=1
	   end if
           a=s(i)-s(j)
           dmsum=dmsum+dabs(a) 
           enddo
           dmsum=dmsum/2.d0
	   
c	Write Time-averaged mod of O.P. and its S.D. to get fluctuations in time, SI,DM
	open (1,file='squid_phase-bg',Access = 'append',STATUS='OLD')
       write (1,*) beta,gam,avz,dsqrt(avz2-avz**2.d0),strength,dmsum
       	close(1)
       write (*,*) beta,gam

c	End of Phase Diagram Loop
	end do
	end do

	end
c	End of Main Program




c------------------------------------------------------------------------------------------------------------------------------------
	subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
	implicit none 
	integer n
	double precision a(n,n), c(n,n)
	double precision L(n,n), U(n,n), b(n), d(n), x(n)
	double precision coeff
	integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
	L=0.0
	U=0.0
	b=0.0

! step 1: forward elimination
	do k=1, n-1
   	do i=k+1,n
      	coeff=a(i,k)/a(k,k)
      	L(i,k) = coeff
      	do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      	end do
   	end do
	end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
	do i=1,n
  	L(i,i) = 1.0
	end do
! U matrix is the upper triangular part of A
	do j=1,n
  	do i=1,j
    	U(i,j) = a(i,j)
  	end do
	end do

! Step 3: compute columns of the inverse matrix C
	do k=1,n
  	b(k)=1.0
  	d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  	do i=2,n
    	d(i)=b(i)
    	do j=1,i-1
      	d(i) = d(i) - L(i,j)*d(j)
    	end do
  	end do
! Step 3b: Solve Ux=d using the back substitution
  	x(n)=d(n)/U(n,n)
  	do i = n-1,1,-1
    	x(i) = d(i)
    	do j=n,i+1,-1
      	x(i)=x(i)-U(i,j)*x(j)
    	end do
    	x(i) = x(i)/u(i,i)
  	end do
! Step 3c: fill the solutions x(n) into column k of C
  	do i=1,n
    	c(i,k) = x(i)
  	end do
  	b(k)=0.0
	end do
	end subroutine inverse	

C-----------------------------------------------------


c	Subroutine for evaluation of the function of the DEs
c	f(L) is the total output

	subroutine rkfn(bet,gam,linv,phiac,omega,phiold,thetaold,f,time)
	implicit none
	integer i,j,L
	parameter (L=64)
	double precision bet,gam,phiac,omega,time
	double precision linv(L,L), phiold(L),thetaold(L),f(L)
	
	do i=1,L
	f(i)=-gam*thetaold(i)-bet*dsin(2.d0*dacos(-1.d0)*phiold(i))
	do j=1,L
	f(i)=f(i)+linv(i,j)*(phiac*dcos(omega*time)-phiold(j))
	end do
	end do


	return
	end









