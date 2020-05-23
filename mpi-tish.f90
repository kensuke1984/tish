program mpitish
!  ************** mpi-tish.f ****************
! Computation of SH synthetic seismograms
! in transversely isotropic media for anisotropic PREM
! using modified DSM operators & modified source representation.
! Synthetics for shallow events can be computed.
!
!                                                 2002.10 K.Kawai
!     2009. ?  Kensuke Konishi
! v0.1.0
! ----------------------------<<constants>>----------------------------
    use parameters
    implicit none
    character(len=160) :: parameter_file
    !c ----------------------------<<variables>>----------------------------
    !c variable for the trial function
    integer:: nnlayer,nlayer(maxnzone)
    integer:: l,m
    double precision:: ra(maxnlay+maxnzone+1),gra(3),plm(3,0:3,maxnr)
    complex(dp):: bvec(3,-2:2,maxnr)
    ! variable for the structure
    integer:: nzone
    integer:: ndc,vnp
    double precision:: rmin,rmax
    double precision:: vrmin(maxnzone),vrmax(maxnzone)
    double precision,dimension(4,maxnzone)::rrho,vsv,vsh
    double precision:: qmu(maxnzone)
    double precision:: vra(maxnlay+2*maxnzone+1)
    double precision:: rho(maxnlay+2*maxnzone+1)
    double precision:: ecL(maxnlay+2*maxnzone+1)
    double precision:: ecN(maxnlay+2*maxnzone+1)
    double precision,dimension(3)::gvra,grho,gecL,gecN
    complex(dp):: coef(maxnzone)
    ! variable for the periodic range
    integer:: np,imin,imax
    double precision:: tlen,omega,omegai
    complex(dp)::comega2
    complex(dp):: u(3,maxnr)
    ! variable for the source
    integer:: spn,ns
    double precision::r0,mt(3,3),spo,mu0,eqlat,eqlon
    ! variable for the station
    integer:: nr,ir
    double precision,dimension(maxnr)::theta,phi,lat,lon
    ! variable for the matrix elements
    complex(dp),dimension(2,maxnlay+1)::a,a0,a2
    double precision,dimension(4*maxnlay)::t,h1,h2,h3,h4
    double precision,dimension(8):: gt,gh1,gh2,gh3,gh4
    complex(dp):: aa(4),ga(8),ga2(2,3),gdr(3)
    complex(dp):: g( maxnlay+1 )
    ! variable for the file
    character(80)::output(maxnr)
    ! variable for grid spacing
    double precision:: tmpr(maxnlay+1)
    double precision:: gridpar(maxnzone),dzpar(maxnzone),vmin(maxnzone)
    double precision:: re,ratc,ratl,maxamp
    integer:: kc,lsuf,ismall,llog
    ! variable for the stack point
    integer:: isp(maxnzone),jsp(maxnzone),ins
    ! other variables
    integer:: i,j,ii,nn,lda,ier
    double precision:: eps,work( 4*maxnlay ),lsq
    complex(dp):: dr(maxnlay+1),z(maxnlay+1)
    complex(dp):: cwork( 4*maxnlay )
    integer:: ltmp(2),iimax

    data lda/ 2 /
    data eps/ -1.d0 /
    !**************MPI***********************************
	include 'mpif.h'
    integer::petot,my_rank,ierr,ista
    complex(dp), allocatable, dimension(:,:,:) :: outputu
    integer, allocatable, dimension (:) :: mpimin, mpimax

    call mpi_init (ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr)
    call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr)

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    allocate(mpimin(PETOT), mpimax(PETOT))
    !cccccccccccccccccccccccccccccccccccccccccccccccc

    !c *************** Inputting and computing the parameters ***************
   ! read input parameters
    call get_command_argument(1, parameter_file)
    !c --- inputting parameter ---
    if (my_rank==0) then
        call pinput_tish( parameter_file,re,ratc,ratl,&
            tlen,np,omegai,imin,imax,nzone,vrmin,vrmax,rrho,vsv,vsh,qmu,&
            r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_BCAST(re,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ratc,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ratl,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(tlen,1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(np,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(omegai, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(imin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(imax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nzone, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(vrmin, maxnzone, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(vrmax, maxnzone, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(rrho, 4*maxnzone, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(vsv, 4*maxnzone, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(vsh, 4*maxnzone, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(qmu, maxnzone, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(r0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(eqlat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(eqlon, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(mt, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(theta, maxnr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(phi, maxnr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(lat(1), maxnr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(lon(1), maxnr, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(output, 80*maxnr, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

    allocate (outputu(3,nr,0:np))

    !c --- computing the required parameters ---
    !c computing and checking the parameters
    rmin = vrmin(1)
    rmax = vrmax(nzone)
    ndc = nzone - 1

    if ( r0<rmin .or. rmax<r0 ) stop 'Location of the source is improper.'

    iimax = imax
    if( (rmax-r0)<shallowdepth) then ! option for shallow events
        !c computing of the number and the location of grid points
        iimax = int(tlen * 2)
        call calgrid( nzone,vrmin,vrmax,vsv,rmin,rmax,iimax,1,tlen,vmin,gridpar,dzpar )
        call calra (nnlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nlayer,ra,re )
        !c --- checking the parameter
        if ( nnlayer>maxnlay ) stop 'The number of grid points is too large.'
        !c computing the stack points
        call calsp( ndc,nlayer,isp,jsp )
        !c computing the source location
        call calspo( ndc,vrmax,nnlayer,r0,rmin,rmax,ra,isp,spo,spn )
        !c computing grids for source computations
        call calgra( isp,ra,r0,spn,spo,gra )
        !c ******************* Computing the matrix elements *******************
        !c computing the structure grid points
        call calstg( nzone,rrho,vsv,vsh,nnlayer,nlayer,ra,rmax,vnp,vra,rho,ecL,ecN)
        call calgstg( spn,rrho,vsv,vsh,gra,gvra,rmax,grho,gecL,gecN,r0,mu0 )
        do i=1,ndc+1
            call calmatc( nlayer(i),vnp,vra,rho,2,0,0,ra( isp(i) ),t( jsp(i) ),work( jsp(i) ) )
            call calmatc( nlayer(i),vnp,vra,ecL,2,1,1,ra( isp(i) ),h1( jsp(i) ),work( jsp(i) ) )
            call calmatc( nlayer(i),vnp,vra,ecL,1,1,0,ra( isp(i) ),h2( jsp(i) ),work( jsp(i) ) )
            call calmatc( nlayer(i),vnp,vra,ecL,0,0,0,ra( isp(i) ),h3( jsp(i) ),work( jsp(i) ) )
            call calmatc( nlayer(i),vnp,vra,ecN,0,0,0,ra( isp(i) ),h4( jsp(i) ),work( jsp(i) ) )
            call caltl( nlayer(i),vnp,vra,rho,ra( isp(i) ),work( jsp(i) ) )
            t(jsp(i):jsp(i)+4*nlayer(i)-1)=(t(jsp(i):jsp(i)+4*nlayer(i)-1)&
             +work(jsp(i):jsp(i)+4*nlayer(i)-1))/2
            call calhl( nlayer(i),vnp,vra,ecL,ra( isp(i) ),work( jsp(i) ) )
            h3(jsp(i):jsp(i)+4*nlayer(i)-1)=(h3(jsp(i):jsp(i)+4*nlayer(i)-1)&
             +work(jsp(i):jsp(i)+4*nlayer(i)-1))/2
            call calhl( nlayer(i),vnp,vra,ecN,ra( isp(i) ),work( jsp(i) ) )
            h4(jsp(i):jsp(i)+4*nlayer(i)-1)=(h4(jsp(i):jsp(i)+4*nlayer(i)-1)&
             +work(jsp(i):jsp(i)+4*nlayer(i)-1))/2
        enddo
        call calmatc( 2,3,gvra,grho,2,0,0,gra,gt, work )
        call calmatc( 2,3,gvra,gecL,2,1,1,gra,gh1,work )
        call calmatc( 2,3,gvra,gecL,1,1,0,gra,gh2,work )
        call calmatc( 2,3,gvra,gecL,0,0,0,gra,gh3,work )
        call calmatc( 2,3,gvra,gecN,0,0,0,gra,gh4,work )
        call caltl( 2,3,gvra,grho,gra,work )
        gt(1:8)=(gt(1:8)+work(1:8))/2
        call calhl( 2,3,gvra,gecL,gra,work )
        gh3(1:8)=(gh3(1:8)+work(1:8))/2
        call calhl( 2,3,gvra,gecN,gra,work )
        gh4(1:8)=(gh4(1:8)+work(1:8))/2

        nn = nnlayer + 1
        ns = isp(spn) + dint(spo)
        ins = 4 * ns - 3

        llog = 0
        do ii=1,2 ! omega-loop
            if(ii==1) then
                if(imin==0) then
                    i=1
                else
                    i=imin
                endif
            endif
            if(ii==2) i=imax
            omega = 2 * pi * i / tlen
            comega2 = dcmplx( omega, -omegai ) * dcmplx( omega, -omegai )
            call callsuf(omega,nzone,vrmax,vsv,lsuf)
            call calcoef(nzone,omega,qmu,coef)

            a0(:,1:nn)=0
            a2(:,1:nn)=0
            do j=1,ndc+1
             !c Computing the coefficient matrix 'a' in the solid part. cala0
                cwork(jsp(j):jsp(j)+4*nlayer(j)-1)= comega2*(t(jsp(j):jsp(j)+4*nlayer(j)-1))&
                -coef(j)*(h1(jsp(j):jsp(j)+4*nlayer(j)-1)-h2(jsp(j):jsp(j)+4*nlayer(j)-1)&
                +h3(jsp(j):jsp(j)+4*nlayer(j)-1)-2*h4(jsp(j):jsp(j)+4*nlayer(j)-1))

                call overlap( nlayer(j),cwork(jsp(j)),a0( 1,isp(j) ) )
                cwork(jsp(j):jsp(j)+4*nlayer(j)-1)=-coef(j)*( h4(jsp(j):jsp(j)+4*nlayer(j)-1) )
                call overlap( nlayer(j),cwork(jsp(j)),a2( 1,isp(j) ) )
            enddo

            kc = 1
            ismall = 0
            maxamp = -1
            ltmp(ii) = maxlmax
            do l=0,maxlmax ! l-loop
                if( 20<ismall  ) then
                    if(ltmp(ii)>l) ltmp(ii) = l
                    exit
                endif
                !c
                tmpr(1:maxnlay+1) = 0
                lsq = dsqrt( l*(l+1d0) )
                !c computing the coefficient matrix elements
                !c --- Initializing the matrix elements
                a(:,1:nn)=0
                ga2=0
                !c Computing the coefficient matrix 'a' in the solid part. cala
                a(1:2,1:nn) = a0(1:2,1:nn) + l*(l+1) * a2(1:2,1:nn)

                !c Computing the coefficient matrix 'a' in the solid part. calga
                aa(1:4) = comega2 * t(ins:ins+3)&
                - coef(spn) * ( h1(ins:ins+3)-h2(ins:ins+3)+h3(ins:ins+3)+(l*(l+1)-2)*h4(ins:ins+3) )
                ga(1:8)=comega2*gt(1:8)&
                - coef(spn)*( gh1(1:8)-gh2(1:8)+gh3(1:8)+(l*(l+1)-2)*gh4(1:8) )


                call overlap( 2,ga,ga2 )
                if( mod(l,100)==0) then
                    call dclisb0_pretreatment( a,nn,1,lda,g,eps,dr,z,ier)
                else
                    call dclisb_pretreatment( a,nn,1,lda,g,eps,dr,z,ier)
                endif

                do m=-2,2 ! m-loop
                    if(m==0) cycle
                    if(l<abs(m)) cycle
                    g(1:nn)=0
                    call calg2( l,m,spo,r0,mt,mu0,coef(spn),ga,aa,ga2,gdr,g( isp(spn) ) )
                    if( mod(l,100)==0) then
                        call dclisb0_kenja( a,nn,1,lda,g,eps,dr,z,ier)
                        tmpr(1:nn) = tmpr(1:nn) + cdabs(g(1:nn))
                    else
                        call dclisb_kenja( a(1,kc),nn-kc+1,1,lda,ns-kc+1,g(kc),eps,dr,z,ier)
                    endif

                    if( mod(l,100)==0) call calcutd(nzone,nlayer,tmpr,ratc,nn,ra,kc)

                    call calamp(g(nn),l,lsuf,maxamp,ismall,ratl)
                enddo ! m-loop
            enddo ! l-loop
        enddo ! omega-loop
        iimax = max(ltmp(1),ltmp(2)) * tlen / lmaxdivf
    endif ! option for shallow events

    !c computing of the number and the location of grid points
    call calgrid( nzone,vrmin,vrmax,vsv,rmin,rmax,iimax,1,tlen, vmin,gridpar,dzpar )

    call calra ( nnlayer, gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nlayer,ra,re )
    !c --- checking the parameter
    if ( nnlayer>maxnlay ) stop 'The number of grid points is too large.'
    !c computing the stack points
    call calsp( ndc,nlayer,isp,jsp )
    !c computing the source location
    call calspo( ndc,vrmax,nnlayer,r0,rmin,rmax,ra,isp,spo,spn )
    !c computing grids for source computations
    call calgra( isp,ra,r0,spn,spo,gra )
    !c ******************* Computing the matrix elements *******************
    !c computing the structure grid points
    call calstg( nzone,rrho,vsv,vsh,nnlayer,nlayer,ra,rmax,vnp,vra,rho,ecL,ecN)
    call calgstg( spn,rrho,vsv,vsh,gra,gvra,rmax,grho,gecL,gecN,r0,mu0 )
    do i=1,ndc+1
        call calmatc( nlayer(i),vnp,vra,rho,2,0,0,ra( isp(i) ),t( jsp(i) ),work( jsp(i) ) )
        call calmatc( nlayer(i),vnp,vra,ecL,2,1,1,ra( isp(i) ),h1( jsp(i) ),work( jsp(i) ) )
        call calmatc( nlayer(i),vnp,vra,ecL,1,1,0,ra( isp(i) ),h2( jsp(i) ),work( jsp(i) ) )
        call calmatc( nlayer(i),vnp,vra,ecL,0,0,0,ra( isp(i) ),h3( jsp(i) ),work( jsp(i) ) )
        call calmatc( nlayer(i),vnp,vra,ecN,0,0,0,ra(isp(i)),h4(jsp(i)),work(jsp(i)))
        call caltl(nlayer(i),vnp,vra,rho,ra(isp(i)),work(jsp(i)))
        t(jsp(i):jsp(i)+4*nlayer(i)-1)=(t(jsp(i):jsp(i)+4*nlayer(i)-1)&
        +work(jsp(i):jsp(i)+4*nlayer(i)-1))/2
        call calhl(nlayer(i),vnp,vra,ecL,ra(isp(i)),work(jsp(i)))
        h3(jsp(i):jsp(i)+4*nlayer(i)-1)=(h3(jsp(i):jsp(i)+4*nlayer(i)-1)&
        +work(jsp(i):jsp(i)+4*nlayer(i)-1))/2
        call calhl(nlayer(i),vnp,vra,ecN,ra(isp(i)),work(jsp(i)))
        h4(jsp(i):jsp(i)+4*nlayer(i)-1)=(h4(jsp(i):jsp(i)+4*nlayer(i)-1)&
        +work(jsp(i):jsp(i)+4*nlayer(i)-1))/2
    enddo
    call calmatc( 2,3,gvra,grho,2,0,0,gra,gt, work )
    call calmatc( 2,3,gvra,gecL,2,1,1,gra,gh1,work )
    call calmatc( 2,3,gvra,gecL,1,1,0,gra,gh2,work )
    call calmatc( 2,3,gvra,gecL,0,0,0,gra,gh3,work )
    call calmatc( 2,3,gvra,gecN,0,0,0,gra,gh4,work )
    call caltl( 2,3,gvra,grho,gra,work )
    gt(1:8)=(gt(1:8)+work(1:8))/2
    call calhl( 2,3,gvra,gecL, gra,work )
    gh3(1:8)=(gh3(1:8)+work(1:8))/2
    call calhl( 2,3,gvra,gecN, gra,work )
    gh4(1:8)=(gh4(1:8)+work(1:8))/2

    ! ******************** Computing the displacement *********************
    nn = nnlayer + 1
    ns = isp(spn) + dint(spo)
    ins = 4 * ns - 3
    llog = 0

    !ccccccccc MPIccccccccccccccccccccccccccc

    ! call simplesplit (imin, imax, PETOT, mpimin, mpimax)
    call trianglesplit (imin, imax, PETOT, mpimin, mpimax)

    do i= mpimin(my_rank+1), mpimax(my_rank+1)! omega-loop
        u(:,1:nr)=0
        if ( i/=0 ) then
            omega = 2 * pi * i / tlen
            comega2 = dcmplx( omega, -omegai ) * dcmplx( omega, -omegai )
            call callsuf(omega,nzone,vrmax,vsv,lsuf)
            plm(:,:,1:nr)=0
            call calcoef( nzone,omega,qmu,coef )

            a0(:,1:nn)=0
            a2(:,1:nn)=0
            do j=1,ndc+1
                !c Computing the coefficient matrix 'a' in the solid part. cala0
                cwork(jsp(j):jsp(j)+4*nlayer(j)-1)= comega2*(t(jsp(j):jsp(j)+4*nlayer(j)-1))&
                -coef(j)*(h1(jsp(j):jsp(j)+4*nlayer(j)-1)-h2(jsp(j):jsp(j)+4*nlayer(j)-1)&
                +h3(jsp(j):jsp(j)+4*nlayer(j)-1)-2*h4(jsp(j):jsp(j)+4*nlayer(j)-1))
                call overlap( nlayer(j),cwork(jsp(j)),a0( 1,isp(j) ) )
                cwork(jsp(j):jsp(j)+4*nlayer(j)-1)=-coef(j)*( h4(jsp(j):jsp(j)+4*nlayer(j)-1) )
                call overlap( nlayer(j),cwork(jsp(j)), a2( 1,isp(j) ) )
            enddo

            kc = 1
            ismall = 0
            maxamp = -1
            llog = maxlmax
            do l=0,maxlmax ! l-loop
                if( ismall>20 ) then
                    if(llog>l) llog = l
                    cycle
                endif

                tmpr(1:maxnlay+1) = 0
                lsq = dsqrt( l*(l+1d0) )
                !c ***** Computing the trial function *****
                do ir=1,nr
                    call calbvec( l,theta(ir),phi(ir),plm(1,0,ir),bvec(1,-2,ir) )
                enddo
                !c computing the coefficient matrix elements
                !c --- Initializing the matrix elements
                a(1:2,1:nn)=0
                ga2=0

                !c Computing the coefficient matrix 'a' in the solid part. cala
                a(1:2,1:nn) = a0(1:2,1:nn) + l*(l+1)* a2(1:2,1:nn)

                !c Computing the coefficient matrix 'a' in the solid part. calga
                aa(1:4)=comega2*t(ins:ins+3)-coef(spn)*( h1(ins:ins+3)-h2(ins:ins+3)+h3(ins:ins+3)+(l*(l+1)-2)*h4(ins:ins+3))
                ga(1:8) = comega2 * gt(1:8)- coef(spn)* ( gh1(1:8)-gh2(1:8)+gh3(1:8)+dble(l*(l+1)-2)*gh4(1:8) )

                call overlap( 2,ga,ga2 )

                if( mod(l,100)==0) then
                    call dclisb0_pretreatment( a,nn,1,lda,g,eps,dr,z,ier)
                else
                    call dclisb_pretreatment( a(1,kc),nn-kc+1,1,lda,ns-kc+1,g(kc),eps,dr,z,ier)
                endif

                do m=-2,2 ! m-loop
                    if(m==0) cycle
                    if(l<abs(m)) cycle
                    g(1:nn)=0
                    call calg2( l,m,spo,r0,mt,mu0,coef(spn), ga,aa,ga2,gdr,g( isp(spn) ) )
                    if( mod(l,100)==0) then
                        call dclisb0_kenja( a,nn,1,lda,g,eps,dr,z,ier)
                        !sum up c of the same l
                        tmpr(1:nn) = tmpr(1:nn) + cdabs(g(1:nn))
                        call calcutd(nzone,nlayer,tmpr,ratc,nn,ra,kc)
                    else
                        call dclisb_kenja( a(1,kc),nn-kc+1,1,lda,ns-kc+1,g(kc),eps,dr,z,ier)
                    endif

                    call calamp(g(nn),l,lsuf,maxamp,ismall,ratl)
                    u(2:3,1:nr) = u(2:3,1:nr) + g(nn) * bvec(2:3,m,1:nr) / lsq
                enddo          ! m-loop
            enddo             ! l-loop
        endif
        outputu (:,1:nr,i) =u (:,1:nr)

    enddo                   ! omega-loop
    if(my_rank /=0) &
        call mpi_send(outputu(1,1,mpimin(my_rank+1)), 3*nr*(mpimax(my_rank+1)-mpimin(my_rank+1)+1),&
            MPI_DOUBLE_complex,0,my_rank, MPI_COMM_WORLD, ierr)

    if(my_rank==0)then
        do i = 1,petot-1
            call mpi_recv(outputu(1,1,mpimin(i+1)),3*nr*(mpimax(i+1)-mpimin(i+1)+1),&
                MPI_DOUBLE_complex,i,i,MPI_COMM_WORLD, ista,ierr)
        enddo
        !c ************************** Files Handling **************************
        do ir = 1 ,nr
            open(1,file=trim(output(ir)),status='replace',form='unformatted',&
           access='stream',convert='big_endian')
            write(1) tlen,np,1,3,omegai,lat(ir),lon(ir),eqlat,eqlon,r0
            do i= imin, imax
                write(1) i,dble(outputu(1,ir,i)),dimag(outputu(1,ir,i))
                write(1) dble(outputu(2,ir,i)),dimag(outputu(2,ir,i))
                write(1) dble(outputu(3,ir,i)),dimag(outputu(3,ir,i))
            enddo
            close(1)
        enddo
    endif
    write(*,*) my_rank, "Ivalice looks to the horizon!"
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call mpi_finalize (ierr)
    stop
end program mpitish
