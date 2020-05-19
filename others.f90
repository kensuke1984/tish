!v0.1.0 Kensuke Konishi
subroutine pinput_tish(parameter_file, &
    re,ratc,ratl,tlen,np,omegai,imin,imax, &
    nzone,vrmin,vrmax,rho,vsv,vsh,qmu, &
    r0,eqlat,eqlon,mt,nr,theta,phi,lat,lon,output)
!------------------------------------------------------------------------
! Parameter Input
!------------------------------------------------------------------------
    use parameters
    implicit none
    character*160, intent(in) :: parameter_file
    integer, intent(out):: np,imin,imax,nzone,nr
    real(dp),intent(out) :: tlen,omegai,re,ratc,ratl
    real(dp),dimension(maxnzone), intent(out):: vrmin,vrmax,qmu
    real(dp),dimension(4,maxnzone), intent(out):: rho,vsv,vsh
    real(dp),dimension(maxnr), intent(out) :: theta,phi,lat,lon
    real(dp),intent(out) :: eqlat,eqlon,r0,mt(3,3)
    character*80,dimension(maxnr),intent(out) :: output
    real(dp) :: stlat,stlon,eqlattmp
    integer i,linenum,io
    logical:: file_exists
    character*80::buffer
    character*80,dimension(1000) :: lines

    inquire(file=parameter_file,exist=file_exists)
    if (.not. file_exists) stop 'parameter file does not exist.'

    linenum=0
    open(unit=1,file=parameter_file,status='old',action='read')
    do
        read(1, '(a)', iostat=io) buffer
        buffer = adjustl(buffer)
        if(buffer(1:1)=='c'.or.buffer(1:1)=='c'.or.buffer(1:1)=='!') cycle
        if(io/=0) exit
        linenum=linenum+1
        lines(linenum) = buffer
    enddo
    close(1)

    read(lines(1),*) tlen,np
    read(lines(2),*) re ! relative error (vertical grid)
    read(lines(3),*) ratc ! ampratio (vertical grid cut-off)
    read(lines(4),*) ratl ! ampratio (for l-cutoff)
    read(lines(5),*) omegai ! omegai
    omegai = -dlog(omegai)/tlen
    read(lines(6),*) imin,imax
    read(lines(7),*) nzone
    if (nzone > maxnzone) stop 'nzone is too large. (pinput)'

! structure
    do i=1,nzone
        read(lines(7+3*(i-1)+1),*) vrmin(i),vrmax(i),rho(1:4,i)
        read(lines(7+3*(i-1)+2),*) vsv(1:4,i)
        read(lines(7+3*(i-1)+3),*) vsh(1:4,i),qmu(i)
    enddo
! source parameter
    read(lines(3*nzone+8),*) r0,eqlat,eqlon
    eqlattmp=eqlat
    call translat(eqlattmp,eqlattmp)
    read(lines(3*nzone+9),*) mt(1,1:3), mt(2,2:3), mt(3,3)
    read(lines(3*nzone+10),*) nr
! station
    if (nr > maxnr) stop 'nr is too large. (pinput)'
    do i=1,nr
        read(lines(3*nzone+10+i),*) lat(i),lon(i)
        stlat = lat(i)
        stlon = lon(i)
        call translat(stlat,stlat)
        call calthetaphi(eqlattmp,eqlon,stlat,stlon,theta(i),phi(i))
    enddo
    theta(1:nr) = theta(1:nr) / 1.8d2 * pi
    phi(1:nr) = phi(1:nr) / 1.8d2 * pi
    do i=1,nr
        read(lines(3*nzone+10+nr+i),'(a)') output(i)
        output(i)=trim(output(i))
    enddo
    return
    end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calthetaphi(ievla,ievlo,istla,istlo,theta,phi)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    double precision,parameter:: pi=3.1415926535897932d0
    !c
    double precision:: ievla,ievlo,istla,istlo
    double precision:: evla,evlo,stla,stlo
    double precision:: theta,phi
    double precision:: gcarc,az
    double precision:: tc,ts
    !c
    !c transformation to spherical coordinates
    !c
    evla = 90.d0 - ievla
    stla = 90.d0 - istla
    !c
    evla = evla / 1.8d2 * pi
    evlo = ievlo / 1.8d2 * pi
    stla = stla / 1.8d2 * pi
    stlo = istlo / 1.8d2 * pi
    !c
    gcarc = dacos( dcos(evla) * dcos(stla)&
        + dsin(evla) * dsin(stla) * dcos(evlo - stlo) )
    !c
    tc = ( dcos(stla) * dsin(evla) &
        - dsin(stla) * dcos(evla) * dcos(stlo - evlo) )&
        / dsin(gcarc)
    ts = dsin(stla) * dsin(stlo - evlo) / dsin(gcarc)
    !c
    az = dacos(tc)
    if( ts < 0.d0 ) az = -1.d0 * az
    !c
    az = az * 1.8d2 / pi

    gcarc = gcarc * 1.8d2 / pi
    !c
    theta = gcarc
    phi   = 180.d0 - az
    return
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine translat(geodetic,geocentric)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    double precision,parameter::flattening = 1.d0 / 298.25d0
    double precision,parameter:: pi=3.1415926535897932d0
    double precision:: geocentric, geodetic

    integer:: flag
    !c      read(5,*) geodetic
    flag = 0
    if(geodetic > 90.d0) then
        geodetic = 1.8d2 - geodetic
        flag = 1
    endif
    !c
    geodetic = geodetic / 1.8d2 * pi
    geocentric = datan( (1.d0 - flattening) * (1.d0 - flattening)&
        * dtan(geodetic) )
    geocentric = geocentric * 1.8d2 / pi
    !c      if(geocentric < 0.d0 ) geocentric = 1.8d2 + geocentric
    if(flag == 1) geocentric = 1.8d2 - geocentric

    !c      write(6,*) 'geocentric latitude', geocentric
    return
end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calgrid( nzone,vrmin,vrmax,vs,rmin,rmax,&
    imax,lmin,tlen,vmin,gridpar,dzpar )
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    double precision,parameter:: pi=3.1415926535897932d0
      !c
    integer:: nzone,imax,lmin
    double precision:: vrmin(*),vrmax(*),vs(4,*)
    double precision:: rmin,rmax,tlen,vmin(*),gridpar(*),dzpar(*)
    integer:: izone,j
    double precision:: coef1,coef2,v(4),vs1,vs2,rh,omega,amax,gtmp
    !c
    do izone=1,nzone
        !c computing the S-velocity at each zone
        v(:) = vs(:,izone)
        vs1 = 0.d0
        vs2 = 0.d0
        do j=1,4
            if ( j==1 ) then
                coef1 = 1.d0
            else
                coef1 = coef1 * ( vrmin(izone) / rmax )
            endif
            if ( j==1 ) then
                coef2 = 1.d0
            else
                coef2 = coef2 * ( vrmax(izone) / rmax )
            endif
            vs1 = vs1 + v(j) * coef1
            vs2 = vs2 + v(j) * coef2
        enddo
        !c computing rh
        rh = vrmax(izone) - vrmin(izone)
        !c computing omega,amax
        omega = 2.d0 * pi * dble(imax) / tlen
        if ( vs1>=vs2 ) then
            vmin(izone) = vs2
        else
            vmin(izone) = vs1
        endif
        amax = vrmax(izone)
        gtmp = ( omega * omega ) / ( vmin(izone) * vmin(izone) ) &
            - ( (dble(lmin)+0.5d0) * (dble(lmin)+0.5d0) )&
            / ( amax * amax )
        if ( gtmp>0.d0 ) then
            dzpar(izone)   = dsqrt( 1.d0/gtmp )
            gridpar(izone) = rh / dzpar(izone)
        else
            dzpar(izone)   = 0.d0
            gridpar(izone) = 0.d0
        endif
    enddo
    !c rearangement of gridpar
    gtmp = sum(gridpar(1:nzone))
    do izone=1,nzone
        if ( gridpar(izone)>0.d0 ) then
            gridpar(izone) = gridpar(izone) / gtmp
        else
            rh = vrmax(izone) - vrmin(izone)
            gridpar(izone) = rh / ( rmax - rmin ) * 0.1d0
        endif
    enddo
    !c re-rearangement of gridpar
    gtmp = sum(gridpar(1:nzone))
    gridpar(1:nzone)=gridpar(1:nzone)/gtmp
    !c
    return
end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calra( maxnlay,maxnzone,&
    nlayer,gridpar,dzpar,nzone,vrmin,vrmax,&
    rmin,rmax,nnl,ra,re )
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Computing the number and the location of grid points.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    double precision,parameter:: pi=3.1415926535897932d0
    !c
    integer:: maxnlay,maxnzone
    integer:: nlayer
    integer:: nzone,nnl(maxnzone)
    double precision:: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax
    double precision:: ra(maxnlay+maxnzone+1)
    integer:: izone,itmp,i,ntmp
    double precision:: rh,re
    !c
    !c Initializing the data
    ra(1:maxnlay+maxnzone+1) = 0.d0
    nnl(1:nzone) = 0
    !c
    !c computing the number and the location of the grid points
    ra(1) = rmin
    itmp = 1
    do  izone=1,nzone
        rh = vrmax(izone) - vrmin(izone)
        if(dzpar(izone)==0.d0) then
            ntmp = 1
        else
            ntmp = int( sqrt(3.3d0 / re ) * rh / dzpar(izone) &
                / 2.d0 / pi  / 7.d-1 + 1 )
        endif
        !c                             ! ntmp (see Geller & Takeuchi 1995 6.2)
        nnl(izone) = ntmp
        if ( nnl(izone)<5 ) nnl(izone)=5
        do i=1,nnl(izone)
            itmp = itmp + 1
            ra(itmp) = vrmin(izone)&
                + rh * dble(i) / dble( nnl(izone) )
        enddo
    enddo
    !c
    !c recouting the total number of grid points
    nlayer =sum(nnl(1:nzone))
    !c
    return
end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calra0( nlayer,nzone,vrmin,vrmax,nnl,ra )
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c Computing the number and the location of grid points.
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    integer:: nlayer,nzone,nnl(nzone)
    double precision:: vrmin(*),vrmax(*),ra(*)
    integer:: izone,itmp,i
    double precision:: rmin,rmax,rh
    !c
    !c computing the number and the location of the grid points
    rmin = vrmin(1)
    rmax = vrmax(nzone)
    ra(1) = rmin
    itmp = 1
    do izone=1,nzone
        rh = vrmax(izone) - vrmin(izone)
        nnl(izone) = int( dble(nlayer) * rh / ( rmax - rmin ) ) + 1
        do i=1,nnl(izone)
            itmp = itmp + 1
            ra(itmp) = vrmin(izone) + rh * dble(i) / dble( nnl(izone) )
        enddo
    enddo
    !c recouting the total number of grid points
    nlayer = sum(nnl)
    !c
    return
end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calgra( isp,ra,r0,spn,spo,gra )
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    integer:: isp(*),spn,itmp
    double precision:: ra(*),r0,spo,gra(*)
    !c
    itmp = isp(spn) + dint( spo )
    gra(1) = ra(itmp)
    gra(2) = r0
    gra(3) = ra(itmp+1)
    !c
    return
end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calsp( ndc,nlayer,isp,jsp )
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c Computing the stack points.
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    integer:: ndc,nlayer(*)
    integer:: isp(*),jsp(*)
    integer:: i
    !c
    !c computation of isp,jsp,ksp,lsp
    isp(1) = 1
    jsp(1) = 1
    do i=1,ndc
        isp(i+1) = isp(i) + nlayer(i)
        jsp(i+1) = jsp(i) + 4 * nlayer(i)
    enddo
    !c
    return
end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calspo( ndc,rdc,nlayer,r0,rmin,rmax,ra,&
    isp,spo,spn )
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Computing the source location.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    integer:: ndc,nlayer,isp(*),spn
    double precision:: rdc(*),r0,rmin,rmax,ra(*),spo
    integer:: itmp
    !c
    !c checking the parameter
    if ( (r0<rmin).or.(r0>rmax) ) stop 'The source location is improper.(calspo)'
    !c computing 'spo'
    if ( r0==rmax ) then
        spo = dble(nlayer) - 0.01d0
        r0 = ra(nlayer) + (spo-dble(nlayer-1)) * ( ra(nlayer+1)-ra(nlayer) )
    !c	  write(6,*) 'r0 is changed to ',r0,spo
    else
        itmp = 2
110     continue
        if ( r0<ra(itmp) ) then
            continue
        else
            itmp = itmp + 1
            goto 110
        endif
        spo = dble(itmp-2) + ( r0-ra(itmp-1) ) / ( ra(itmp)-ra(itmp-1) )
        !c temporal handling
        if ( (spo-dble(itmp-2))<0.01d0 ) then
            spo = dble(itmp-2) + 0.01d0
            r0 = ra(itmp-1) + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
        !c	    write(6,*) 'r0 is changed to ',r0,spo
        endif
        if ( (spo-dble(itmp-2))>0.99d0 ) then
            spo = dble(itmp-2) + 0.99d0
            r0 = ra(itmp-1) + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
        !c	    write(6,*) 'r0 is changed to ',r0,spo
        endif
!c
    endif
!c computing 'spn'
    spn = 0
    itmp = 1
120 continue
    spn = spn + 1
    if ( r0<=rdc(itmp) ) then
        continue
    else
        itmp = itmp + 1
        goto 120
    endif
    !c changing 'spo'
    spo = spo - dble( isp(spn) - 1 )
    !c
    return
end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calstg( nzone,rrho,vsv,vsh,&
    nlayer,nnl,ra,rmax,vnp,vra,rho,&
    ecL,ecN )
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c Computing the structure grid points.
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nzone,nlayer,nnl(*),vnp
    double precision:: rrho(4,*),vsv(4,*),vsh(4,*),ra(*),rmax
    double precision:: vra(*),rho(*),ecL(*),ecN(*)
    double precision:: trho,tvsv,tvsh,coef
    integer:: izone,i,j,itmp,jtmp
    !c
    !c initializing the data
    do  i=1,nlayer+nzone+1
        vra(i) = 0.d0
        rho(i) = 0.d0
        ecL(i) = 0.d0
        ecN(i) = 0.d0
    enddo
    !c computing the structure grid points
    itmp = 0
    jtmp = 0
    do   izone=1,nzone
        do  i=1,nnl(izone)+1
            itmp = itmp + 1
            jtmp = jtmp + 1
            vra(itmp) = ra(jtmp)
            !c --- evaluating the density and elastic constants at this point
            trho = 0.d0
            tvsv = 0.d0
            tvsh = 0.d0
            do  j=1,4
                if ( j==1 ) then
                    coef = 1.d0
                else
                    coef = coef * ( vra(itmp) / rmax )
                endif
                trho = trho + rrho(j,izone) * coef
                tvsv  = tvsv  + vsv(j,izone)   * coef
                tvsh  = tvsh  + vsh(j,izone)   * coef
            enddo
            rho(itmp) = trho
            ecL(itmp)  = rho(itmp) * tvsv * tvsv
            ecN(itmp)  = rho(itmp) * tvsh * tvsh
        enddo
        jtmp = jtmp - 1
    enddo
    vnp = itmp
    !c
    return
end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calgstg( spn,rrho,vsv,vsh,&
    ra,vra,rmax,rho,ecL,ecN,r0,mu0 )
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Computing the structure grid points.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: spn
    double precision:: rrho(4,*),vsv(4,*),vsh(4,*)
    double precision:: ra(*),rmax
    double precision:: vra(*),rho(*),ecL(*),ecN(*),r0,mu0
    double precision:: trho,tvsv,tvsh,coef
    integer:: i,j
    !c
    !c initializing the data
    do  i=1,3
        vra(i) = 0.d0
        rho(i) = 0.d0
        ecL(i) = 0.d0
        ecN(i) = 0.d0
    enddo
    !c computing the structure grid points
    do i=1,3
        vra(i) = ra(i)
        !c --- evaluating the density and elastic constants at this point
        trho = 0.d0
        tvsv = 0.d0
        tvsh = 0.d0
        do  j=1,4
            if ( j==1 ) then
                coef = 1.d0
            else
                coef = coef * ( vra(i) / rmax )
            endif
            trho = trho + rrho(j,spn) * coef
            tvsv  = tvsv  + vsv(j,spn)   * coef
            tvsh  = tvsh  + vsh(j,spn)   * coef
        enddo
        rho(i) = trho
        ecL(i)  = rho(i) * tvsv * tvsv
        ecN(i)  = rho(i) * tvsh * tvsh
    enddo
    !c
    mu0 = ecL(2)
    !c
    return
end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calcoef( nzone,omega,q,coef )
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    double precision,parameter:: pi=3.1415926535897932d0
    !c
    integer:: izone,nzone
    double precision:: omega,q(*)
    complex(kind(0d0)):: coef(*)
    double precision:: aa,bb
    !c
    do izone=1,nzone
        if ( omega==0.d0 ) then
            aa = 1.d0
        else
            aa = 1.d0 + dlog( omega / ( 2.d0 * pi ) ) / ( pi * Q(izone) )
        endif
        bb = 1.d0 / ( 2.d0 * Q(izone) )
        coef(izone) = dcmplx( aa, bb ) * dcmplx( aa, bb )
    enddo
    !c
    return
end
!c

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calcutd(nzone,nnl,tmpr,rat,nn,ra,kc)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: nzone,nn,kc,nnl(*)
    !c	complex(kind(0d0)) tmpc(*)
    complex(kind(0d0)):: tmpr(*)
    double precision:: rat,ra(*)
    integer:: nc
    !c
    double precision:: cU(nn),rc
    double precision:: maxamp,amp(nn)
    integer:: iz,jz,jj,i,ml(nzone),tzone
    !c
    cU(1:nn) = 0
    !c
    iz = 2
    jz = 1
    cU(1:nn) = tmpr(1:nn)
    !c
    maxamp = -1.d0
    do i=1,nn
        amp(i) = cU(i)
        if(maxamp<amp(i)) maxamp = amp(i)
    enddo
    !c
    maxamp = maxamp * rat ! threshold value
    if(maxamp==0.d0) then
        kc = 1
        return
    endif
    !c
    do i=1,nn
        if(amp(i)>maxamp) then
            nc = i
            goto 140
        endif
    enddo
140 continue

    !c
    i = 1
    do jj=1,nzone
        i = i + nnl(jj)
        ml(jj) = i
    enddo
    !c
    do  jj=nzone,1,-1
        if(ml(jj)>nc) tzone = jj
    enddo
    !c
    rc = ra(nc)
    kc = nc
    !c

    return
end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine callsuf(omega,nzone,vrmax,vsv,lsuf)
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c
    implicit none
    integer:: nzone,lsuf
    double precision:: omega,vrmax(*),vsv(4,*)
    !c
    double precision:: tvs,coef
    integer:: i
    !c
    tvs = 0.d0
    do i=1,4
        if(i==1) then
            coef = 1.d0
        else
            coef = coef
        endif
        tvs = tvs + ( vsv(i,nzone) ) * coef
    enddo
    !c
    lsuf = int(omega * vrmax(nzone) / tvs - 0.5d0) + 1
    return
end
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calamp(g,l,lsuf,maxamp,ismall,ratl)
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer:: l,lsuf,ismall
    double precision:: maxamp,ratl
    complex(kind(0d0)):: g
    !c
    double precision:: amp,ampratio
    !c
    ampratio = 0.d0
    amp = cdabs(g)
    if( amp>maxamp ) maxamp = amp
    if ( amp/=0.d0.and.maxamp/=0.d0 ) then
        ampratio = amp / maxamp
    endif
    if( ampratio<ratl.and.l>lsuf ) then
        ismall = ismall + 1
    else
        ismall = 0
    endif
    !c
    return
end
!c
