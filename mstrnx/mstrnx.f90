module paras
!     Parameters for MSTRN X
!--History
! 02.11.18 Created
! 04. 4. 7 Add for CFC absorption.
!--
  implicit none

  INCLUDE 'param.inc'  !! read kln, kbnd, kch

  integer,parameter:: kln1=kln+1 !! maximum number of boundarys
  integer,parameter:: kmol=7     !! # of gaserous spicies
  integer,parameter:: kpcl=11    !! # of aerosol spicies
  integer,parameter:: kpclc=3    !! aerosol distribution data
  integer,parameter:: kprg=5     !! surface parameter
  integer,parameter:: kvar=6     !! optical variables
  integer,parameter:: kp=26      !! # of pressure
  integer,parameter:: kt=3       !! # of temperature
  integer,parameter:: ksfc=7     !! # of surface spicies
  integer,parameter:: kmax2=6    !! # of optical valuables
  integer,parameter:: kaplk=5    !! # of dimension for planck function
  integer,parameter:: kplnk=3    !! # of order for plank function
  integer,parameter:: kabs=3     !! max # of absorption in each band
  integer,parameter:: kcfc=28    !! # of CFC species
  integer,parameter:: kflg=7     !! # of flags about bands
  integer,parameter:: kbnd1=kbnd+1
  real(8),parameter:: pstd=1013.25d0, tstd=273.15d0, avg=6.0221367d23
  real(8),parameter:: rgas=8.31451,grav=9.80665,airm=29.0

end module paras


subroutine ck(nln,zb,pb,tb,lyrflag,h2o,co2,o3,n2o,co,ch4,o2,wvn,wgt,dtau)
  use paras
  implicit none
  integer,intent(in)::nln,lyrflag
  real(8),intent(in)::zb(kln1),pb(kln1),tb(kln1),&
       h2o(kln),co2(kln),o3(kln),n2o(kln),co(kln),ch4(kln),o2(kln)
  real(8),intent(out)::wvn(kbnd1),wgt(kch,kbnd),dtau(kln,kch,kbnd)

  integer::l,nln1
  integer::iug=11
  character::err*64
  integer::nch(kbnd)
  
! atmospheric condition
  real(8)::cgas(kln,kmol),ccfc(kcfc),zl(kln),&
       pl(kln),tl(kln)
  real(8)::zbt(kln1),pbt(kln1),tbt(kln1)
! transfer
  character::cbnd*2

!--- exec
! open read files

  nln1 = nln+1

  if (nln>kln) then
     write(*,*) 'Error: increase kln in param.inc and param.h to at least ', nln, ' and recompile!'
     stop
  endif

! linear interpolation for center of layers
  if (lyrflag /= 0) then
     do l=1,nln
        zl(l)=sum(zb(l:l+1))*0.5d0
        pl(l)=pb(l)
        tl(l)=tb(l)
     enddo
  else
     do l=1,nln1
        zbt(l)=zb(nln1-l+1)
        tbt(l)=tb(nln1-l+1)
        pbt(l)=log(pb(nln1-l+1))
     enddo
     do l=1,nln
        zl(l)=sum(zb(l:l+1))*0.5d0
     enddo
     call lint(nln1,zbt,tbt,nln,zl,tl,err)
     if(err/='') return
     call lint(nln1,zbt,pbt,nln,zl,pl,err)
     if(err/='') return
     pl(1:nln)=exp(pl(1:nln))
  end if

! copy trace gas concentrations to final destination
  do l = 1, nln
     cgas(l,1) = h2o(l)
     cgas(l,2) = co2(l)
     cgas(l,3) = o3 (l)
     cgas(l,4) = n2o(l)
     cgas(l,5) = co (l)
     cgas(l,6) = ch4(l)
     cgas(l,7) = o2 (l)
  enddo

! CFCs concentration
!   1: CFC-11  ,  2: CFC-12   ,  3: CFC-13   ,  4: CFC-14  ,    5: CFC-113  , 
!   6: CFC-114 ,  7: CFC-115  ,  8: HCFC-21  ,  9: HCFC-22 ,   10: HCFC-123 ,
!  11: HCFC-124, 12: HCFC-141b, 13: HCFC-142b, 14: HCFC-225ca, 15: HCFC-225cb, 
!  16: HFC-32  , 17: HFC-125  , 18: HFC-134  , 19: HFC-134a  , 20: HFC-143a  ,
!  21: HFC-152a, 22: SF6      , 23: ClONO2   , 24: CCl4    ,   25: N2O5     , 
!  26: C2F6    , 27: HNO4       
!  ccfc(1:kcfc)=1.d-10   *** default by Nakajima et al.
  ccfc(1:kcfc)=0.


  write(cbnd,'(i2)') kbnd
  OPEN(IUG,FILE='mstrnx.data/PARAG.'//cbnd,status='old')

!! radiative transfer
  call dtrn3(iug,nln,zb,pb,tl,cgas,ccfc,wvn,wgt,dtau,nch,err)
  if(err/='') then
     write(*,*) 'Error: ',err; stop
  endif

! close files
  close(iug);
return
end subroutine ck

subroutine sdp2(nln,p,dp,np)
! derive absorption coefficient from pt table and fitting.
!--history
! 02.02.27  modified knupt in ckdng.f (Zhang Hua)
! 02.11.25  modified from knpt1.      (Miho Sekiguchi)
! 03.09.09  modified from sdp for free source form.(Miho Sekiguchi)
!--
  implicit none
  integer,intent(in)::nln
  real(8),intent(in)::p(nln)
  integer,intent(out)::np(nln)
  real(8),intent(out)::dp(nln)
! works
  integer,parameter::kp=26
  real(8)::prs(1:kp)=(/ &
       1.0130d+3,  0.63096d+3, 0.39811d+3,  0.25119d+3,  0.15849d+3, &
       1.0000d+2,  0.63096d+2, 0.39811d+2,  0.25119d+2,  0.15849d+2, &
       1.0000d+1,  0.63096d+1, 0.39811d+1,  0.25119d+1,  0.15849d+1, &
       1.0000d+0,  0.63096d+0, 0.39811d+0,  0.25119d+0,  0.15849d+0, &
       1.0000d-1,  0.63096d-1, 0.39811d-1,  0.25119d-1,  0.15849d-1,&
       1.0000d-2/)
  integer::l,ip
!--exec
  do l=1,nln
     do ip = 1, kp
        if (p(l) >= prs(ip)) exit
     enddo
     np(l)=ip
     if(np(l) > kp) np(l)=kp
     if(np(l) <= 1) np(l)=2
     dp(l)=log10(p(l)/prs(np(l)-1))/log10(prs(np(l))/prs(np(l)-1))
  enddo
end subroutine sdp2

subroutine tdok2(nln,akt,t,knu)
!----------------------------------------------------------------------c
!                using shi's formula: (t/to)**(a+bt)
!----------------------------------------------------------------------c
  implicit none
  integer,intent(in)::nln
  real(8),intent(inout)::akt(3,nln)
  real(8),intent(in)::t(nln)
  real(8),intent(out)::knu(nln)
  integer::l
  real(8)::a,b,at1,at2,bt1,bt2
!--exec
  at1 = log10(200.d0/260.d0); at2 = log10(320.d0/260.d0)
  
  do l=1,nln
     if(akt(2,l)<=-20.d0.and.akt(1,l)>-20.d0.and.akt(3,l)>-20.d0) &
          akt(2,l)=(akt(1,l)+akt(3,l))*0.5d0
     bt1 = akt(1,l)-akt(2,l)
     bt2 = akt(3,l)-akt(2,l)
     if(bt1==0.d0 .and. bt2==0.d0) then
        knu(l)=10**akt(2,l); cycle
     endif
     b   = (bt2/at2-bt1/at1)/120.d0
     a   = bt2/at2 - b*320.d0
     
     knu(l)=10**(akt(2,l))*(t(l)/260.d0)**(a+b*t(l))
  enddo
  return
end subroutine tdok2
subroutine dtrn3(iug,nln,zb,pl,tl,cgas,ccfc,wv,wgt,dtau,nch,err)
!subroutine dtrn3(iug,iup,iuv,ams,nln,pl,tb,tl,gtmp,cpcl,gdcfrc,cgas,ccfc, &
!     prg,fd,fu,err)
! flux calculation
!--- history
! 93. 3.18 (quarante) hitran 1992  k-distributi
!     4.14 (teruyuki) registered
! 94. 6.11 change fitting yy to sqrt(sqrt(b))
!          add layer temperature for planck 2-dim. expansion term
! 94. 6.15 use layer and layer boundary p,t and avoid z
!     6.23 revise cplk
!    11.25 change h2o continuum treatment
!    11.29 change planck function treatment
! 94.12. 2 (yoko)     partial cloudiness approximation by semi-random method.
!    12. 5 polish up h2o continuum treatment
!    12.11 kd treatment; pt-interpolation -> function fitting
!    12.14 add cfcs. number of optical properties ifbnd changes 7 -> 8
! 95. 3. 8 revise o2 continuum optical thickness calculation
! 97.04.28  modified from mstr7e kpcl=8 -> kpcl=11
!				 iup -> iupp,iupg,iupv
!				 rmode interpolation(linear) by Y.Tsushima
! 02.11.21 modidfied from dtrn22 by M.Sekiguchi
! 04. 4. 7 Add CFCs for mstrnX by M.Sekiguchi.
!--- input
! iug      I           read unit number of gas line absorption para file
! iup      I           read unit number of other extinction para file
! iuv      I           read unit number of source file of iup
! nln       I        number of layers.
! pl      R(nln)     atmospheric pressure at the layer (mb) (94. 6.15)
! tb      R(nln1)    temperature of layer boundaries(k)
! tl      R(nln)     temperature of layers(k) (94. 6.11)
! gtmp     R         ground temperature (k)
! gdcfrc  R(knln)     cloud cover rate of sublayer.
! cpcl R(nln,npcl,3) parameter packet for particulates in sublayers. (l,i,j)
!                    l: layer number (from top to bottom)
!                       1: water cloud               2: ice cloud
!                       3: dust-like                 4: soot
!                       5: volcanic-ash              6: h2so4
!                       7: rural                     8: sea salt
!                       9: urban                    10: tropo.
!                      11: yellow dust
!                    j=1: concentration in ppmv
!                         = aerosol volume(cm3/cm3)*t*p0/t0/p*1.0e6
!                         with t0=273, p0=1atm.
!                         type 7,8,9,10 = dry concentration(parapc.voldp)
!                                       = total concentration(parapc.volp)
!                    j=2  mode radius(cm)
!                    j=3 undecided for future use.
!
! cgas  R(nln,ngas)   gas concentration in the sublayer (ppmv)
!  1: h2o, 2: co2, 3: o3, 4: n2o, 5: co, 6: ch4, 7: o2
!
! ccfc   R(kcfc)     CFCs concentrations
!   1: CFC-11  ,  2: CFC-12   ,  3: CFC-13   ,  4: CFC-14  ,    5: CFC-113  , 
!   6: CFC-114 ,  7: CFC-115  ,  8: HCFC-21  ,  9: HCFC-22 ,   10: HCFC-123 ,
!  11: HCFC-124, 12: HCFC-141b, 13: HCFC-142b, 14: HCFC-225ca, 15: HCFC-225cb, 
!  16: HFC-32  , 17: HFC-125  , 18: HFC-134  , 19: HFC-134a  , 20: HFC-143a  ,
!  21: HFC-152a, 22: SF6      , 23: ClONO2   , 24: CCl4    ,   25: N2O5     , 
!  26: C2F6    , 27: HNO4     , 28: SF5CF3  
!
! prg     r(5)       parameter packet for surface reflection.
!                      (1)             (2)           (3)
!                      1 ocean         -ams-         -tr- (automatically set
!                      2 wet land      water content      inside this loutine)
!                      3 dry land       -              -
!                      4 low plants     -              -
!                      5 forest         -              -
!                      6 snow           -              -
!                      7 ice            -              -
!--- output
! fd      r(nln1,2)  downward flux at the sublayer interface (wm-2)
!                     (from top to bottom)
!                    j=1 shorter than 4 micron, j=2 longer than 4 micron
! fu      r(nln1,2)   upward   flux at the sublayer interface.
! err     c*64       error code.  if ' ' then no error.
!--
  use paras
  implicit none
  integer,intent(in)::iug,nln
  real(8),intent(in)::zb(kln1),pl(kln),tl(kln), &
       cgas(kln,kmol),ccfc(kcfc)
  integer,intent(out)::nch(kbnd)
  real(8),intent(out)::wv(kbnd1)
  real(8),intent(out)::wgt(kch,kbnd),dtau(kln,kch,kbnd)
  character::err*64

! gttbl
  integer::nbnd,nda,np,nt,nabs(kbnd),iabs(kabs,kbnd),nflg, &
       iflgb(kflg,kbnd)
  real(8)::prs(kp),temp(kt), &
       akd(kch,kp,kt,kabs,kbnd),skd(kch,kp,kt,kbnd),acfc(kcfc,kbnd), &
       sr(kbnd,ksfc)

! gas concentration for absorption
  real(8)::amtpb(kmol,kln),wbrd(kln)

! p-t table
  integer::mp(nln)
  real(8)::akt(kt,kln),knu(kln),sku(kln),tkd(kln),ddp(nln)

! work 
  integer::iw,l,ia,ich
  real(8)::dz(kln)
!--- exec

! read optical parameters
     call gttbl(iug,nbnd,nda,np,nt,nabs,iabs,nflg,iflgb,nch, &
          wv,prs,temp,wgt,akd,skd,acfc,sr)

  do l=1,nln
!! get p, t, dp, dz from grid data (94. 6.15)
!  dz is actually the air column (km), not the layer width!
     dz(l)=(zb(l)-zb(l+1))/1000.0*pl(l)/pstd/tl(l)*tstd
     amtpb(1:kmol,l)=cgas(l,1:kmol)*dz(l)*1.d-1
     wbrd(l)=dz(l)*1.d5
  enddo

! for interpolation of layers (p,T)
  call sdp2(nln,pl,ddp,mp)

! loop for bands
  do iw = 1,nbnd

! loop for channels
     
     do ich = 1, nch(iw)
        tkd(1:nln)=0.d0
        if(iflgb(1,iw)>0) then
           if(iflgb(7,iw)>0) then
              do l=1,nln
                 tkd(l)=sum(10**acfc(1:kcfc,iw)*ccfc(1:kcfc))*wbrd(l)
              enddo
           endif
           do ia=1,nabs(iw)
              do l=1,nln
                 akt(1:kt,l)=akd(ich,mp(l)-1,1:kt,ia,iw) &
                      +(akd(ich,mp(l),1:kt,ia,iw)-akd(ich,mp(l)-1,1:kt,ia,iw))&
                      *ddp(l)
              enddo
              call tdok2(nln,akt(1:kt,1:nln),tl,knu(1:nln))
              tkd(1:nln)=tkd(1:nln)+knu(1:nln)*amtpb(iabs(ia,iw),1:nln)
              
           enddo

           if(iflgb(5,iw)>0) then
              do l=1,nln
                 akt(1:kt,l)=skd(ich,mp(l)-1,1:kt,iw)+(skd(ich,mp(l),1:kt,iw)&
                      -skd(ich,mp(l)-1,1:kt,iw))*ddp(l)                     
              enddo
              call tdok2(nln,akt(1:3,1:nln),tl,sku(1:nln))
              tkd(1:nln)=tkd(1:nln)+sku(1:nln)*amtpb(1,1:nln)**2/ &
                   (amtpb(1,1:nln)+wbrd(1:nln))
           endif
        endif

        do l=1,nln
           dtau(l,ich,iw) = tkd(l)
        enddo
     enddo
  enddo
  err = ''
  return
end subroutine dtrn3
subroutine gttbl(iug,nbnd,nda,np,nt,nabs,iabs,nflg,iflgb,nch,wv, &
     prs,temp,wgt,akd,skd,acfc)
!subroutine gttbl(iug,nbnd,nda,np,nt,nabs,iabs,nflg,iflgb,nch,wv, &
!     prs,temp,wgt,akd,skd,acfc,nsfc,nplk,naplk,aplnk,fsol,sr,ry,qmol,q,npcl, &
!     nvar,r0,err)
! read optical parameters
!--- history
! 93. 3.28 created by quarante
!     4.26 registered by T. Nakajima
! 94.11.25 change h2o continuum treatment
!    11.29 change planck function treatment
!    12.14 add cfcs. number of optical properties ifbnd changes 7 -> 8
! 97.04.28  modified from mstr7e kpcl=8 -> kpcl=11
!				 iup -> iup,iupg,iuv
!				 rmode interpolation(linear) by Y.Tsushima
! 02.11.18 modified from mstrn8 by M. Sekiguchi 
! 04. 4. 7 Add CFCs for mstrnX by M. Sekiguchi
!--- input
! iug       I        read unit number of parameters of gas line absorption
! iup       I        read unit number of parameters of other extinction
! iuv       I        read unit number of source file of iup
!--- output
! nbnd      I        # of bands
! nda       I        # of streams
! np        I        # of pressures for P-T table
! nt        I        # of temperatures for P-T table
! nabs   I(kbnd)     # of absorbers in each band
! iabs  I(kabs,kbnd) # of absorbers spicies in each band
! nflg  I(kbnd)      # of flags
! iflgb I(kflg,kbnd) optical flag
! 1:line abs, 2:scat, 3:plnk, 4:solar, 5:h2o self, 6:CFC
! wv     R(kbnd1)    boundary wavenumbers [cm-1]
! prs    R(kp)       pressures for P-T table
! temp   R(kt)       temperatures for P-T table
! nch    I(kbnd)     # of subintervals in each band
! wgt   R(kch,kbnd)   weights of subintervals in each band
! akd  R(kch,kp,kt,kabs,kbnd)  absoption coefficients
! skd  R(kch,kp,kt,kbnd)  absoption coefficients for H2O self broadening
! acfc  R(kcfc,kbnd) ansorption coefficients for CFC
! npcl      I        # of particle spicies
! nplk      I        order of planck function fitting
! naplk     I        # of planck function fitting
! aplnk R(kbnd,kaplk) coefficients of planck function fitting
! fsol   R(kbnd)     solar constants in each bands [W/m**2]
! sr    R(kbnd,ksfc) surface condition in each bands
! ry     R(kbnd)     rayleigh scattering in each bands
! qmol   R(kmax2)    
! q   R(kbnd,kpcl,kvar,kmax2) coefficients of particles
! nvar  I(kpcl,2)    1: effective radius change or humidity growth
!                    2: # of values
! r0   R(kpcl,kvar)  values in each particles
!--
  use paras, only: kbnd,kbnd1,kabs,kp,kt,kch,kcfc,kaplk,ksfc,kmax2, &
       kpcl,kvar,kflg
  implicit none

! file number
  integer,intent(in)::iug

! gas parameters
  integer,intent(out)::nbnd,nda,np,nt,nabs(kbnd),iabs(kabs,kbnd),nflg, &
       iflgb(kflg,kbnd),nch(kbnd)
  real(8),intent(out)::wv(kbnd1),prs(kp),temp(kt),wgt(kch,kbnd), &
       akd(kch,kp,kt,kabs,kbnd),skd(kch,kp,kt,kbnd),acfc(kcfc,kbnd)
  integer::iw,ip,it,ia,ncfc

!--- exec
! read gas parameters
  read(iug,*) nbnd,nda,np,nt,nflg,ncfc

! band boundaries
  read(iug,*)
  read(iug,*) wv(1:nbnd+1)

! log(pressure) grids
  read(iug,*)
  read(iug,*) prs(1:np)

! temperature grids
  read(iug,*)
  read(iug,*) temp(1:nt)

! quantities for each band
  do iw=1,nbnd

!! optical properties flag
     read(iug,*) 
     read(iug,*) iflgb(1:nflg,iw)

!! number of subintervals
     read(iug,*) 
     read(iug,*) nch(iw)

!! weights for channels
     read(iug,*)
     read(iug,*) wgt(1:nch(iw),iw)

!! molecules
     read(iug,*)
     read(iug,*) nabs(iw)
!!! major gas absorption
     if(nabs(iw)>0) then
        do ia=1,nabs(iw)
           read(iug,*) iabs(ia,iw)
           do it=1,nt; do ip=1,np
              read(iug,*) akd(1:nch(iw),ip,it,ia,iw)
           enddo; enddo
        enddo
     endif
!!! H2O continuum
     if(iflgb(5,iw)>0) then
        read(iug,*) 
        do it=1,nt; do ip=1,np
           read(iug,*) skd(1:nch(iw),ip,it,iw)
        enddo; enddo
     endif
!!! CFC absorption
     if(iflgb(7,iw)>0) then
        read(iug,*) 
        read(iug,*) acfc(1:kcfc,iw)
     endif

   enddo
   return
end subroutine gttbl
SUBROUTINE LINT(K,X,Y,N,A,B,ERR)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::K,N
  REAL(8),INTENT(IN)::X(K),Y(K),A(N)
  REAL(8),INTENT(OUT)::B(N)
  INTEGER::I,J,MIN,MAX
  REAL(8)::WGT
  CHARACTER::err*64
!--exec
  ERR=''
  MIN=COUNT(A(:)<X(1))
  MAX=COUNT(A(:)>X(K))
  DO I=1,N
     IF(I<=MIN) THEN
        WGT=(X(1)-A(I))/(X(2)-X(1))
        B(I)=Y(1)-(Y(2)-Y(1))*WGT
     ELSEIF(I<=N-MAX) THEN
        DO J=1,K-1
           IF(X(J)<=A(I) .and. X(J+1)>A(I)) then
              WGT=(A(I)-X(J))/(X(J+1)-X(J))
              B(I)=Y(J)*(1-WGT)+Y(J+1)*WGT
              EXIT
           ENDIF
        ENDDO
     ELSE
        WGT=(A(I)-X(K))/(X(K)-X(K-1))
        B(I)=Y(K)+(Y(K)-Y(K-1))*WGT
     ENDIF
  ENDDO
END SUBROUTINE LINT

