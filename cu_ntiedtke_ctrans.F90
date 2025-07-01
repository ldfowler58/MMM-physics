!==================================================================================================================
 module cu_ntiedtke_ctrans
 use ccpp_kind_types,only: kind_phys
 use cu_ntiedtke_common,only: cmfcmin


 implicit none
 private
 public:: cu_ntiedtke_ctrans_run,     &
          cu_ntiedtke_ctrans_init,    &
          cu_ntiedtke_ctrans_finalize


 contains


!==================================================================================================================
!>\section arg_table_cu_ntiedtke_ctrans_init
!!\html\include cu_ntiedtke_ctrans_init.html
!!
 subroutine cu_ntiedtke_ctrans_init(errmsg,errflg)
!==================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!------------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine cu_ntiedtke_ctrans_init

!==================================================================================================================
!>\section arg_table_cu_ntiedtke_ctrans_finalize
!!\html\include cu_ntiedtke_ctrans_finalize.html
!!
 subroutine cu_ntiedtke_ctrans_finalize(errmsg,errflg)
!==================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!------------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine cu_ntiedtke_ctrans_finalize

!==================================================================================================================
 subroutine cu_ntiedtke_ctrans_run(klon,klev,nchem,ldcum,lddraf,kctype,kcbot,kctop,kdtop,grav,ztmst,do_scav, &
                                   fscav,chem,ptenc,ghti,paph,pmfu,pmfd,pmfude_rate,pmfdde_rate,errmsg,errflg)
!==================================================================================================================

!--- input arguments:
 logical,intent(in):: do_scav

 integer,intent(in):: klon,klev
 integer,intent(in):: nchem
 integer,intent(in),dimension(klon):: kctype,kcbot,kctop,kdtop

 logical,intent(in),dimension(klon):: ldcum,lddraf

 real(kind=kind_phys),intent(in):: grav,ztmst
 real(kind=kind_phys),intent(in),dimension(klon,klev):: pmfu,pmfude_rate
 real(kind=kind_phys),intent(in),dimension(klon,klev):: pmfd,pmfdde_rate
 real(kind=kind_phys),intent(in),dimension(klon,klev+1):: paph,ghti

 real(kind=kind_phys),intent(in),dimension(nchem):: fscav
 real(kind=kind_phys),intent(in),dimension(klon,klev,nchem):: chem

!--- inout arguments:
 real(kind=kind_phys),dimension(klon,klev,nchem):: ptenc

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).
 integer,intent(out):: &
    errflg      ! output error flag (-).

!--- local variables and arrays:
 logical,dimension(klon):: lldcum,llddraf3

 integer:: ic,jl,jk

 real(kind=kind_phys):: zcons,zmfmax
 real(kind=kind_phys),dimension(klon):: zmfs
 real(kind=kind_phys),dimension(klon,klev):: zmfuus,zmfudr
 real(kind=kind_phys),dimension(klon,klev):: zmfdus,zmfddr

!------------------------------------------------------------------------------------------------------------------

 zcons=1./(grav*ztmst)

 do ic = 1,nchem
    ptenc(:,:,ic) = 0._kind_phys
 enddo


!--- no convective transport for mid-level convection:
 do jl = 1,klon
    if(ldcum(jl) .and. kctype(jl) /= 3 .and. kcbot(jl)-kctop(jl) >= 1 ) then
       lldcum(jl)   = .true.
       llddraf3(jl) = lddraf(jl)
    else
       lldcum(jl) = .false.
       llddraf3(jl) = .false.
    endif
 enddo


!--- check and correct mass fluxes for CFL criterium:
 zmfs(:) = 1.
 do jk = 2,klev
    do jl = 1,klon
       if(lldcum(jl) .and. jk >= kctop(jl) ) then
          zmfmax = (paph(jl,jk)-paph(jl,jk-1))*0.8*zcons

          if(pmfu(jl,jk) > zmfmax) then
             zmfs(jl) = min(zmfs(jl),zmfmax/pmfu(jl,jk))
          endif
       endif
    enddo
 enddo

 do jk = 1,klev
    do jl = 1,klon
       if(lldcum(jl) .and. jk >= kctop(jl)-1) then
          zmfuus(jl,jk) = pmfu(jl,jk)*zmfs(jl)
          zmfudr(jl,jk) = pmfude_rate(jl,jk)*zmfs(jl)
       else
          zmfuus(jl,jk) = 0.
          zmfudr(jl,jk) = 0.
       endif

       if(llddraf3(jl) .and. jk >= kdtop(jl)-1) then
          zmfdus(jl,jk) = pmfd(jl,jk)*zmfs(jl)
          zmfddr(jl,jk) = pmfdde_rate(jl,jk)*zmfs(jl)
       else
          zmfdus(jl,jk) = 0.
          zmfddr(jl,jk) = 0.
       endif
    enddo
 enddo


!--- call to subroutine that computes the convective transport of chemical species:
 call cuctracer(klon,klev,nchem,kctop,kdtop,lldcum,llddraf3,grav,ztmst,do_scav,ghti,paph,zmfuus,zmfdus, &
                zmfudr,zmfddr,chem,ptenc)


 end subroutine cu_ntiedtke_ctrans_run

!==================================================================================================================
 subroutine cuctracer(klon,klev,ktrac,kctop,kdtop,ldcum,lddraf,grav,ztmst,do_scav,ght,paph,pmfu,pmfd, &
                      pudrate,pddrate,pcen,ptenc)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: klon,klev,ktrac
 integer,intent(in),dimension(klon):: kctop,kdtop

 logical,intent(in):: do_scav
 logical,intent(in),dimension(klon):: ldcum,lddraf

 real(kind=kind_phys),intent(in):: grav,ztmst
 real(kind=kind_phys),intent(in),dimension(klon,klev):: pmfu,pudrate
 real(kind=kind_phys),intent(in),dimension(klon,klev):: pmfd,pddrate
 real(kind=kind_phys),intent(in),dimension(klon,klev+1):: ght,paph
 real(kind=kind_phys),intent(in),dimension(klon,klev,ktrac):: pcen

!--- inout arguments:
 real(kind=kind_phys),intent(inout),dimension(klon,klev,ktrac):: ptenc

!--- variables and arrays:
 integer:: ik,jk,jl,jn

 logical,dimension(klon,klev):: llcumask,llcumbas

 real(kind=kind_phys):: zzp,zmfa,zerate,zposi
 real(kind=kind_phys),dimension(klon,klev):: zdp
 real(kind=kind_phys),dimension(klon,klev,ktrac):: zcen,zcu,zcd,zmfc,ztenc


!------------------------------------------------------------------------------------------------------------------

!--- initialization:
 do jk = 2,klev
    do jl = 1,klon
        llcumask(jl,jk) = .false.
        if(ldcum(jl)) then
           zdp(jl,jk) = grav/(paph(jl,jk+1)-paph(jl,jk))
           if( jk >= kctop(jl)-1) llcumask(jl,jk) = .true.
        endif
    enddo
 enddo


!--- loop over all chemical species:
 do jn = 1,ktrac
    !define chemical species at half levels:
    do jk = 2,klev
       ik = jk-1
       do jl = 1,klon
          zcen(jl,jk,jn) = pcen(jl,jk,jn)
          zcd(jl,jk,jn)  = pcen(jl,ik,jn)
          zcu(jl,jk,jn)  = pcen(jl,ik,jn)
          zmfc(jl,jk,jn) = 0.
          ztenc(jl,jk,jn)= 0.
       enddo
    enddo
    do jl = 1,klon
       zcu(jl,klev,jn) = pcen(jl,klev,jn)
    enddo

    !compute updraft values:
    do jk = klev-1,3,-1
       ik = jk + 1
       do jl = 1,klon
          if(llcumask(jl,jk)) then
             zerate = pmfu(jl,jk) - pmfu(jl,ik) + pudrate(jl,jk)
             zmfa   = 1./max(cmfcmin,pmfu(jl,jk))
             if(jk >= kctop(jl)) then
                zcu(jl,jk,jn) = (pmfu(jl,ik)*zcu(jl,ik,jn) + &
                zerate*pcen(jl,jk,jn)-pudrate(jl,jk)*zcu(jl,ik,jn))*zmfa
             endif
          endif
       enddo
    enddo

    !compute downdraft values:
    do jk = 3,klev
       ik = jk - 1
       do jl = 1,klon
          if(lddraf(jl) .and. jk == kdtop(jl)) then
             !note: in order to avoid final negative tracer values at LFS
             !the allowed value of ZCD depends on the jump in mass flux
             !at the LFS:
             zcd(jl,jk,jn) = 0.1*zcu(jl,jk,jn) + 0.9*pcen(jl,ik,jn)
          elseif(lddraf(jl).and.jk>kdtop(jl)) then
             zerate = -pmfd(jl,jk) + pmfd(jl,ik) + pddrate(jl,jk)
             zmfa = 1./min(-cmfcmin,pmfd(jl,jk))
             zcd(jl,jk,jn) = (pmfd(jl,ik)*zcd(jl,ik,jn) - &
             zerate*pcen(jl,ik,jn)+pddrate(jl,jk)*zcd(jl,ik,jn))*zmfa
          endif
       enddo
    enddo

    !in order to avoid negative tracer at KLEV, then adjust ZCD:
    jk = klev
    ik = jk - 1
    do jl = 1,klon
       if(lddraf(jl)) then
          zposi = -zdp(jl,jk) *(pmfu(jl,jk)*zcu(jl,jk,jn) + &
                  pmfd(jl,jk)*zcd(jl,jk,jn)-(pmfu(jl,jk)+pmfd(jl,jk))*pcen(jl,ik,jn))
          if(pcen(jl,jk,jn)+zposi*ztmst < 0.) then
             zmfa = 1./min(-cmfcmin,pmfd(jl,jk))
             zcd(jl,jk,jn) = ((pmfu(jl,jk)+pmfd(jl,jk))*pcen(jl,ik,jn) - &
                  pmfu(jl,jk)*zcu(jl,jk,jn)+pcen(jl,jk,jn) / &
                  (ztmst*zdp(jl,jk)))*zmfa
          endif
       endif
    enddo
 enddo


 do jn = 1,ktrac
    !compute fluxes:
    do jk = 2,klev
       ik = jk - 1
       do jl = 1,klon
          if(llcumask(jl,jk)) then
             zmfa = pmfu(jl,jk) + pmfd(jl,jk)
             zmfc(jl,jk,jn) = pmfu(jl,jk)*zcu(jl,jk,jn) + &
                              pmfd(jl,jk)*zcd(jl,jk,jn) - zmfa*zcen(jl,ik,jn)
          endif
       enddo
    enddo

    !compute tendencies:
    do jk = 2 , klev - 1
       ik = jk + 1
       do jl = 1,klon
          if(llcumask(jl,jk)) then
             ztenc(jl,jk,jn) = zdp(jl,jk)*(zmfc(jl,ik,jn)-zmfc(jl,jk,jn))
          endif
       enddo
    enddo
    jk = klev
    do jl = 1,klon
       if(ldcum(jl)) ztenc(jl,jk,jn) = -zdp(jl,jk)*zmfc(jl,jk,jn)
    enddo
 enddo


 do jn = 1,ktrac
    !update tendencies:
    do jk = 2,klev
       do jl = 1,klon
          if(llcumask(jl,jk)) then
             ptenc(jl,jk,jn) = ptenc(jl,jk,jn)+ztenc(jl,jk,jn)
          endif
       enddo
    enddo
 enddo

 end subroutine cuctracer

!==================================================================================================================
 end module cu_ntiedtke_ctrans
!==================================================================================================================
