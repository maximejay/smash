!%      (MD) Module Differentiated.
!%
!%      Subroutine
!%      ----------
!%
!%      - upstream_discharge
!%      - linear_routing
!%      - kinematic_wave1d
!%      - lag0_time_step
!%      - lr_time_step
!%      - kw_time_step

module md_routing_operator

    use md_constant !% only : sp
    use mwd_setup !% only: SetupDT
    use mwd_mesh !% only: MeshDT
    use mwd_options !% only: OptionsDT
    use mwd_returns !% only: ReturnsDT
    use mwd_discontinuities

    implicit none

contains

    subroutine upstream_discharge(mesh, row, col, ac_q, qup)

        implicit none

        type(MeshDT), intent(in) :: mesh
        integer, intent(in) :: row, col
        real(sp), dimension(mesh%nac), intent(in) :: ac_q
        real(sp), intent(out) :: qup

        integer :: i, row_imd, col_imd, k
        integer, dimension(8) :: drow = (/1, 1, 0, -1, -1, -1, 0, 1/)
        integer, dimension(8) :: dcol = (/0, -1, -1, -1, 0, 1, 1, 1/)

        qup = 0._sp

        do i = 1, 8

            row_imd = row + drow(i)
            col_imd = col + dcol(i)

            if (row_imd .lt. 1 .or. row_imd .gt. mesh%nrow .or. col_imd .lt. 1 .or. col_imd .gt. mesh%ncol) cycle
            k = mesh%rowcol_to_ind_ac(row_imd, col_imd)

            if (mesh%flwdir(row_imd, col_imd) .eq. i) qup = qup + ac_q(k)

        end do

    end subroutine upstream_discharge
    
    subroutine hydraulics_discontinuities(setup, mesh, time_step, row, col, k, ac_hd, ac_qz)
        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        integer, intent(in) :: row
        integer, intent(in) :: col
        integer, intent(in) :: k
        integer, intent(in) :: time_step
        real(sp), dimension(mesh%nac), intent(inout) :: ac_hd
        real(sp), dimension(mesh%nac, setup%nqz), intent(inout) :: ac_qz
        
        integer ::  index_dam, index_input_q
        
        !Hydraulics discontiuities
        select case (mesh%hydraulics_discontinuities%discontinuities_code(row, col))
        
        case(1)
            index_dam=mesh%hydraulics_discontinuities%discontinuities_rank(row,col)
            
            call laminage(setup%dt, mesh%hydraulics_discontinuities%dam_hv(index_dam,:,:),&
            &mesh%hydraulics_discontinuities%dam_hq(index_dam,:,:), ac_qz(k, setup%nqz), ac_hd(k), ac_qz(k, setup%nqz))
        
        case(2)
            index_input_q=mesh%hydraulics_discontinuities%discontinuities_rank(row,col)
            ac_qz(k, setup%nqz)=ac_qz(k, setup%nqz) + mesh%hydraulics_discontinuities%input_q(index_input_q,time_step)
            
        end select
        
    end subroutine hydraulics_discontinuities

    subroutine linear_routing(dx, dy, dt, flwacc, llr, hlr, qup, q)

        implicit none

        real(sp), intent(in) :: dx, dy, dt, flwacc
        real(sp), intent(in) :: llr
        real(sp), intent(inout) :: hlr, qup, q

        real(sp) :: hlr_imd

        qup = (qup*dt)/(1e-3_sp*(flwacc - dx*dy))

        hlr_imd = hlr + qup

        hlr = hlr_imd*exp(-dt/(llr*60._sp))

        q = q + (hlr_imd - hlr)*1e-3_sp*(flwacc - dx*dy)/dt

    end subroutine linear_routing

    subroutine kinematic_wave1d(dx, dy, dt, akw, bkw, qlijm1, qlij, qim1j, qijm1, qij)

        implicit none

        real(sp), intent(in) :: dx, dy, dt
        real(sp), intent(in) :: akw, bkw
        real(sp), intent(in) :: qlijm1, qlij, qim1j, qijm1
        real(sp), intent(inout) :: qij

        real(sp) :: wqlijm1, wqlij, wqim1j, wqijm1
        real(sp) :: dtddx, n1, n2, n3, d1, d2, rhs, rsd, rsd_d
        integer :: iter, maxiter

        !% Avoid numerical issues
        wqlijm1 = max(1e-6_sp, qlijm1)
        wqlij = max(1e-6_sp, qlij)
        wqim1j = max(1e-6_sp, qim1j)
        wqijm1 = max(1e-6_sp, qijm1)

        dtddx = dt/dx

        d1 = dtddx
        d2 = akw*bkw*((wqijm1 + wqim1j)/2._sp)**(bkw - 1._sp)

        n1 = dtddx*wqim1j
        n2 = wqijm1*d2
        n3 = dtddx*(wqlijm1 + wqlij)/2._sp

        !% Linearized solution
        qij = (n1 + n2 + n3)/(d1 + d2)

        !% Non-Linear solution solved with Newton-Raphson
        !% Commented while testing Linearized solution

!~         rhs = n1 + akw*wqijm1**bkw + n3

!~         iter = 0
!~         maxiter = 2
!~         rsd = 1._sp

!~         do while (abs(rsd) > 1e-6 .and. iter < maxiter)

!~             rsd = dtddx*qij + akw*qij**bkw - rhs
!~             rsd_d = dtddx + akw*bkw*qij**(bkw - 1._sp)

!~             qij = qij - rsd/rsd_d

!~             qij = max(qij, 0._sp)

!~             iter = iter + 1

!~         end do

    end subroutine kinematic_wave1d

    subroutine lag0_time_step(setup, mesh, options, returns, time_step, ac_qtz, ac_qz)

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(OptionsDT), intent(in) :: options
        type(ReturnsDT), intent(inout) :: returns
        integer, intent(in) :: time_step
        real(sp), dimension(mesh%nac, setup%nqz), intent(in) :: ac_qtz
        real(sp), dimension(mesh%nac, setup%nqz), intent(inout) :: ac_qz

        integer :: i, j, row, col, k, time_step_returns
        real(sp) :: qup

        ac_qz(:, setup%nqz) = ac_qtz(:, setup%nqz)

        ! Skip the first partition because boundary cells are not routed
        do i = 2, mesh%npar

            ! Tapenade does not accept 'IF' condition within OMP directive. Therefore, the routing loop
            ! is duplicated ... Maybe there is another way to do it.
            if (mesh%ncpar(i) .ge. options%comm%ncpu) then
#ifdef _OPENMP
                !$OMP parallel do schedule(static) num_threads(options%comm%ncpu) &
                !$OMP& shared(setup, mesh, ac_qtz, ac_qz, i) &
                !$OMP& private(j, row, col, k, qup)
#endif
                do j = 1, mesh%ncpar(i)

                    row = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 1)
                    col = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 2)
                    k = mesh%rowcol_to_ind_ac(row, col)

                    if (mesh%active_cell(row, col) .eq. 0 .or. mesh%local_active_cell(row, col) .eq. 0) cycle

                    call upstream_discharge(mesh, row, col, ac_qz(:, setup%nqz), qup)

                    ac_qz(k, setup%nqz) = ac_qz(k, setup%nqz) + qup

                    !$AD start-exclude
                    !internal fluxes
                    if (returns%internal_fluxes_flag) then
                        if (allocated(returns%mask_time_step)) then
                            if (returns%mask_time_step(time_step)) then
                                time_step_returns = returns%time_step_to_returns_time_step(time_step)
                                returns%internal_fluxes( &
                                    row, &
                                    col, &
                                    time_step_returns, &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + 1: &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + setup%n_routing_fluxes &
                                    ) = (/qup/)
                            end if
                        end if
                    end if
                    !$AD end-exclude

                end do
#ifdef _OPENMP
                !$OMP end parallel do
#endif
            else

                do j = 1, mesh%ncpar(i)

                    row = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 1)
                    col = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 2)
                    k = mesh%rowcol_to_ind_ac(row, col)

                    if (mesh%active_cell(row, col) .eq. 0 .or. mesh%local_active_cell(row, col) .eq. 0) cycle

                    call upstream_discharge(mesh, row, col, ac_qz(:, setup%nqz), qup)

                    ac_qz(k, setup%nqz) = ac_qz(k, setup%nqz) + qup

                    !$AD start-exclude
                    !internal fluxes
                    if (returns%internal_fluxes_flag) then
                        if (allocated(returns%mask_time_step)) then
                            if (returns%mask_time_step(time_step)) then
                                time_step_returns = returns%time_step_to_returns_time_step(time_step)
                                returns%internal_fluxes( &
                                    row, &
                                    col, &
                                    time_step_returns, &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + 1: &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + setup%n_routing_fluxes &
                                    ) = (/qup/)
                            end if
                        end if
                    end if
                    !$AD end-exclude

                end do

            end if

        end do

    end subroutine lag0_time_step
    
    subroutine lr_time_step(setup, mesh, options, returns, time_step, ac_qtz, ac_llr, ac_hlr, ac_hd, ac_qz)

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(OptionsDT), intent(in) :: options
        type(ReturnsDT), intent(inout) :: returns
        integer, intent(in) :: time_step
        real(sp), dimension(mesh%nac, setup%nqz), intent(in) :: ac_qtz
        real(sp), dimension(mesh%nac), intent(in) :: ac_llr
        real(sp), dimension(mesh%nac), intent(inout) :: ac_hlr
        real(sp), dimension(mesh%nac), intent(inout) :: ac_hd
        real(sp), dimension(mesh%nac, setup%nqz), intent(inout) :: ac_qz

        integer :: i, j, row, col, k, time_step_returns, index_dam, index_input_q
        real(sp) :: qup, qr

        ac_qz(:, setup%nqz) = ac_qtz(:, setup%nqz)
        
        ! Skip the first partition because boundary cells are not routed
        do i = 2, mesh%npar

            ! Tapenade does not accept 'IF' condition within OMP directive. Therefore, the routing loop
            ! is duplicated ... Maybe there is another way to do it.
            if (mesh%ncpar(i) .ge. options%comm%ncpu) then
#ifdef _OPENMP
                !$OMP parallel do schedule(static) num_threads(options%comm%ncpu) &
                !$OMP& shared(setup, mesh, ac_qtz, ac_llr, ac_hlr, ac_hd, ac_qz, i) &
                !$OMP& private(j, row, col, k, qup)
#endif
                do j = 1, mesh%ncpar(i)

                    row = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 1)
                    col = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 2)
                    k = mesh%rowcol_to_ind_ac(row, col)

                    if (mesh%active_cell(row, col) .eq. 0 .or. mesh%local_active_cell(row, col) .eq. 0) cycle

                    call upstream_discharge(mesh, row, col, ac_qz(:, setup%nqz), qup)

                    call linear_routing(mesh%dx(row, col), mesh%dy(row, col), setup%dt, mesh%flwacc(row, col), &
                    & ac_llr(k), ac_hlr(k), qup, ac_qz(k, setup%nqz))
                    
                    !Hydraulics discontiuities
                    call hydraulics_discontinuities(setup, mesh, time_step, row, col, k, ac_hd, ac_qz)
                    
                    !$AD start-exclude
                    !internal fluxes
                    if (returns%internal_fluxes_flag) then
                        if (allocated(returns%mask_time_step)) then
                            if (returns%mask_time_step(time_step)) then
                                time_step_returns = returns%time_step_to_returns_time_step(time_step)
                                returns%internal_fluxes( &
                                    row, &
                                    col, &
                                    time_step_returns, &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + 1: &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + setup%n_routing_fluxes &
                                    ) = (/qup/)
                            end if
                        end if
                    end if
                    !$AD end-exclude

                end do
#ifdef _OPENMP
                !$OMP end parallel do
#endif
            else

                do j = 1, mesh%ncpar(i)

                    row = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 1)
                    col = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 2)
                    k = mesh%rowcol_to_ind_ac(row, col)

                    if (mesh%active_cell(row, col) .eq. 0 .or. mesh%local_active_cell(row, col) .eq. 0) cycle

                    call upstream_discharge(mesh, row, col, ac_qz(:, setup%nqz), qup)

                    call linear_routing(mesh%dx(row, col), mesh%dy(row, col), setup%dt, mesh%flwacc(row, col), &
                    & ac_llr(k), ac_hlr(k), qup, ac_qz(k, setup%nqz))
                    
                    !Hydraulics discontiuities
                    call hydraulics_discontinuities(setup, mesh, time_step, row, col, k, ac_hd, ac_qz)

                    !$AD start-exclude
                    !internal fluxes
                    if (returns%internal_fluxes_flag) then
                        if (allocated(returns%mask_time_step)) then
                            if (returns%mask_time_step(time_step)) then
                                time_step_returns = returns%time_step_to_returns_time_step(time_step)
                                returns%internal_fluxes( &
                                    row, &
                                    col, &
                                    time_step_returns, &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + 1: &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + setup%n_routing_fluxes &
                                    ) = (/qup/)
                            end if
                        end if
                    end if
                    !$AD end-exclude

                end do

            end if

        end do

    end subroutine lr_time_step

    subroutine kw_time_step(setup, mesh, options, returns, t, ac_qtz, ac_akw, ac_bkw, ac_qz)
        ! Bug fixed on the parallel adjoint linked to variable names by changing of variable names :
        ! Added local variable t instead of time_step as input to kw_time_step -> no change
        ! Added local variable t_returns instead of time_step_returns to kw_time_step
        ! -> modified openMP, unmodified non-openMP

        implicit none

        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        type(OptionsDT), intent(in) :: options
        type(ReturnsDT), intent(inout) :: returns
        integer, intent(in) :: t
        real(sp), dimension(mesh%nac, setup%nqz), intent(in) :: ac_qtz
        real(sp), dimension(mesh%nac), intent(in) :: ac_akw, ac_bkw
        real(sp), dimension(mesh%nac, setup%nqz), intent(inout) :: ac_qz

        integer :: i, j, row, col, k, t_returns
        real(sp) :: qlijm1, qlij, qim1j, qijm1

        ac_qz(:, setup%nqz) = ac_qtz(:, setup%nqz)

        ! Skip the first partition because boundary cells are not routed
        do i = 2, mesh%npar

            ! Tapenade does not accept 'IF' condition within OMP directive. Therefore, the routing loop
            ! is duplicated ... Maybe there is another way to do it.
            if (mesh%ncpar(i) .ge. options%comm%ncpu) then
#ifdef _OPENMP
                !$OMP parallel do schedule(static) num_threads(options%comm%ncpu) &
                !$OMP& shared(setup, mesh, ac_qtz, ac_akw, ac_bkw, ac_qz, i) &
                !$OMP& private(j, row, col, k, qlijm1, qlij, qim1j, qijm1)
#endif
                do j = 1, mesh%ncpar(i)

                    row = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 1)
                    col = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 2)
                    k = mesh%rowcol_to_ind_ac(row, col)

                    if (mesh%active_cell(row, col) .eq. 0 .or. mesh%local_active_cell(row, col) .eq. 0) cycle

                    qlijm1 = ac_qtz(k, setup%nqz - 1)
                    qlij = ac_qtz(k, setup%nqz)
                    qijm1 = ac_qz(k, setup%nqz - 1)

                    call upstream_discharge(mesh, row, col, ac_qz(:, setup%nqz), qim1j)

                    call kinematic_wave1d(mesh%dx(row, col), mesh%dy(row, col), setup%dt, &
                    & ac_akw(k), ac_bkw(k), qlijm1, qlij, qim1j, qijm1, ac_qz(k, setup%nqz))

                    !$AD start-exclude
                    !internal fluxes
                    if (returns%internal_fluxes_flag) then
                        if (allocated(returns%mask_time_step)) then
                            if (returns%mask_time_step(t)) then
                                t_returns = returns%time_step_to_returns_time_step(t)
                                returns%internal_fluxes( &
                                    row, &
                                    col, &
                                    t_returns, &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + 1: &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + setup%n_routing_fluxes &
                                    ) = (/qim1j/)
                            end if
                        end if
                    end if
                    !$AD end-exclude
                end do
#ifdef _OPENMP
                !$OMP end parallel do
#endif
            else

                do j = 1, mesh%ncpar(i)

                    row = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 1)
                    col = mesh%cpar_to_rowcol(mesh%cscpar(i) + j, 2)
                    k = mesh%rowcol_to_ind_ac(row, col)

                    if (mesh%active_cell(row, col) .eq. 0 .or. mesh%local_active_cell(row, col) .eq. 0) cycle

                    qlijm1 = ac_qtz(k, setup%nqz - 1)
                    qlij = ac_qtz(k, setup%nqz)
                    qijm1 = ac_qz(k, setup%nqz - 1)

                    call upstream_discharge(mesh, row, col, ac_qz(:, setup%nqz), qim1j)

                    call kinematic_wave1d(mesh%dx(row, col), mesh%dy(row, col), setup%dt, &
                    & ac_akw(k), ac_bkw(k), qlijm1, qlij, qim1j, qijm1, ac_qz(k, setup%nqz))

                    !$AD start-exclude
                    !internal fluxes
                    if (returns%internal_fluxes_flag) then
                        if (allocated(returns%mask_time_step)) then
                            if (returns%mask_time_step(t)) then
                                t_returns = returns%time_step_to_returns_time_step(t)
                                returns%internal_fluxes( &
                                    row, &
                                    col, &
                                    t_returns, &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + 1: &
                                    setup%n_snow_fluxes + setup%n_hydro_fluxes + setup%n_routing_fluxes &
                                    ) = (/qim1j/)
                            end if
                        end if
                    end if
                    !$AD end-exclude
                end do

            end if

        end do

    end subroutine kw_time_step
    
    !*************************************************************************
    !      PROCEDURE DE LAMINAGE
    !      npdt = nombre de pas de temps de la chronique de débit entrant
    !      dt = pas de temps de la chronique de débit entrant (secondes)
    !      n_HV = nombre de points de la relation hauteur/volume (m - Mm3)
    !      n_HQ = nombre de points de la relation hauteur/débit (m - m3/s)
    !      rel_HV = relation hauteur/volume
    !      rel_HQ = relation hauteur/débit
    !      x1 = debit entrant (m3/s)
    !      x2 = cote du plan d'eau (m)
    !      x3 = debit sortant (m3/s)
    !*************************************************************************
    subroutine laminage(dt, rel_hv, rel_hq, x1, x2, x3)

        implicit none

        real, intent(in)    :: dt
        real, dimension(:,:), intent(in) :: rel_hv
        real, dimension(:,:), intent(in) :: rel_hq
        real, intent(inout)    :: x1, x2, x3
        integer :: i, n_hv,n_hq
        real zini, z, qbid, volout, volin, volt, zoutmax, qoutmoy
        logical :: dt_min
        
        n_hv=size(rel_hv,dim=2)
        n_hq=size(rel_hq,dim=2)
        z=x2
        qbid=0.
        qoutmoy = 0.
        
        !cdl hv ou hq ?? Peut-on avoir des débits nulles dans smash ?
        if (z .lt. rel_hv(1,1)) then
            z=rel_hv(1,1)
        end if

        do i=1,int(dt/60.)             ! boucle sur les pas de temps en minutes
        
            call h_to_v(n_hv, rel_hv, z,volt)
            
            volin = x1* 60. / 1000000.  ! transformation des m3/s en mm3 par min
            volt = volt + volin
            
            call v_to_h(n_hv, rel_hv, volt, z)
            call h_to_q(n_hq, rel_hq, z, qbid)
            
            volout = qbid* 60. / 1000000.
            volt = volt - volout
            
            call v_to_h(n_hv, rel_hv, volt, z)
            
            qoutmoy = qoutmoy + qbid
            
        enddo ! fin de la boucle sur l'heure
        x2 = z
        x3 = qoutmoy/(dt/60.)
        
!~         else ! on n'est pas dans un evenement pluvieux donc on reste au pas de temps horaire
        
!~             VOLin = x1 / 1000000. * dt ! Transformation des m3/s de l'hydrogramme en Mm3 par pas de temps
!~             VOLt = VOLt + VOLin
           
!~             CALL V_to_H(n_HV, rel_HV, VOLt, Z)
!~             CALL H_to_Q(n_HQ, rel_HQ, Z, Qbid)
           
!~             VOLout= Qbid / 1000000. * dt
!~             VOLt = VOLt - VOLout
           
!~             CALL V_to_H(n_HV, rel_HV, VOLt, Z)
           
!~             x2 = Z
!~             x3 = Qbid
           
!~         endif

!~         Zini = x2 
        
        contains
        
        subroutine v_to_h(nv, rel_hv, volume, hauteur)

            implicit none

            integer, intent(in) :: nv
            real, intent(in)  :: volume
            real, intent(out) :: hauteur
            real, dimension(2,nv) :: rel_hv
            integer :: ii

            !.... rel_hv(1,) == hauteur et rel_hv(2,) == volume      
            do ii = 1,nv-1
            if ((volume.ge.rel_hv(2,ii)).and.(volume.le.rel_hv(2,ii+1))) then
                hauteur = rel_hv(1,ii) + (volume-rel_hv(2,ii))*(rel_hv(1,ii+1)-rel_hv(1,ii))/(rel_hv(2,ii+1)-rel_hv(2,ii))
                exit
             endif
            enddo

        end subroutine v_to_h
          
        !**************************************************************************
        !.....Calcul de V(Mm3) pour un H(m) donné
        !**************************************************************************
        subroutine h_to_v(nv, rel_hv, hauteur,volume)

            implicit none

            integer, intent(in) :: nv
            real, intent(in)  :: hauteur
            real, intent(out) :: volume
            real, dimension(2,nv) :: rel_hv
            integer :: ii

            do ii = 1,nv-1
             if ((hauteur.ge.rel_hv(1,ii)).and.(hauteur.le.rel_hv(1,ii+1))) then
                volume = rel_hv(2,ii) + (hauteur-rel_hv(1,ii))*(rel_hv(2,ii+1)-rel_hv(2,ii))/(rel_hv(1,ii+1)-rel_hv(1,ii))
                exit
             endif
            enddo

        end subroutine h_to_v

        !**************************************************************************
        !.....Calcul de Q(m3/s) sortant pour un H(m) donné 
        !**************************************************************************
        subroutine h_to_q(nq, rel_hq, hauteur, qsortant)

            implicit none

            integer, intent(in) :: nq
            real, intent(in)  :: hauteur
            real, intent(out) :: qsortant
            real, dimension(2,nq) :: rel_hq
            integer :: ii

            do ii = 1,nq-1
             if ((hauteur.ge.rel_hq(1,ii)).and.(hauteur.le.rel_hq(1,ii+1))) then
                qsortant = rel_hq(2,ii) + (hauteur-rel_hq(1,ii))*(rel_hq(2,ii+1)-rel_hq(2,ii))/(rel_hq(1,ii+1)-rel_hq(1,ii))
                exit
             endif
            enddo

        end subroutine h_to_q

      
    END SUBROUTINE laminage

end module md_routing_operator
