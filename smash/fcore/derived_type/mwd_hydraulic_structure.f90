module mwd_hydraulic_structure
    
    use md_constant !% only: sp, lchar
    use mwd_setup !% only: SetupDT
    use mwd_mesh !% only: MeshDT
    
    type Hydraulic_StructureDT
        
        real(sp), dimension(:,:,:), allocatable :: dam_hv
        real(sp), dimension(:,:,:), allocatable :: dam_hq
        real(sp), dimension(:,:), allocatable :: inflow
    
    end type Hydraulic_StructureDT
    
    contains
    
    subroutine Hydraulic_StructureDT_initialise(this, setup, mesh)

        !% Notes
        !% ----- 
        !% discontinuitiesDT initialisation subroutine

        implicit none

        type(Hydraulic_StructureDT), intent(inout) :: this
        type(SetupDT), intent(in) :: setup
        type(MeshDT), intent(in) :: mesh
        
        integer :: nmax_val
        
        nmax_val=100
        
        allocate(this%dam_hv(mesh%ndam,2,nmax_val))
        allocate(this%dam_hq(mesh%ndam,2,nmax_val))
        allocate(this%inflow(mesh%ninflow,setup%ntime_step))
        
        this%dam_hv=-99.
        this%dam_hq=-99.
        this%inflow=-99.

    end subroutine Hydraulic_StructureDT_initialise
    
    subroutine Hydraulic_StructureDT_reallocate(this, key, ndam, ninflow, ntime_step, nmax_val)

        implicit none

        type(Hydraulic_StructureDT), intent(inout) :: this
        character(lchar), intent(in) :: key
        integer, intent(in) :: ndam, ninflow, ntime_step, nmax_val

        select case(key)
        
        case("dam_hv")
            if (allocated(this%dam_hv)) then
                deallocate(this%dam_hv)
                allocate(this%dam_hv(ndam, 2, nmax_val))
            end if
            this%dam_hv=-99.
        
        case("dam_hq")
            if (allocated(this%dam_hq)) then
                deallocate(this%dam_hq)
                allocate(this%dam_hq(ndam, 2, nmax_val))
            end if
            this%dam_hq=-99.
        
        end select

    end subroutine Hydraulic_StructureDT_reallocate

    subroutine Hydraulic_StructureDT_copy(this, this_copy)

        implicit none

        type(Hydraulic_StructureDT), intent(in) :: this
        type(Hydraulic_StructureDT), intent(out) :: this_copy

        this_copy = this

    end subroutine Hydraulic_StructureDT_copy


end module mwd_hydraulic_structure
