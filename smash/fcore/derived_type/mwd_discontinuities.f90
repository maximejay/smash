module mwd_discontinuities
    
    use md_constant !% only: lchar
    
    type discontinuitiesDT
        
        character(lchar), dimension(:), allocatable :: discontinuities_name !$F90W char-array
        character(lchar), dimension(:), allocatable :: discontinuities_type !$F90W char-array
        integer, dimension(:,:), allocatable :: discontinuities_pos
        integer, dimension(:,:), allocatable :: discontinuities_rank
        integer, dimension(:,:), allocatable :: discontinuities_code
        integer :: n_dam
        integer :: n_input_q
        
        real(sp), dimension(:,:,:), allocatable :: dam_hv
        real(sp), dimension(:,:,:), allocatable :: dam_hq
        real(sp), dimension(:,:), allocatable :: input_q
    
    end type discontinuitiesDT
    
    
    contains
    
    subroutine discontinuitiesDT_initialise(this, nd, nrow, ncol, ntime_step, ndam, ninput_q, nmax_val)

        !% Notes
        !% ----- 
        !% discontinuitiesDT initialisation subroutine

        implicit none

        type(discontinuitiesDT), intent(inout) :: this
        integer, intent(in) :: nd
        integer, intent(in) :: nrow
        integer, intent(in) :: ncol
        integer, intent(in) :: ntime_step
        integer, intent(in) :: ndam
        integer, intent(in) :: ninput_q
        integer, intent(in) :: nmax_val
        
        allocate(this%discontinuities_name(nd))
        allocate(this%discontinuities_type(nd))
        allocate(this%discontinuities_pos(nd,2))
        
        allocate(this%discontinuities_rank(nrow,ncol))
        allocate(this%discontinuities_code(nrow,ncol))
        
        this%discontinuities_name="..."
        this%discontinuities_type="zero"
        
        this%discontinuities_rank=0
        this%discontinuities_code=0
        
        this%n_dam=ndam
        this%n_input_q=ninput_q
        
        allocate(this%dam_hv(ndam,2,nmax_val))
        allocate(this%dam_hq(ndam,2,nmax_val))
        allocate(this%input_q(ninput_q,ntime_step))
        
        this%dam_hv=-99.
        this%dam_hq=-99.
        this%input_q=-99.

    end subroutine discontinuitiesDT_initialise

    subroutine discontinuitiesDT_copy(this, this_copy)

        implicit none

        type(discontinuitiesDT), intent(in) :: this
        type(discontinuitiesDT), intent(out) :: this_copy

        this_copy = this

    end subroutine discontinuitiesDT_copy

end module mwd_discontinuities
