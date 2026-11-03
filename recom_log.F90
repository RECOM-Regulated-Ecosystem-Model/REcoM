!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  REcoM  |  recom_log.F90
!  Signed, quietly coloured terminal output for REcoM, so REcoM lines stand out
!  from FESOM's own (blue/cyan) output.
!
!  [REcoM] ── REcoM init ──────────────────────────────────────────────
!  [REcoM]    BGC tracers               33
!  [REcoM]    DIN                    0.010 ..      43.200   mmol m-3
!  [REcoM]    ✓ initialisation done
!  [REcoM] ── step 42 ─────────────────────────────────────────────────
!  [REcoM] REcoM_Forcing
!  [REcoM] ├─ REcoM_sms
!  [REcoM] └─ ciso
!  [REcoM] ▲ WARNING  DIC close to upper bound   n=10342
!
!  Usage:
!     use recom_log, only: recom_section, recom_trace, recom_kv, recom_range, &
!                          recom_ok, recom_info, recom_warn, recom_error
!     if (mype == 0) call recom_section('REcoM init')
!     if (mype == 0) call recom_kv('BGC tracers', bgc_num)
!     if (mype == 0) call recom_kv('coccos', enable_coccos)
!     if (mype == 0) call recom_range('DIN', glo_min, glo_max, 'mmol m-3')
!     if (mype == 0) call recom_ok('initialisation done')
!     if (recom_debug .and. mype == 0) call recom_trace('REcoM_Forcing')
!     if (recom_debug .and. mype == 0) call recom_trace('REcoM_sms', 2)
!     if (recom_debug .and. mype == 0) call recom_trace('ciso', 2, last=.true.)
!     call recom_warn('DIC close to upper bound', 'n=10342')
!
!  Switches (e.g. set from a namelist):
!     recom_use_color   = .false.  -> no ANSI codes (clean log files)
!     recom_use_unicode = .false.  -> ASCII only:  --  |-  `-  |  ok  !  x
!
!  Note: Fortran counts bytes, not characters. The box-drawing and status
!  glyphs are multi-byte UTF-8, so they only appear in fixed prefixes and in
!  repeated rules, never in strings that get padded or measured with len().
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module recom_log
    ! Deliberately independent of recom_modules.F90: recom_config uses this module, so
    ! taking wp from recom_declarations (same source file) would make the two source
    ! files depend on each other and break the build order.
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    private

    logical, public :: recom_use_color   = .true.
    logical, public :: recom_use_unicode = .true.

    !---------------------------------------------------------------------------
    ! Palette (256-colour ANSI codes, without the leading ESC)
    !---------------------------------------------------------------------------
    character(len=*), parameter :: C_RESET = '[0m'
    character(len=*), parameter :: C_TAG   = '[38;5;71m'   ! [REcoM] signature    (mid green)
    character(len=*), parameter :: C_L1    = '[38;5;108m'  ! level-1 calls, labels (sage)
    character(len=*), parameter :: C_L2    = '[38;5;71m'   ! level-2 calls        (mid green)
    character(len=*), parameter :: C_L3    = '[38;5;65m'   ! level-3+ calls       (dark green)
    character(len=*), parameter :: C_DIM   = '[38;5;244m'  ! structure: rules, units, details
    character(len=*), parameter :: C_WARN  = '[38;5;179m'  ! muted amber
    character(len=*), parameter :: C_ERR   = '[38;5;167m'  ! muted red

    !---------------------------------------------------------------------------
    ! Layout
    !---------------------------------------------------------------------------
    integer, parameter :: LINE_WIDTH  = 62  ! visible width of a section line after the tag
    integer, parameter :: LABEL_WIDTH = 18  ! key column of recom_kv / recom_range
    integer, parameter :: VALUE_WIDTH = 12  ! one number

    !---------------------------------------------------------------------------
    ! recom_kv: aligned "label   value   unit" for real, integer, logical, text
    !---------------------------------------------------------------------------
    interface recom_kv
        module procedure kv_real, kv_int, kv_ints, kv_logical, kv_text
    end interface recom_kv

    public :: recom_tag, recom_section, recom_trace, recom_kv, recom_range
    public :: recom_ok, recom_info, recom_warn, recom_error, recom_itoa

contains

    !===========================================================================
    ! Helpers
    !===========================================================================

    ! Colour a piece of text (plain text when colours are off)
    function paint(code, text) result(s)
        character(len=*), intent(in) :: code, text
        character(len=:), allocatable :: s
        if (recom_use_color) then
            s = achar(27) // code // text // achar(27) // C_RESET
        else
            s = text
        end if
    end function paint

    ! Unicode glyph, or its ASCII fallback
    function glyph(uni, ascii) result(s)
        character(len=*), intent(in) :: uni, ascii
        character(len=:), allocatable :: s
        if (recom_use_unicode) then
            s = uni
        else
            s = ascii
        end if
    end function glyph

    ! One real, right-aligned in VALUE_WIDTH columns; exponent notation for
    ! very large or very small magnitudes
    function fmt_real(val) result(s)
        real(kind=wp), intent(in) :: val
        character(len=:), allocatable :: s
        character(len=VALUE_WIDTH) :: buf
        if (val /= 0.0_wp .and. (abs(val) >= 1.0e6_wp .or. abs(val) < 1.0e-3_wp)) then
            write(buf, '(es12.4)') val
        else
            write(buf, '(f12.3)') val
        end if
        s = buf
    end function fmt_real

    ! Left-aligned label, padded (or cut) to LABEL_WIDTH
    function pad_label(label) result(s)
        character(len=*), intent(in) :: label
        character(len=LABEL_WIDTH) :: s
        s = adjustl(label)
    end function pad_label

    ! Integer as a trimmed string, for building detail text: 'n=' // recom_itoa(n)
    function recom_itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=12) :: buf
        write(buf, '(i0)') i
        s = trim(buf)
    end function recom_itoa

    !===========================================================================
    ! The signature
    !===========================================================================
    function recom_tag() result(s)
        character(len=:), allocatable :: s
        s = paint(C_TAG, '[REcoM]')
    end function recom_tag

    !===========================================================================
    ! Section header:  [REcoM] ── title ───────────────────────────────
    !===========================================================================
    subroutine recom_section(title)
        character(len=*), intent(in) :: title
        character(len=:), allocatable :: rule_char
        integer :: n_rule

        rule_char = glyph('─', '-')
        ! visible width = 2 (lead) + 1 + title + 1 + rule = LINE_WIDTH
        n_rule = max(3, LINE_WIDTH - len_trim(title) - 4)

        write(*, '(a)') recom_tag() // ' ' // &
                paint(C_DIM, repeat(rule_char, 2)) // ' ' // &
                paint(C_L1, trim(title)) // ' ' // &
                paint(C_DIM, repeat(rule_char, n_rule))
    end subroutine recom_section

    !===========================================================================
    ! Call tree
    !   level 1 :  REcoM_Forcing
    !   level 2 :  ├─ REcoM_sms         (└─ with last = .true.)
    !   level 3 :  │  ├─ sinking
    !===========================================================================
    subroutine recom_trace(msg, level, last)
        character(len=*), intent(in) :: msg
        integer, intent(in), optional :: level
        logical, intent(in), optional :: last

        integer :: lev
        logical :: is_last
        character(len=:), allocatable :: branch, prefix

        lev = 1
        if (present(level)) lev = max(1, level)
        is_last = .false.
        if (present(last)) is_last = last

        if (is_last) then
            branch = glyph('└─ ', '`- ')
        else
            branch = glyph('├─ ', '|- ')
        end if

        select case (lev)
        case (1)
            write(*, '(a)') recom_tag() // ' ' // paint(C_L1, msg)
        case (2)
            prefix = paint(C_DIM, branch)
            write(*, '(a)') recom_tag() // ' ' // prefix // paint(C_L2, msg)
        case default
            prefix = paint(C_DIM, repeat(glyph('│  ', '|  '), lev - 2) // branch)
            write(*, '(a)') recom_tag() // ' ' // prefix // paint(C_L3, msg)
        end select
    end subroutine recom_trace

    !===========================================================================
    ! Aligned key / value lines
    !   [REcoM]    BGC tracers               33
    !   [REcoM]    DIN                    0.010 ..      43.200   mmol m-3
    !===========================================================================
    subroutine write_kv(label, value, unit)
        character(len=*), intent(in) :: label, value
        character(len=*), intent(in), optional :: unit
        character(len=:), allocatable :: tail

        tail = ''
        if (present(unit)) tail = '   ' // paint(C_DIM, trim(unit))

        write(*, '(a)') recom_tag() // '    ' // paint(C_L1, pad_label(label)) // ' ' // &
                value // tail
    end subroutine write_kv

    subroutine kv_real(label, val, unit)
        character(len=*), intent(in) :: label
        real(kind=wp), intent(in) :: val
        character(len=*), intent(in), optional :: unit
        call write_kv(label, fmt_real(val), unit)
    end subroutine kv_real

    subroutine kv_int(label, val, unit)
        character(len=*), intent(in) :: label
        integer, intent(in) :: val
        character(len=*), intent(in), optional :: unit
        character(len=VALUE_WIDTH) :: buf
        write(buf, '(i12)') val
        call write_kv(label, buf, unit)
    end subroutine kv_int

    ! Integer list (e.g. tracer IDs), eight per line
    subroutine kv_ints(label, vals)
        character(len=*), intent(in) :: label
        integer, intent(in) :: vals(:)
        integer, parameter :: PER_LINE = 8
        character(len=:), allocatable :: line
        character(len=12) :: buf
        integer :: i

        if (size(vals) == 0) then
            call write_kv(label, '  -')
            return
        end if

        line = ''
        do i = 1, size(vals)
            write(buf, '(i0)') vals(i)
            line = line // '  ' // trim(buf)
            if (mod(i, PER_LINE) == 0 .or. i == size(vals)) then
                if (i <= PER_LINE) then
                    call write_kv(label, line)
                else
                    call write_kv('', line)
                end if
                line = ''
            end if
        end do
    end subroutine kv_ints

    subroutine kv_logical(label, val)
        character(len=*), intent(in) :: label
        logical, intent(in) :: val
        character(len=VALUE_WIDTH) :: buf
        if (val) then
            buf = 'on'
            buf = adjustr(buf)
            call write_kv(label, paint(C_L2, buf))
        else
            buf = 'off'
            buf = adjustr(buf)
            call write_kv(label, paint(C_DIM, buf))
        end if
    end subroutine kv_logical

    subroutine kv_text(label, val)
        character(len=*), intent(in) :: label, val
        character(len=VALUE_WIDTH) :: buf
        buf = val
        if (len_trim(val) <= VALUE_WIDTH) then
            buf = adjustr(buf)
            call write_kv(label, buf)
        else
            call write_kv(label, val)
        end if
    end subroutine kv_text

    ! min .. max of a field, e.g. the sanity check after recom_init
    subroutine recom_range(label, vmin, vmax, unit)
        character(len=*), intent(in) :: label
        real(kind=wp), intent(in) :: vmin, vmax
        character(len=*), intent(in), optional :: unit
        call write_kv(label, fmt_real(vmin) // paint(C_DIM, ' ..') // fmt_real(vmax), unit)
    end subroutine recom_range

    !===========================================================================
    ! Status lines
    !   [REcoM]    ✓ initialisation done
    !   [REcoM] ▲ WARNING  message   detail
    !   [REcoM] ✗ ERROR    message   detail
    !===========================================================================
    subroutine recom_ok(msg)
        character(len=*), intent(in) :: msg
        write(*, '(a)') recom_tag() // '    ' // paint(C_L2, glyph('✓', 'ok')) // ' ' // &
                paint(C_L1, msg)
    end subroutine recom_ok

    subroutine recom_info(msg)
        character(len=*), intent(in) :: msg
        write(*, '(a)') recom_tag() // '    ' // msg
    end subroutine recom_info

    subroutine recom_warn(msg, detail)
        character(len=*), intent(in) :: msg
        character(len=*), intent(in), optional :: detail
        character(len=:), allocatable :: tail
        tail = ''
        if (present(detail)) tail = '   ' // paint(C_DIM, detail)
        write(*, '(a)') recom_tag() // ' ' // &
                paint(C_WARN, glyph('▲', '!') // ' WARNING') // '  ' // msg // tail
    end subroutine recom_warn

    subroutine recom_error(msg, detail)
        character(len=*), intent(in) :: msg
        character(len=*), intent(in), optional :: detail
        character(len=:), allocatable :: tail
        tail = ''
        if (present(detail)) tail = '   ' // paint(C_DIM, detail)
        write(*, '(a)') recom_tag() // ' ' // &
                paint(C_ERR, glyph('✗', 'x') // ' ERROR  ') // '  ' // msg // tail
    end subroutine recom_error

end module recom_log
