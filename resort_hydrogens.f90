! name:    resort_hydrogens.f90
! author:  nbehrnd@yahoo.com
! date:    [2024-12-04 Wed]
! edit:    [2024-12-05 Thu]
! license: GPL version 2.

! This Fortran program edits the sequence of atoms in xyz files such
! that hydrogen atoms follow directly their heavy/non-hydrogen atoms
! they might bind to.  For this, it is presumed a distance of less
! than 1.7 Angstrom between H and non-H atom alone is a sufficient
! criterion.  For comparison: a typical C(sp3)-H bond in an alkane
! is about 1.1 A long, the one in hydrogen iodide about 1.6 A.
!
! Internally, atoms are organized as either `hydrogens`, or
! `non-hydrogens`.  Because subsequent tests probe atoms of one set
! with atoms of the other, this implementation is not suitable for
! data with molecular hydrogen/H-H bonds.
!
! Compiles without warning with gfortran 14.2.0 by e.g., either
!
! ```bash
! gfortran resort_hydrogens.f90 -o exe -Wall --std=f2018
! gfortran resort_hydrogens.f90 -o exe -Wall --std=f2023
! ```

program resort
   use iso_fortran_env, only: dp => real64, wp => int32, &
      stdin => input_unit, stdout => output_unit, stderr => error_unit
   implicit none

   character(len=200) :: file_in, file_out
   character(len=3)   :: save_option
   logical :: save_file = .false.

   integer(wp) :: i, j, k, new_unit, error
   integer(wp) :: number_of_atoms  ! limit of int32(huge): 2147483647
   character(len=200) :: title_line
   character(len=3) :: atom_label
   real(dp) :: coordinates(3)

   type atom
      character(len=3) :: atom_label
      real(dp) :: coordinates(3)
   end type atom

   type(atom), allocatable :: hydrogens(:), non_hydrogens(:), output(:)

   ! consider X-H distances shorter than 1.7 A as a bond
   real(dp), parameter :: threshold = 1.7_dp

   ! get in touch with the user
   if (command_argument_count() == 0) then
      ! data stream piped from the CLI
      new_unit = stdin
   else if (command_argument_count() == 1) then
      ! read a file as sole argument
      new_unit = 10
      call get_command_argument(1, file_in)

      if ((trim(file_in) == "-h") .or. (trim(file_in) == "--help")) then
         print *, "Either read the xyz file as an argument, e.g."
         print *, ""
         print *, "    ./exe input.xyz [-s]"
         print *, ""
         print *, "with an optional save (`-s`) as `input.xyz_resort.xyz`, or"
         print *, "pipe a structure via std to the executable (`exe` below) like"
         print *, ""
         print *, "    cat input.xyz | ./exe | obabel -ixyz -omol2"
         stop

      else
         open (new_unit, file=file_in, status="old", action="read", iostat=error)
         if (error /= 0) stop "The input xyz file '" // trim(file_in) // "' is inaccessible."
      end if
   else if (command_argument_count() == 2) then
      ! read a file as an argument with an automatic save by -s
      new_unit = 10
      call get_command_argument(1, file_in)
      call get_command_argument(2, save_option)
      open (new_unit, file=file_in, status="old", action="read", iostat=error)
      if (error /= 0) stop "The input xyz file '" // trim(file_in) // "' is inaccessible."
      if (save_option == "-s") then
          file_out = trim(file_in)//"_resort.xyz"
          save_file = .true.
      end if
   end if

   do i = 1, 2
      if (i == 1) then
         read (new_unit, *, iostat=error) number_of_atoms
         if (error /= 0) stop "Error reading number of atoms."
      else if (i == 2) then
         read (new_unit, "(A)", iostat=error) title_line
         if (error /= 0) stop "Error reading title line."
      end if
   end do

   allocate (non_hydrogens(number_of_atoms))
   allocate (hydrogens(number_of_atoms))
   allocate (output(number_of_atoms))

   non_hydrogens%atom_label = "X"
   hydrogens%atom_label = "X"

   i = 1
   do
      read (new_unit, *, iostat=error) atom_label, coordinates
      if (error /= 0) exit  ! e.g., at the end of the input file read

      if (atom_label == "H") then
         hydrogens(i)%atom_label = atom_label
         hydrogens(i)%coordinates = coordinates
      else
         non_hydrogens(i)%atom_label = atom_label
         non_hydrogens(i)%coordinates = coordinates
      end if
      i = i + 1
   end do

   close (new_unit)

   ! resort the hydrogens to follow suite their corresponding non-H atom
   !
   ! Every record not labeled by `X` participates in the checks below.
   ! `i` and `j` iterate over the `non_hydrogens` or `hydrogens` -- both
   ! arrays still indicate the place holder atoms of the initialization
   ! by `X` -- while `k` helps to manage the array of rearranged atoms.

   k = 1
   outer: do i = 1, size(non_hydrogens)
      if (non_hydrogens(i)%atom_label == "X") then
         cycle outer
      end if
      output(k)%atom_label = non_hydrogens(i)%atom_label
      output(k)%coordinates = non_hydrogens(i)%coordinates
      k = k + 1

      ! for this non-H atom `i`, identify suitable `j` hydrogen atoms
      inner: do j = 1, size(hydrogens)
         if (hydrogens(j)%atom_label == "X") then
            cycle inner
         end if

         if (distance(non_hydrogens(i)%coordinates, hydrogens(j)%coordinates) < threshold) then
            output(k)%atom_label = hydrogens(j)%atom_label
            output(k)%coordinates = hydrogens(j)%coordinates
            k = k + 1

            ! in subsequent cycles, do not test this very hydrogen atom
            hydrogens(j)%atom_label = "X"

            ! print *, non_hydrogens(i)%atom_label, hydrogens(j)%atom_label, &
            ! distance(non_hydrogens(i)%coordinates, hydrogens(j)%coordinates

         end if
      end do inner
   end do outer

   deallocate (non_hydrogens); deallocate (hydrogens)


   ! report the results to either file, or CLI
   !
   ! Because OpenBabel explicitly states the `10.5` format as the one adopted
   ! by them, explicit white space is used to provide "the same look".
   ! <https://github.com/openbabel/openbabel>
   ! <https://open-babel.readthedocs.io/en/latest/FileFormats/XYZ_cartesian_coordinates_format.html>

   if (save_file .eqv. .true.) then
      open (20, file=file_out, action="write", status="replace", iostat=error)
      if (error /= 0) stop "Unable to provide a permanent record."
      write (20, "(I0)") number_of_atoms
      write (20, "(A)") trim(title_line)
      do i = 1, size(output)
         write (20, "(A, 3(F10.5, 5x))") output(i)
      end do
      close (20)
   else
      write (*, "(I0)") number_of_atoms
      write (*, "(A)") trim(title_line)
      do i = 1, size(output)
         write (*, "(A, 3(F10.5, 5x))") output(i)
      end do
   end if

contains

   pure function distance(a, b)
      ! calculate the distance between points `a(x,y,z)` and `b(x,y,z)`
      real(dp), intent(in) :: a(3), b(3)
      real(dp) :: distance
      distance = sqrt( &
         (a(1) - b(1))**2 + &
         (a(2) - b(2))**2 + &
         (a(3) - b(3))**2)
   end function distance

end program resort
