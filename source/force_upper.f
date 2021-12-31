      subroutine force_upper(alf)
      dimension alf(80)
      character*4 alfchar(80)
      character*1 cchar, char_new
c      equivalence (alf(1), alfchar(1))
c
c     force alf to uppercase -
      do jkl = 1,80
         l = jkl
         write(cchar,'(a1)') alf(l)
         ialph = ichar('a')
         i=ichar(cchar)-ialph+1
c
c        force uppercase
c
         if ( i .ge. 33 .and. i .lt. 59 ) then
             i = i - 32
            letter=i
            char_new = char(letter+ialph-1)
            read (char_new,'(a4)' ) alf(l)
         endif
      enddo
      return
      end
