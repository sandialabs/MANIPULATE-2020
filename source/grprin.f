c***********************************************************************
c* this routine reads the energy grid, and cross sections from a
c* groupr type output tape.
c***********************************************************************
c* author        : jw vandenburg
c* creation date : 05 march 1992
c* last modified : 27 march 1992
c* organization  : snl (sandia national laboratories)
c***********************************************************************
c* data is read into the following common arrays.
c*
c* nenergy (scalar) : number of energy bins.
c* energy  (vector) : energy bounds of interval i,i+1.
c* emid    (vector) : midpoint energy of interval.
c* array   (array ) : cross section data set (i,j)
c*                      i - energy index corresponding to emid(i).
c*                      j - location to load reaction (1-harwired).
c***********************************************************************
      subroutine grprin(infile,irxn)
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "location.cmn"
      character infile*(*),fdate*24
      character*106 optical_old
      character*256 temp_file
      character*256 njoy2016_dir
      character*1 cdummy
      common /guide/ icon(40)
      common /io/ nt5, nt6, nfile, nplot, npun
      common /datain/ nenergy, energy(1001), array(1000,15),
     1                emid(1001)
      dimension enghst(2001), xschst(2001)
c
c open the input tape.
c
      if (icon(9) < 1) then
        write (6,7823)
 7823   format (1x, '*** Enter GRPRIN ')
      end if

      !Default to NJOY2016 by default
      idir = '..'//slash//'NJOY2016'//slash//'groupr'//slash

      call getenv('opt', optical)
      optical_old = optical

      istyle = 0
      if (optical == '') then
        optical='/odsk1'
        write (*,9013) optical
 9013   format (1x, 'blank optical default filled', 1x, a6)
      end if

c     redirection for NJOY99
      if (icon(16) == 4) then
        idir = '..'//slash//'njoy99l'//slash//'output'//slash//
     &         'groupr'//slash
      end if

c     redirection for NJOY2016 in nonstandard directory
      if (icon(7) == -1) then
        optical_old = optical
        optical = ''
        njoy2016_dir = '/mnt/e/sync_sandialabs/NJOY2016/'
        idir = trim(njoy2016_dir)//'groupr/'
      end if
      if (icon(16) == 2) then
        idir = '/application/njoy94/output/groupr/'
      else if (icon(16) == 1) then
c       special re-direction for pka groupr files
        idir = '/application/njoy-91.118/output/groupr/'
      end if

      temp_file = trim(optical)//trim(idir)//trim(infile)
      write (*,9012) trim(optical)//trim(idir) //trim(infile),
     1               trim(temp_file),trim(optical),trim(idir),
     1               trim(infile)
9012  format (1x, 'opening file: =', a,'=',/,
     &        1x, '  temp_file: =', a,'=',/,
     &        1x, '     optical: =', a,'=',/,
     &        1x, '        idir: =', a,'=',/,
     &        1x, '      infile: =', a,'=',/)
      open(unit=7,file=trim(temp_file),status='unknown')
c
c read in the zaid of the material along with the atomic weight.
c
      read(7,9023) cdummy
 9023 format (a1)
      read(7,100)zaid,aweight
      izaid = zaid
c
c read in the temperature and the number of groups.
c
      read(7,200)temp,nenergy
c
c read in the group structure.
c first begin by reading the first line.
c
      if (istyle == 0) then
        ird = 4
        read(7,300)dummy1,dummy2,(energy(i),i=1,4)
      else
        ird = 5
        read(7,300)dummy1,(energy(i),i=1,5)
      end if
c
c calculate how many more rows of full data are left and
c read them in accordingly.
c
      togo = nenergy-(ird-1)
      nlines= int(togo/6.0)
      ncount = ird+1
      do i = 1,nlines
        read(7,400)(energy(j),j=ncount,ncount+5)
        ncount = ncount + 6
      end do
c
c calculate the number of energy bounds left
c in the remaining line and read in the data.
c
      left = nenergy+1-(nlines*6+ird)
      read(7,400)(energy(i),i=ncount,j+left)
c
c calculate the midpoint energy an fill that array.
c
      do i = 1,nenergy
        emid(i) = (energy(i)+energy(i+1))*0.5
      end do
c
c now scan the tape for the reaction identifier.
c
      irxntemp = 0
   10 read(7,600,end=20)irxntemp
      if (irxntemp==irxn) then
        go to 30
      end if
      go to 10
c
c requested reaction not found on tape.
c die with error
c
   20 write(*,*)
      write(*,*)'error : invalid reaction number :',irxn
      write(*,*)'error : program termination (grprin)'
      stop 'invalid-reaction'
      write(6,*)'error : invalid reaction number :',irxn
c
c continue processing the data.
c
c clear the xs array
c
   30 do i = 1,1000
        array(i,1) = 1.0e-35
      end do
c
c read in the cross section data.
c
      flag = 0.0
      backspace(unit=7)
      read(7,650)nend
      read(7,650)nstart
      backspace(unit=7)
      do i = nstart,nend
        read(7,110)index,irxn2
        if (irxn2/=irxn) go to 40
        read(7,1900)dummy,array(index,1)
        if (array(index,1)<1e-35) array(index,1) = 1.e-35
      end do
   40 continue
c
c close and rewind the input tape.
c
      rewind(unit=7)
      close(unit=7)
c
c write the material and problem dependent properties to
c stdio.
c
      write(*,700)izaid,aweight,temp,nenergy,fdate()
c
c write stdio the lower, and upper energy bounds, the midpoint of the`
c energy interval and the cross section data.
c
      write(*,800)irxn
      do i = 1,nenergy
        write(*,900)energy(i),energy(i+1),emid(i),array(i,1)
      end do
c
c   If point interrogation is requested - do interpolation
c
      if (icon(17) == 1) then
        read (nt5, *) npoints
c       build histogram for interpolation
        nlast = 0
        do jl = 1,nenergy
          nlast = nlast + 1
          enghst(nlast) = energy(jl)
          xschst(nlast) = array(jl,1)
          nlast = nlast + 1
          enghst(nlast) = energy(jl+1)
          xschst(nlast) = array(jl,1)
        end do
        if (npoints > 0) then
c         set constant histogram interpoaltion
          imethod = 1
          do jk=1,npoints
            read (nt5,*) engpt, vary_eng
            xscpt  = fitmd(engpt, nlast, enghst(1), xschst(1),
     1                     imethod)
            xscpt1 = fitmd(engpt + vary_eng, nlast, enghst(1),
     1                     xschst(1), imethod )
            xscpt2 = fitmd(engpt - vary_eng, nlast, enghst(1),
     1                     xschst(1), imethod )
            vary1 = abs(xscpt1-xscpt)/xscpt*100.0
            vary2 = abs(xscpt2-xscpt)/xscpt*100.0
            vary = max(vary1, vary2)
 876        format (1x, 'debug interog ', 8g14.7)
            write (nt6, 987) engpt, xscpt, vary
 987        format (1x, 'interrogation point = ', g14.7,
     1      ' MeV, xsec = ', g14.7, ' barn, variation = ',
     2      g14.7, ' %')
          end do
        end if
      end if
c
c return to the calling routine.
c
      optical = optical_old

      return
  100 format(2(e11.6))
 1900 format(2(e11.6))
  110 format(63x,i3,6x,i3)
  200 format(e11.6,12x,i10)
  300 format(6(e11.6))
  400 format(6(e11.6),i4,2(i5))
  500 format(6(e12.5))
  600 format(72x,i3)
  650 format(63x,i3)
  700 format(//,1x,'material zaid              : ',i7,/,
     11x,'material atomic weight     : ',f7.3,/,
     21x,'cross section temperature  : ',f7.3,/,
     31x,'number of energy groups    : ',i7,/,
     41x,'time and date of execution : ',a24,//)
  800 format(7x,'upper',8x,'  lower ',8x,'midpoint',9x,'reaction',/,
     1       7x,'(ev) ',8x,'  (ev)  ',8x,'  (ev)  ',9x,' (',i3,')  ')
  900 format(3(3x,1pe12.5),5x,1pe12.5)

      end subroutine grprin