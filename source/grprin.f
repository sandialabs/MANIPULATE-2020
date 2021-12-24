      subroutine grprin(infile,irxn)
c***********************************************************************
c* this routine reads the energy grid, and cross sections from a 
c* groupr type output tape.  
c***********************************************************************
c* author        : jw vandenburg
c* creation date : 05 march 1992
c* last modified : 27 march 1992
c* organization  : snl ( sandia national laboratories )
c***********************************************************************
c* data is read into the following common arrays.
c* 
c* nenergy ( scalar ) : number of energy bins.
c* energy  ( vector ) : energy bounds of interval i,i+1.
c* emid    ( vector ) : midpoint energy of interval.
c* array   ( array  ) : cross section data set (i,j)
c*                      i - energy index corresponding to emid(i).
c*                      j - location to load reaction ( 1-harwired ).
c***********************************************************************
c234567890
      implicit real (a-h,o-z)
      implicit integer (i-n)
      character infile*(*),fdate*24
      character*145 idir, jdir, kdir
      character*106 optical, optical_old
      character*256  temp_file1, temp_file2
      character*256 njoy2016_dir
      integer le2016
      character*1 cdummy
      common /guide/ icon(40)
      character*15 ename
      common /location / optical, idir, jdir, kdir
      common /io/ nt5, nt6, nfile, nplot, npun
c
c common blocks.
c
      common /datain/ nenergy, energy(1001), array(1000,15)
     1         , emid(1001)
      dimension enghst(2001), xschst(2001)
c
c open the input tape.
c
      if ( icon(9) .lt. 1) then 
         write (6,7823)
 7823    format (1x, '*** Enter GRPRIN ')
      endif
      ename = 'opt'
      lename = lnblnk(ename)
      call getenv(ename(1:lename), optical)
      optical_old = optical
      istyle = 0
      if (optical .eq. '') then 
           optical='/odsk1'
           write (*,9013) optical
 9013      format (1x, 'blank optical default filled', 1x, a6)
      endif
      leopt = lnblnk(optical)
      idir = '/application/njoy-91.118/output/groupr/'
c     redirection for NJOY99
      if ( icon(16) .eq. 4) then 
           optical_old = optical
           optical = ''
           istyle = 0
           leopt = lnblnk(optical)
           idir = '../njoy99l/output/groupr/'
      endif
c     redirection for NJOY97
      if ( icon(16) .eq. 3) then 
           optical_old = optical
           optical = ''
           istyle = 0
           leopt = lnblnk(optical_old)
           idir = optical_old(1:leopt)//'../NJOY2016/groupr/'
           if ( icon(7) .eq. 1) then
               njoy2016_dir = '/mnt/e/sync_sandialabs/NJOY2016/' 
               le2016 = lnblnk(njoy2016_dir)
               idir = njoy2016_dir(1:le2016)//'../NJOY2016/groupr/'
           endif
      endif
      if ( icon(9) .lt. 1) then 
         iposition = 1
         write (6,5823) iposition
 5823    format (1x, '*** GRPRIN check ', i5)
      endif
      if ( icon(16) .eq. 2) then 
           idir = '/application/njoy94/output/groupr/'
      endif
      if ( icon(16) .eq. 1) then 
c         special re-direction for pka groupr files
          idir = '/application/njoy-91.118/output/groupr/'
      endif
      iblank2 = lnblnk(idir)
      iblank = lnblnk(infile)

c      temp_file1 = optical(1:leopt)//idir(1:iblank2)
c     1   //infile(1:iblank)
      temp_file1 = idir(1:iblank2)
     1   //infile(1:iblank)
      temp_file2 = "/esata/codes/app/njoy99l/output/groupr/" //
     &            "ASTM_b10_ASTM_1018_107_xsec_tpl"
      iend1 = lnblnk(temp_file1)
      iend2 = lnblnk(temp_file2)
      write (*,9012) optical(1:leopt)//idir(1:iblank2)
     1   //infile(1:iblank), temp_file1(1:iend1),
     1   temp_file2(1:iend2),
     &   optical(1:leopt), idir(1:iblank2), 
     &   infile(1:iblank)
9012  format (1x, 'opening file: =', a,'=',/,
     &        1x, '  temp_file1: =', a,'=',/,
     &        1x, '  temp_file2: =', a,'=',/,
     &        1x, '     optical: =', a,'=',/,
     &        1x, '        idir: =', a,'=',/,
     &        1x, '      infile: =', a,'=',/)
c      write (6,7846)
7846  format (/,1x, "temporarily override grprin.f line 84 file ")
      open(unit=7,file=temp_file1(1:iend1),
     1   status='unknown')
c      open(unit=7,file=optical(1:leopt)//idir(1:iblank2)
c     1   //infile(1:iblank),
c     1   status='unknown')
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
      if ( istyle .eq. 0) then 
         ird = 4
         read(7,300)dummy1,dummy2,(energy(i),i=1,4)
      else
         ird = 5
         read(7,300)dummy1,(energy(i),i=1,5)
      endif
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
      read(7,400)(energy(i),i=j,j+left)
c
c calculate the midpoint energy an fill that array.
c
      do i = 1,nenergy
        emid(i) = (energy(i)+energy(i+1))/2.0
      end do
c
c now scan the tape for the reaction identifier.
c	
      irxntemp = 0
   10 read(7,600,end=20)irxntemp
      if ( irxntemp.eq.irxn ) then
        go to 30
      endif
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
        if ( irxn2.ne.irxn ) go to 40
        read(7,1900)dummy,array(index,1)
        if (array(index,1).lt.1e-35) array(index,1) = 1.e-35
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
      if ( icon(17) .eq. 1) then 
         read (nt5, *) npoints
c        build histogram for interpolation
         nlast = 0
         do jl = 1,nenergy
           nlast = nlast + 1
           enghst(nlast) = energy(jl)
           xschst(nlast) = array(jl,1)
           nlast = nlast + 1
           enghst(nlast) = energy(jl+1)
           xschst(nlast) = array(jl,1)
         enddo
         if ( npoints .gt. 0) then 
c           set constant histogram interpoaltion
            imethod = 1
            do jk=1,npoints
               read (nt5,*) engpt, vary_eng
               xscpt  = fitmd(engpt, nlast, enghst(1), xschst(1),
     1                 imethod )
               xscpt1 = fitmd(engpt + vary_eng, nlast, enghst(1), 
     1                  xschst(1), imethod  )
               xscpt2 = fitmd(engpt - vary_eng, nlast, enghst(1), 
     1                  xschst(1), imethod  )
               vary1 = abs(xscpt1-xscpt)/xscpt*100.
               vary2 = abs(xscpt2-xscpt)/xscpt*100.
               vary = max(vary1, vary2)
c               write (6,876) engpt, vary_eng, xscpt, xscpt1, 
c     1          xscpt2, vary1, vary2, vary
 876           format (1x, 'debug interog ', 8g14.7)
               write (nt6, 987) engpt, xscpt, vary
 987           format (1x, 'interrogation point = ', g14.7, 
     1         ' MeV, xsec = ', g14.7, ' barn, variation = ',
     2         g14.7, ' %')
            enddo
         endif
      endif
c
c return to the calling routine.
c
      optical = optical_old
      return
c  100 format(2(1pe11.6))
  100 format(2(e11.6))
 1900 format(2(e11.6))
  110 format(63x,i3,6x,i3)
c  200 format(1pe11.6,12x,i10)
c  300 format(6(1pe11.6))
c  400 format(6(1pe11.6),i4,2(i5))
c  500 format(6(1pe12.5))
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
      end
