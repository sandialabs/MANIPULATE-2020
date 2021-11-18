      subroutine integral_pka_grprin(infile,irxn)
c***********************************************************************
c* this routine reads the energy grid, and recoil spectra from a 
c* groupr type pka output tape.  
c***********************************************************************
c* author        : PJ Griffin
c* creation date : 23 November 2013
c* last modified : 23 November 2013
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
      character*256  temp_file1, temp_file2, pkadir
      character*350  temp_file_out
      character*6 reaction_label
      character*15 energy_label
      character*80 line
      character*1 cdummy
      common /guide/ icon(40)
      character*15 ename
      common /location / optical, idir, jdir, kdir
      common /io/ nt5, nt6, nfile, nplot, npun
c
      common /PKA/ ipka_low_energy(1000), ipka_high_energy(1000), 
     &             pka_spectrum(1000,1000), 
     &             pka_spectrum_dede(1000,1000),
     &             robinson_DE(1000)
c
c common blocks.
c
      common /datain/ nenergy, energy(1001), array(1000,15)
     1         , emid(1001)
      dimension enghst(2001), xschst(2001)
c
c open the input tape.
c
c     zero out arrays
      ipka_low_energy(:) = 0
      ipka_high_energy(:) = 0 
      pka_spectrum(:,:) = 0.0 
      pka_spectrum_dede(:,:) = 0.0
      robinson_DE(:) = 0.0
      do i = 1,1000
        array(i,1) = 0.0
      end do
c
      if ( icon(9) .lt. 1) then 
         write (6,7823)
 7823    format (1x, '*** Enter INTEGRAL_PKA_GRPRIN ')
      endif
      ename = 'opt'
      lename = lnblnk(ename)
      call getenv(ename(1:lename), optical)
      optical_old = optical
      istyle = 0
      if (optical .eq. '') then 
           optical='/sync_Project-Git'
           write (*,9013) optical
 9013      format (1x, 'blank optical default filled', 1x, a6)
      endif
      leopt = lnblnk(optical)
      idir = '/NJOY-2012/pka/'
      if ( icon(7) .eq. 1) then 
           idir = '/NJOY-2016/pka/'
      endif
c
c     redirection for NJOY99
      if ( icon(16) .eq. 4) then 
           optical_old = optical
           optical = ''
           istyle = 0
           leopt = lnblnk(optical)
           idir = '../njoy99l/output/groupr/'
      endif
c     redirection for NJOY-2012/NJOY-2016
      if ( icon(16) .eq. 3) then 
           optical_old = optical
           optical = ''
           istyle = 0
           leopt = lnblnk(optical_old)
           idir = optical_old(1:leopt)//'../NJOY-2012/pka/'
           if ( icon(7) .eq. 1) then 
            idir = optical_old(1:leopt)//'../NJOY-2016/pka/'
           endif
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
      temp_file1 = idir(1:iblank2)
     1   //infile(1:iblank)
      iend1 = lnblnk(temp_file1)
      write (*,9012) optical(1:leopt)//idir(1:iblank2)
     1   //infile(1:iblank), temp_file1(1:iend1),
     &   optical(1:leopt), idir(1:iblank2), 
     &   infile(1:iblank)
9012  format (1x, 'opening file: =', a,'=',/,
     &        1x, '  temp_file1: =', a,'=',/,
     &        1x, '     optical: =', a,'=',/,
     &        1x, '        idir: =', a,'=',/,
     &        1x, '      infile: =', a,'=',/)
      open(unit=7,file=temp_file1(1:iend1),
     1   status='unknown')
c
c read in the zaid of the material along with the atomic weight.
c
      read(7,9023) cdummy
 9023 format (a1)
      read(7,100) zaid,aweight
      izaid = zaid
c
c read in the temperature and the number of groups.
c
      read(7,200) temp,nenergy
      if ( nenergy .gt.1000) then 
          write (6,5612) nenergy
 5612     format (1x, 'ERROR in pka_grprin energy dimension ', i5)
          stop 'integral_pka_grprin dimension'
      endif
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
   10 read(7,600,end=20) irxntemp
c      if ( irxntemp.eq.irxn ) then
        go to 30
c      endif
      go to 10
c
c requested reaction not found on tape. 
c die with error
c
   20 write(*,*)
      write(*,*)'error : invalid reaction number :',irxn
      write(*,*)'error : program termination (grprin)'
      stop 'rsn error'
c
c continue processing the data.
c
c clear the xs array
c
   30 continue
c
c read in the pka data.
c
      flag = 0.0
      backspace(unit=7)
      read(7,650) nend
      read(7,650) nstart
      backspace(unit=7)
      do i = nstart,nend
        read (7,9111) line
 9111   format (a)
        if ( icon(9) .lt. 0) write (6,9112) line
 9112   format (1x, 'integral echo: ', a)
        read(line, 9110) ilow, inumber, index, irxn2
        ihigh = ilow + inumber - 1
c        read(7,9110) ilow, ihigh, index,irxn2
         ipka_low_energy(index) = ilow
         ipka_high_energy(index) = ihigh-1
 9110 format(11x, 11x, 11x, 8x, i3, 8x, i3, 8x, i3, 6x, i3)
c        if ( irxn2.ne.irxn ) go to 40
        read(7,2900) dummy, (pka_spectrum(index,jk), 
     &                 jk = ilow, ihigh-1)
        if ( icon(9) .lt. 0) then 
           write (6, 1971) ilow, ihigh-1, index, irxn2, 
     &        irxn, pka_spectrum(index, ilow) 
 1971      format (1x, 'integral echo-2: ', 5i5, g14.7)
        endif
        if ( ilow .gt. 1) then 
           do jk=1,ilow-1
                pka_spectrum(index,jk) = 0.0
                pka_spectrum_dede(index,jk) = 0.0
           enddo
        endif
        if ( ihigh-1 .lt. nend) then 
           do jk=ihigh,nend
                pka_spectrum(index,jk) = 0.0
                pka_spectrum_dede(index,jk) = 0.0
           enddo
        endif
c        do jk=ilow, ihigh-1
c           if (pka_spectrum(index,jk).lt.1e-35) 
c     &                pka_spectrum(index,jk) = 1.e-35
c        enddo
 2900 format(6(e11.6))
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
      write(*,700) izaid,aweight,temp,nenergy,fdate()
c
c write stdio the lower, and upper energy bounds, the midpoint of the`
c energy interval and the cross section data.
c
      write(*,800) irxn
      do i = 1,nenergy
        write(*,900) energy(i),energy(i+1),emid(i)
      end do
c
c     Renormalize PKA spectrum to unity
c          ??? should this be normalized like a number fraction or like a differential fluence
c
      do ie = 1, nenergy
         if ( icon(9) .lt. 0) then 
             write (6,924) ie, ipka_low_energy(ie), 
     &            ipka_high_energy(ie), 
     &            pka_spectrum(ie,ipka_low_energy(ie))
 924         format (1x, 'integral renorm: ', 3i5, g14.7)
         endif
         sum = 0.0
         do je = ipka_low_energy(ie), ipka_high_energy(ie)
            sum = sum + pka_spectrum(ie,je)
         enddo
         if ( sum .le. 0.0) sum = -1.0
         do je = ipka_low_energy(ie), ipka_high_energy(ie)
            pka_spectrum(ie,je) = pka_spectrum(ie,je)/sum
         enddo
      enddo
c
c    Calculate Robinson Damage Energy - normalized to one recoil particle
c
      zr = zaid
      ar = aweight
      zl = zr
      al = aweight
      displ_th = 25.
      write (6,5613) infile(1:iblank), 
     &         irxn, zr, ar, zl, al, displ_th
 5613 format (/,/, 5x, ' Robinson Effective Damage energy',
     &    ' computed using: ',/,
     &             5x, ' File                           = ', a,/,
     &             5x, ' Reaction                       = ', i3,/,
     &             5x, ' Z/A for recoil ion             = ', 2g14.7,/,
     &             5x, ' Z/A for lattice atoms          = ', 2g14.7,/,
     &             5x, ' Displacement threshold energy  = ', g14.5, 
     &             2x, ' eV',/)
      do ie = 1, nenergy
         sum = 0.0
         do je = ipka_low_energy(ie), ipka_high_energy(ie)
c            DE = 1.0
            DE = df(emid(je), zr, ar, zl, al, displ_th)
            sum = sum + pka_spectrum(ie,je)*DE
         enddo
         robinson_DE(ie) = sum
      enddo
c
c     Compute the dE/dE integral PKA spectrum
c
      do ie = 1, nenergy
         do je = ipka_low_energy(ie), ipka_high_energy(ie)
            pka_spectrum_dede(ie,je) = pka_spectrum(ie,je)*
     &          emid(je)/(energy(je+1) - energy(je))
         enddo
      enddo
c      do ie = 1, nenergy
c         tum = 0.0
c         do je = ipka_low_energy(ie), ipka_high_energy(ie)
c            tum = tum + pka_spectrum_dede(ie,je) * 
c     &            (energy(je+1) - energy(je))*1.e-6
c         enddo
c         if ( tum .le. 0.0) tum = 1.0
c         do je = ipka_low_energy(ie), ipka_high_energy(ie)
c            pka_spectrum_dede(ie,je) = pka_spectrum_dede(ie,je)/tum
c         enddo
c      enddo
c
c     Write a number fraction histogram representation of the integral PKA spectra
c
      iblank = lnblnk(infile)
      leopt = lnblnk(optical_old)
      pkadir = optical_old(1:leopt)//'output/pka/'
      kblank2 = lnblnk(pkadir)
      temp_file_out = pkadir(1:kblank2)
     1   //infile(1:iblank)
      iend1 = lnblnk(temp_file_out)
      write (reaction_label, 672) irxn
 672  format ('rxn',i3.3)
c
      do ie = 1, nenergy
         write (energy_label, 671) emid(ie)
 671     format ('ENG',1pe12.6)
         if ( icon(9) .lt. 0) then 
             write (6,814) 
     &          temp_file_out(1:iend1)//'_FF_'//reaction_label(1:6)
     &          //'_'//energy_label(1:15)
 814          format (1x, 'PKA spectrum write: ', a)
             write (6,914) ie, ipka_low_energy(ie), 
     &            ipka_high_energy(ie), 
     &            pka_spectrum(ie,ipka_low_energy(ie))
 914         format (1x, 'rehash: ', 3i5, g14.7)
         endif
         sintg = 0.0
         do je = ipka_low_energy(ie), ipka_high_energy(ie)
            sintg = sintg + pka_spectrum(ie,je)
         enddo
         if ((ipka_low_energy(ie)  .gt. 0) .and. 
     &       (ipka_high_energy(ie) .ge. ipka_low_energy(ie)) 
     &       .and. sintg .gt. 0.0 ) then 
c            open(unit=47,file=temp_file_out(1:iend1)//'_FF_'
c     &      //reaction_label(1:6)
c     &          //'_'//energy_label(1:15), status='unknown')
c            write (47, 815) energy(ipka_low_energy(ie)), 
c     &                   pka_spectrum(ie,ipka_low_energy(ie))*1.e-6
c            do je = ipka_low_energy(ie), ipka_high_energy(ie)
c               write (47, 815) energy(je),   pka_spectrum(ie,je)
c               write (47, 815) energy(je+1), pka_spectrum(ie,je)
c 815           format (1x, 2g15.7)
c            enddo
c            write (47, 815) energy(ipka_high_energy(ie)+1), 
c     &                   pka_spectrum(ie,ipka_high_energy(ie))*1.e-6
c            close (unit=47)
         endif
      enddo
c
c     Write a dE/dE histogram representation of the PKA spectra
c
c
c     Write a number fraction histogram representation of the PKA dE/dE spectra
c
      iblank = lnblnk(infile)
      leopt = lnblnk(optical_old)
      pkadir = optical_old(1:leopt)//'output/pka/'
      kblank2 = lnblnk(pkadir)
      temp_file_out = pkadir(1:kblank2)
     1   //infile(1:iblank)
      iend1 = lnblnk(temp_file_out)
      write (reaction_label,2672) irxn
2672  format ('rxn',i3.3)
c
      do ie = 1, nenergy
         write (energy_label, 2671) emid(ie)
2671     format ('ENG',1pe12.6)
         if ( icon(9) .lt. 0) write (6,2814) 
     &          temp_file_out(1:iend1)//'_dEdE_'//reaction_label(1:6)
     &          //'_'//energy_label(1:15)
2814     format (1x, 'PKA dEdE spectrum write: ', a)
         sintg = 0.0
         do je = ipka_low_energy(ie), ipka_high_energy(ie)
            sintg = sintg + pka_spectrum_dede(ie,je)
         enddo
         if ((ipka_low_energy(ie)  .gt. 0) .and. 
     &       (ipka_high_energy(ie) .ge. ipka_low_energy(ie))
     &        .and. sintg .gt.0.0 ) then 
c            open(unit=47,file=temp_file_out(1:iend1)//'_dede_'
c     &      //reaction_label(1:6)
c     &          //'_'//energy_label(1:15), status='unknown')
c            write (47, 2815) energy(ipka_low_energy(ie)), 
c     &                   pka_spectrum_dede(ie,ipka_low_energy(ie))
c     &                   *1.e-6
c            do je = ipka_low_energy(ie), ipka_high_energy(ie)
c               write (47, 2815) energy(je),   pka_spectrum_dede(ie,je)
c               write (47, 2815) energy(je+1), pka_spectrum_dede(ie,je)
c2815           format (1x, 2g15.7)
c            enddo
c            write (47, 2815) energy(ipka_high_energy(ie)+1), 
c     &                   pka_spectrum_dede(ie,ipka_high_energy(ie))
c     &                   *1.e-6
c            close (unit=47)
         endif
c
      enddo
c
c    Output Robinson Damage Energy summary
c
c      open(unit=57,file=temp_file_out(1:iend1)//'.DE',
c     &    status='unknown')
c      write (57, 451) (emid(ie), robinson_DE(ie),ie=1,nenergy)
c 451  format (1x, 2g15.7)
c      close (unit=57)
c
c     Printout of selected the recoil for selective energies
c
      if ( nenergy .eq. 770) then 
          do ie = 1, nenergy
             if ( ie .eq. 316 .or. ie .eq. 361 .or. 
     &            ie .eq. 406 .or. ie. eq. 451 .or. 
     &            ie .eq. 541 .or. ie .eq. 640 .or. 
     &            ie. eq. 641. or. ie .eq. 721) then 
                if ((ipka_low_energy(ie)  .gt. 0) .and. 
     &              (ipka_high_energy(ie) .ge. 
     &               ipka_low_energy(ie))) then 
c                    write (6, 7812) emid(ie)*1.e-6
 7812               format (1x,'Recoil spectrum for neutron energy = ',
     &              g14.7, ' MeV',/,/,6x,
     &              'Bin', 3x, '  Lower PKA   ',4x, ' Upper PKA  ',4x,
     &                   ' Number ',9x,'dE/dE',/,6x,
     &              '   ', 3x, ' Energy (MeV) ',4x, 'Energy (MeV)',4x,
     &                   'Fraction',4x,'     ')
                    icount = 0
c                    do je = ipka_low_energy(ie), ipka_high_energy(ie)
c                       icount = icount + 1
c                       write (6,7813) icount, energy(je)*1.e-6, 
c     &                    energy(je+1)*1.e-6, 
c     &                    pka_spectrum(ie,je), pka_spectrum_dede(ie,je)
c                    enddo
 7813               format (3x,i5,3x,1pg14.7,3x,1pg14.7,2x,g14.7, 
     &                  2x, g14.7, 2x)
c                    write (6,491) robinson_DE(ie)*1.e-6, 
c     &                    robinson_DE(ie)/emid(ie)*100.
 491                format (/, 4x, 'Robinson Damage Energy = ', 1pg14.7,
     &              '  MeV',/, 4x, 'Damage Energy Percent  = ', 2pg12.4,
     &               ' %',/)
                endif
             endif
          enddo
      endif
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
  800 format(7x,'upper',8x,'  lower ',8x,'midpoint',9x,'        ',/,
     1       7x,'(ev) ',8x,'  (ev)  ',8x,'  (ev)  ',9x,' (',i3,')  ')
  900 format(3(3x,1pe12.5),5x,1pe12.5)
      end
