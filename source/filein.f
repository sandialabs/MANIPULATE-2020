      subroutine filein (imode, file_name, ierr, iflag)
      include "location.cmn"
      character*22 lab
      common /misc/ lab(20), imaterial, noption, moption
      character*80 title, dummy
      character*250 file_name
      common /datain/ nenergy, energy(1001), array(1000,15),
     1                emid(1001)
      common /datahld/ nenergy_hld(15), energy_hld(1001,15),
     1         array_hld(1000,15), emid_hld(1001,15)
      dimension energy_old(1001), array_old(1000,15),
     &          diff(1001), diffeng(1002)
      character*80 comment
      common /labx/ icomment, comment(20)
      common /guide/ icon(40)
      common /io/ nt5, nt6, nfile, nplot
      character*1 char
      character*80 id
      character*145 ldir
      character*250 name
      character*250 outfile
      common /whatever/ outfile
      external grprin,grspin,spectra_out,spectra_format
c
c      ****************************************************
c
c      fill default names for disk retrieves
c
      if (icon(9) .lt. 0) then
        write (6,7823) nenergy, imode, ierr, iflag, trim(file_name)
 7823   format (1x, '*** Enter FILEIN ', 4i5,/,
     &          1x, '                 ', a)
      endif
      nfile = 34
      call getenv('opt', optical)
c      jdir = '/app/manipulate-2/response/'
      jdir = 'response/'
c      kdir = '/app/manipulate-2/spectrum/'
      kdir = 'spectrum/'
c      ldir = '/app/manipulate-2/'
      ldir = ''
c      write (*,9012) trim(optical)//trim(jdir)
c9012  format (1x, 'opening directory location: ', a41)

       do jk=1,1000
          do kl = 1,10
             array(jk,kl) = 0.0
          enddo
          emid(jk) = 0.0
          energy(jk) = 0.0
       enddo
       emid(1001) = 0.0
       energy(1001) = 0.0
       nenergy = 0
       if ( imode .eq. 0 .or. imode .eq. 12) then
c
c         njoy matxs format - sandii energy structure
c         if imode = 12, convert from differential number to number fraction
c
cje
cje override this type of input and use the grprin
cje routine to readin the groupr format.
cje
cje       open(unit=nfile, form='unformatted',
cje  1         file='rsp_neu:'//file_name//'_njoy.damage_sandii_mtx'
cje  2         , status='old', iostat=ilook, err=909)
cje       call matxsin(nfile,nt6)
          if ( iflag .gt. 1) then
               call grprin(file_name,iflag)
          elseif (abs(iflag) .eq. 1) then
               call grspin(file_name, iflag)
          endif
          if ( icon(9) .lt. 0) then
           write (nt6, 6672)
 6672      format (1x, '*** filein spectra ***')
           write (nt6, 4592) (energy(jk), jk=1,nenergy+1)
 4592      format (1x, 'energy bounds: ', 10g14.7)
          endif
c         check data validity
c
          if ( nenergy .ne. 620 .and. nenergy .ne. 640 .and.
     &         nenergy .ne. 770 .and. nenergy .ne. 89  .and.
     &         nenergy .ne. 48  .and. nenergy .ne. 725) then
              write (nt6,2034) nenergy
2034          format (1x, '*** warning *** number of energy',
     1        ' points = ', i5, /,
     2        1x, 15x, 'should be /770/725/640/8948 structure',
     &        ' for this filein option')
          endif
c
c          convert energy grid from ev to mev
c
          if ( energy(1) .gt. 10.e+6 .or.
     1         energy(nenergy+1) .gt. 10.e+6) then
            do ie = 1,nenergy+1
                energy(ie) = energy(ie)*1.e-6
            enddo
          endif
          do ie=1,nenergy
             emid(ie) = 0.5*(energy(ie) + energy(ie+1))
          enddo
c         convert from diff. number to number fraction
          if (imode .eq. 12) then
             xsum = 0.0
             do jk=1,nenergy
                diff(jk) = array(jk,1)
                zap  = array(jk,1)*(energy(jk) - energy(jk+1) )
                array(jk,1) = abs(zap)
                xsum = xsum + array(jk,1)
             enddo
             if ( xsum .le. 0.0) xsum = 1.0
             if ( abs( xsum - 1.00) .gt. 1.e-3) then
                 write (6,9253) xsum
 9253            format (1x, 'renormalized filein imode=12 spectrum',
     &           ' number fraction ',
     1           'from ', 1pg14.7, ' to 1.0 ')
             endif
             do jk=1,nenergy
                array(jk,1) = array(jk,1)/xsum
             enddo
          endif
c
c         check number fraction normalization
c
          if ( imode .eq. 0 .and. icon(6) .eq. 0) then
             xsum = 0.0
             do jk=1,nenergy
                xsum = xsum + array(jk,1)
             enddo
             if ( xsum .le. 0.0) xsum = 1.0
             if ( abs( xsum - 1.00) .gt. 1.e-3) then
                 write (6,2253) xsum
 2253            format (1x, 'renormalized filein imode=0 spectrum',
     &           ' number fraction ',
     1           'from ', 1pg14.7, ' to 1.0 ')
             endif
             do jk=1,nenergy
               array(jk,1) = array(jk,1)/xsum
             enddo
          endif
c
c         convert damage data from ev-b to mev-mb
c
cje       do imat=1,3
c             imat = 1
c             do ie = 1,nenergy
c                 array(ie,imat) = array(ie,imat)*1.e-3
c             enddo
cje       enddo
c
c      following three lines produce a formated output listing
c
c       exclude this output for known response function forms
c           mod 3/8/2008 by PJG to inhibit this for folding with
c           response function
c
         if ( icon(1) .eq. 5 .or. icon(1) .eq. -7 .or.
     1        icon(1) .eq. 7) then
         else
             if ( icon(9) .lt. 1) then
                call spectra_out
             endif
         endif
c
c      OK from here on
c
       elseif ( imode .eq. 10) then
c
c         response format - sandii energy structure
c         arranged in energy low to high energy,
c         read opposite inverted
c
          if ( icon(13) .eq. 0) then
             if ( icon(5) .eq. 0) then
                open(unit=nfile, form='formatted',
     1          file=trim(optical)//trim(jdir)
     1          //'sand641.nrg', status='old', iostat=ilook, err=909)
                ilmt = 641
             elseif ( icon(5) .eq. 1) then
                open(unit=nfile, form='formatted',
     1          file=optical//trim(jdir)//'sand771.nrg'
     2          , status='old', iostat=ilook, err=909)
                ilmt = 771
             endif
          else
             open(unit=nfile, form='formatted',
     1         file=trim(optical)//trim(jdir)
     1          //'sand621.nrg'
     2         , status='old', iostat=ilook, err=909)
             ilmt = 621
          endif
          read(nfile, 781) title
          read (nfile,*) ienergy
          if ( icon(12) .ne. 1) then
             read (nfile,*) (energy(jk), jk=ilmt,1,-1)
          else
             read (nfile,*) (energy(jk), jk=1,ilmt)
          endif
          close (unit=nfile)
          open(unit=nfile, form='formatted',
     1         file=trim(optical)//trim(jdir)
     1         //trim(file_name)
     2         , status='old', iostat=ilook, err=909)
          read (nfile, 781) title
 781      format (a80)
          read (nfile, *) renorm
          read(nfile, *) nenergy
c
c         check data validity
c
          if ( nenergy .ne. ilmt-1 ) then
              write (nt6,5034) nenergy
5034          format (1x, '*** warning *** number of energy',
     1        ' points = ', i5, /,
     2        1x, 15x, 'should be 640 structure for this option')
          endif
          read (nfile, *) (array(jk,1), jk=nenergy,1,-1)
c
c         define energy grid midpoints
c
          do ie=1,nenergy
             emid(ie) = 0.5*(energy(ie) + energy(ie+1))
          enddo
c
c         convert damage data according to normalization
c
          nmat = 1
          do imat=1,1
             do ie = 1,nenergy
                 array(ie,imat) = array(ie,imat)*renorm
             enddo
          enddo
       elseif ( imode .eq. -10) then
c
c         response format - sandii energy structure
c         arranged in energy high to low energy,
c         read opposite inverted
c
          if ( icon(13) .eq. 0) then
             if ( icon(5) .eq. 0) then
                open(unit=nfile, form='formatted',
     1          file=trim(optical)//trim(jdir)
     1          //'sand641.nrg'
     2          , status='old', iostat=ilook, err=909)
                ilmt = 641
             elseif ( icon(5) .eq. 1) then
                open(unit=nfile, form='formatted',
     1          file=trim(optical)//trim(jdir)
     1          //'sand771.nrg'
     2          , status='old', iostat=ilook, err=909)
                ilmt = 771
             endif
          else
             open(unit=nfile, form='formatted',
     1         file=trim(optical)//trim(jdir)
     1         //'sand621.nrg'
     2         , status='old', iostat=ilook, err=909)
             ilmt = 621
          endif
          read(nfile, 781) title
          read (nfile,*) ienergy
          if ( icon(12) .ne. 1) then
             read (nfile,*) (energy(jk), jk=1,ilmt)
          else
             read (nfile,*) (energy(jk), jk=ilmt,1,-1)
          endif
          close (unit=nfile)
          open(unit=nfile, form='formatted',
     1         file=trim(optical)//trim(jdir)
     1         //trim(file_name)
     2         , status='old', iostat=ilook, err=909)
          read (nfile, 781) title
          read (nfile, *) renorm
          read(nfile, *) nenergy
c
c         check data validity
c
          if ( nenergy .ne. ilmt-1 ) then
              write (nt6,5034) nenergy
          endif
          read (nfile, *) (array(jk,1), jk=1,nenergy)
c
c         define energy grid midpoints
c
          do ie=1,nenergy
             emid(ie) = 0.5*(energy(ie) + energy(ie+1))
          enddo
c
c         convert damage data according to normalization
c
          nmat = 1
          do imat=1,1
             do ie = 1,nenergy
                 array(ie,imat) = array(ie,imat)*renorm
             enddo
          enddo
       elseif ( imode .eq. 11) then
c
c         response format - arbitrary energy structure
c         arranged in energy low to high energy,
c         read opposite inverted
c
          open(unit=nfile, form='formatted',
     1         file=trim(optical)//trim(jdir)
     1         //trim(file_name)//'.nrg'
     2         , status='old', iostat=ilook, err=909)
          read(nfile, 781) title
          read (nfile,*) nenergy
          nenergy1 = nenergy + 1
          read (nfile,*) (energy(jk), jk=1,nenergy1)
          close (unit=nfile)
          open(unit=nfile, form='formatted',
     1         file=trim(optical)//trim(jdir)
     1         //trim(file_name)
     2         , status='old', iostat=ilook, err=909)
          read (nfile, 781) title
          read (nfile, *) renorm
          read(nfile, *) menergy
c
c         check data validity
c
          if ( nenergy .ne. menergy ) then
              write (nt6,6034) nenergy, menergy
6034          format (1x, '*** warning *** number of energy',
     1        ' points = ', i5, /,
     2        1x, 15x, 'should be same structure for this option')
          endif
          read (nfile, *) (array(jk,1), jk=1,nenergy)
c
c         define energy grid midpoints
c
          do ie=1,nenergy
             emid(ie) = 0.5*(energy(ie) + energy(ie+1))
          enddo
c
c         convert damage data according to normalization
c
          nmat = 1
          do imat=1,1
             do ie = 1,nenergy
                 array(ie,imat) = array(ie,imat)*renorm
             enddo
          enddo
       elseif ( imode .eq. 5) then
c
c         njoy matxs format - sandii energy structure
c
          noption = 1
          open(unit=nfile, form='unformatted',
     1         file=trim(optical)//trim(jdir)
     1         //trim(file_name)//'.pka'
     2         , status='old', iostat=ilook, err=909)
cje       call matxsin(nfile,nt6)
c
c         check data validity
c
          if ( nenergy .ne. 620 .and. nenergy .ne. 640 .and.
     &         nenergy .ne. 770  .and. nenergy .ne. 89 .and.
     &         nenergy .ne. 48 . and. nenergy .ne. 725) then
              write (nt6,2034) nenergy
          endif
c
c          convert energy grid from ev to mev
c
          do ie = 1,nenergy+1
              energy(ie) = energy(ie)*1.e-6
          enddo
          do ie=1,nenergy
             emid(ie) = 0.5*(energy(ie) + energy(ie+1))
          enddo
c
c         convert damage data from ev to kev and
c         kinematic kerma limit from ev-b to (mev-mb) and
c         divide recoil energy by total cross section
c
          do ie = 1,nenergy
              array(ie,2) = array(ie,2)*1.e-3
              array(ie,3) = array(ie,3)*1.e-3
          enddo
          do ie = 1,nenergy
              value = array(ie,1)
              if ( value .le. 0.0) then
                  write (nt6, 8723) ie, value
 8723             format (1x, '*** warning *** total cross sec',
     1            'tion zero ',
     1            5x, i5, 3x, g14.7)
                  value = 1.
              endif
              array(ie,3) = array(ie,3)/value
          enddo
       elseif ( imode .eq. 4 ) then
          open(unit=nfile, form='formatted',
     1         file=trim(file_name)
     2         , status='old', iostat= ilook, err=909)
           read(nfile,*) nenergy, (emid(jk), array(jk,1),
     1        jk=1,nenergy)
       elseif ( imode .eq. -4 ) then
          open(unit=nfile, form='formatted',
     1         file=trim(file_name)
     2         , status='old', iostat= ilook, err=909)
           read(nfile,*) nenergy, (emid(jk), array(jk,1),
     1        jk=nenergy,1,-1)
       elseif ( imode .eq. 2 .or. imode .eq. 1) then
          name = trim(file_name)
          if (name(1:6) .ne. '/esata') then
             write (nt6, 4512) trim(optical),
     1           trim(ldir), trim(file_name)
 4512            format (1x, 'filein name components (optical) = ',
     &                    a190,/,
     1                   1x, '                       (ldir)    = ',
     &                    a190,/,
     2                   1x, '                     (file_name) = ',
     &                    a190,/)
             name = trim(optical)//trim(ldir)//trim(file_name)
          endif
          lend = lnblnk(name)
          write (nt6,8726) name(1:lend)
          open(unit=nfile, form='formatted',
     1         file=name,
     2         status='old', iostat=ilook, err=909)
          read(nfile,*) nenergy
          read(nfile,1012) char
          read(nfile,1012) char
 1012     format (a1)
          read(nfile,*) (energy(jk), array(jk,1), jk=1,nenergy+1)
c
c         if energy is in ev, convert to mev
c
          emaxg = 0.0
          do jkl=1,nenergy+1
              emaxg = max(emaxg,energy(jkl))
          enddo
          if ( emaxg .gt. 10.e+6) then
               do jkl=1,nenergy+1
                       energy(jkl) = energy(jkl)*1.e-6
               enddo
          endif
c
          icomment = 1
          comment(1) = trim(file_name)
c          write (6,8419) (energy(jkl), array(jkl,1), jkl=1,nenergy+1)
c8419      format (1x, '** debug echo energy/flux   ', g14.7, 4x, g14.7)
          if ( iflag .ne. 3 .and. outfile(1:8) .ne.
     1         'spectrum' ) write (nt6, 8261) iflag
          call spectra_format
        elseif ( imode .eq. 22 ) then
 8261      format (1x, '*** probable error ***, option only',
     1     ' defined for differential number spectrum, 3 /= ', i5,/,
     2     5x, ' spectra_format should not be called ')
c
c          special mode to read from endf/b-vi cf-252 spectrum file
c          assumed format is differential number (??), unit 1/ev
c
c          interpolate on curve to sand-ii energy grid
c          input with energy in ev
c
          name = trim(optical)//trim(ldir)//trim(file_name)
          lend = lnblnk(name)
c          write (nt6,8726) name(1:lend)
 8726     format (1x, 'file = ', a80,'=eof')
          open(unit=nfile, form='formatted',
     1         file=name(1:lend),
     2         status='old', iostat=ilook, err=909)
c          write (nt6,8726) name(1:lend)
          read(nfile,*) nenergy_old
          read(nfile,1012) char
          read(nfile,1012) char
          read(nfile,*) (energy_old(jk), array_old(jk,1),
     1                  jk=1,nenergy_old+1)
          icomment = 1
          comment(1) = trim(file_name)
c
c         convert data
c
          if ( icon(5) .eq. 0) then
              nenergy = 640
              open(unit=nfile, form='formatted',
     1         file=trim(optical)//trim(jdir)//'sand641.nrg'
     2         , status='old', iostat=ilook, err=909)
          elseif ( icon(5) .eq. 1) then
              nenergy = 770
              open(unit=nfile, form='formatted',
     1         file=trim(optical)//trim(jdir)//'sand771.nrg'
     2         , status='old', iostat=ilook, err=909)
          elseif ( icon(5) .eq. 2) then
              nenergy = 89
              open(unit=nfile, form='formatted',
     1         file=trim(optical)//trim(jdir)//'esata.nrg'
     2         , status='old', iostat=ilook, err=909)
          elseif ( icon(5) .eq. 3) then
              nenergy = 48
              open(unit=nfile, form='formatted',
     1         file=trim(optical)//trim(jdir)//'nuget49.nrg'
     2         , status='old', iostat=ilook, err=909)
          elseif ( icon(5) .eq. 6) then
              nenergy = 725
              open(unit=nfile, form='formatted',
     1         file=trim(optical)//trim(jdir)//'IAEA725.nrg'
     2         , status='old', iostat=ilook, err=909)
          else
              write (6,7245) nenergy, icon(5)
 7245         format (1x, 'icon(5) error ', 2i6)
              stop 'icon(5)-3'
          endif
          read(nfile, 781) title
          read (nfile,*) ienergy
          if ( icon(5) .eq. 0) then
             ilmt = 641
          elseif (icon(5) .eq. 1) then
             ilmt = 771
          endif
          if ( icon(12) .ne. 1) then
             read (nfile,*) (energy(jk), jk=ilmt,1,-1)
          else
             read (nfile,*) (energy(jk), jk=1,ilmt)
          endif
          close (unit=nfile)
c
c         fill interpolated array
c
c         use linear interpolation
c
          interp_mode = 2
c
          do jkl=1,nenergy
            x = 0.5*(energy(jkl)+energy(jkl+1))*1.e+6
            y = fitmd(x, nenergy_old, energy_old(1),
     1                array_old(1,1),interp_mode)
            array(jkl,1) = y
c            write (nt6,8792) x, y
c8792        format (1x, 'energy/fluence pair = ', g14.7, 4x, g14.7)
          enddo
           if ( iflag .ne. 3) write (nt6, 8261) iflag
          call spectra_format
      elseif ( imode .eq. 9) then
c
c          user defined input format xxx.flux file
c
          open(unit=nfile, form='formatted',
     1      file=trim(optical)//trim(kdir)//trim(file_name)//'.flux',
     2      status='old', iostat= ilook, err=909)
           read(nfile,*) icomment, nenergy, id
           read(nfile,8901) (comment(jk),jk=1,icomment)
 8901      format (a80)
           read(nfile,*) (energy(jk), array(jk,1),
     1        jk=1,nenergy)
           energy(nenergy+1) = energy(nenergy)*1.01
           if ( iflag .ne. 3) write (nt6, 8261) iflag
           call spectra_format
       elseif ( imode .eq. 6) then
c
c         njoy matxs format - arbitrary energy structure
c
          open(unit=nfile, form='unformatted',
     1         file=trim(file_name)
     2         , status='old', iostat=ilook, err=909)
cje       call matxsin(nfile,nt6)
c
c          convert energy grid from ev to mev
c
          do ie = 1,nenergy+1
              energy(ie) = energy(ie)*1.e-6
          enddo
          do ie=1,nenergy
             emid(ie) = 0.5*(energy(ie) + energy(ie+1))
          enddo
c
c         convert damage data from ev-b to mev-mb
c
          do imat=1,3
             do ie = 1,nenergy
                 array(ie,imat) = array(ie,imat)*1.e-3
             enddo
          enddo
       elseif (imode .eq. 23) then
c
c         njoy matxs format - sandii energy structure
c         if imode = 23, convert from differential number to number fraction
c          built upon imode=12 option
c
cje
cje override this type of input and use the grprin
cje routine to readin the groupr format.
cje
cje       open(unit=nfile, form='unformatted',
cje  1         file='rsp_neu:'//file_name//'_njoy.damage_sandii_mtx'
cje  2         , status='old', iostat=ilook, err=909)
cje       call matxsin(nfile,nt6)
          if ( iflag .gt. 1) then
               call grprin(file_name,iflag)
          elseif (abs(iflag) .eq. 1) then
               call grspin(file_name, iflag)
          endif
c         check data validity
c
          if ( nenergy .ne. 620 .and. nenergy .ne. 640  .and.
     &         nenergy .ne. 770 .and. nenergy .ne. 89   .and.
     &         nenergy .ne. 48 .and. nenergy .ne. 725) then
              write (nt6,2734) nenergy
2734          format (1x, '*** warning *** number of energy',
     1        ' points = ', i5, /,
     2        1x, 15x,
     &        'should be 620/640/770/89/48/725 structure for this',
     3        ' option')
          endif
c
c          convert energy grid from ev to mev
c
          if ( energy(1) .gt. 10.e+6 .or.
     1         energy(nenergy+1) .gt. 10.e+6) then
            do ie = 1,nenergy+1
                energy(ie) = energy(ie)*1.e-6
            enddo
          endif
          do ie=1,nenergy
             emid(ie) = 0.5*(energy(ie) + energy(ie+1))
          enddo
c         convert from diff. number to number fraction and check normalization
             xsum = 0.0
             diffeng(nenergy+1) = energy(nenergy+1)
             do jk=1,nenergy
                diffeng(jk) = energy(jk)
                diff(jk) = array(jk,1)
                zap  = array(jk,1)*(energy(jk) - energy(jk+1) )*1.e+6
                array(jk,1) = abs(zap)
                xsum = xsum + array(jk,1)
             enddo
             if ( xsum .le. 0.0) xsum = 1.0
             if ( abs( xsum - 1.00) .gt. 1.e-3) then
                 write (6,9753) xsum
 9753            format (1x, 'normalized spectrum number fraction ',
     1           'is ', 1pg14.7, ' not 1.0 ')
             endif
             do jk=1,nenergy
                array(jk,1) = array(jk,1)/xsum
             enddo
c
c         refill differential number data to interface with spectra_format
c
c
          do jk=1,nenergy
               energy(jk) = diffeng(jk)
               array(jk,1) = diff(jk)
          enddo
          energy(nenergy+1) = diffeng(nenergy+1)
c
c         provide spectrum definition tabulation
c
         call spectra_format
c
       elseif ( imode .eq. 24) then
c
c         NJOY-2012 MF=5 spectrum format
c
          open(unit=nfile, form='formatted',
     1         file=trim(file_name)
     2         , status='old', iostat=ilook, err=909)
c
c          space down to energy grid
c
           read (nfile, 89) dummy
           read (nfile, 89) dummy
 89        format (a80)
c          read energy grid
           read (nfile, *) dum1, dum2, nenergy
           read (nfile, 891) dum1, dum2, (energy(jk),jk=1,4)
 891       format (6e11.6)
           read (nfile,891) (energy(jk), jk=5,nenergy+1)
c
c          find MF=5 MT=18 entry
c
           read (nfile, 89) dummy
           read (nfile, 89) dummy
           read (nfile, 89) dummy
           if ( dummy(71:75) .ne. " 5 18") then
              write (6,3651) dummy(71:75)
 3651         format (1x, 'ERROR in filein for imode=24 option: ',
     &          2x, ' label = ', a5)
           endif
c
c          read spectrum number fractions
c
           read (nfile,891) (array(jk,1),jk=1,nenergy)
c
c          convert energy grid from ev to mev
c
          do ie = 1,nenergy+1
              energy(ie) = energy(ie)*1.e-6
          enddo
          do ie=1,nenergy
             emid(ie) = 0.5*(energy(ie) + energy(ie+1))
          enddo
c
c         check number fraction normalization
c
          xsum = 0.0
          do jk=1,nenergy
             xsum = xsum + array(jk,1)
          enddo
          if ( xsum .le. 0.0) xsum = 1.0
          if ( abs( xsum - 1.00) .gt. 1.e-3) then
              write (6,2253) xsum
          endif
          do jk=1,nenergy
            array(jk,1) = array(jk,1)/xsum
          enddo
          if ( icon(1) .eq. 5 .or. icon(1) .eq. -7 .or.
     1        icon(1) .eq. 7) then
          else
             if ( icon(9) .lt. 1) then
                call spectra_out
             endif
          endif
       elseif ( imode .eq. 25) then
c
c         LSL format used for spectrum
c
          open(unit=nfile, form='formatted',
     1         file=trim(file_name)
     2         , status='old', iostat=ilook, err=909)
c
c          space down to energy grid
c
           read (nfile, 89) dummy
           read (nfile, 89) dummy
c          read energy grid
           read (nfile, *) nenergy1
           nenergy1 = abs(nenergy1)
           nenergy = nenergy1 - 1
           if ( icon(9) .lt. 3) then
              write (6, 6712) nenergy
 6712         format (1x, 'filein imode=23 nenergy = ', i5)
           endif
           read (nfile, 89) dummy
           read (nfile, *) (energy(jk), jk=1,nenergy+1)
           read (nfile, 89) dummy
c
c          read spectrum number fractions
c
           read (nfile,*) (array(jk,1),jk=1,nenergy)
c
c          convert energy grid from ev to mev
c
          do ie = 1,nenergy+1
              energy(ie) = energy(ie)*1.e-6
          enddo
          do ie=1,nenergy
             emid(ie) = 0.5*(energy(ie) + energy(ie+1))
          enddo
c
c         check number fraction normalization
c
          xsum = 0.0
          do jk=1,nenergy
             xsum = xsum + array(jk,1)
          enddo
          if ( xsum .le. 0.0) xsum = 1.0
          if ( abs( xsum - 1.00) .gt. 1.e-3) then
              write (6,2253) xsum
          endif
          do jk=1,nenergy
            array(jk,1) = array(jk,1)/xsum
          enddo
          if ( icon(1) .eq. 5 .or. icon(1) .eq. -7 .or.
     1        icon(1) .eq. 7) then
          else
             if ( icon(9) .lt. 1) then
                call spectra_out
             endif
          endif
       else
           write (nt6,2033) icon(2)
2033       format (1x, '*** input format not defined *** ',
     1     i2)
           ierr = 1
           stop 'format'
       endif
       close (unit=nfile,err=2910)
2910   continue
       if ( icon(9) .lt. 0) then
         write (6,9823) nenergy
 9823    format (1x, '*** Exit FILEIN ', i5)
       endif
       return
 909   continue
       write (nt6,910) ilook, trim(file_name)
910    format (1x, '*** error ', i5, ' opening filein ', a250)
       return
       end
