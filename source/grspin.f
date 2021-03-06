      subroutine grspin(infile,irxn)
c***********************************************************************
c* this routine reads the energy grid, and cross sections from a 
c* tabular type output tape.  
c***********************************************************************
c***********************************************************************
c234567890
      implicit real (a-h,o-z)
      implicit integer (i-n)
c      character infile*(*)
      character*(*) infile
      character*145 idir, jdir, kdir, ldir, idir_src
      character*106 optical
      character*15 ename
      character*250 new_file_name, name
      common /location / optical, idir, jdir, kdir
      common /guide/ icon(40)
      dimension elocal(1001), alocal(1001), dlocal(1001)
      character*250 outfile
      common /whatever/ outfile
c
c common blocks.
c
      common /datain/ nenergy, energy(1001), array(1000,15)
     1         , emid(1001)
c
c get the name of the input tape.
c
      iblank = lnblnk(infile)
      new_file_name = infile(1:iblank)
      if ( icon(9) .lt. 0) then 
         write (6,7823) irxn, nenergy, icon(12)
 7823    format (1x, '*** Enter GRSPIN ', 3i5)
         iend = lnblnk(new_file_name)
         write (6,17843) iend, new_file_name(1:iend)
17843    format (1x, '*** GRSPIN filename in locb', i5, 2x, a)
      endif
c  
      if ( new_file_name(1:1) .ne. '/') then 
         ename = 'opt'
         lename = lnblnk(ename)
         call getenv(ename(1:lename), optical)
c         if (optical .eq. '') then 
c             optical=''
c             write (*,8013) optical
c 8013        format (1x, 'blank optical default filled at abc', 1x, a)
c         endif
         leopt = lnblnk(optical)
c         idir = '/app/manipulate-2/response/'
c         idir_src = '/app/manipulate-2/spectrum/'
         idir = 'response/'
         idir_src = 'spectrum/'
         iblank2 = lnblnk(idir)
         if ( icon(9) < 0) then 
           write (6, 4523) optical(1:leopt), idir(1:iblank2),
     &          infile(1:iblank)
 4523      format (1x, 'name components: optical = ', a,/,
     &             1x, '                 idir    = ', a,/,
     &             1x, '                 infile  = ', a,/)
         endif
         new_file_name =  infile(1:iblank)
         if (new_file_name(3:16) .ne. "\sync_sandialabs"      ) then 
            new_file_name = optical(1:leopt) // idir(1:iblank2) 
     &         // infile(1:iblank)
         endif
         write (6,9012) infile(1:iblank)
9012     format (1x, 'opening file: ', a)
      endif
c
c open the input tape
c
      ientry = 0
      iend = lnblnk(new_file_name)
      if ( icon(9) .lt. 0) then 
         write (6,7843) iend, new_file_name(1:iend)
 7843    format (1x, '*** GRSPIN filename in loca ', i5, 2x, a)
      endif
      open(unit=7,file=new_file_name(1:iend),status='unknown')
c
c read in the data from the input tape.
c
 101  continue
      do i=1,1000
         read (7,*, end=9013) emid(i), array(i,1)
         if ( icon(9) .lt. -1) then 
            write (6,7812) i, emid(i), array(i,1)
 7812       format (1x, 'Raw GRSPIN read: ', i5, 1x, 2g14.7)
         endif
         nenergy = i
      enddo
      write (*, 9014) 
 9014 format (1x, '*** error ***, array too long in grspin ', 
     1     2x)
      stop 'grspin-long'
 9013 continue
      close (unit=7)
      if ( nenergy .eq. 0 .and. ientry .eq. 0) then 
c         error may be due to needing file in source directory
c         rather than response directory
          ientry = ientry + 1
          write (6,7845) nenergy, new_file_name
 7845     format (1x, 'GRSPIN file read error ', i5, 3x, a)
c          close (unit=7)
c          iblank2 = lnblnk(idir_src)
c          new_file_name = optical(1:leopt)//idir_src(1:iblank2)
c     1       //infile(1:iblank)
c          write (*,9012) optical(1:leopt)//idir_src(1:iblank2)
c     1       //infile(1:iblank)
c9112      format (1x, 'opening alternate file: ', a57)
c          open(unit=7,file=new_file_name,status='unknown')
c          go to 101
      endif
c
c Look at the energy order - we expect a low to high order
c     set iflip = -1 if we need to reverse the order
c
      iflip = 1
      if ( emid(1) .gt. emid(2) ) then 
           iflip = -1
      endif
c
c if a reversed energy grid due to one part of 
c the spectrum, reverse the grid now.
c
      if (irxn.eq.-1 .or. iflip .eq. -1) then 
        do jk=1,nenergy
            dlocal(jk) = emid(jk)
            alocal(jk) = array(jk,1)
        enddo
        do jk=1,nenergy
            emid(jk)   =  dlocal(nenergy+1-jk)
            array(jk,1) = alocal(nenergy+1-jk)
        enddo
      end if
c
c calculate the energy bounds from midpoints an fill that array.
c
      if ( emid(1) .lt. emid(2) ) then
         if ( emid(1) .gt. 1.e-4) then  
             energy(1) = 1.e-4
         else
             energy(1) = 1.e-10
         endif
         do i = 2,nenergy+1
           energy(i) = 2.0*emid(i-1) - energy(i-1)
         end do
      else
         if ( icon(9) .lt. 0) then
            write (6,7892) nenergy, emid(1)
 7892       format (1x, 'GRSPIN emid upper check ', i5, g14.7)
         endif
         if ( nenergy .eq. 640) then
            if ( emid(1) .gt. 15.0e+6 ) then
               energy(1) = 20.0e+6
            else
               energy(1) = 20.0
            endif
         endif
c         if ( nenergy .eq. 89) then
c            if ( emid(1) .gt. 15.0e+6 ) then
c               energy(1) = 20.0e+6
c            else
c               energy(1) = 20.0
c            endif
c         endif
         if ( nenergy .eq. 770) then
            if ( emid(1) .gt. 140.0e+6 ) then
               energy(1) = 150.0e+6
            else
               energy(1) = 150.0
            endif
         endif
         do i=2,nenergy+1
           energy(i) = 2.0*emid(i-1) - energy(i-1)
         enddo
         if ( icon(9) .lt. 0) then
            write (6,7893) nenergy, emid(1), energy(1)
 7893       format (1x, 'GRSPIN emid after check ', i5, 2g14.7)
         endif
      endif
c
c     fix exact energy bounds if 640 sand-ii energy grid
c
      if ( icon(9) .lt. 0) then 
         write (6, 2781) nenergy, energy(1), energy(nenergy+1)
 2781    format (1x, 'GRPSIN energy grid check ', i7, 2g14.7)
      endif
      if ( (nenergy .eq. 640 .and. energy(1) .eq. 20.0)  
     &      ) then 
         write (6,9361)
 9361    format (/,1x, 'sand-ii energy grid is imposed ',/) 
c
c       fill default names for disk retrieves
c
         ename = 'opt'
         lename = lnblnk(ename)
         call getenv(ename(1:lename), optical)
c         if (optical .eq. '') then 
c             optical=''
c         endif
         mblank2= lnblnk(optical)
c         jdir = '/app/manipulate-2/response/'
         jdir = 'response/'
         jblank2 = lnblnk(jdir)
c         kdir = '/app/manipulate-2/spectrum/'
         kdir = 'spectrum/'
         kblank2 = lnblnk(kdir)
c         ldir = '/app/manipulate-2/'
         ldir = '/'
         lblank2 = lnblnk(ldir)
c
         nfile = 34
         name = optical(1:mblank2)//jdir(1:jblank2)//'sand641.nrg'
         if ( icon(5) .eq. 0) then 
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 641
          elseif ( icon(5) .eq. 1) then 
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 771
          elseif ( icon(5) .eq. 2) then 
            name = optical(1:mblank2)//jdir(1:jblank2)//'nuget90.nrg'
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 90
          endif
          read(nfile, 781) title
 781      format (a80)
          read (nfile,*) isand_energy
          read (nfile,*) (energy(jk), jk=ilmt,1,-1)
          close (unit=nfile)
         do jk=1,ilmt-1
             emid(jk) = 0.5*(energy(jk) + energy(jk+1))
         enddo
      elseif ( (nenergy .eq. 770 .and. energy(1) .eq. 150.0)  
     &      ) then 
         write (6,9362)
 9362    format (/,1x, 'extended sand-ii energy grid is imposed ',/) 
c
c       fill default names for disk retrieves
c
         ename = 'opt'
         lename = lnblnk(ename)
         call getenv(ename(1:lename), optical)
         if (optical .eq. '') then 
             optical=''
         endif
         mblank2= lnblnk(optical)
c         jdir = '/app/manipulate-2/response/'
         jdir = 'response/'
         jblank2 = lnblnk(jdir)
c         kdir = '/app/manipulate-2/spectrum/'
         kdir = 'spectrum/'
         kblank2 = lnblnk(kdir)
c         ldir = '/app/manipulate-2/'
         ldir = ''
         lblank2 = lnblnk(ldir)
c
         nfile = 34
         name = optical(1:mblank2)//jdir(1:jblank2)//'sand771.nrg'
         if ( icon(5) .ne. 1 .and. nenergy .eq. 770) then 
            write (6,1823) icon(5), nenergy, emid(1)
1823        format (1x, 'WARNING: icon(5) not set for 770 SAND-II',
     &       'structure ', 2i5, g14.7)
         endif
         if ( icon(5) .eq. 0) then 
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 771
          elseif ( icon(5) .eq. 1) then 
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 771
          endif
          read(nfile, 781) title
          read (nfile,*) isand_energy
          read (nfile,*) (energy(jk), jk=ilmt,1,-1)
          close (unit=nfile)
         do jk=1,ilmt-1
             emid(jk) = 0.5*(energy(jk) + energy(jk+1))
         enddo
      elseif ( (nenergy .eq. 89 )  
     &      ) then 
         write (6,9462)
 9462    format (/,1x, 'NuGET energy grid is imposed ',/) 
c
c       fill default names for disk retrieves
c
         ename = 'opt'
         lename = lnblnk(ename)
         call getenv(ename(1:lename), optical)
         if (optical .eq. '') then 
             optical=''
         endif
         mblank2= lnblnk(optical)
c         jdir = '/app/manipulate-2/response/'
         jdir = 'response/'
         jblank2 = lnblnk(jdir)
c         kdir = '/app/manipulate-2/spectrum/'
         kdir = 'spectrum/'
         kblank2 = lnblnk(kdir)
c         ldir = '/app/manipulate-2/'
         ldir = ''
         lblank2 = lnblnk(ldir)
c
         nfile = 34
         name = optical(1:mblank2)//jdir(1:jblank2)//'nuget90.nrg'
         if ( icon(5) .ne. 2 .and. nenergy .eq. 89) then 
            write (6,823) icon(5), nenergy, emid(1)
 823        format (1x, 'WARNING: icon(5) not set for 89 NuGET',
     &       'structure ', 2i5, g14.7)
         endif
         if ( icon(5) .eq. 2) then 
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 90
          else 
            write (6, 3491) icon(5), nenergy
 3491       format (1x, 'GRSPIN open error NuGET expected: ', 2i5)
          endif
          read(nfile, 781) title
          read (nfile,*) isand_energy
          if (emid(1) .gt. 15.) then 
             read (nfile,*) (energy(jk), jk=ilmt,1,-1)
          else
             read (nfile,*) (energy(jk), jk=1,ilmt)
          endif
          close (unit=nfile)
         do jk=1,ilmt-1
             emid(jk) = 0.5*(energy(jk) + energy(jk+1))
         enddo
         if ( icon(9) .lt. 0) then 
            write (6,4512) (energy(jk), jk=1,ilmt+1)
 4512       format (1x, 'GRSPIN energy bounds: ', 5g14.7)
            write (6,4513) (emid(jk),jk=1,ilmt)
 4513       format (1x, 'GRSPIN emid           ', 5g14.7)
         endif
      elseif ( (nenergy .eq. 725 .and. energy(726) .eq. 60.0)  
     &      ) then 
         write (6,7362)
 7362    format (/,1x, 'extended IAEA energy grid is imposed ',/) 
c
c       fill default names for disk retrieves
c
         ename = 'opt'
         lename = lnblnk(ename)
         call getenv(ename(1:lename), optical)
         if (optical .eq. '') then 
             optical=''
         endif
         mblank2= lnblnk(optical)
c         jdir = '/app/manipulate-2/response/'
         jdir = 'response/'
         jblank2 = lnblnk(jdir)
c         kdir = '/app/manipulate-2/spectrum/'
         kdir = 'spectrum/'
         kblank2 = lnblnk(kdir)
c         ldir = '/app/manipulate-2/'
         ldir = ''
         lblank2 = lnblnk(ldir)
c
         nfile = 34
         name = optical(1:mblank2)//jdir(1:jblank2)//'IAEA725.nrg'
         if ( icon(5) .ne. 6 .and. nenergy .eq. 725) then 
            write (6,5823) icon(5), nenergy, emid(1)
5823        format (1x, 'WARNING: icon(5) not set for 725-grp IAEA ',
     &       'structure ', 2i5, g14.7)
         endif
         if ( icon(9) .lt. 0) then 
            write (6,2891) name, optical(1:mblanks), jdir(1:jblank2)
 2891       format (1x, 'File = ', a)
         endif
         if ( icon(5) .eq. 0) then 
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 726
          elseif ( icon(5) .eq. 6) then 
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 726
          endif
          read(nfile, 781) title
          read (nfile,*) isand_energy
          read (nfile,*) (energy(jk), jk=1, ilmt)
          close (unit=nfile)
         do jk=1,ilmt-1
             emid(jk) = 0.5*(energy(jk) + energy(jk+1))
         enddo
      endif
c
c reverse order of elements if required
c
      if (icon(12).eq.1) then 
        elocal(nenergy+1) = energy(nenergy+1)
        do jk=1,nenergy
            elocal(jk) = energy(jk)
            dlocal(jk) = emid(jk)
            alocal(jk) = array(jk,1)
        enddo
        do jk=1,nenergy
            energy(jk) =  elocal(nenergy+2-jk)
            emid(jk)   =  dlocal(nenergy+1-jk)
            array(jk,1) = alocal(nenergy+1-jk)
        enddo
        energy(nenergy+1) = elocal(1)
      end if
c
c
c
c return to the calling routine.
c
      if ( icon(9) .lt. 0) then 
         write (6,9823) irxn, nenergy
 9823    format (1x, '*** Exit GRSPIN ', 2i5)
      endif
      return
 909  continue
      write (6,9023) name
 9023 format (1x, 'error opening energy file in grspin ', /,
     1       5x, 'file = ', a70)
      stop 'grspin'
      end
