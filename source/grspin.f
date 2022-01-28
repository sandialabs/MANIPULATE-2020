c***********************************************************************
c* this routine reads the energy grid, and cross sections from a
c* tabular type output tape.
c***********************************************************************
      subroutine grspin(infile,irxn)
      implicit real (a-h,o-z)
      implicit integer (i-n)
      include "location.cmn"
      character*(*) infile
      character*145 idir_src
      character*250 new_file_name, name
      common /guide/ icon(40)
      dimension elocal(1001), alocal(1001), dlocal(1001)
      character*250 outfile
      common /whatever/ outfile
      common /datain/ nenergy, energy(1001), array(1000,15),
     1                emid(1001)
c
c get the name of the input tape.
c
      new_file_name = trim(infile)
      if (icon(9) < 0) then
        write (6,7823) irxn, nenergy, icon(12)
 7823   format (1x, '*** Enter GRSPIN ', 3i5)
        write (6,17843) trim(new_file_name)
17843   format (1x, '*** GRSPIN filename in locb', a)
      end if

      if (new_file_name(1:1) /= trim(slash) .and.
     &    new_file_name(2:2) /= ":") then
        call getenv('opt', optical)
        idir = 'response'//trim(slash)
        idir_src = 'spectrum'//trim(slash)
        if (icon(9) < 0) then
          write (6, 4523) trim(optical), trim(idir), trim(infile)
 4523     format (1x, 'name components: optical = ', a,/,
     &            1x, '                 idir    = ', a,/,
     &            1x, '                 infile  = ', a,/)
      end if
      if (index(new_file_name,"sync_sandialabs") == 0) then
        new_file_name = trim(optical) // trim(idir) // trim(infile)
      end if
      write (6,9012) trim(new_file_name)
9012    format (1x, 'opening file: ', a)
      end if
c
c open the input tape
c
      ientry = 0
      if (icon(9) < 0) then
        write (6,7843) trim(new_file_name)
 7843   format (1x, '*** GRSPIN filename in loca ', a)
      end if
      open(unit=7,file=trim(new_file_name),status='unknown')
c
c read in the data from the input tape.
c
 101  continue
      do i=1,1000
        read (7,*, end=9013) emid(i), array(i,1)
        if (icon(9) < -1) then
          write (6,7812) i, emid(i), array(i,1)
 7812     format (1x, 'Raw GRSPIN read: ', i5, 1x, 2g14.7)
        end if
        nenergy = i
      end do
      write (*, 9014)
 9014 format (1x, '*** error ***, array too long in grspin ',2x)
      stop 'grspin-long'
 9013 continue
      close (unit=7)
      if (nenergy == 0 .and. ientry == 0) then
c         error may be due to needing file in source directory
c         rather than response directory
          ientry = ientry + 1
          write (6,7845) nenergy, new_file_name
 7845     format (1x, 'GRSPIN file read error ', i5, 3x, a)
c          close (unit=7)
c          iblank2 = len_trim(idir_src)
c          new_file_name = optical(1:leopt)//idir_src(1:iblank2)
c     1       //infile(1:iblank)
c          write (*,9012) optical(1:leopt)//idir_src(1:iblank2)
c     1       //infile(1:iblank)
c9112      format (1x, 'opening alternate file: ', a57)
c          open(unit=7,file=new_file_name,status='unknown')
c          go to 101
      end if
c
c Look at the energy order - we expect a low to high order
c     set iflip = -1 if we need to reverse the order
c
      iflip = 1
      if (emid(1) > emid(2)) then
           iflip = -1
      end if
c
c if a reversed energy grid due to one part of
c the spectrum, reverse the grid now.
c
      if (irxn==-1 .or. iflip == -1) then
        do jk=1,nenergy
            dlocal(jk) = emid(jk)
            alocal(jk) = array(jk,1)
        end do
        do jk=1,nenergy
            emid(jk)   =  dlocal(nenergy+1-jk)
            array(jk,1) = alocal(nenergy+1-jk)
        end do
      end if
c
c calculate the energy bounds from midpoints an fill that array.
c
      if (emid(1) < emid(2)) then
         if (emid(1) > 1.e-4) then
             energy(1) = 1.e-4
         else
             energy(1) = 1.e-10
         end if
         do i = 2,nenergy+1
           energy(i) = 2.0*emid(i-1) - energy(i-1)
         end do
      else
         if (icon(9) < 0) then
            write (6,7892) nenergy, emid(1)
 7892       format (1x, 'GRSPIN emid upper check ', i5, g14.7)
         end if
         if (nenergy == 640) then
            if (emid(1) > 15.0e+6) then
               energy(1) = 20.0e+6
            else
               energy(1) = 20.0
            end if
         end if
c         if (nenergy == 89) then
c            if (emid(1) > 15.0e+6) then
c               energy(1) = 20.0e+6
c            else
c               energy(1) = 20.0
c            end if
c         end if
         if (nenergy == 770) then
            if (emid(1) > 140.0e+6) then
               energy(1) = 150.0e+6
            else
               energy(1) = 150.0
            end if
         end if
         do i=2,nenergy+1
           energy(i) = 2.0*emid(i-1) - energy(i-1)
         end do
         if (icon(9) < 0) then
            write (6,7893) nenergy, emid(1), energy(1)
 7893       format (1x, 'GRSPIN emid after check ', i5, 2g14.7)
         end if
      end if
c
c     fix exact energy bounds if 640 sand-ii energy grid
c
      if (icon(9) < 0) then
         write (6, 2781) nenergy, energy(1), energy(nenergy+1)
 2781    format (1x, 'GRPSIN energy grid check ', i7, 2g14.7)
      end if

      !Setup file location/naming variables
      call getenv('opt', optical)
      jdir = 'response'//slash
      kdir = 'spectrum'//slash
      nfile = 34
      if (nenergy == 640 .and. energy(1) == 20.0) then
         write (6,9361)
 9361    format (/,1x, 'sand-ii energy grid is imposed ',/)
c
c       fill default names for disk retrieves
c
         name = trim(optical)//trim(jdir)//'sand641.nrg'
         if (icon(5) == 0) then
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 641
          else if (icon(5) == 1) then
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 771
          else if (icon(5) == 2) then
            name = trim(optical)//trim(jdir)//'nuget90.nrg'
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 90
          end if
          read(nfile, 781) title
 781      format (a80)
          read (nfile,*) isand_energy
          read (nfile,*) (energy(jk), jk=ilmt,1,-1)
          close (unit=nfile)
         do jk=1,ilmt-1
             emid(jk) = 0.5*(energy(jk) + energy(jk+1))
         end do
      else if (nenergy == 770 .and. energy(1) == 150.0) then
         write (6,9362)
 9362    format (/,1x, 'extended sand-ii energy grid is imposed ',/)
c
c       fill default names for disk retrieves
c
         name = trim(optical)//trim(jdir)//'sand771.nrg'
         if (icon(5) /= 1 .and. nenergy == 770) then
            write (6,1823) icon(5), nenergy, emid(1)
1823        format (1x, 'WARNING: icon(5) not set for 770 SAND-II',
     &       'structure ', 2i5, g14.7)
         end if
         if (icon(5) == 0) then
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 771
          else if (icon(5) == 1) then
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 771
          end if
          read(nfile, 781) title
          read (nfile,*) isand_energy
          read (nfile,*) (energy(jk), jk=ilmt,1,-1)
          close (unit=nfile)
         do jk=1,ilmt-1
             emid(jk) = 0.5*(energy(jk) + energy(jk+1))
         end do
      else if (nenergy == 89) then
         write (6,9462)
 9462    format (/,1x, 'NuGET energy grid is imposed ',/)
c
c       fill default names for disk retrieves
c
         name = trim(optical)//trim(jdir)//'nuget90.nrg'
         if (icon(5) /= 2 .and. nenergy == 89) then
            write (6,823) icon(5), nenergy, emid(1)
 823        format (1x, 'WARNING: icon(5) not set for 89 NuGET',
     &       'structure ', 2i5, g14.7)
         end if
         if (icon(5) == 2) then
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 90
          else
            write (6, 3491) icon(5), nenergy
 3491       format (1x, 'GRSPIN open error NuGET expected: ', 2i5)
          end if
          read(nfile, 781) title
          read (nfile,*) isand_energy
          if (emid(1) > 15.) then
             read (nfile,*) (energy(jk), jk=ilmt,1,-1)
          else
             read (nfile,*) (energy(jk), jk=1,ilmt)
          end if
          close (unit=nfile)
         do jk=1,ilmt-1
             emid(jk) = 0.5*(energy(jk) + energy(jk+1))
         end do
         if (icon(9) < 0) then
            write (6,4512) (energy(jk), jk=1,ilmt+1)
 4512       format (1x, 'GRSPIN energy bounds: ', 5g14.7)
            write (6,4513) (emid(jk),jk=1,ilmt)
 4513       format (1x, 'GRSPIN emid           ', 5g14.7)
         end if
      else if (nenergy == 725 .and. energy(726) == 60.0) then
         write (6,7362)
 7362    format (/,1x, 'extended IAEA energy grid is imposed ',/)
c
c       fill default names for disk retrieves
c
         name = trim(optical)//trim(jdir)//'IAEA725.nrg'
         if (icon(5) /= 6 .and. nenergy == 725) then
            write (6,5823) icon(5), nenergy, emid(1)
5823        format (1x, 'WARNING: icon(5) not set for 725-grp IAEA ',
     &       'structure ', 2i5, g14.7)
         end if
         if (icon(9) < 0) then
            write (6,2891) name, optical(1:mblanks), trim(jdir)
 2891       format (1x, 'File = ', a)
         end if
         if (icon(5) == 0) then
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 726
          else if (icon(5) == 6) then
            open(unit=nfile, form='formatted',
     1        file=name,
     2        status='old', iostat=ilook, err=909)
             ilmt = 726
          end if
          read(nfile, 781) title
          read (nfile,*) isand_energy
          read (nfile,*) (energy(jk), jk=1, ilmt)
          close (unit=nfile)
         do jk=1,ilmt-1
             emid(jk) = 0.5*(energy(jk) + energy(jk+1))
         end do
      end if
c
c reverse order of elements if required
c
      if (icon(12)==1) then
        elocal(nenergy+1) = energy(nenergy+1)
        do jk=1,nenergy
            elocal(jk) = energy(jk)
            dlocal(jk) = emid(jk)
            alocal(jk) = array(jk,1)
        end do
        do jk=1,nenergy
            energy(jk) =  elocal(nenergy+2-jk)
            emid(jk)   =  dlocal(nenergy+1-jk)
            array(jk,1) = alocal(nenergy+1-jk)
        end do
        energy(nenergy+1) = elocal(1)
      end if
c
c return to the calling routine.
c
      if (icon(9) < 0) then
         write (6,9823) irxn, nenergy
 9823    format (1x, '*** Exit GRSPIN ', 2i5)
      end if
      return
 909  continue
      write (6,9023) name
 9023 format (1x, 'error opening energy file in grspin ', /,
     1        5x, 'file = ', a70)
      stop 'grspin'
      end subroutine grspin