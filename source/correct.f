       subroutine correct (resp, file_name)
       common /pltlab/ lab_lead(15), lab_trail(15), lab_nchar(15),
     1 lab_len(15), icurve(15), lab_file(15)
       character*80 lab_file, mcnploc
       character*250  file_name
       character*20 lab
       character*100 file_corr, file_id
       common /sself/ file_corr
       common /misc/ lab(20), imaterial
       common /datain/ nenergy, energy(1001), array(1000,15)
     1         , emid(1001)
       common /datahld/ nenergy_hld(15), energy_hld(1001,15), 
     1          array_hld(1000,15)
     1         , emid_hld(1001,15)
       common /guide/ icon(40)
       common /io/ nt5, nt6, nfile, nplot
       common /fpnew/ alpha
       dimension  resp(1000)
       dimension corx(1000)
               iselect_material = 1
              isum = 3
              mcnploc = '/mnt/e/sync_sandialabs/SNL-LSL/self-shield'
c              mcnploc = '/esata/codes/app/lsl-2/self-shield'
              lnloc = lnblnk(mcnploc)
              iendx = lnblnk(file_name)
              call prune(file_name, ilead_1, itrail_1, 
     1           length_1, nchar_1) 
              file_id = mcnploc(1:lnloc)//'/'
     1                 //file_name(ilead_1+1:nchar_1-itrail_1)//
     1                 '.summary'
c 8923         format (1x,'apply special self-shielding/cover correction',
c     1         /, 15x, ' file = ', a100)
                if ( file_name(ilead_1+1:nchar_1-itrail_1) .eq. 
     1               'null#-void-bare') then   
                   write (nt6, 8923)
 8923              format (1x, 'Self-shielding non-existence',
     1            ' bypass used',/,/)
                   do jk=1,770
                       corx(jk) = 1.0
                   enddo
                else
                   open (unit=29, 
     1                 file=file_id, status='old', err=9099)
                   read (29, 2871) header
 2871              format (a100)
                   read (29,*) (corx(ikl), ikl=1,640)
c                   write (6,3321) (corx(ikl),ikl=1,640)
 3321              format ('ss-fct   ', 5(1pe13.6,1x),/, 
     1             (   '         ', 5(1pe13.6,1x) )    )
                   close (unit=29)
c                  do not use any attenution above 20 MeV
                   do jk=641,770
                       corx(jk) = 1.0
                   enddo
                endif
c
c             save self-shielding correction function
c
              lab_file(isum) = file_name
              lab_lead(isum) = ilead_1
              lab_trail(isum) = itrail_1
              lab_nchar(isum) = nchar_1
              lab_len(isum) = length_1
              nenergy_hld(isum) = nenergy
c
c             apply self-shielding correction factor to response
c
              do ip=1,nenergy
                  fact2 =corx(ip)
                  resp(ip) = resp(ip)*fact2
              enddo
       if ( icon(9) < 0) then 
           write (6,672) nenergy
 672       format (1x, 'CORRECT EXIT: ', i5)
       endif
       return
 9099  continue
       jstatus = 1
       write (nt6, 6723) jstatus, file_id
 6723  format (1x, 'Error in correct file open ', i5, 2x, a90)
       stop 'Open'
       end

