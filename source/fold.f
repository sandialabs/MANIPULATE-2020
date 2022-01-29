       subroutine fold
       common /pltlab/ lab_lead(15), lab_trail(15), lab_nchar(15),
     1 lab_len(15), icurve(15), lab_file(15)
       character*250 lab_file
       character*100 file_corr
       integer spc_tag_end, xsec_tag_end
       character*50 spc_tag, xsec_tag
       common /tag/ spc_tag, xsec_tag, spc_tag_end, xsec_tag_end
       common /sself/ file_corr
       common /datain/ nenergy, energy(1001), array(1000,15)
     1         , emid(1001)
       common /datahld/ nenergy_hld(15), energy_hld(1001,15),
     1          array_hld(1000,15)
     1         , emid_hld(1001,15)
       common /guide/ icon(40)
       common /io/ nt5, nt6, nfile, nplot
       common /fpnew/ alpha
       dimension resp_frac(641), resp_eng(641)
       character*2256 FMT9011!, FMT8011
       sum = 0.0
       sum2 = 0.0
       sum3 = 0.0
       tot1 = 0.0
       tot2 = 0.0
       flu1 = 0.0
       flu2 = 0.0
       ibegsrc = lab_lead(1) + 1
       iendsrc = lab_nchar(1) - lab_trail(1)
       itrailsrc = lab_trail(1)
       lensrc = lab_len(1)
       ncharsrc = lab_nchar(1)
       ibegrsp = lab_lead(2) + 1
       iendrsp = lab_nchar(2) - lab_trail(2)
       itrailrsp = lab_trail(2)
       lenrsp = lab_len(2)
       ncharrsp = lab_nchar(2)
       lencor = lnblnk(file_corr)
          if ( icon(13) .eq. 1) then
           write (nt6,7823)
 7823      format (15x, 'attenuation factor applied')
          endif
       write (FMT9011, '("/,/,/, 2x, 24Hspectrum folding option ,/,
     1  15x, 16Hsource file:    , a",I3.3,",/,15x, 16Hresponse file:  ,
     2  a",I3.3,",/,/,10x, 5H mid , 10x, 6Hsource, 9x, 8Hresponse, 9x,
     3  7Hfolding,7x, 9Hincrement, /,9x, 6Henergy, 10x, 6Hnumber, 9x,
     4  8Hfunction, 2x,/,9x, 5H(mev), 10x, 8Hfraction, /" )') lensrc,
     5  lenrsp
c       write (nt6, 3409 ) FMT9011
c3409  format (1x, 'FMT9011: ', a)
c      write (FMT8011, '("/,/,/, 2x, 24Hspectrum folding option ,/,
c    1  15x, 16Hsource file:    , a",I3.3,",/,15x, 16Hresponse file:  ,
c    2  a",I3.3,",/, 15x, 18Hcorrection file:   , a",I3.3,",/,/,10x,
c    3  5H mid , 10x, 6Hsource, 9x, 8Hresponse, 9x, 7Hfolding,7x,
c    4  9Hincrement, /,9x, 6Henergy, 10x, 6Hnumber, 9x, 8Hfunction,2x,
c    5  /,9x, 5H(mev), 10x, 8Hfraction, /" )') lensrc, lenrsp, lencor
c       write (nt6, 3408 ) FMT8011
c3408   format (1x, 'FMT8011: ', a)
       if ( icon(9) .lt. 2 .and. icon(13) .ne. 2 ) then
c
c      ??? some format issue here with the Linux version - need to uncomment and fix ???
c
c          write (nt6,FMT9011) lab_file(1)(ibegsrc:iendsrc),
c     1                     lab_file(2)(ibegrsp:iendrsp)
c          write (30,FMT9011) lab_file(1)(ibegsrc:iendsrc),
c     1                     lab_file(2)(ibegrsp:iendrsp)
c          write (31,FMT9011) lab_file(1)(ibegsrc:iendsrc),
c     1                     lab_file(2)(ibegrsp:iendrsp)
       elseif ( icon(9) .lt. 2 .and. icon(13) .eq. 2 ) then
c
c      ??? some format issue here with the Linux version - need to uncomment and fix ???
c
c          write (nt6,FMT8011) lab_file(1)(ibegsrc:iendsrc),
c     1                     lab_file(2)(ibegrsp:iendrsp),
c     2                     file_corr(1:lencor)
c          write (30,FMT8011) lab_file(1)(ibegsrc:iendsrc),
c     1                     lab_file(2)(ibegrsp:iendrsp),
c     2                     file_corr(1:lencor)
c          write (31,FMT8011) lab_file(1)(ibegsrc:iendsrc),
c     1                     lab_file(2)(ibegrsp:iendrsp),
c     2                     file_corr(1:lencor)
c
       elseif (icon(13) .eq. 2) then
c
c      ??? some format issue here with the Linux version - need to uncomment and fix ???
c
c          write (nt6,FMT5311) lab_file(1)(ibegsrc:iendsrc),
c     1                     lab_file(2)(ibegrsp:iendrsp),
c     2                     file_corr(1:lencor)
       else
c
c      ??? some format issue here with the Linux version - need to uncomment and fix ???
c
c          write (nt6,FMT9311) lab_file(1)(ibegsrc:iendsrc),
c     1                     lab_file(2)(ibegrsp:iendrsp)
       endif
c8311   format (/,/,/, 2x, 'spectrum folding option ',/,
c     1        15x, 'source file:      ', a<lensrc>,/,
c     2        15x, 'response file:    ', a<lenrsp>,/,
c     3        15x, 'correction file:  ', a<lencor>,/)
c8011   format (/,/,/, 2x, 'spectrum folding option ',/,
c     1        15x, 'source file:      ', a<lensrc>,/,
c     2        15x, 'response file:    ', a<lenrsp>,/,
c     3        15x, 'correction file:  ', a<lencor>,/,/,
c     3        10x, ' mid ', 10x, 'source', 9x,
c     4        'response', 9x, 'folding',7x, 'increment', /,
c     5        9x, 'energy', 10x, 'number', 9x,
c     6        'function', 2x,/,
c     7        9x, '(mev)', 10x, 'fraction', /)
c5311   format (/,/,/, 2x, 'spectrum folding option ',/,
c     1        15x, 'source file:    ', a<lensrc>,/,
c     2        15x, 'response file:  ', a<lenrsp>,/,
c     3        15x, 'correction file:  ', a<lencor>,/)
c9311   format (/,/,/, 2x, 'spectrum folding option ',/,
c     1        15x, 'source file:    ', a<lensrc>,/,
c     2        15x, 'response file:  ', a<lenrsp>,/)
c9011   format (/,/,/, 2x, 'spectrum folding option ',/,
c     1        15x, 'source file:    ', a<lensrc>,/,
c     2        15x, 'response file:  ', a<lenrsp>,/,/,
c     3        10x, ' mid ', 10x, 'source', 9x,
c     4        'response', 9x, 'folding',7x, 'increment', /,
c     5        9x, 'energy', 10x, 'number', 9x,
c     6        'function', 2x,/,
c     7        9x, '(mev)', 10x, 'fraction', /)
       do ie = 1,nenergy_hld(1)
          xdeflt = emid_hld(ie,1)
          if ( xdeflt .le. 0.0) xdeflt = 1.0
          dif = (emid_hld(ie,1) - emid_hld(ie,2))/ xdeflt
          dif = abs(dif)
          if ( dif .gt. 1.e-2) then
              write (nt6, 9012) ie, emid_hld(ie,1), emid_hld(ie,2)
9012          format (1x, '*** error in fold energy grid *** ',/,
     1        15x, i5, 2g14.7)
              stop 'e-grid'
          endif
          sum = sum + array_hld(ie,1)*array_hld(ie,2)
          xincrement = array_hld(ie,1)*array_hld(ie,2)
          resp_frac(ie) = sum
          tot1 = tot1 + array_hld(ie,1)
          tot2 = tot2 + array_hld(ie,2)
          if ( emid_hld(ie,1) .gt. 10.e-3) then
               sum2 = sum2 + array_hld(ie,1)*array_hld(ie,2)
               flu1 = flu1 + array_hld(ie,1)
               flu2 = flu2 + array_hld(ie,2)
          endif
          if ( emid_hld(ie,1) .gt. 3.e+0) then
               sum3 = sum3 + array_hld(ie,1)*array_hld(ie,2)
               flu3_1 = flu3_1 + array_hld(ie,1)
               flu3_2 = flu3_2 + array_hld(ie,2)
          endif
          if ( icon(9) .lt. 2) then
             write (nt6,9014) ie, emid_hld(ie,1),
     1            array_hld(ie,1),
     1            array_hld(ie,2), sum, xincrement
9014          format (1x, i4, 5(g14.7,2x))
             write (30,9714) emid_hld(ie,1),
     1            xincrement
             write (31,9714) emid_hld(ie,1),
     1            sum
9714          format (1x, 5(g14.7,2x))
          endif
       enddo
       sum2_norm = sum2/flu1
c
c      calculate e05, e95,and e50 (mode) of response
c
c      re-normalize incremental response
       do ie = 1,nenergy_hld(1)
          resp_frac(ie) = resp_frac(ie)/sum
          resp_eng(ie) = emid_hld(ie,1)
       enddo
c      interpolate
c
c         use linear interpolation
c
       interp_mode = 2
       x = 0.05
       e05 = fitmd(x, nenergy_hld(1), resp_frac(1),
     1      resp_eng(1), interp_mode)
       x = 0.10
       e10 = fitmd(x, nenergy_hld(1), resp_frac(1),
     1      resp_eng(1), interp_mode)
       x = 0.25
       e25 = fitmd(x, nenergy_hld(1), resp_frac(1),
     1      resp_eng(1), interp_mode)
       x = 0.75
       e75 = fitmd(x, nenergy_hld(1), resp_frac(1),
     1      resp_eng(1), interp_mode)
       x = 0.90
       e90 = fitmd(x, nenergy_hld(1), resp_frac(1),
     1      resp_eng(1), interp_mode)
       x = 0.95
       e95 = fitmd(x, nenergy_hld(1), resp_frac(1),
     1      resp_eng(1), interp_mode)
       x = 0.50
       e50 = fitmd(x, nenergy_hld(1), resp_frac(1),
     1      resp_eng(1), interp_mode)
       write (nt6, 9015) sum
 9015  format (/,/,/, 5x, 'final folded summation: ',
     1         g14.7 )
          if ( icon(13) .eq. 1) then
           write (nt6,7863) alpha
 7863      format (5x, 'attenuation factor applied', 5x, g14.7)
          endif
       write (nt6, 9335) sum3, sum2, sum2_norm, flu1, tot1, flu2,
     1      tot2, e05, e10, e25, e75, e90, e95, e50
 9335  format ( 5x, 'final folded summation above 3 mev: ',
     1         g14.7 ,/,
     1          5x, 'final folded summation above 10 kev: ',
     1         g14.7 ,/,
     1          5x, ' ---- if normed to fluence > 10 kev: ',
     2         g14.7,/,
     1          5x, 'first function = ', g14.7,
     2              ' of a total response = ', g14.7,/,
     1          5x, 'second function = ', g14.7,
     2              ' of a total response = ', g14.7,/,
     1          5x, 'E05 energy = ', g14.7,/,
     1          5x, 'E10 energy = ', g14.7,/,
     1          5x, 'E25 energy = ', g14.7,/,
     1          5x, 'E75 energy = ', g14.7,/,
     1          5x, 'E90 energy = ', g14.7,/,
     1          5x, 'E95 energy = ', g14.7,/,
     1          5x, 'E50 energy = ', g14.7,/,
     1           1x)
       write (39, 2919) spc_tag(1:spc_tag_end),
     1      xsec_tag(1:xsec_tag_end), sum,
     1      e05, e10, e25, e50, e75, e90, e95
 2919  format ('EXTRACT ENG    ', a, 1x, a, 1x, 8g14.7)
       return
       end
