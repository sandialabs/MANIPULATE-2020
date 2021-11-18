      subroutine covlate(iresp_call)
c
c     Special manipulation of covariance quantities
c
c        icov(1) = 1    propogate standard deviation of an integral
c                       cross section uncertainty component
c                       (first quantity has a covariance, second is fixed)
c                  2    input covariance file, output tecplot file
c                  3    combine covariance components
c                  4    propogate standard deviation of an integral
c                       spectrum uncertainty component
c                       (second quantity has a covariance, first is fixed)
c                  5    combine weighted covariance matrices (TBD)
c
c        icov(2)         override standard response (only for 640/770/89/48 generated grp)
c                = 0        no action, use data from covariance array
c                = 1        override response with more accurate info
c            
c
c        icov(3)         override standard cov (only for generated 640/770/89/48/175 grp)
c                = 0        no action, use data from covariance array
c                = 1        override cov/cor with fully uncorrelated data
c                                in native group structure
c                = 2        override cov/cor with fully uncorrelated data
c                                in 640 group structure frame-of-ref.
c                = 3        override cov/cor with fully correlated data
c                                in 640 group structure frame-of-ref.
c            
c
c        icovfmt        format for covariance read
c                = 0       known function - no covariance
c                = 1       covfils output summary
c                = 2       SNLRML .lsl dosimetry xsec format
c                = 3       snlcov manipulate format
c                = 4       known func. plus std. dev. - no correlations
c                = 5       SNLRML .lsl spectrum format 
c                          (self-covariance only)
c                          (correlation, Std. Dev.,and number frac. data)
c                = 6       ENDF/B-VI absolute covariance format (e.g. Cf252)
c                          (abs. cov., energy grid, number frac, and 
c                           rel. std. dev.)
c                = 7       PDF (Probability Distribution Function) for spectrum
c                          uncertainty - used for mono-energetic sources
c                = 8       fcov input format (no trial, only covariance)
c
c        icode          action for covread
c                = 0       use default energy grid
c                = 1       expand covariance to 640 SAND-II energy grid
c                = 2       expand covariance to 770 extended SAND-II energy grid
c                = 3       expand in NuGET 89 group neutron structure
c                = 4       expand in NuGET 48 group photon structure
c                = 5       expand in Vitamin-J 175 group neutron structure
c
c
cje
cje addition character assignments.
cje
       integer iresp_call
       character*250 outfile, file_name
       common /whatever/ outfile
cje
       common /datain/ nenergy, energy(1001), array(1000,15)
     1         , emid(1001)
       dimension range(1000)
       dimension percent(1000)
       character*100 file_corr
       integer spc_tag_end, xsec_tag_end
       character*50 spc_tag, xsec_tag
       common /tag/ spc_tag, xsec_tag, spc_tag_end, xsec_tag_end
       common /sself/ file_corr
       common /fpnew/ alpha
       common /datahld/ nenergy_hld(15), energy_hld(1001,15), 
     1          array_hld(1000,15)
     1         , emid_hld(1001,15)
       common /guide/ icon(40)
       common /io/ nt5, nt6, nfile, nplot, npun
       common /pltlab/ lab_lead(15), lab_trail(15), lab_nchar(15),
     1 lab_len(15), icurve(15), lab_file(15)
       character*80 lab_file
       character*20 lab
      character*250 ovr
      character*145 idir, jdir, kdir, job, ename
      character*106 optical, xoptical
c      character*55 outfile
      common /location / optical, idir, jdir, kdir
      character*80 name
       common /misc/ lab(20), imaterial, noption, moption
      common /statsum/ arr_mean(1000), arr_sq(1000), 
     1  arr_std(1000),
     1  arr_var(1000)
      common /statinfo/ efactor, number_of_files, xoptical
      common /acov1/ icov(40), iself, ovr
      common /bcov1/ icoveng1, coveng1(1001), covrsp1a(1000), 
     1  covrsp1b(1000), cov1(1001,1001), cor1(1001,1001), 
     1  icoveng2, coveng2(1001), covrsp2a(1000), covrsp2b(1000), 
     1  stddev1a(1000), stddev1b(1000), stddev2a(1000), stddev2b(1000),
     1  cov2(1000,1000),
     2  cor2(1000,1000)
      common /covcmb/ icmbeng1(5), cmbeng1(5,1001), cmbrsp1a(5,1000),
     1    cmbv1(5,5,1000,1000), cmbr1(5,5,1000,1000),
     2    cmbstd1a(5,1000)
      common /covcmb1/ icmbeng, cmbeng(1001), cmbrsp(1000),
     1    cmbv(1000,1000), cmbr(1000,1000), cmbstd(1000)
      character*250 fmt1, fmt2, new_file_name, fmt3, fmt4
      common /mone/ pdf(1000)
      dimension scale(5,1000)
      dimension coveng1_mid(1000)
      dimension abs_cov(1001,1001), row_norm(1001)
      dimension abs2_cov(1001,1001), rel_corr(1001, 1001)
      character*256 epath
      common /new_io/ epath
      kdir = 'spectrum/'
      kblank2 = lnblnk(kdir)
      if (icon(9) < 0) then 
          write (6,6712) iresp_call, nenergy, kdir(1:kblank2)
 6712     format (1x, 'COVLATE entered ', 2i5, 2x, a)
      endif
      kstatus = 0
      if ( icov(1) .eq. 1) then 
c
c        propogate std. dev. of an integral
c
         read (nt5,*) icovfmt1, fmt1, icode
         read (nt5,*) icovfmt2, fmt2, 
     1      imode, iselect_material
         if ( icovfmt1 .ge. 0 .and. icovfmt1 .le. 3) then 
             if ( icode .ne. 1 .and. icode .ne. 2 .and.
     &            icode .ne. 3) then 
                icode = 1
             endif
             istatus = 1
             if ( icon(9) < 0) then 
               write (6, 7821) istatus,icode,kdir
 7821          format (1x, 'covlate covread call: ', 2i5, 2x, a)
             endif
             call covread(icovfmt1, fmt1, icoveng1, coveng1, covrsp1a, 
     1       stddev1a, covrsp1b, stddev1b, cov1, cor1, icode)
             if ( nenergy .eq. 770 .and. icoveng1 .eq. 640) then 
               write (nt6, 2891) nenergy, icoveng1
 2891          format (1x, 'covlate energy group change: 770 reduced', 
     &              ' to 640 for fold ', 2i5)
               nenergy = 640
             else
               write (nt6, 2892) nenergy, icoveng1
 2892          format (1x, 'covlate ', 
     &              ' bypass check a ', 2i5)
             endif
             if ( icon(9) .le. 0) then 
                write (6,5629) (covrsp1a(jk), jk=1,icoveng1)
 5629           format (1x, 'covread spectrum input ', 5g14.7)
                write (6,2629) (coveng1(jk), jk=1,icoveng1)
 2629           format (1x, 'covread energy input ', 5g14.7)
             endif
         else
             write (nt6, 7801) icov(1), icovfmt1
 7801        format (1x, '*** ERROR on covariance format option, ',
     1       'covfmt = ', 5i5)
             stop 'icovfmt-a'
         endif
         if ( icovfmt2 .eq. -1 ) then 
            icode = 0
            fraction = 1.0
            file_name = fmt2
            leopt = lnblnk(xoptical)
            jblank2 = lnblnk(jdir)
            kdir = 'spectrum/'
            kblank2 = lnblnk(kdir)
            ireverse= iselect_material
            if ( iselect_material.le.0 ) iselect_material = 1
            call prune (file_name, ilead, itrail, length, nchar)
            ierr = 0
            iflast = lnblnk(file_name)
            new_file_name = xoptical(1:leopt)//kdir(1:kblank2)//
c     1       "spectrum/" //
     1       file_name(1:iflast)
            if ( icon(9) < 0) then 
              write (6,4528) iresp_call, new_file_name, 
     &            file_name(iflast-3:iflast)
 4528         format (1x, 'covlate filein calla: ',i5, a,/,
     &                1x, '                last: ', '=',a,'=')
              write (6,6891) xoptical(1:leopt), kdir(1:kblank2), 
     &            file_name(1:iflast)
 6891         format (1x, 'portions: xoptical = ', a,/,
     &                1x, '              kdir = ', a,/,
     &                1x, '         file_name = ',a)
            endif
            if (file_name(iflast-3:iflast) .eq. "frac" ) then 
               call filein(imode, new_file_name, ierr, ireverse)     
            elseif ( iresp_call == 1) then 
               call filein_response(imode, new_file_name, 
     &              ierr, ireverse)
            else
               call filein(imode, new_file_name, ierr, ireverse)
            endif
c
c           save source data
c
            isum = 1
            lab_file(isum) = file_name
            lab_lead(isum) = ilead
            lab_trail(isum) = itrail
            lab_nchar(isum) = nchar
            lab_len(isum) = length
            nenergy_hld(isum) = nenergy
            do ip=1,nenergy
               emid_hld(ip,isum) = emid(ip)
               array_hld(ip,isum) = array(ip,iselect_material)
            enddo
            do jk=1,icoveng1
             coveng1_mid(jk) = (coveng1(jk) + coveng1(jk+1))/2.0
            enddo
            zap = abs(coveng1_mid(icoveng1)*1.e-6 - 
     &              emid_hld(nenergy,isum))/
     1             (coveng1_mid(icoveng1)*1.e-6)
            zap_reverse = abs(coveng1_mid(icoveng1)*1.e-6 - 
     &                emid_hld(1,isum))/
     1             (coveng1_mid(icoveng1)*1.e-6)
            if ( nenergy .ne. icoveng1 .or. 
     &         (zap .gt. 0.1 .and. zap_reverse .gt. 0.1) ) then 
               jstatus = 1
               write (nt6, 3612) jstatus, nenergy, icoveng1, 
     &          emid_hld(1,isum),
     1          emid_hld(nenergy,isum), coveng1_mid(1)*1.e-6,
     1          coveng1_mid(icoveng1)*1.e-6, zap, zap_reverse
 3612          format (1x, 'Covariance fold energy mis-match ', 
     1         3i8, 6g14.7)
               zap_part = abs(coveng1_mid(icoveng1)*1.e-6 - 
     &                emid_hld(1,isum))
               write (6,9390)  coveng1_mid(icoveng1)*1.e-6, 
     &           emid_hld(1,isum), zap_part, zap_reverse
 9390          format (1x, 'zap_reverse details ', 5g14.7)
               if (icon(9) .le. 2) then 
                  write (6,2390) iresp_call
 2390             format (1x, 'covlate iresp_call = ', i5)
                  write (6,2391) (emid_hld(jk,isum), jk=1,nenergy)
 2391             format (1x, 'emid_hld:   ', 5g14.7)
                  write (6,2392) (coveng1_mid(jk)*1.e-6, jk=1,icoveng1)
 2392             format (1x, 'coveng1_mid:', 5g14.7)
               endif
            endif
         else
             write (nt6, 7801) icov(1), icovfmt2
             stop 'icovfmt-c'
         endif
c
c        Check for self-covariance - otherwise folding doesn't mean anything
c
         if ( iself .ne. 0) then 
              write (nt6, 6804) iself
 6804         format (1x, 'ISELF = ', i4, /, 1x, 
     1         '*** WARNING *** Folding uncertainty may not have any',
     2          ' meaning',/,
     2        '                Response function 1 and associated',
     3         ' std. dev. used ')
         endif
c
c
c        perform integral - folding
c
         if ( icon(9) .le. -2) then
           write (6,6739) icoveng1, nenergy
 6739      format (1x, 'Full uncertainty fold: ', 2i4)
         endif
         sum = 0.0
         var = 0.0
         do jc1=1,icoveng1
           comp1 = array_hld(jc1,isum)*covrsp1a(jc1)
           sum = sum + comp1
           do jc2=1,nenergy
              comp2 = array_hld(jc2,isum)*covrsp1a(jc2)
              var = var + comp1*stddev1a(jc1)*comp2*
     1             stddev1a(jc2)*cor1(jc1,jc2)
              if ( icon(9) .le. 0) then
                write (6,6729) jc1,jc2, array_hld(jc1,isum), 
     &            covrsp1a(jc1), stddev1a(jc1), array_hld(jc2,isum), 
     &            covrsp1a(jc2), stddev1a(jc2)
 6729           format (1x, 'Fold partials: ', 2i4, 6g14.7)
              endif
           enddo
         enddo
         var2 = sqrt(var)
         var3 = var2/sum
c
c        find average std. dev.
c
         tum = 0.0
         avg = 0.0
         jstart = 0
         if ( icon(9) < -4 ) then 
            write (6,5198)
 5198       format (1x, 'Profile of avg std dev: ')
         endif
         do jc1=1,icoveng1
             avg = avg + covrsp1a(jc1)*stddev1a(jc1)*array_hld(jc1,isum)
             tum = tum + covrsp1a(jc1)*array_hld(jc1,isum)
             if ( jstart .eq.0 .and. avg .gt. 0) jstart = jc1
             if ( icon(9) < -4) then 
               write (6, 4579) jc1, avg, tum, covrsp1a(jc1), 
     &             stddev1a(jc1), 
     &             array_hld(jc1,isum)
 4579          format (1x,  i5, 5g14.7   )
             endif
         enddo
         avg2 = avg/tum
c
c        output uncertainty on integral
c
          if ( icon(13) .ne. 2) then 
             if ( icov(2) .ne. 1) then 
               ln1 = lnblnk(fmt1)
               ln2 = lnblnk(fmt2)
               write (nt6, 6801) fmt1(1:ln1), fmt2(1:ln2), sum, 
     &            var3, avg2
 6801          format (1x, 'Covariance Folding:',/,
     1               1x, '       Func. 1:    = ', a,/,
     2               1x, '       Func. 2:    = ', a,/,
     3               1x, '       Fold        = ', g14.7,/,
     4               1x, '       Uncertainty = ', g14.7, ' %',/,
     5               1x, '    Avg. Std. Dev. = ', g14.7, ' %') 
               write (39, 2919) spc_tag(1:spc_tag_end), 
     1              xsec_tag(1:xsec_tag_end), 
     1              sum, var3, avg2
 2919          format ('EXTRACT COV1   ', a, 1x, a, 1x, 7g14.7,/)
             else
               ln1 = lnblnk(fmt1)
               ln2 = lnblnk(fmt2)
               ln3 = lnblnk(ovr)
               write (nt6, 6871) fmt1(1:ln1), ovr(1:ln3), fmt2(1:ln2), 
     &                 sum, var3, avg2
 6871          format (1x, 'Covariance Folding:',/,
     1               1x, '       Func. 1:    = ', a,/,
     1               1x, '          override response = ', a,/,
     2               1x, '       Func. 2:    = ', a,/,
     3               1x, '       Fold        = ', g14.7,/,
     4               1x, '       Uncertainty = ', g14.7, ' %',/,
     5               1x, '    Avg. Std. Dev. = ', g14.7, ' %') 
               write (39, 2929) spc_tag(1:spc_tag_end), 
     1              xsec_tag(1:xsec_tag_end), 
     1              sum, var3, avg2
 2929          format ('EXTRACT COV2   ', a, 1x, a, 1x, 7g14.7,/)
             endif
          else
             if ( icov(2) .ne. 1) then 
               ln1 = lnblnk(fmt1)
               ln2 = lnblnk(fmt2)
               ln3 = lnblnk(file_corr)
               write (nt6, 6807) fmt1(1:ln1), fmt2(1:ln2), 
     &            file_corr(1:ln3), 
     1            sum, var3, avg2
 6807          format (1x, 'Covariance Folding:',/,
     1               1x, '       Func. 1:       = ', a100,/,
     2               1x, '       Func. 2:       = ', a100,/,
     2               1x, '       Self-shielding = ', a100,/,
     3               1x, '       Fold           = ', g14.7,/,
     4               1x, '       Uncertainty    = ', g14.7, ' %',/,
     5               1x, '    Avg. Std. Dev.    = ', g14.7, ' %') 
               write (39, 2939) spc_tag(1:spc_tag_end), 
     1              xsec_tag(1:xsec_tag_end), 
     1              sum, var3, avg2
 2939          format ('EXTRACT COV3   ', a, 1x, a, 1x, 7g14.7,/)
             else
               ln1 = lnblnk(fmt1)
               ln2 = lnblnk(fmt2)
               ln3 = lnblnk(ovr)
               ln4 = lnblnk(file_corr)
               write (nt6, 6877) fmt1(1:ln1), ovr(1:ln3), 
     &            file_corr(1:ln4), 
     1            fmt2(1:ln2), sum, var3, avg2
 6877          format (1x, 'Covariance Folding:',/,
     1               1x, '       Func. 1:       = ', a,/,
     1               1x, '          override response = ', a,/,
     2               1x, '       Self-shielding = ', a,/,
     2               1x, '       Func. 2:       = ', a,/,
     3               1x, '       Fold           = ', g14.7,/,
     4               1x, '       Uncertainty    = ', g14.7, ' %',/,
     5               1x, '    Avg. Std. Dev.    = ', g14.7, ' %') 
               write (39, 2949) spc_tag(1:spc_tag_end), 
     1              xsec_tag(1:xsec_tag_end), 
     1              sum, var3, avg2
 2949          format ('EXTRACT COVXSC ', a, 1x, a, 1x, 7g14.7,/)
             endif
          endif
c
c      do a profile of the response
c
          if ( icon(9) .lt. 1) then 
            write (nt6, 7804)
 7804       format (/,/,/, 5x, 
     1      'Profile of folded response and uncertainty',
     1      /, 2x, 'Bin  ', 5x, ' Eng. ', 11x, 'Frac.', 9x, 'Std Dev',
     1         11x, 'response', 10x, 'source', 9x, 'Sens. Func', 
     2         7x, 'Sens. Intg.'/,
     2         2x, 'Numb.',5x, '(MeV) ', 11x, 'Resp.', 9x, '       ',/) 
            tot = 0.0
            do jc1 = 1,icoveng1
                tot = tot + covrsp1a(jc1)*array_hld(jc1,isum)
            enddo       
            cum = 0.0
            jstart = jstart - 2
            do jc1 = 1,icoveng1
                cum = cum + covrsp1a(jc1)*array_hld(jc1,isum)
                if ( jc1 .gt. jstart .or. jc1 .eq. 1) then 
                  write (nt6, 6803) jc1, coveng1(jc1), cum/tum,
     1            stddev1a(jc1), covrsp1a(jc1), array_hld(jc1,isum),
     2            covrsp1a(jc1)*array_hld(jc1,isum)/tot, cum/tot
                endif
            enddo       
6803        format (3x, i5, g14.7, 2x, g14.7, 3x, g14.7, 
     &              3x, g14.7, 3x, g14.7, 3x, g14.7, 3x, g14.7)
          endif
      endif
      if ( icov(1) .eq. 4) then 
c
c        propogate std. dev. of an integral
c          covariance for spectrum
c
         read (nt5,*) icovfmt2, fmt2, 
     1      imode, iselect_material
         read (nt5,*) icovfmt1, fmt1, icode
         if ( (icovfmt1 .ge. 0 .and. icovfmt1 .le. 3) .or. 
     1         icovfmt1 .eq. 5 .or.icovfmt1 .eq. 6 .or. 
     2         icovfmt1 .eq. 7 ) then 
              if ( icode .ne. 1 .and. icode .ne. 2 .and. 
     &             icode .ne. 3 .and. icode .ne. 4 .and. 
     &             icode .ne. 5) then 
                icode = 1
              endif
c             if ( icovfmt1 .eq. 7) icode = 0
             istatus = 2
             if ( icon(9) < 0) then 
                write (6, 7821) istatus,icode, kdir
             endif
             call covread(icovfmt1, fmt1, icoveng1, coveng1, covrsp1a, 
     1       stddev1a, covrsp1b, stddev1b, cov1, cor1, icode)
             if ( icon(9) .le. 0) then 
                write (6,5639) (covrsp1a(jk), jk=1,icoveng1)
 5639           format (1x, 'covread(2) spectrum input ', 5g14.7)
             endif
             if ( nenergy .eq. 770 .and. icoveng1 .eq. 640) then 
               write (nt6, 2891) nenergy, icoveng1
               nenergy = 640
             else
               write (nt6, 2893) nenergy, icoveng1
 2893          format (1x, 'covlate ', 
     &              ' bypass check b ', 2i5)
             endif
         else
             write (nt6, 7801) icov(1), icovfmt1
             stop 'icovfmt-d'
         endif
         if ( icovfmt2 .eq. -1 ) then 
            icode = 0
            fraction = 1.0
            file_name = fmt2
            leopt = lnblnk(xoptical)
            jblank2 = lnblnk(jdir)
            kblank2 = lnblnk(kdir)
            ireverse= iselect_material
            if ( iselect_material.le.0 ) iselect_material = 1
            call prune (file_name, ilead, itrail, length, nchar)
            ierr = 0
            iflast = lnblnk(file_name)
            new_file_name = xoptical(1:leopt)//jdir(1:jblank2)//
     1       file_name(1:iflast)
c            call filein(imode, new_file_name, ierr, ireverse)
            if ( icon(9) < 0) then 
              write (6,4527) iresp_call, new_file_name
 4527         format (1x, 'covlate filein callb: ',i5, a)
            endif
            if ( iresp_call == 1) then 
               call filein_response(imode, new_file_name, 
     &              ierr, ireverse)
            else
               call filein(imode, new_file_name, ierr, ireverse)
            endif
            if ( nenergy .eq. 770 .and. icoveng1 .eq. 640) then 
               write (nt6, 3891) nenergy, icoveng1
 3891          format (1x, 'covlate filein energy group change: ',
     &              '770 reduced', 
     &              ' to 640 for fold ', 2i5)
               nenergy = 640
            else
               write (nt6, 3892) nenergy, icoveng1
 3892          format (1x, 'covlate ', 
     &              ' bypass check c ', 2i5)
            endif
c
c           save source data
c
            isum = 1
            lab_file(isum) = file_name
            lab_lead(isum) = ilead
            lab_trail(isum) = itrail
            lab_nchar(isum) = nchar
            lab_len(isum) = length
            nenergy_hld(isum) = nenergy
            do ip=1,nenergy
               emid_hld(ip,isum) = emid(ip)
               array_hld(ip,isum) = array(ip,iselect_material)
            enddo
c
c          Apply self-shielding correction if needed
c
            if ( icon(13) .eq. 2) then 
                read(nt5,*) file_name
                file_corr = file_name
               kstatus = kstatus + 1
c               write (6,8712) kstatus, new_file_name
c 8712          format (1x, 'COVLATE correct call', 4x, i5, a)
                call correct(array_hld(1,isum), file_name)
            endif
            do jk=1,icoveng1
             coveng1_mid(jk) = (coveng1(jk) + coveng1(jk+1))/2.0
            enddo
            zap = abs(coveng1_mid(icoveng1)*1.e-6 - 
     &                emid_hld(nenergy,isum))/
     1             (coveng1_mid(icoveng1)*1.e-6)
            zap_reverse = abs(coveng1_mid(icoveng1)*1.e-6 - 
     &                emid_hld(1,isum))/
     1             (coveng1_mid(icoveng1)*1.e-6)
            if ( nenergy .ne. icoveng1 .or. 
     &         (zap .gt. 0.01 .and. zap_reverse .gt. 0.1) ) then 
               jstatus = 2
               write (nt6, 3612) jstatus, nenergy, icoveng1, 
     &          emid_hld(1,isum),
     1          emid_hld(nenergy,isum), coveng1_mid(1)*1.e-6,
     1          coveng1_mid(icoveng1)*1.e-6, zap, zap_reverse
               zap_part = abs(coveng1_mid(icoveng1)*1.e-6 - 
     &                emid_hld(1,isum))
               write (6,9390)  coveng1_mid(icoveng1)*1.e-6, 
     &           emid_hld(1,isum), zap_part, zap_reverse
               if (icon(9) .le. 2) then 
                  write (6,3390) iresp_call
 3390             format (1x, 'covlate-2 iresp_call = ', i5)
                  write (6,2391) (emid_hld(jk,isum), jk=1,nenergy)
                  write (6,2392) (coveng1_mid(jk)*1.e-6, 
     &                  jk=1,icoveng1)
               endif
            endif
         else
             write (nt6, 7801) icov(1), icovfmt2
             stop 'icovfmt-e'
         endif
c
c        Check for self-covariance - otherwise folding doesn't mean anything
c
         if ( iself .ne. 0) then 
              write (nt6, 6804) iself
         endif
c
c
c        perform integral - folding
c
         if ( icovfmt1 .ne. 7) then 
            sum = 0.0
            var = 0.0
            do jc1=1,icoveng1
              comp1 = array_hld(jc1,isum)*covrsp1a(jc1)
              if ( icon(9) .le. 0) then 
                   write (nt6, 8926) jc1, coveng1(jc1), 
     1             array_hld(jc1,isum), 
     1             covrsp1a(jc1), 
     1             comp1, sum
              endif
 8926         format (1x, 'Folding partials ', i5, 9g14.7)
              sum = sum + comp1
              do jc2=1,nenergy
                 comp2 = array_hld(jc2,isum)*covrsp1a(jc2)
                 var = var + comp1*stddev1a(jc1)*comp2*
     1                stddev1a(jc2)*cor1(jc1,jc2)
              enddo
            enddo
            var2 = sqrt(var)
            var3 = var2/sum
c
c           find average std. dev.
c
            tum = 0.0
            avg = 0.0
            jstart = 0
            do jc1=1,icoveng1
                avg = avg + covrsp1a(jc1)*stddev1a(jc1)*
     1          array_hld(jc1,isum)
                tum = tum + covrsp1a(jc1)*array_hld(jc1,isum)
                if ( jstart .eq.0 .and. avg .gt. 0) jstart = jc1
            enddo
            avg2 = avg/tum
         endif
         if ( icovfmt1 .eq. 7) then 
c
c           For PDF uncertainty of monoenergetic spectrum
c               exact unc. given by pdf std. dev.
c               avg unc. = exact unc.
c
c            calculate mean and square moment
c
             tum = 0.0
             sum = 0.0
             sum2 = 0.0
             do jc = 1,icoveng1
                zam = array_hld(jc,isum)
                sum  = sum  + pdf(jc)*zam
                if ( icon(9) .le. 0) then 
                   write (6,349) jc, sum, zap, pdf(jc)
 349               format (1x, 'covlate_sum: ', i5, 2x, 3g14.7)
                endif
             enddo
             do jc = 1,icoveng1
c                zam = covrsp1a(jc)*array_hld(jc,isum)
                zam = array_hld(jc,isum)
                zap = (zam - sum)**2
                tum = zap*pdf(jc)
                sum2 = sum2 + tum
                zap = sum*sum - sum2
              if ( icon(9) .lt. 2) then 
                   write (nt6, 8926) jc, coveng1(jc), 
     1             array_hld(jc,isum), 
     1             covrsp1a(jc), pdf(jc),
     1             sum, sum2, tum, zap
              endif
             enddo
             var = sum2
             if ( var .lt. 0.0) then 
                    write (6,8268) var
 8268               format (1x, 'Warning of negative variance ', g14.7)
                    var = abs(var)
                    var2 = sqrt(var)
             elseif ( var .eq. 0.0) then 
                    var2 = 0.0
             else
                    var2 = sqrt(var)
             endif
             var3 = var2/sum*100.
             avg2 = var3
             if ( icon(9) .lt. 2) then 
                  write (nt6, 6734) sum, sum2, var, var2, 
     1            var3, avg2, tum
             endif
 6734        format (1x, 'spec partials ', 9g14.7)
         endif
c
c        output uncertainty on integral
c
          if ( icov(1) .eq. 1) then 
          if ( icon(13) .ne. 2) then 
             if ( icov(2) .ne. 1) then 
               ln1 = lnblnk(fmt1)
               ln2 = lnblnk(fmt2)
               write (nt6, 6301) fmt1(1:ln1), fmt2(1:ln2), 
     &           sum, var3, avg2
 6301          format (1x, 'Covariance  Spectrum Folding (1):',/,
     1               1x, '       Func. 1:    = ', a100,/,
     2               1x, '       Func. 2:    = ', a100,/,
     3               1x, '       Fold        = ', g14.7,/,
     4               1x, '       Uncertainty = ', g14.7, ' %',/,
     5               1x, '    Avg. Std. Dev. = ', g14.7, ' %') 
               write (39, 2959) spc_tag(1:spc_tag_end), 
     1              xsec_tag(1:xsec_tag_end), 
     1              sum, var3, avg2
 2959          format ('EXTRACT COV5   ', a, 1x, a, 1x, 7g14.7,/)
             else
               ln1 = lnblnk(fmt1)
               ln2 = lnblnk(fmt2)
               ln3 = lnblnk(ovr)
               write (nt6, 6371) fmt1(1:ln1), ovr(1:ln3), fmt2(1:ln2), 
     &           sum, var3, avg2
 6371          format (1x, 'Covariance Spectrum Folding (2):',/,
     1               1x, '       Func. 1:    = ', a,/,
     1               1x, '          override response = ', a,/,
     2               1x, '       Func. 2:    = ', a,/,
     3               1x, '       Fold        = ', g14.7,/,
     4               1x, '       Uncertainty = ', g14.7, ' %',/,
     5               1x, '    Avg. Std. Dev. = ', g14.7, ' %') 
               write (39, 2969) spc_tag(1:spc_tag_end), 
     1              xsec_tag(1:xsec_tag_end), 
     1              sum, var3, avg2
 2969          format ('EXTRACT COV6   ', a, 1x, a, 1x, 7g14.7,/)
             endif
          else
             if ( icov(2) .ne. 1) then
               ln1 = lnblnk(fmt1)
               ln2 = lnblnk(fmt2)
               ln3 = lnblnk(file_corr) 
               write (nt6, 6308) fmt1(1:ln1), fmt2(1:ln2), 
     &           file_corr(1:ln3), 
     1           sum, var3, avg2
 6308          format (1x, 'Covariance  Spectrum Folding (3):',/,
     1               1x, '       Func. 1:       = ', a,/,
     2               1x, '       Func. 2:       = ', a,/,
     2               1x, '       Self-shielding = ', a,/,
     3               1x, '       Fold           = ', g14.7,/,
     4               1x, '       Uncertainty    = ', g14.7, ' %',/,
     5               1x, '    Avg. Std. Dev.    = ', g14.7, ' %') 
               write (39, 2979) spc_tag(1:spc_tag_end), 
     1              xsec_tag(1:xsec_tag_end), 
     1              sum, var3, avg2
 2979          format ('EXTRACT COV7   ', a, 1x, a, 1x, 7g14.7,/)
             else
               ln1 = lnblnk(fmt1)
               ln2 = lnblnk(fmt2)
               ln3 = lnblnk(ovr)
               ln4 = lnblnk(file_corr)
               write (nt6, 6378) fmt1(1:ln1), ovr(1:ln3), 
     &           file_corr(1:ln4), 
     1           fmt2(1:ln2), sum, var3, avg2
 6378          format (1x, 'Covariance Spectrum Folding (4):',/,
     1               1x, '       Func. 1:       = ', a,/,
     1               1x, '          override response = ', a,/,
     2               1x, '       Self-shielding = ', a,/,
     2               1x, '       Func. 2:       = ', a,/,
     3               1x, '       Fold           = ', g14.7,/,
     4               1x, '       Uncertainty    = ', g14.7, ' %',/,
     5               1x, '    Avg. Std. Dev.    = ', g14.7, ' %') 
               write (39, 2989) spc_tag(1:spc_tag_end), 
     1              xsec_tag(1:xsec_tag_end), 
     1              sum, var3, avg2
 2989          format ('EXTRACT COV8   ', a, 1x, a, 1x, 7g14.7,/)
             endif
          endif
          elseif ( icov(1) .eq. 4) then 
          if ( icon(13) .ne. 2) then 
             if ( icov(2) .ne. 1) then 
               ln1 = lnblnk(fmt1)
               ln2 = lnblnk(fmt2)
               write (nt6, 5301) fmt1(1:ln1), fmt2(1:ln2), 
     &               sum, var3, avg2
 5301          format (1x, 'Covariance  Spectrum Folding (5):',/,
     1               1x, '       Spc. Func. 1:    = ', a,/,
     2               1x, '       Spc. Func. 2:    = ', a,/,
     3               1x, '       Spc. Fold        = ', g14.7,/,
     4               1x, '       Spc. Uncertainty = ', g14.7, ' %',/,
     5               1x, '    Spc. Avg. Std. Dev. = ', g14.7, ' %') 
               write (39, 2999) spc_tag(1:spc_tag_end), 
     1              xsec_tag(1:xsec_tag_end), 
     1              sum, var3, avg2
 2999          format ('EXTRACT COV9   ', a, 1x, a, 1x, 7g14.7,/)
             else
               ln1 = lnblnk(fmt1)
               ln2 = lnblnk(fmt2)
               ln3 = lnblnk(ovr)
               write (nt6, 5371) fmt1(1:ln1), ovr(1:ln3), fmt2(1:ln2),
     &              sum, var3, avg2
 5371          format (1x, 'Covariance Spectrum Folding (6):',/,
     1               1x, '       Spc. Func. 1:    = ', a,/,
     1               1x, '          Spc. override response = ', a,/,
     2               1x, '       Spc. Func. 2:    = ', a,/,
     3               1x, '       Spc. Fold        = ', g14.7,/,
     4               1x, '       Spc. Uncertainty = ', g14.7, ' %',/,
     5               1x, '    Spc. Avg. Std. Dev. = ', g14.7, ' %') 
               write (39, 3919) spc_tag(1:spc_tag_end), 
     1              xsec_tag(1:xsec_tag_end), 
     1              sum, var3, avg2
 3919          format ('EXTRACT COV11  ', a, 1x, a, 1x, 7g14.7,/)
             endif
          else
             if ( icov(2) .ne. 1) then 
               ln1 = lnblnk(fmt1)
               ln2 = lnblnk(fmt2)
               ln3 = lnblnk(file_corr)
               write (nt6, 5308) fmt1(1:ln1), fmt2(1:ln2), 
     &           file_corr(1:ln3), 
     1           sum, var3, avg2
 5308          format (1x, 'Covariance  Spectrum Folding (7):',/,
     1               1x, '       Spc. Func. 1:       = ', a,/,
     2               1x, '       Spc. Func. 2:       = ', a,/,
     2               1x, '       Self-shielding = ', a,/,
     3               1x, '       Spc. Fold           = ', g14.7,/,
     4               1x, '       Spc. Uncertainty    = ', g14.7, ' %',/,
     5               1x, '    Spc. Avg. Std. Dev.    = ', g14.7, ' %') 
               write (39, 3929) spc_tag(1:spc_tag_end), 
     1              xsec_tag(1:xsec_tag_end), 
     1              sum, var3, avg2
 3929          format ('EXTRACT COVSPC ', a, 1x, a, 1x, 7g14.7,/)
             else
               ln1 = lnblnk(fmt1)
               ln2 = lnblnk(fmt2)
               ln3 = lnblnk(file_corr)
               write (nt6, 5378) fmt1(1:ln1), fmt2(1:ln2), 
     &           file_corr(1:ln3), 
     1           fmt2, sum, var3, avg2
 5378          format (1x, 'Covariance Spectrum Folding (8):',/,
     1               1x, '       Spc. Func. 1:       = ', a,/,
     1               1x, '          Spc. override response = ', a,/,
     2               1x, '       Spc. Self-shielding = ', a,/,
     2               1x, '       Spc. Func. 2:       = ', a,/,
     3               1x, '       Spc. Fold           = ', g14.7,/,
     4               1x, '       Spc. Uncertainty    = ', g14.7, ' %',/,
     5               1x, '    Spc. Avg. Std. Dev.    = ', g14.7, ' %') 
               write (39, 3939) spc_tag(1:spc_tag_end), 
     1              xsec_tag(1:xsec_tag_end), 
     1              sum, var3, avg2
 3939          format ('EXTRACT COV13  ', a, 1x, a, 1x, 7g14.7,/)
             endif
          endif
          else
             write (nt6, 2671) icov(1)
2671         format (1x, 'No output format for this option ', i5)
          endif
c
c      do a profile of the response
c
          if ( icon(9) .lt. 1) then 
            write (nt6, 7804)
            tot = 0.0
            do jc1 = 1,icoveng1
                tot = tot + covrsp1a(jc1)*array_hld(jc1,isum)
            enddo     
            cum = 0.0
            jstart = jstart - 2
            do jc1 = 1,icoveng1
                cum = cum + covrsp1a(jc1)*array_hld(jc1,isum)
                if ( jc1 .gt. jstart .or. jc1 .eq. 1) then 
                  write (nt6, 6803) jc1, coveng1(jc1), cum/tum,
     1            stddev1a(jc1), covrsp1a(jc1), array_hld(jc1,isum),
     2            covrsp1a(jc1)*array_hld(jc1,isum)/tot, cum/tot
                endif
            enddo       
          endif
      endif
      if ( icov(1) .eq. 3) then 
c
c        combine covariance components
c
         read (nt5, *) ndiag, noff
         if ( ndiag .gt. 5 .or. ndiag .lt. 1) then 
           write (nt6, 9056) ndiag
 9056      format (1x, '*** Too many/few components to combine: ',
     1     ' ndiag = ', i5 )
           stop 'ndiag'
         endif
         do jreact=1,ndiag
            do ireact=1,ndiag
               do jk1=1,1000
                 do jk2=1,1000
                   cmbv1(jreact,ireact,jk1,jk2) = 0.0
                   cmbr1(jreact,ireact,jk1,jk2) = 0.0
                 enddo
               enddo
            enddo
         enddo
         do jk=1,ndiag
           read (nt5,*) icovfmt1, fmt1, icode
           if ( icovfmt1 .ge. 0 .and. icovfmt1 .le. 3) then 
             icode = 1
             istatus = 3
             if ( icon(9) < 0) then 
               write (6, 7821) istatus, icode, kdir
             endif
               call covread(icovfmt1, fmt1, icoveng1, coveng1, covrsp1a, 
     1         stddev1a, covrsp1b, stddev1b, cov1, cor1, icode)
               if ( icon(9) .le. 0) then 
                  write (6,5649) (covrsp1a(jr), jr=1,icoveng1)
 5649             format (1x, 'covread (3) spectrum input ', 5g14.7)
               endif
           else
               write (nt6, 5801) icov(1), icovfmt1
 5801          format (1x, '*** ERROR on covariance format option, ',
     1         'covfmt = ', 2i5)
               stop 'icovfmt-b'
           endif
           if ( iself .ne. 0) then 
               write (nt6, 3529) jk, iself
 3529          format (1x, 'iself invalid for jk= ', i5, 'iself=', i5)
               stop 'iself-2'
           endif
           do ip=1,1000
               scale(jk,ip) = 0.0
           enddo
c          read scaling factor
           read (nt5,*) fmt3, 
     1        imode, iselect_material
           fraction = 1.0
           file_name = fmt3
           leopt = lnblnk(xoptical)
           jblank2 = lnblnk(jdir)
           ireverse= iselect_material
           if ( iselect_material.le.0 ) iselect_material = 1
           call prune (file_name, ilead, itrail, length, nchar)
           ierr = 0
           iflast = lnblnk(file_name)
           new_file_name = xoptical(1:leopt)//jdir(1:jblank2)//
     1      file_name(1:iflast)
c           call filein(imode, new_file_name, ierr, ireverse)
            if ( icon(9) < 0) then 
              write (6,4523) iresp_call, new_file_name
 4523         format (1x, 'covlate filein callc: ',i5, a)
            endif
            if ( iresp_call == 1) then 
               call filein_response(imode, new_file_name, 
     &              ierr, ireverse) 
            else
               call filein(imode, new_file_name, ierr, ireverse)
            endif
c
c          save response data
c
           isum = 3
           lab_file(isum) = file_name
           lab_lead(isum) = ilead
           lab_trail(isum) = itrail
           lab_nchar(isum) = nchar
           lab_len(isum) = length
           nenergy_hld(isum) = nenergy
           do ip=1,nenergy
              emid_hld(ip,isum) = emid(ip)*1.e6
              array_hld(ip,isum) = array(ip,iselect_material)
              scale(jk,ip) = array_hld(ip,isum)
           enddo
c
c          Not implemented for self-shielding correction
c
           if ( icon(13) .eq. 2) then 
              write (nt6, 3761)
 3761         format (1x, 'Self-shielding not permitted in this case !')
              stop 'Self-Shielding'
           endif
c          place results in holding array
           icmbeng1(jk) = icoveng1
           do jk1=1,641
             cmbeng1(jk,jk1) = coveng1(jk1)
           enddo
c          check for disagreement between covariance 
c               energy grid and response/scale energy grid
           if ( nenergy .ne. icoveng1) then 
              write (nt6, 3521) nenergy, icoveng1
 3521         format (1x, 'covlate energy grid error ', 2i5)
              stop 'covlate e grid 1'
           endif
           do jk1=1,nenergy
              engx = 0.5*(coveng1(jk1)+coveng1(jk1+1))
              zap = abs(engx-emid_hld(jk1,isum))
              zap2 = zap/engx
              if ( zap2 .gt. 0.0001) then 
                 write (nt6, 3523) jk, jk1, engx, 
     1              emid_hld(jk1,isum), zap2
 3523            format (1x, 'covlate energy table error ', 2i5, 3g14.7)
                 stop 'covlate e table error 1'
              endif
           enddo
           do jk1=1,640
             cmbstd1a(jk,jk1) = stddev1a(jk1)
             cmbrsp1a(jk,jk1) = covrsp1a(jk1)*scale(jk,jk1)
           enddo
           do jk1=1,640
             do jk2=1,640
               cmbv1(jk,jk,jk1,jk2) = cov1(jk1,jk2)*scale(jk,jk1)
     1         *scale(jk,jk2)
             enddo
           enddo
         enddo
c        read off-diagonal couplings,where they exist
         if ( noff .gt. 0) then 
         do jk=1,noff
           read (nt5, *) idif1,idif2
           if ( idif1 .le. 0 .or. idif1 .gt. ndiag) then 
              write (nt6, 7603) idif1, ndiag
 7603         format (1x, 'idif = ', i5, 'max = ', i4)
              stop 'idif1'
           endif
           if ( idif2 .le. 0 .or. idif2 .gt. ndiag) then 
              write (nt6, 7603) idif2, ndiag
              stop 'idif2'
           endif
           read (nt5,*) icovfmt1, fmt1, icode
           if ( icovfmt1 .ge. 0 .and. icovfmt1 .le. 3) then 
             icode = 1
             istatus = 4
             if(  icon(9) < 0) then 
                write (6, 7821) istatus
             endif
               call covread(icovfmt1, fmt1, icoveng1, coveng1, covrsp1a, 
     1         stddev1a, covrsp1b, stddev1b, cov1, cor1, icode)
               if ( icon(9) .le. 0) then 
                  write (6,5659) (covrsp1a(js), js=1,icoveng1)
 5659             format (1x, 'covread (5) spectrum input ', 5g14.7)
               endif
           else
               write (nt6, 7801) icov(1), icovfmt1
               stop 'icovfmt-f'
           endif
           if ( iself .eq. 0) then 
               write (nt6, 4529) jk, iself
 4529          format (1x, 'warning: iself one for off-diagonal',
     1          ' jk= ', i5, '  iself=', i5)
c               stop 'iself-3'
           endif
           do jk1=1,640
             do jk2=1,640
               cmbv1(idif1,idif2,jk1,jk2) = cov1(jk1,jk2)*
     1              scale(idif1,jk1)*scale(idif2,jk2)
               cmbv1(idif2,idif1,jk1,jk2) = cov1(jk1,jk2)*
     1              scale(idif2,jk1)*scale(idif1,jk2)
             enddo
           enddo
         enddo
         endif
c        combine covariance and stddev terms
c            - stddev assumes uncorrolated input parts
         do jk1=1,640
           zap = 0.0
           zap2 = 0.0
           do ireact=1,ndiag
               xap = cmbstd1a(ireact,jk1)*cmbrsp1a(ireact,jk1)
               zap = zap + xap**2
               zap2 = zap2 + cmbrsp1a(ireact,jk1)
           enddo
           if ( zap .gt. 0.0) then 
               zapx = sqrt(zap)
           elseif (zap .eq. 0.0) then 
               zapx = 0.0
           else
               write (nt6, 7462) jk1, zap
 7462          format (1x, 'negative square of CMB at ', i5, g14.7)
               stop 'CMB-SQR'
           endif
           cmbrsp(jk1) = zap2
           if ( zap2 .ne. 0.0) then
              cmbstd(jk1) = zapx/zap2
           elseif ( zap .eq. 0.0) then 
              cmbstd(jk1) = 0.0
           else 
               write (nt6, 3823) jk1, zap, zap2
 3823          format (1x, 'covlate combo divide by zero ', 
     1         i5, 2g14.7)
               stop 'covlate div zap2'
           endif
         enddo
c        generate covariance 
         do jk1=1,640
            do jk2=1,640
               zap = 0.0
               do ireact = 1,ndiag
                  do jreact = 1,ndiag
                     zap = zap + cmbv1(ireact,jreact,jk1,jk2)
                  enddo
               enddo
               cmbv(jk1,jk2) = zap
            enddo
         enddo
c        fix std. dev. to reflect covariance
         do jk1=1,640
            if ( cmbrsp(jk1) .ne. 0.0) then 
               zay = sqrt(cmbv(jk1,jk1))/cmbrsp(jk1)
            else
               zay = 0.0
            endif
            cmbstd(jk1) = zay*100.0
         enddo
c        regeneate correlation matrix
         do jk1=1,640
             dif1 = cmbrsp(jk1)*cmbstd(jk1)*0.01
            do jk2=1,640
               dif2 = cmbrsp(jk2)*cmbstd(jk2)*0.01
               dif = dif1*dif2
               if ( dif .ne. 0.0) then 
                 cmbr(jk1,jk2) = cmbv(jk1,jk2)/dif
               else
                 cmbr(jk1,jk2) = 0.0
               endif
c              fix diagonal matrix elements
               if ( jk1 .eq. jk2 .and. cmbr(jk1,jk2) .ne. 0.0) then 
                   zap = abs(cmbr(jk1,jk2) - 1.0)
                   if ( zap .gt. 1.e-2) then 
                      write (nt6,7931) jk1,jk2, cmbr(jk1,jk2), 
     1                cmbv(jk1,jk2), dif1, dif2, zap
 7931                 format (1x, 'CMBR .NE. 1.0 check value ',
     1                2i5, 5g14.7)
c                      stop 'CMBR ERROR 1'
                   endif 
                   cmbr(jk1,jk2) = 1.0
               endif
               if ( cmbr(jk1,jk2) .gt. 1.0 ) then 
                  zap = abs(cmbr(jk1,jk2)-1.0)
                  if ( zap .gt. 0.001) then 
                      write (nt6, 4825) jk1,jk2, cmbr(jk1,jk2)
 4825                 format (1x, 'Override correlation gt 1.0 ', 
     1                 2i5, g14.7) 
                  endif
                  cmbr(jk1,jk2) = 1.0
               endif
               if ( cmbr(jk1,jk2) .lt. -1.0) then 
                  zap = abs(cmbr(jk1,jk2)+1.0)
                  if ( zap .gt. 0.001) then 
                      write (nt6, 4225) jk1,jk2, cmbr(jk1,jk2)
 4225                 format (1x, 'Override correlation lt -1.0 ', 
     1                 2i5, g14.7) 
                  endif
                  cmbr(jk1,jk2) = -1.0
               endif
             enddo
         enddo
c        store new covariance data
         read (nt5,*) fmt4
         ifmt4 = lnblnk(fmt4)
         nfile = 38
         name = xoptical(1:leopt)//jdir(1:jblank2)//
     1          'covar/'//fmt4(1:ifmt4)//'.snlcov'
         lend = lnblnk(name)
         file_name = name(1:lend)
         if ( icon(9) < 0) then 
           write (6,5624)
 5624      format (1x, 'covlate open: 5624 ')
         endif
          open(unit=nfile, form='formatted',
     1         file=file_name,
     2         status='unknown', iostat=ilook, err=909)
c        write snlcov type file
         icoveng1 = icmbeng1(1)
c        check for all the same energy bins
         do ireact=1,ndiag
            if ( icmbeng1(1) .ne. icmbeng1(ireact)) then 
                  write (nt6, 6681) ireact
 6681             format (1x, 'Covlate combo energy error ', 
     1            i5)
            endif
         enddo
         do jk1=1,icoveng1+1
            do ireact=1,ndiag
               if ( cmbeng1(1,jk1) .ne. cmbeng1(ireact,jk1) ) then
                  write (nt6, 5681) ireact, jk1, cmbeng1(1,jk1), 
     1            cmbeng1(ireact,jk1)
 5681             format (1x, 'Covlate combo energy error ', 
     1            2i5, 2g14.7)
               endif
            enddo
            cmbeng(jk1) = cmbeng1(1,jk1)
         enddo
         write (nfile,8011) icoveng1
 8011    format (i10)
         write (nfile,8012) (cmbeng(jk),jk=1,icoveng1+1)
 8012    format (1x, 8(g14.7,1x))
         write (nfile,8012) (cmbrsp(jk),jk=1,icoveng1)
         write (nfile,8012) (cmbstd(jk),jk=1,icoveng1)
         write (nfile,8012) (cmbrsp(jk),jk=1,icoveng1)
         write (nfile,8012) (cmbstd(jk),jk=1,icoveng1)
         write (nfile,8012) ((cmbv(jk1,jk2), jk1=1,icoveng1), 
     1      jk2=1,icoveng1)
         write (nfile,8012) ((cmbr(jk1,jk2), jk1=1,icoveng1), 
     1      jk2=1,icoveng1)
         close (unit=nfile)
c         write (nt6, 8801) icov(1)
c 8801    format (1x, 'ICOV option not yet implemented ', i5)
c         stop 'cov combo tbd'
      endif
      if ( icov(1) .eq. 2) then
c
c         input covariance, produce tecplot output
c 
         read (nt5,*) icovfmt1, fmt1, icode
         if ( icovfmt1 .eq. 3 .or. icovfmt1 .eq. 2 
     1   .or. icovfmt1 .eq. 1 .or. icovfmt1 .eq. 6 
     1   .or. icovfmt1 .eq. 8 .or. icovfmt1 .eq. 5 ) then 
c            icovfmt1 = 7 (PDF) is not a valid option here
             istatus = 5
c             write (6, 7821) istatus
             call covread(icovfmt1, fmt1, icoveng1, coveng1, covrsp1a, 
     1       stddev1a, covrsp1b, stddev1b, cov1, cor1, icode)
c             do jl1 = 1, icoveng1
c             do jl2 = 1, icoveng1
c                write (6,17829) jl1, jl2, cor1(jl1,jl2)
c17829           format (1x, 'COVLATE receive ', 2i5, g14.7)
c             enddo
c             enddo
         else
             write (nt6, 7801) icov(1), icovfmt1, icode
             stop 'icovfmt-g'
         endif
c
c        plot tecplot output of correlation matrix
c          
         ename = 'job'
         lename = lnblnk(ename)
         call getenv(ename(1:lename), job)
         ijob = lnblnk(job)
         nfile = 34
c         file_name = '/app/manipulate-2/tecplot/data/'
c     1     //job(1:ijob)//'-cor.data'
         ipath = lnblnk(epath)
         indx = lnblnk(fmt1)
c         file_name = epath(1:ipath) // '/tecplot/data/'
c     1     // fmt1(1:indx)//'-tecplot.data'
         file_name ='/sync_Linux/Manipulate-2020/tecplot/data/'
     1     //fmt1(1:indx)//'-tecplot.data'
         if ( icon(9) < 0) then 
           write (6,5623)
 5623      format (1x, 'covlate open: 5623 ')
         endif
         open(unit=nfile, form='formatted',
     1     file=file_name,
     2     status='unknown', iostat=ilook, err=909)
c        locate reaction threshold
         ithresh = 1
         inonzero = 0
         if ( icon(9) < 0) then 
           write (6,5625) icoveng1
 5625      format (1x, 'covlate: 5625  ', i6)
         endif
         do jk1=1,icoveng1+1
             if ( cor1(jk1,jk1) .eq. 0.0 .and .inonzero .eq. 0) then
                 ithresh = jk1 + 1
             else
                 inonzero = inonzero + 1
             endif
         enddo
          write (nfile,5632) icoveng1+2-ithresh, icoveng1+2-ithresh
 5632    format (1x, 'TITLE = "LSL-format Correlation" ',/,
     1           1x, 'VARIABLES = "Eng-1", "Eng-2", "Rel. Cov." '
     2           ,/,
     2           1x, 'ZONE T="TEST-1", I=',i3, ', J=', i3, 
     3           ', F=POINT ')
c         write (6,8921) ithresh, icoveng1
 8921    format (1x, 'DEBUG covlate dimension ', 2i5)
         do jk1=ithresh,icoveng1+1
             cor1(jk1,icoveng1+1) = 0.0
             cor1(icoveng1+1,jk1) = 0.0
             do jk2=ithresh,icoveng1+1
                write(nfile,5671) alog10(coveng1(jk1)),
     1           alog10(coveng1(jk2)), 
     1            cor1(jk1,jk2)
             enddo
         enddo
         zap = 0.0
c         do jk1=1,icoveng1
c            write (nfile, 5671) alog10(coveng1(jk1)), 
c     1       alog10(coveng1(icoveng1+1)), 
c     1        zap
c            write (nfile, 5671) alog10(coveng1(icoveng1+1)), 
c     1      alog10(coveng1(jk1)), 
c     1        zap
c         enddo
c         write (nfile, 5671) alog10(coveng1(icoveng1+1)),
c     1     alog10(coveng1(icoveng1+1)), 
c     1     zap
 5671    format (2x, g14.7, 3x, g14.7, 3x, g14.7)
         close (unit=nfile)
c
c       form absolute covariance
c
         inc = 0
         do jk1=1,icoveng1
            do jk2=jk1, icoveng1
               inc = inc + 1
c
c              Convert relative covariance to absolute correlation
c
               dummx = covrsp1a(jk1)*covrsp1a(jk2)*stddev1a(jk1)*
     1                 stddev1a(jk2)*0.01*0.01
               abs_cov(jk1,jk2) = cor1(jk1,jk2)*dummx
c               write (6,345) jk1, jk2, abs_cov(jk1,jk2), dummx, 
c     &                 covrsp1a(jk1), covrsp1a(jk2), stddev1a(jk1), 
c     &                 stddev1a(jk2)
c 345           format (1x, 'abs_cov debug: ', 2i4, 6g14.7)
               abs_cov(jk2,jk1) = abs_cov(jk1,jk2)
            enddo
         enddo
c     
c       check normalization condition for covariance values
c        e.g. sum of elements inany row/column is zero 
c        ENDF requires this to a precision of 1.E-5
c
c        only check for a spectrum - not a response function
c
         ispc_cov = 0
c
         if ( ispc_cov == 1) then 
            flunorm = 0.0
            do jp = 1, icoveng1
              flunorm = flunorm + covrsp1a(jp)
            enddo
            write (6,1592) icoveng1, flunorm
 1592       format (1x, 'Overall fluence normalization = ', i5, g14.7)
c
            do jk1 = 1,icoveng1
               row_norm(jk1) = 0.0
               do jk2 = 1,icoveng1
                 row_norm(jk1) = row_norm(jk1) + abs_cov(jk1,jk2)
               enddo
               if ( abs(row_norm(jk1)) .gt. 1.e-5) then 
c                   inhibit write since renorm is done
c                   write (6,7827) jk1, row_norm(jk1)
c 7827              format (1x, 'initial row norm error ', i5, 2x, g14.7)
               endif
            enddo  
c
c           In case of a normalization error on covariance constraint,
c               e.g. sum over any row or column is not zero,
c               treat these as the covariance for an unnormalized
c               spectrum and apply the normalization condition
c               using derived formula from Smith, pg. 140
c
           do jk1 = 1, icoveng1
              do jk2 = 1, icoveng1
                 sum = 0.0
                 do il1 = 1, icoveng1
                    do il2 = 1, icoveng1
                        offset1 = 0.0
                        offset2 = 0.0
                        if ( jk1 .eq. il1) offset1 = flunorm
                        if ( jk2 .eq. il2) offset2 = flunorm
                        sum = sum + (offset1-covrsp1a(jk1))
     &                             *abs_cov(il1,il2)
     &                             *(offset2-covrsp1a(jk2))
                    enddo
                 enddo
                 abs2_cov(jk1, jk2) = sum/flunorm**4
                 if ( icon(9) < 0) then 
                   write (6,645) jk1, jk2, abs2_cov(jk1,jk2)
 645               format (1x, 'abs2_cov debug: ', 2i4, 6g14.7)
                 endif
             enddo
           enddo 
c
c         also fix the fluence normalization - just in case
c
           do jp = 1, icoveng1
             covrsp1a(jp) = covrsp1a(jp)/flunorm
           enddo
c
c         When the covariance matrix was modified, the standard deviations
c         may have been slightly modified.  Thus, update this data to 
c         be consistent.  Note this is the relative std. dev.
c
          do jp = 1,icoveng1
              stddev1a(jp) = sqrt(abs2_cov(jp,jp))*100./covrsp1a(jp)
               if ( icon(9) < 0) then
                 write (6,2982) jp, stddev1a(jp), abs_cov(jp,jp), 
     &                covrsp1a(jp)
 2982            format (1x, 'renormed stddev: ', i4, 6g14.7)
               endif
          enddo
c
c     
c          check normalization condition for modified covariance values
c           e.g. sum of elements in any row/column is zero 
c           ENDF requires this to a precision of 1.E-5
c
            iflag = 0
            error_max = 0.0
            do jk1 = 1,icoveng1
               row_norm(jk1) = 0.0
               do jk2 = 1,icoveng1
                 row_norm(jk1) = row_norm(jk1) + abs2_cov(jk1,jk2)
               enddo
               error_max = max(error_max, row_norm(jk1))
               if ( abs(row_norm(jk1)) .gt. 1.e-5) then 
                   iflag = iflag + 1
                   write (6,6827) jk1, row_norm(jk1)
 6827             format (1x, 'renormalized row norm error ', 
     1            i5, 2x, g14.7)
               endif
            enddo   
            write (6,7290) iflag, error_max
 7290       format (1x, 'Absolute covariance renormalization had ',
     &            i4, ' normalization problems ',/, 1x,
     &            'Largest row normalization = ', g14.7) 
c
c           end of normalization logic segment
c
         endif  
c
c        write ENDF format interface of renormalized 
c              absolute covariance values
c
c         TBD
c
c        derive new relative covariance matrix 
c        (rel_corr) and output renormalized covariance 
c        files in LSL format
c
         inc = 0
         do jk1=1,icoveng1
            do jk2=jk1, icoveng1
               inc = inc + 1
c
c              Convert absolute covariance to relative correlation
c
               dummx = covrsp1a(jk1)*covrsp1a(jk2)*stddev1a(jk1)*
     1                 stddev1a(jk2)*0.01*0.01
               rel_corr(jk1,jk2) = abs2_cov(jk1,jk2)/dummx
              if ( icon(9) < 0) then 
               write (6,445) jk1, jk2, rel_corr(jk1,jk2), dummx,
     &                 abs_cov(jk1,jk2),  
     &                 covrsp1a(jk1), covrsp1a(jk2), stddev1a(jk1), 
     &                 stddev1a(jk2)
 445           format (1x, 'rel_cov covlate debug: ', 2i4, 7g14.7)
               endif
               rel_corr(jk2,jk1) = rel_corr(jk1,jk2)
            enddo
         enddo
c
c         TBD 
c
c        print LSL format covariance in original energy grid
c
         write (nt6, 7832) 
 7832    format (/,/,1x, 
     1    'Output renormalized covariance data in LSL - ',
     1     'format ',/,/)
c                                                                    
c Title cards                                                       
c                                                                  
         open (unit=78, file='renorm_covariance.lsl',
     &         status='unknown')
         write(78,8562)                                                
 8562    format(    '*COR    (LIBRARY)    (MAT.#)    (TEMP)K')         
c                                                                    
c Energy grid                                                         
c                                                                   
         write(78,8996)                                                 
 8996    format('*Number of Energies plus 1')                      
         write(78,8995) icoveng1+1                                          
         write(78,8993)                                                 
 8993    format('*Energy Grid ( eV )')          
 8995    format(i5)                             
         write(78,8990) (coveng1(jk),jk=1,icoveng1+1)  
c                                                  
c Cross section                                                       
c                                                                    
         write(78,8998)                                                 
 8998    format('*Number Fractions ')                           
         write(78,7454) (covrsp1a(jk), jk=1,icoveng1)                                                     
c                                                                   
c Standard deviation                                               
c                                                                   
         write(78,8991)                                            
 8991    format('*% Standard Deviation')                            
         write(78,8990) (stddev1a(i1),i1=1,icoveng1)                           
c                                                                  
c Correlation coefficients                                        
c                                                                  
         write(78,8992)                                              
 8992    format('*Correlation Coefficient -- Upper Triangular')       
         do i1 = 1,icoveng1       
            write(78,8990) ((100.0*rel_corr(i1,i2)),i2=i1,icoveng1)              
         end do                                                       
 8990    format((1x,8(1pe10.3,1x)))                                    
 7454    format((1x,1p8e10.3))                                         
 7990    format((1x,15(i4,1x)))                                   
         close (unit=78)
c
c       end of LSL re-write logic                                
c
      endif
c
      if ( icov(1) .eq. 5) then 
c
c        Combine covariance matrices
c
         if ( icon(9) <= 0) then 
            write (6,672)
 672        format (1x, 'covlate covariance combination option',
     &        ' is in development ')
         endif
      endif
c
      if ( icon(9) < 0) then 
         write (6,3649) 
 3649    format (1x, 'EXIT COVLATE: return ')
c         stop 'covlate 5649'
      endif
      return
 909  continue
      length = len(file_name)
      write (nt6,910) ilook, file_name(1:length)
910   format (1x, '*** covlate error ', i5, ' opening file ', a70)
      stop 'covlate read error'
      end
