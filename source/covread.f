      subroutine covread(icovfmt, fmt, icoveng, coveng, 
     1           covrsp1, stddev1, covrsp2, stddev2, 
     1           cov, cor, icode)
       common /io/ nt5, nt6, nfile, nplot, npun
       common /datain/ nenergy, energy(1001), array(1000,15)
     1         , emid(1001)
       common /pltlab/ lab_lead(15), lab_trail(15), lab_nchar(15),
     1 lab_len(15), icurve(15), lab_file(15)
       character*80 lab_file
      dimension alf(80),ilf(80),
     1 ytem(41),item(41)
      character*250 name, ovr, over_ride 
       common /guide/ icon(40)
      common /acov1/ icov(40), iself, ovr
      character*250 fmt, fmt3
      character*80 dummy
       character*100 file_corr
       common /sself/ file_corr
      common /mone/ pdf(1000)
      dimension engpdf(1001), outpdf(1001)
      dimension cov(1001,1001), coveng(1001), covrsp1(1000), 
     1   covrsp2(1000), cor(1001,1001)
      real cov_cov
      dimension covr(1001,1001), covreng(1001), corr(1001,1001),
     1    covrrsp1(1000), covrrsp2(1000),
     1    stdrdev1(1000),
     1    stdrdev2(1000), cov_cov(1001,1001)
      dimension idummy(1001,1001)
      dimension xstdrdev1(1000)
      dimension xcov_cov(1001,1001)
      dimension xhv(38226)
      dimension stddev1(1000), stddev2(1000)
      dimension engsnd(1001), emidsnd(1000) 
      real, dimension(1001) :: eigen_save
      character*145 idir, jdir, kdir, ldir2, ldir3, ldir4, job
      character*145 ldir5
      character*106 optical, xoptical
      character*15 ename
      character*250 file_name, new_file_name
      common /statinfo/ efactor, number_of_files, xoptical
      common /location / optical, idir, jdir, kdir
      character*250 outfile
      common /whatever/ outfile
      dimension xarray(1000), yarray(1000),itran(1000)
      common /datahld/ nenergy_hld(15), energy_hld(1001,15), 
     1          array_hld(1000,15)
     1         , emid_hld(1001,15)
      equivalence (alf(1),ilf(1)),(item(1),ytem(1))
      integer*4 icov_1, icov_2, icov_3, icov_4, ilocate
      integer*4 irenorm_flag
      equivalence (icov_1, lcov_1), (icov_2, lcov_2), 
     1     (icov_3, lcov_3),   (icov_4, lcov_4)
c      logical ok
c     character*4 inull
      character*4 lcov_1, lcov_2, lcov_3, lcov_4
      character*80 label
      logical pos_def
      real, allocatable    :: input_replace(:,:)
c
c     Control Flags
c
c        icov(1) = 1    propogate standard deviation of an integral
c                       cross section uncertainty component
c                       (first quantity has a covariance, second is fixed)
c                  2    input covariance file, output tecplot file
c                  3    combine covariance components
c                  4    propogate standard deviation of an integral
c                       spectrum uncertainty component
c                       (second quantity has a covariance, first is fixed)
c                  5    input covariance file supporting combination
c
c        icov(2)         override standard response (only for generated 640/770 grp)
c                = 0        no action, use data from covariance array
c                = 1        override response with more accurate info
c            
c
c        icov(3)         override standard cov (only for generated 640/770 grp)
c                = 0        no action, use data from covariance array
c                = 1        override cov/cor with fully uncorrelated data
c                                in native group structure
c                = 2        override cov/cor with fully uncorrelated data
c                                in 640/770 group structure frame-of-ref.
c                = 3        override cov/cor with fully correlated data
c                                in 640/770 group structure frame-of-ref.
c            
c
c        icovfmt        format for covariance read
c                = 0       known manipulate function - no covariance
c                = 1       covfils output summary (Rel. Stddev., Rel. Cov.)
c                = 2       SNLRML .lsl dosimetry xsec format 
c                          (self-covariance only)
c                          (correlation and Std. Dev. data)
c                = 3       manipulate snlcov format
c                = 4       known func. plus std. dev. - no correlations
c                = 5       SNLRML .lsl spectrum format 
c                          (self-covariance only)
c                          (correlation and Std. Dev. data)
c                = 6       ENDF/B-VI absolute covariance format (e.g. Cf252)
c                          (abs. cov., energy grid, number frac, and 
c                           rel. std. dev.)
c                = 7       PDF (Probability Distribution Function) for spectrum
c                          uncertainty - used for mono-energetic sources
c                = 8       fcov input format (no trial, only covariance)
c
c
c        icode          action for covread
c                = 0       leave in native energy group structure
c                = 1       expand covariance to 640 SAND-II energy grid
c                = 2       expand covariance to 770 extended SAND-II energy grid
c                = 3       use 89 NuGET neutron energy grid
c                = 4       use 48 NuGET gamma energy grid
c                = 5       use 175 Vitamin-J neutron energy grid
c
c
c
      if (icon(9) < 0) then 
          write (6,6712) icode, icovfmt
 6712     format (1x, 'COVREAD entered ', 2i5)
      endif
      iscale = 0
      ename = 'opt'
      lename = lnblnk(ename)
      call getenv(ename(1:lename), optical)
c      if (optical .eq. '') then 
c           optical=''
c           write (*,9013) optical
c 9013      format (1x, 'blank optical default filled', 1x, a6)
c      endif
      mblank2 = lnblnk(optical)
      ename = 'job'
      lename = lnblnk(ename)
      call getenv(ename(1:lename), job)
      jblank3= lnblnk(job)
c      jdir = '/app/manipulate-2/response/'
      jdir = 'response/'
      jblank2 = lnblnk(jdir)
c      kdir = '/app/manipulate-2/spectrum/'
      kdir = 'spectrum/'
      kblank2 = lnblnk(kdir)
c      ldir2 = '/app/manipulate-2/'
      ldir2 = ''
      lblank2 = lnblnk(ldir2)
c      ldir3 = '/app/NJOY-2012/correlation/'
      ldir3 = '..\NJOY-2012\correlation\'
      if ( icon(7) .eq. 1) then 
         ldir3 = '../NJOY2016/correlation/'
      endif
      lblank3 = lnblnk(ldir3)
c      ldir4 = '/app/covfil/output/'
      ldir4 = '../covfil/output/'
      lblank4 = lnblnk(ldir4)
c      ldir5 = '/app/lsl/library/'
      ldir5 = '../lsl/library/'
      lblank5 = lnblnk(ldir5)
      lfmt = lnblnk(fmt)
c
c     assume known manipulate response function with no uncertainty
c     assume SAND-II 640/770 group structure
c
      nfile = 34
      name = optical(1:mblank2)//jdir(1:jblank2)//'sand641.nrg'
      ilmt = 641
      if (icode .eq. 1) then
         name = optical(1:mblank2)//jdir(1:jblank2)//'sand641.nrg'
         ilmt = 641
      elseif ( icode .eq. 2) then 
         name = optical(1:mblank2)//jdir(1:jblank2)//'sand771.nrg'
         ilmt = 771
      elseif ( icode .eq. 3) then 
         name = optical(1:mblank2)//jdir(1:jblank2)//'nuget90.nrg'
         ilmt = 90
      elseif ( icode .eq. 4) then 
         name = optical(1:mblank2)//jdir(1:jblank2)//'nuget49.nrg'
         ilmt = 49
      elseif ( icode .eq. 5) then 
         name = optical(1:mblank2)//jdir(1:jblank2)//'vit176.nrg'
         ilmt = 176
      endif
      lend = lnblnk(name)
      open(unit=nfile, form='formatted',
     1      file=name(1:lend)
     2      , status='old', iostat=ilook, err=1909)
      read(nfile, 781) title
 781  format (a1)
      read (nfile,*) iengsnd
      if ( icon(12) .ne. 1) then 
c        read high to low in units of MeV
         read (nfile,*) (engsnd(jk), jk=ilmt,1,-1)
      else
c        read low to high energy in units of MeV
         read (nfile,*) (engsnd(jk), jk=1,ilmt)
      endif
      close (unit=nfile)
      if ( icon(9) < 0) then 
        write (6,8519) (engsnd(jk), jk=1,ilmt)
 8519   format (1x, 'engsnd energies ', 5g14.7)
      endif
      do ie=1,640
         emidsnd(ie) = 0.5*(engsnd(ie) + engsnd(ie+1))
      enddo
      if ( icon(9) < 0) then 
          write (6,568) icovfmt
 568      format (1x, 'COVREAD flags: ', 5i6)
      endif
      if ( icovfmt .eq. 0) then 
          iself = 0
          ifmt = lnblnk(fmt)
          name = optical(1:mblank2)//ldir2(1:lblank2)//fmt(1:ifmt)
          lend = lnblnk(name)
          write (nt6,8726) name(1:lend)
 8726     format (1x, 'covread open file = ', a90)
          open(unit=nfile, form='formatted',
     1         file=name(1:lend),
     2         status='old', iostat=ilook, err=1909)
c
c         assume low to high energy in MeV
c
          read(nfile,*) (covreng(jk), covrrsp1(jk), jk=1,640)
c
c         check energy grid
          if ( covreng(1) .ne. 0.10250e-09 .or. 
     1         covreng(2) .ne. 0.10750E-09 .or.
     2         covreng(640) .ne. 0.19950E-02 ) then
                write (nt6,9023) covreng(1), covreng(2), covreng(640)
 9023           format (1x, 'Energy grid error in covread ', 3g14.7)
                stop 'CORREAD eng grid'
          endif
c
c         move data to output arrays
c
          icoveng = 640
          if ( nenergy == 89) then 
            icoveng = 89
          endif
          if ( nenergy == 175) then 
            icoveng = 175
          endif
          do jk=1,640
             coveng(jk) = covreng(jk)
             covrsp1(jk) = covrrsp1(jk)
             covrsp2(jk) = covrrsp1(jk)
          enddo
c
c         fill-in default covariance array
c
          do jk1=1,icoveng
             do jk2 = 1,icoveng
                 corr(jk1,jk2) = 0.0
             enddo
          enddo
          do jk1 = 1,icoveng
             corr(jk1,jk1) = 1.0
          enddo
c
c         fill-in default std. dev. array      
c
          do jk=1,770
             stddev1(jk) = 0.0
             stddev2(jk) = 0.0
          enddo
c
c         create covariance array
c
          do jk1=1,icoveng
             do jk2 = 1,icoveng
                 covr(jk1,jk2) = 0.0
             enddo
          enddo
c
      elseif ( icovfmt .eq. 1) then 
c
c        read COVFIL output file
c
c        check for covariance data
         nfile = 38
         name = optical(1:mblank2)//ldir4(1:lblank4)//
     1          fmt(1:lfmt)//'.out'
         lend = lnblnk(name)
          open(unit=nfile, form='formatted',
     1         file=name(1:lend),
     2         status='old', iostat=ilook, err=1909)
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
 1923    format (a80)
         read (nfile,1924) icovreng1, idum2 ,idum3
 1924    format (54x, i4, 3x, i4, 3x, 4x, i4)
         if ( icon(9) < 0) then 
           write (6,671) icovreng1
 671       format (1x, 'covread icovreng1 = ', i5)
         endif
         if ( icovreng1 .gt. 641 .or. icovreng1 .le. 0) then
              write (6,7934) icovreng1
 7934         format (1x, '**** number of covariance points exceeds',
     1        ' dimension **** ', i10) 
              stop 'cov-pts'
         endif
         icovreng = icovreng1 - 1
c        energy grid
c        data stored high to low, we need low to high, so reverse
         read (nfile,7933) (covreng(jk),jk=icovreng1,1,-1)
         if ( covreng(1) .gt. 1.e-4 .or. 
     1        covreng(icovreng1) .lt. 19.6e6) then 
             write (nt6, 7231) covreng(icovreng1), covreng(1)
 7231        format (1x, 'COVFILs energy bounds not large enough ',
     1       2g14.7)
             stop 'covfil-45'
         endif
 7933    format (8e10.3)
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile,1924) idum1, idum2 ,idum3
         if ( idum3 .ne. icovreng) then 
              write (6,7932) idum3, icovreng
 7932         format (1x, 'Apparent conflict in covariance',
     1        ' energy groups',
     1        /, 10x, 'idum3 = ', i5, 'icovreng = ', i10)
              stop 'COVENG-1a'
         endif
c        cross section (b)
         read (nfile,7933) (covrrsp1(jk), jk=icovreng,1,-1)
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile,1924) idum1, idum2 ,idum3
         if ( idum3 .ne. icovreng) then 
              write (6,7932) idum3, icovreng
              stop 'COVENG-1'
         endif
c        relative standard deviation
         read (nfile,7933) (stdrdev1(jk), jk=icovreng,1,-1)
c        convert stddev into a percentage basis from a fractional basis
         do jk=1,icovreng
              stdrdev1(jk) = stdrdev1(jk)*100.0
         enddo
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile,1924) idum1, idum2 ,idum3
         if ( idum3 .ne. icovreng) then 
              write (6,7932) idum3, icovreng
              stop 'COVENG-1b'
         endif
c        cross section (b)
         read (nfile,7933) (covrrsp2(jk), jk=icovreng,1,-1)
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile,1924) idum1, idum2 ,idum3
         if ( idum3 .ne. icovreng) then 
              write (6,7932) idum3, icovreng
              stop 'COVENG-1'
         endif
c        relative standard deviation
         read (nfile,7933) (stdrdev2(jk), jk=icovreng,1,-1)
c        convert stddev into a percentage basis from a fractional basis
         do jk=1,icovreng
              stdrdev2(jk) = stdrdev2(jk)*100.0
         enddo
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile,1924) idum1, idum2 ,idum3
         if ( idum3 .ne. icovreng) then 
              write (6,7932) idum3, icovreng
              stop 'COVENG-1'
         endif
c
c        check for same function covariance
c
         iself = 0
         do jk=1,icovreng
           if ( covrrsp1(jk) .ne. covrrsp2(jk) ) iself = iself + 1
         enddo
c        relative covariance matrix
         do i1=icovreng,1,-1
            read (nfile,7933) (covr(i1,i2), i2=icovreng,1,-1)
         enddo
c        convert relative to absolute covariance matrix
         do jk1=1,icovreng
             do jk2=1,icovreng
                covr(jk1,jk2) = covr(jk1,jk2)*covrrsp1(jk1)*
     1          covrrsp2(jk2)
             enddo
         enddo
c        construct correlation matrix 
         do jk1=icovreng,1,-1
            dif1 = stdrdev1(jk1)*0.01*covrrsp1(jk1)
            do jk2=1,icovreng
               dif2 = stdrdev2(jk2)*0.01*covrrsp2(jk2)
               dif = dif1*dif2
               if ( dif .ne. 0.0) then
                  corr(jk1,jk2) = covr(jk1,jk2)/dif
               else
                   corr(jk1,jk2) = 0.0
               endif
               if ( jk1 .eq. jk2 .and. covr(jk1,jk2) .ne. 0.0 
     1             .and .iself .eq. 0 ) then 
                   zap = abs(corr(jk1,jk2) - 1.0)
                   if ( zap .gt. 1.e-3) then 
                      write (6,7931) jk1,jk2, corr(jk1,jk2), 
     1                covr(jk1,jk2), dif1, dif2, zap
 7931                 format (1x, 'CORR .NE. 1.0 check value ',
     1                2i5, 5g14.7)
c                      stop 'CORR ERROR 1'
                   endif 
                   corr(jk1,jk2) = 1.0
               endif
            enddo
         enddo
         close (unit=nfile)
      elseif ( icovfmt .eq. 2) then 
c
c        read in correlation coefficients from lsl format and convert
c        to covariance format.
c
c         nin = 641
c         mpl1 = nin
c        check for covariance data
         nfile = 38
         iself = 0
         name = optical(1:mblank2)//ldir3(1:lblank3)//
     1          fmt(1:lfmt)//'.lsl'
         lend = lnblnk(name)
         if ( icon(9) < -1) then 
            write (6,3651) mblank2, lblank3, lfmt, 
     1          optical(1:mblank2), ldir3(1:lblank3),
     1          fmt(1:lfmt)
 3651       format (1x, 'icovfmt=2 check ', 3i5,/, (5x, a))
         endif
c
c         over_ride = "G:\Projects\Manipulate-2010\snl-work\..\" //
c     &          "..\NJOY-2012\correlation\" //
c     &          "ni58_IRDF2002_final_Irdf2002p.dat_103_cov_tpl.lsl"
c         name = over_ride
         lend = lnblnk(name)
c
          open(unit=nfile, form='formatted',
     1         file=name(1:lend),
     2         status='old', iostat=ilook, err=1909)
         if ( icon(9) <= -1 ) then 
            write (6,5449) lend, name(1:lend)
 5449       format (1x, 'LSL file: ', i5, 1x, a   )
         endif
         read (nfile,5) alf
         read (nfile,5) alf
    5    format(80a1)
         read (nfile,*) icovreng1
         icovreng = icovreng - 1
         if ( icon(9) <= 0 ) then 
            write (6,5639) icovreng1
 5639       format (1x, 'ICOVRENG1: ', i9,/,
     &             (1x, '        ', 6g14.7)   )
         endif
c         if ( icode == 1 .and. icovreng1 == 90) then 
c           over-ride icode for this case from 640 to 89 group structure
c            icode = 3
c         endif
         icov_reverse = 0
         if ( icovreng1 < 0) then 
c          set flag to read in covariance matrix from high to low energy
           icov_reverse = -1
           icovreng1 = iabs(icovreng1)
           icovreng = icovreng - 1
         endif
         if ( icovreng1 .gt. 661 .or. icovreng1 .le. 0) then
              write (6,7904) icovreng1, name
 7904         format (1x, '**** number of covariance points exceeds',
     1        ' dimension (b) **** ', i10, 1x, a150) 
              stop 'cov-pts-1'
         endif
         icovreng = icovreng1 - 1
         read (nfile,5) alf

         if ( icov_reverse == 0) then 
            read (nfile,*) (covreng(jk),jk=1,icovreng1)
         elseif ( icov_reverse == -1) then 
            read (nfile,*) (covreng(jk),jk=icovreng1, 1, -1)
         else
            write (6,672) icov_reverse
 672        format ( 1x, 'icov_reverse invalid ', i5)
         endif
         if ( icon(9) <= -1 ) then 
            write (6,5623) (covreng(jk), jk=1,icovreng1)
 5623       format (1x, 'COVRENG: ', 6g14.7,/,
     &             (1x, '         ', 6g14.7)   )
         endif
         if ( icon(5) .eq. 3) then 
            if ( covreng(1) .gt. 1000. .or. 
     1           covreng(icovreng1).lt. 50.0e6) then 
               write (nt6, 5231) covreng(icovreng1), covreng(1)
               stop 'covfil-46a1'
            endif
         else
            if ( covreng(1) .gt. 1.e-4 .or. 
     1           covreng(icovreng1).lt. 19.6e6) then 
               write (nt6, 5231) covreng(icovreng1), covreng(1)
 5231          format (1x, 'LSL energy bounds not large enough ',
     1         2g14.7)
               stop 'covfil-46a'
            endif
         endif
         read (nfile,5) alf
         if ( icon(9) .lt. -1) then 
           write (6,4812) alf
         endif
         if ( icov_reverse == 0) then 
            read (nfile,*) (covrrsp1(jk), jk=1,icovreng)
         elseif ( icov_reverse == -1) then 
            read (nfile,*) (covrrsp1(jk), jk=icovreng, 1, -1)
         else
            write (6,672) icov_reverse
         endif
         if ( icon(9) <= 0 ) then 
            write (6,5679) (covrrsp1(jk), jk=1,icovreng)
 5679       format (1x, 'values: ', 6g14.7,/,
     &             (1x, '        ', 6g14.7)   )
         endif
c
c        If LSL-format spectrum - these are number fractions
c        ensure they are normalized to unity
c
         if ( icon(9) <= 0) then 
            write (6,5619) icov(1), icovreng, nenergy
 5619       format (1x, 'COVREAD icov(1) check: ', 5i6)
         endif
         if ( icov(1) .eq. 4) then 
            sumx = 0.0
            do jk = 1,icovreng
              sumx = sumx + covrrsp1(jk)
            enddo
            gsum = abs(sumx - 1.0)
            if ( gsum .gt. 1.1 .or. gsum .lt. 0.9 ) then
               write (nt6, 7243) sumx
 7243          format (1x, 'Spectrum number fractions renormalized',
     1         ', sum = ', g14.7)
            endif
            do jk = 1,icovreng
              covrrsp1(jk) = covrrsp1(jk)/sumx
            enddo
         endif
         do jk=1,icovreng
            covrrsp2(jk) = covrrsp1(jk)
         enddo
         read (nfile,5) alf
         if ( icon(9) .lt. -1) then 
           write (6,4812) alf
         endif
         if ( icov_reverse == 0) then 
            read (nfile,*) (stdrdev1(jk), jk=1,icovreng)
         elseif ( icov_reverse == -1) then 
            read (nfile,*) (stdrdev1(jk), jk=icovreng, 1, -1)
         else
            write (6,672) icov_reverse
         endif
         do jk=1,icovreng
            stdrdev2(jk) = stdrdev1(jk)
         enddo
         if ( icon(9) <= -1 ) then 
            write (6,5629) (stdrdev1(jk), jk=1,icovreng)
 5629       format (1x, 'STDDEV: ', 6g14.7,/,
     &             (1x, '        ', 6g14.7)   )
         endif
         read (nfile,5) alf
         if ( icon(9) .lt. -1) then 
           write (6,4812) alf
 4812      format (1x, 'alf separator: ', 80a1)
         endif
         if ( icov_reverse == 0) then 
            do i1=1,icovreng
               read (nfile,*) (cov_cov(i1,i2), i2=i1,icovreng)
c              note correlations below threshold are zero
                if ( cov_cov(i1,i1) .eq. 0.) then
                      cov_cov(i1,i1) = 100.
                      write (6,6672) i1
6672                  format (1x, 'Warning: zero relative correlation',
     &                'coefficient modified ', g14.7)
                endif
                if ( cov_cov(i1,i1) .eq. 1.0e4) then
                    if ( iscale < 1) then
                      iscale = iscale + 1
                      write (6,5672) i1
5672                  format (1x,'Warning: 1.0E4 relative correlation',
     &                ' coefficient modified ', g14.7)
                    endif
                    do i2=i1, icovreng
                      cov_cov(i1,i2) = cov_cov(i1,i2)*0.01
                    enddo
                endif
            enddo    
            do i1 = 1, icovreng
               do i2 = i1, icovreng
                 cov_cov(i2,i1) = cov_cov(i1,i2)
               enddo
            enddo               
         elseif ( icov_reverse == -1) then 
            if ( icon(9) .lt. 0) then 
              write (6,782)
 782          format (1x, '*** WARNING *** correlation matrix',
     &          ' MUST be ',
     &        'reversed in energy' )
              write (6,7834) icovreng 
 7834         format (1x, 'Debug icovreng: ', i5)
            endif
            do i1=icovreng, 1, -1
               read (nfile,*) (cov_cov(i1,i2), i2=i1, 1, -1)
               if ( icon(9) .lt. -1) then
                 write (6,6723) i1, (cov_cov(i1,i2), i2=i1, 1, -1)
 6723            format (1x, 'cov ', i5, 5g14.7, /, (9x, 5g14.7)  )
               endif
c              note correlations below threshold are zero
c                if ( cov_cov(i1,i1) .eq. 0.) cov_cov(i1,i1) = 1.00
                if ( cov_cov(i1,i1) .eq. 0.) then 
                      cov_cov(i1,i1) = 1.00
                      write (6,4672) i1
4672                  format (1x, 'Warning: zero relative correlation',
     &                'coefficient rev modified ', g14.7)
                endif
                if ( cov_cov(i1,i1) .eq. 1.0e4)  then 
                      write (6,3672) i1
3672                  format (1x, 'Warning: 1.0E4 relative correlation',
     &                ' coefficient rev modified ', g14.7)
                      do i2=i1, 1, -1
                        cov_cov(i1,i2) = cov_cov(i1,i2)*0.01
                      enddo
                endif
            enddo
c           complete the triangular matrix
            do i1 = icovreng, 1, -1
               do i2 = i1, 1, -1
                 cov_cov(i2,i1) = cov_cov(i1,i2)
               enddo
            enddo          
         else
            write (6,672) icov_reverse
         endif
         if ( icon(9) .lt. 0) then 
           write (6,5490) cov_cov(1,1)
 5490      format (1x, 'covread renorm first element: ', g14.7)
         endif
         irenorm_flag = 0
         if ( cov_cov(1,1) .gt. 1.1 .and. cov(1,1) .le. 110.0) then
           irenorm_flag = 1
         endif
         if ( cov_cov(1,1) .gt. 110.0 .and. cov(1,1) .le. 110.0E2) then
           irenorm_flag = 2
         endif
         do jk1=1,icovreng
            if ( icon(9) .lt. 0) then 
              write (6,5491) jk1, irenorm_flag, cov_cov(jk1,jk1)
 5491         format (1x, 'covread renorm first element: ', 
     &                2i5, 1x, g14.7)
            endif
c            if ( cov_cov(jk1,jk1) == 100.) then
            if ( irenorm_flag .eq. 1) then
c              re-normalize correlation matrix back to unity on the diagonal
               do jk2=1,icovreng
                  if ( jk1 .ge. jk2) then
                      corr(jk1,jk2) = 0.01*cov_cov(jk2,jk1)
                  else
                      corr(jk1,jk2) = 0.01*cov_cov(jk1,jk2)
                  endif
               enddo
c            elseif ( cov_cov(jk1,jk1) == 10000.) then
            elseif ( irenorm_flag .eq. 2) then
c              re-normalize correlation matrix back to unity on the diagonal
               do jk2=1,icovreng
                  if ( jk1 .ge. jk2) then
                      corr(jk1,jk2) = 0.0001*cov_cov(jk2,jk1)
                  else
                      corr(jk1,jk2) = 0.0001*cov_cov(jk1,jk2)
                  endif
               enddo
            else
c              make the correlation matrix a real number
               do jk2=1,icovreng
                  if ( jk1 .ge. jk2) then
                      corr(jk1,jk2) = cov_cov(jk2,jk1)
                  else
                      corr(jk1,jk2) = cov_cov(jk1,jk2)
                  endif
               enddo
            endif
         enddo
         close (unit=nfile)
c
c        Check correlation matrix for eigenvalues
c
         allocate( input_replace(1001, 1001),
     $         STAT=iStat)
         if (iStat/=0) then 
           write (6,25000) iStat
           stop 'Memory-allocation error in covread.'
         endif
25000   format ('Subroutine covread memory-(de)allocation error: ',
     &          2i5)
         ipass = 0
 1010    continue
         ipass = ipass + 1
         if ( icon(9) < 0) then 
            write (6,7623) ipass
 7623       format (1x, 'covread ipass loop: ', i5)
         endif
         input_replace = 0
         call eigen_out(corr, icovreng, input_replace, pos_def, ipass, 
     &       eigen_save)
         if ( icon(9) < 0) then 
            write (6,7829) icovreng, nenergy, 
     &           icov(1), ipass, pos_def
 7829       format (1x, 'COVREAD after eigen_out ', 4i5, l)
         endif
         if ( icon(9) .lt. -1) then 
            write (6,3482)
 3482       format (1x, 'post eigen_out covread corr reset ')
            write (6,5481) (corr(1,jk), jk=1,5)
 5481       format (1x, ' post eigen_out corr = ', 7g14.7)
            write (6,3481) (input_replace(1,jk), jk=1,5)
 3481       format (1x, ' post eigen_out input_replace = ', 7g14.7)
         endif
         if ( ipass < 750 .and. .not. pos_def) then
c             enforce fluence unity normalization on the replacement 
c             matrix - but only for the first 50 tries at positive definite
              if ( icon(9) < -1) then 
                 write (6,78291) 
78291            format (1x, 'COVREAD first inner ipass loop ', 5i5)
              endif
              if ( ipass < 50) then 
                if ( icon(9) < -1) then 
                   write (6,78293) ipass, icov(1)
78293              format (1x, 'COVREAD inner-a ipass loop ', 5i5)
                endif
                if ( icov(1) /= 4) then 
                   if ( icon(9) < 0) then 
                     write (6,78294) ipass, icov(1)
78294                format (1x, 'COVREAD inner-b ipass loop ', 5i5)
                   endif
                   if ( ipass < 2 ) write (6,5698)
 5698              format (1x, 'Covariance is for a cross section ', 
     &             'so do not enforce normalization condition ')
                else
                  if ( icon(9) < -1) then 
                     write (6,78295) 
78295                format (1x, 'COVREAD inner-c ipass loop ', 5i5)
                  endif
c                 NOTE: these lines were commented out for responses
                  call cov_norm(icovreng, covrrsp1, stdrdev1, 
     &                 input_replace, ipass)
                endif
              endif
              if ( icon(9) .lt. -1) then 
                 write (6,3402)
 3402            format (1x, 'covread corr reset ')
                 write (6,5401) (corr(1,jk), jk=1,5)
 5401            format (1x, ' corr = ', 7g14.7)
                 write (6,3401) (input_replace(1,jk), jk=1,5)
 3401            format (1x, ' input_replace = ', 7g14.7)
              endif
              corr = input_replace
              go to 1010
         elseif ( ipass == 1 .and. icov(1) == 4) then 
              if ( icon(9) < -1) then 
                 write (6,78292) 
78292            format (1x, 'COVREAD second inner ipass loop ', 5i5)
              endif
c
c             since corr matrix was not changed, fill it for the 
c             normalization test with the original
c
              input_replace = corr
c
              call cov_norm(icovreng, covrrsp1, stdrdev1, 
     &             input_replace, ipass)
              if ( icon(9) < -1) then 
                 write (6,78297) 
78297            format (1x, 'COVREAD third inner ipass loop ', 5i5)
              endif
              corr = input_replace
              go to 1010
         endif
         ineg = 0
         do jk=1,icovreng
           if (eigen_save(jk) <= 0.0 ) then
            ineg = ineg + 1
           endif
         enddo
         if ( icon(9) < 2) then 
            write (6,102) ipass, ineg, ( eigen_save(jk), jk=1,
     $               icovreng )
         endif
 102     format (/,1x, 'Eigenvalue passes = ', i5, ' with ', i5, 
     &                 ' final negative entries ',/,
     &             1x, "Final Eigenvalues: ", 5g14.7,/,
     &       (     1x, "                   ", 5g14.7)   )
         if ( icon(9) <= -1) then 
            write (6,7823) icov_reverse, ipass,icovreng1, icovreng
 7823       format (1x, 'eigen_out return in covread ', 4i5)
            if ( icov_reverse == -1) then 
c             write out reversed version of covariance matrix
              open(unit=68, form='formatted',
     1         file='correlation_reverse.lsl',
     2         status='unknown', iostat=ilook, err=1909)
              label = 'replacement correlation'
              write (68, 901) label 
              label = 'energies covread_1'
              write (68, 901) label 
              write (68, 17901) (covreng(jk),jk=icovreng1, 1, -1)
              label = 'values'
              write (68, 901) label 
              write (68, 17901) (covrrsp1(jk), jk=icovreng, 1, -1)
              label = 'std. dev.'
              write (68, 901) label 
              write (68, 17901) (stdrdev1(jk), jk=icovreng, 1, -1)
              label = 'correlation x 1'
              write (68, 901) label 
              do jk1 = icovreng, 1, -1
                 write (68, 17901) (corr(jk1, jk2)*1., 
     &                  jk2 = jk1, 1, -1)
              enddo
              close (unit=68)
            endif
         endif
c
         if ( .not. pos_def .or. ipass > 1) then
c           write out a copy of a proper lsl-format covariance matrix
            write (6,6729) pos_def, ipass
 6729       format (/,/,1x, 'correlation replacement matrix written ',
     *       l, 2x, i5, /, /)
            open(unit=68, form='formatted',
     1         file='correlation_replacement.lsl',
     2         status='unknown', iostat=ilook, err=1909)
            label = 'energies from covread'
            write (68, 901) label 
            write (68, 17901) (covreng(jk),jk=1,icovreng1)
            label = 'values'
            write (68, 901) label 
            write (68, 17901) (covrrsp1(jk), jk=1,icovreng)
            label = 'std. dev.'
            write (68, 901) label 
            write (68, 16990) (stdrdev1(jk), jk=1,icovreng)
            label = 'replacement correlation'
            write (68, 901) label 
 901        format (a)
            do jk1 = 1, icovreng
               write (68, 6990) (corr(jk1, jk2)*1., 
     &                  jk2 = jk1, icovreng)
 6990          format((1x,8(1pe12.5,1x)))
16990          format((1x,5(1pe12.5,1x)))
 7901          format (2x, g14.7, 2x, g14.7, 2x, g14.7, 2x, 
     &              g14.7, 2x, g14.7, 2x, g14.7, 2x, g14.7)
17901          format (2x, g14.7, 2x, g14.7, 2x, g14.7, 2x, 
     &              g14.7, 2x, g14.7)
            enddo
            close (unit=68)
c
c           write out correlation matrix in SigmaPlot format
c
            open(unit=68, form='formatted',
     1         file='correlation_replacement.sigplt',
     2         status='unknown', iostat=ilook, err=1909)
            do jk1=1, icovreng
                write (68,2295) covreng(jk1), covreng(jk1),
     &                   (corr(jk1,jk2), jk2 = 1, icovreng)
 2295           format (1x, g12.4, 90(', ', g12.4, 2x))
            enddo
            close (unit=68)
c
c           write out percent stdard devaiton for SigmaPlot input
c
c
c           write out correlation matrix in SigmaPlot format
c
            open(unit=68, form='formatted',
     1         file='pct_stddev_replacement.sigplt',
     2         status='unknown', iostat=ilook, err=1909)
            do jk1=1, icovreng
                write (68,3295) stdrdev1(jk1)
 3295           format (1x, g14.7)
            enddo
            close (unit=68)
c
         endif
c
         deallocate(input_replace, STAT=iStat)
         if (iStat/=0) then 
           write (6,25000) iStat
           stop 'Memory-(de)allocation error in eigen_out.'
         endif
         if ( icon(9) <= -1 ) then 
            write (6,9679) (covrrsp1(jk), jk=1,icovreng)
 9679       format (1x, 'after dalocate values: ', 6g14.7,/,
     &             (1x, '                       ', 6g14.7)   )
         endif
          
c               TBD
c
c        Check for over-ride of covariance matrix
c
      elseif ( icovfmt .eq. 3) then 
c
c        read manipulate snlcov  file
c
         nfile = 38
         name = optical(1:mblank2)//jdir(1:jblank2)//
     1          'covar/'//fmt(1:lfmt)//'.snlcov'
         lend = lnblnk(name)
          open(unit=nfile, form='formatted',
     1         file=name(1:lend),
     2         status='unknown', iostat=ilook, err=1909)
c        read snlcov type file
         read (nfile,8011) icovreng
         if ( icovreng .gt. 640 .or. icovreng .le. 0) then
              write (6,7904) icovreng, name
              stop 'cov-pts-2'
         endif
         icovreng1 = icovreng+1
         read (nfile,8012) (covreng(jk),jk=1,icovreng+1)
         if ( covreng(1) .gt. 1.e-4 .or. 
     1        covreng(icovreng1).lt. 19.6e6) then 
             write (nt6, 6231) covreng(icovreng1), covreng(1)
 6231        format (1x, 'SNLCOV energy bounds not large enough ',
     1       2g14.7)
              stop 'covfil-47'
         endif
         read (nfile,8012) (covrrsp1(jk),jk=1,icovreng)
         read (nfile,8012) (stdrdev1(jk),jk=1,icovreng)
         read (nfile,8012) (covrrsp2(jk),jk=1,icovreng)
         read (nfile,8012) (stdrdev2(jk),jk=1,icovreng)
c
c        check for same function covariance
c
         iself = 0
         do jk=1,icovreng
           if ( covrrsp1(jk) .ne. covrrsp2(jk) ) iself = iself + 1
         enddo
c
         read (nfile,8012) ((covr(jk1,jk2), jk1=1,icovreng),
     1          jk2=1,icovreng)
         read (nfile,8012) ((corr(jk1,jk2), jk1=1,icovreng), 
     1         jk2=1,icovreng)
         close (unit=nfile)
 8934        format (1x, 'No covariance data in COVREAD ')
      elseif ( icovfmt .eq. 5) then 
c
c        read in correlation coefficients from lsl format and convert
c        to covariance format.
c          ----- spectrum input format ---
c
c         nin = 641
c         mpl1 = nin
c        check for covariance data
         nfile = 38
         iself = 0
         name = optical(1:mblank2)//ldir3(1:lblank3)//
     1          fmt(1:lfmt)//'.lsl'
         lend = lnblnk(name)
          open(unit=nfile, form='formatted',
     1         file=name(1:lend),
     2         status='old', iostat=ilook, err=1909)
         read (nfile,5) alf
         read (nfile,5) alf
         read (nfile,*) icovreng1
         if ( icovreng1 .gt. 641 .or. icovreng1 .le. 0) then
              write (6,7904) icovreng1, name
              stop 'cov-pts-1'
         endif
         icovreng = icovreng1 - 1
         if ( icon(9) < 0) then 
            write (6,4519) icovreng, icovfmt
 4519       format (1x, 'COVREAD icovfmt=5 portion ', 5i5)
         endif
         read (nfile,5) alf
         read (nfile,*) (covreng(jk),jk=1,icovreng1)
         if ( covreng(1) .gt. 1.e-4 .or. 
     1        covreng(icovreng1).lt. 19.6e6) then 
             write (nt6, 5231) covreng(icovreng1), covreng(1)
              stop 'covfil-46b'
         endif
         read (nfile,5) alf
         read (nfile,*) (covrrsp1(jk), jk=1,icovreng)
         if ( icon(9) .le. 0) then 
           write (6,2614) (covrrsp1(jr), jr=1,icovreng)
 2614      format (1x, 'covread lsl input: ', 5g14.7)
         endif
c
c        If LSL-format spectrum - these are number fractions
c        ensure they are normalized to unity
c
         if ( icov(1) .eq. 4 .or. icov(1) .eq. 2) then 
            sumx = 0.0
            do jk = 1,icovreng
              sumx = sumx + covrrsp1(jk)
            enddo
            gsum = abs(sumx - 1.0)
            if ( gsum .gt. 1.1 .or. gsum .lt. 0.9 ) then
               write (nt6, 7243) sumx
            endif
            do jk = 1,icovreng
              covrrsp1(jk) = covrrsp1(jk)/sumx
            enddo
         endif
         do jk=1,icovreng
            covrrsp2(jk) = covrrsp1(jk)
         enddo
         read (nfile,5) alf
         read (nfile,*) (stdrdev1(jk), jk=1,icovreng)
         do jk=1,icovreng
            stdrdev2(jk) = stdrdev1(jk)
         enddo
         read (nfile,5) alf
         do i1=1,icovreng
            read (nfile,*) (cov_cov(i1,i2), i2=i1,icovreng)
c           note correlations below threshold are zero
c             if ( cov_cov(i1,i1) .eq. 0) cov_cov(i1,i1) = 100
         enddo
         do jk1=1,icovreng
            do jk2=1,icovreng
               if ( jk1 .ge. jk2) then
                   corr(jk1,jk2) = 0.01*cov_cov(jk2,jk1)
               else
                   corr(jk1,jk2) = 0.01*cov_cov(jk1,jk2)
               endif
            enddo
         enddo
         close (unit=nfile)
c
c        Check correlation matrix for eigenvalues
c
c         call eigen_out(corr, icovreng)
         allocate( input_replace(1001, 1001),
     $         STAT=iStat)
         if (iStat/=0) then 
           write (6,25001) iStat
           stop 'Memory-allocation error in covread.'
         endif
25001    format ('Subroutine covread memory-(de)allocation error: ',
     &          2i5)
         ipass = 0
         input_replace = 0.0
         call eigen_out(corr, icovreng, input_replace, pos_def, ipass, 
     &       eigen_save)
c
c        Check for over-ride of covariance matrix
c
         deallocate(input_replace, STAT=iStat)
         if (iStat/=0) then 
           write (6,25000) iStat
           stop 'Memory-(de)allocation error in eigen_out.'
         endif
      elseif ( icovfmt .eq. 6) then 
c
c        read in correlation coefficients from ENDF format and convert
c        to covariance format.
c            energy grid is in eV from low to high energy
c
c          ----- spectrum input format ---
c
c        check for covariance data
         nfile = 38
         iself = 0
         name = optical(1:mblank2)//ldir3(1:lblank3)//
     1          fmt(1:lfmt)//'.endf'
         lend = lnblnk(name)
          open(unit=nfile, form='formatted',
     1         file=name(1:lend),
     2         status='old', iostat=ilook, err=1909)
         read (nfile,5) alf
         read (nfile,*) icovreng1
         if ( icovreng1 .gt. 641 .or. icovreng1 .le. 0) then
              write (6,7904) icovreng1, name
              stop 'cov-pts-1'
         endif
         icovreng = icovreng1 - 1
         read (nfile,*) (idum, dum1, dum2, covrrsp1(jk), 
     1      xstdrdev1(jk), idum2,idum3,idum4, jk=1,icovreng)
         read (nfile,5) alf
         read (nfile,5) alf
         read (nfile,5) alf
         read (nfile,5) alf
         read (nfile,5) alf
         ilist = icovreng*(icovreng+1)/2
         if ( ilist .gt. 38226) then 
            stop 'ilist'
         endif
         read (nfile,5681) (covreng(jk), jk=1,icovreng1), 
     1     (xhv(jk),jk=1,ilist)
c
c        Define new rel. std. dev. that have better norm
c
         inc = 0
         do jk1=1,icovreng
            do jk2=jk1,icovreng
               inc = inc + 1
               if ( jk1 .eq. jk2) then 
                  dummx = covrrsp1(jk1)*covrrsp1(jk2)*xstdrdev1(jk1)*
     1                    xstdrdev1(jk2)*0.01*0.01
                  corx = xhv(inc)/dummx
                  dif = abs(1.0 - corx)
                  if ( corx .gt. 1.1 .or. corx .lt. 0.9) then 
                     write (nt6, 8362) jk1, dif, corx
 8362                format (1x, 'STD RENORM OFF ', i5, g14.7, g14.7)
                  endif
                  stdrdev1(jk1) = xstdrdev1(jk1)*sqrt(corx)
               endif
            enddo
         enddo
         inc = 0
         do jk1=1,icovreng
            do jk2=jk1,icovreng
               inc = inc + 1
c
c              Convert absolute covariance to relative correlation
c
               dummx = covrrsp1(jk1)*covrrsp1(jk2)*stdrdev1(jk1)*
     1                 stdrdev1(jk2)*0.01*0.01
               xcov_cov(jk1,jk2) = xhv(inc)/dummx
               xcov_cov(jk2,jk1) = xcov_cov(jk1,jk2)
            enddo
         enddo
c     1    ((xcov_cov(i1,i2), i1=12,icovreng), i2=1,icovreng)
 5681    format (6e11.6)
         if ( covreng(1) .gt. 1.e-4 .or. 
     1        covreng(icovreng1).lt. 19.6e6) then 
             write (nt6, 5231) covreng(icovreng1), covreng(1)
              stop 'covfil-46c'
         endif
c
c        If LSL-format spectrum - these are number fractions
c        ensure they are normalized to unity
c
         if ( icov(1) .eq. 4) then 
            sumx = 0.0
            do jk = 1,icovreng
              sumx = sumx + covrrsp1(jk)
            enddo
            gsum = abs(sumx - 1.0)
            if ( gsum .gt. 1.1 .or. gsum .lt. 0.9 ) then
               write (nt6, 7243) sumx
            endif
            do jk = 1,icovreng
              covrrsp1(jk) = covrrsp1(jk)/sumx
            enddo
         endif
         do jk=1,icovreng
            covrrsp2(jk) = covrrsp1(jk)
         enddo
         do jk=1,icovreng
            stdrdev2(jk) = stdrdev1(jk)
         enddo
         do jk1=1,icovreng
            do jk2=1,icovreng
               if ( jk1 .ge. jk2) then
                   corr(jk1,jk2) = xcov_cov(jk2,jk1)
               else
                   corr(jk1,jk2) = xcov_cov(jk1,jk2)
               endif
            enddo
         enddo
         do i1=1,icovreng
c           note correlations below threshold are zero
             if ( corr(i1,i1) .eq. 0) corr(i1,i1) = 1.00
         enddo
         close (unit=nfile)
c
c        print LSL format covariance in original energy grid
c
         write (nt6, 7432) 
 7432    format (/,/,1x, 
     1    'Print input covariance data in LSL - ',
     1     'not SNL-LSL format ',/,/)
c                                                                    
c Title cards                                                       
c                                                                  
         open (unit=78, file='covariance.lsl', status='unknown')
         write(78,8562)                                                
 8562    format(    '*COR    (LIBRARY)    (MAT.#)    (TEMP)K')         
c                                                                    
c Energy grid                                                         
c                                                                   
         write(78,8996)                                                 
 8996    format('*Number of Energies plus 1')                      
         write(78,8995) icovreng1                                          
         write(78,8993)                                                 
 8993    format('*Energy Grid ( eV )')          
 8995    format(i5)                             
         write(78,8990) (covreng(jk),jk=1,icovreng1)  
c                                                  
c Cross section                                                       
c                                                                    
         write(78,8998)                                                 
 8998    format('*Number Fractions ')                           
         write(78,7454) (covrrsp1(jk), jk=1,icovreng)                                                     
c                                                                   
c Standard deviation                                               
c                                                                   
         write(78,8991)                                            
 8991    format('*% Standard Deviation')                            
         write(78,8990) (stdrdev1(i1),i1=1,icovreng)                           
c                                                               
c Write standard deviation to another file for input           
c   into a plotting package                                   
c                                                              
         open (unit=77, file='covariance.stdout', status='unknown')
         write(77,9101)                                         
 9101    format(' Standard Deviation',/,                         
     1          ' (input is in _templegraph_ format)',/)        
         do i1 = 1,icovreng                                          
            write(77,9102) covreng(i1), stdrdev1(i1)                        
            write(77,9102) covreng(i1+1), stdrdev1(i1)                    
         end do                                                 
 9102    format(2(5x,1pe10.3))                                     
         close(unit=77)                                           
c                                                                  
c Correlation coefficients                                        
c                                                                  
         write(78,8992)                                              
 8992    format('*Correlation Coefficient -- Upper Triangular')       
         do i1 = 1,icovreng       
            write(78,8990) ((corr(i1,i2)),i2=i1,icovreng)              
         end do                                                       
 8990    format((1x,8(1pe10.3,1x)))                                    
 7454    format((1x,1p8e10.3))                                         
 7990    format((1x,15(i4,1x)))                                   
         close (unit=78)
      elseif ( icovfmt .eq. 7) then 

c
c        read in PDF for uncertainty calculation - ignore
c        true covariance input.
c
c          ----- mono-energetic spectrum input format ---
c
c        read-in normal spectrum file
c
         icovreng = 640
         imode = 12
         iselect_material = -1
         fraction = 1.0
         leopt = lnblnk(xoptical)
         jblank2 = lnblnk(jdir)
         file_name = fmt(1:lfmt)//'.ref-list'
         ireverse= iselect_material
         if ( iselect_material.le.0 ) iselect_material = 1
         call prune (file_name, ilead, itrail, length, nchar)
         ierr = 0
         iflast = lnblnk(file_name)
         new_file_name = xoptical(1:leopt)//jdir(1:jblank2)//
     1           file_name(1:iflast)
         if ( icov(1) .eq. 4) then 
c             spectrum NOT response here
              new_file_name = xoptical(1:leopt)//kdir(1:kblank2)//
     1            file_name(1:iflast)
         else
              stop 'ERR 3, icov invalid'
         endif
         if ( icov(1) .eq. 4) then
            call filein(imode, new_file_name, ierr, ireverse)
         else
            call filein_response(imode, new_file_name, 
     &          ierr, ireverse)
         endif
c
c        Dump array_hld into covrsp1
c           convert to number fraction
c
         ip = iabs(iselect_material)
         do jk1=1,icovreng
            covrrsp1(jk1) = array(jk1,ip)
            stdrdev1(jk1) = 0.0
            covreng(jk1) = energy(jk1)*1.e6
         enddo
         covreng(icovreng+1) = energy(icovreng+1)*1.e6
c
c        read-in PDF data
c
         nfile = 38
         iself = 0
         name = optical(1:mblank2)//ldir3(1:lblank3)//
     1          fmt(1:lfmt)//'.pdf'
         lend = lnblnk(name)
          open(unit=nfile, form='formatted',
     1         file=name(1:lend),
     2         status='old', iostat=ilook, err=1909)
         read (nfile,1923) dummy
         read (nfile,*) ipdf
         if ( ipdf .gt. 1000 .or. ipdf .le. 0) then 
              stop 'pdf dim'
         endif
         read (nfile,*) (engpdf(jk), outpdf(jk), jk=1,ipdf)
         close (unit=nfile)
c        linear interpolation along pdf values
         do jk=1,icovreng
            pdf(jk) = fitmd(covreng(jk), ipdf, engpdf,
     1        outpdf, 2)
         enddo
c        renorm pdf
         sum = 0.0
         do jk=1,icovreng
             sum = sum + pdf(jk)
         enddo
         do jk=1,icovreng
             pdf(jk) = pdf(jk)/sum
         enddo
c
c        spectrum is number fraction - ensure normalized
c
         if ( icov(1) .eq. 4) then 
            sumx = 0.0
            do jk = 1,icovreng
              sumx = sumx + covrrsp1(jk)
            enddo
            gsum = abs(sumx - 1.0)
            if ( gsum .gt. 1.1 .or. gsum .lt. 0.9 ) then
               write (nt6, 7243) sumx
            endif
            do jk = 1,icovreng
              covrrsp1(jk) = covrrsp1(jk)/sumx
            enddo
         endif
         do jk=1,icovreng
            covrrsp2(jk) = covrrsp1(jk)
         enddo
         do jk=1,icovreng
            stdrdev2(jk) = stdrdev1(jk)
         enddo
         do jk1=1,icovreng
            do jk2=1,icovreng
               if ( jk1 .eq. jk2) then
                   corr(jk1,jk2) = 1.0
               else
                   corr(jk1,jk2) = 0.0
               endif
            enddo
         enddo
c
c         move data to output arrays
c
          icoveng = 640
          if ( icon(9) < 0) then 
            write (6,562) icoveng, nenergy
 562        format (1x, 'covread energy bounds check: ', 2i6)
          endif
          if ( nenergy == 89 .or. icovreng == 89) then 
            icoveng = 89
          endif 
          if ( nenergy == 175 .or. icovreng == 175) then 
            icoveng = 175
          endif
          do jk=1,640
             coveng(jk) = covreng(jk)
             covrsp1(jk) = covrrsp1(jk)
             covrsp2(jk) = covrrsp1(jk)
             stddev1(jk) = stdrdev1(jk)
             stddev2(jk) = stdrdev2(jk)
          enddo
c
c         fill-in default covariance array
c
          do jk1=1,icoveng
             do jk2 = 1,icoveng
                 corr(jk1,jk2) = 0.0
             enddo
          enddo
          do jk1 = 1,icoveng
             corr(jk1,jk1) = 1.0
          enddo
      elseif ( icovfmt .eq. 8) then 
c
c        read fcov output file
c
c        check for covariance data
         nfile = 38
         name = optical(1:mblank2)//ldir5(1:lblank5)//
     1          fmt(1:lfmt)//'.lib'
         lend = lnblnk(name)
          open(unit=nfile, form='formatted',
     1         file=name(1:lend),
     2         status='old', iostat=ilook, err=1909)
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, 1923) dummy
         read (nfile, *) id1, id2
         read (nfile, *) (xdummy, ik=id1,id2)
         read (nfile,6924) icovreng1
 6924    format (i10)
         if ( icovreng1 .gt. 641 .or. icovreng1 .le. 0) then
              write (6,7934) icovreng1
              stop 'cov-pts'
         endif
         icovreng = icovreng1 - 1
c        energy grid
c        data stored low to high
 4139    format (5g14.7)
         read (nfile,*) (covreng(jk),jk=1,icovreng1)
c         do jk1=1,icovreng1
c           covreng(jk1) = covreng(jk1)*1.e6
c         enddo
         if ( covreng(1) .gt. 1.e-4 .or. 
     1        covreng(icovreng1) .lt. 19.6e6) then 
             write (nt6, 7231) covreng(icovreng1), covreng(1)
             stop 'covfil-145'
         endif
c        cross section (b) - dummy entries
         do jk1=1,icovreng
           covrrsp1(jk1) = 1.0
         enddo
         if ( icov(1) .eq. 4) then 
            sumx = 0.0
            do jk = 1,icovreng
              sumx = sumx + covrrsp1(jk)
            enddo
            gsum = abs(sumx - 1.0)
            if ( gsum .gt. 1.1 .or. gsum .lt. 0.9 ) then
               write (nt6, 7243) sumx
            endif
            do jk = 1,icovreng
              covrrsp1(jk) = covrrsp1(jk)/sumx
            enddo
         endif
         do jk=1,icovreng
            covrrsp2(jk) = covrrsp1(jk)
         enddo
c        relative standard deviation
         read (nfile,*) (stdrdev1(jk), jk=1,icovreng)
 3379    format (5g14.7)
         do jk=1,icovreng
            stdrdev1(jk) = stdrdev2(jk)
         enddo
c
c        check for same function covariance
c
         iself = 0
         do jk=1,icovreng
           if ( covrrsp1(jk) .ne. covrrsp2(jk) ) iself = iself + 1
         enddo
         do i1=1,icovreng
            read (nfile,*) (cov_cov(i1,i2), i2=i1,icovreng)
c           note correlations belwo threshold are zero
c             if ( cov_cov(i1,i1) .eq. 0) cov_cov(i1,i1) = 1000
         enddo
         do jk1=1,icovreng
            do jk2=1,icovreng
               if ( jk1 .ge. jk2) then
                   corr(jk1,jk2) = 0.001*cov_cov(jk2,jk1)
               else
                   corr(jk1,jk2) = 0.001*cov_cov(jk1,jk2)
               endif
            enddo
         enddo
         close (unit=nfile)
      endif
c
c     Check for over-ride of covariance matrix
c
      if ( icov(3) .eq. 0) then 
            continue
      elseif (icov(3) .eq. 1) then 
c
c        override generated covariance with an uncorrolated matrix
c        on native energy structure
c 
            do jk1=1,icovreng
               dif1 = stddev1(jk1)*0.01*covrsp1(jk1)
               do jk2 = 1,icovreng
                   if ( jk1 .eq. jk2) then 
                       corr(jk1,jk1) = 1.0
                       dif2 = stdrdev2(jk2)*0.01*covrrsp2(jk2)
                       covr(jk1,jk1) = dif1*dif2
                       if ( covrrsp1(jk1) .eq. 0.0 .or. 
     1                      covrrsp2(jk2) .eq. 0.0 ) then
                                    corr(jk1,jk1) = 0.0
                                    covr(jk1,jk1) = 0.0
                       endif
                   else
                       covr(jk1,jk2) = 0.0
                       corr(jk1,jk2) = 0.0
                endif
               enddo
            enddo
       elseif (icov(3) .eq. 2) then
c
c           use default generation procedure, fix regridded cov/cor
c 
            continue
       elseif (icov(3) .eq. 3) then
c
c           use default generation procedure, fix regridded cov/cor
c 
            continue
       else 
            write (6,1457) icov(3)
 1457       format (1x, 'Illegal icov field ', i5)
            stop 'ICOV-3'
       endif
       if ( icon(9) <= 0) then 
            write (6,7869) icode, icoveng, icovreng, nenergy
 7869       format (1x, 'COVREAD before default fill ', 5i10)
       endif
      if ( icode .eq. 1 .or. icode .eq. 2 .or. 
     &     icode .eq. 3 .or. icode .eq. 4 .or.
     &     icode .eq. 5) then 
c
c         fill default energy grid
c
          if ( icon(9) .le. 0) then
             write (6,2829) icode
 2829        format (1x, 'covread icode: ', 3i5)
          endif
          if ( icode .eq. 1) then
             icoveng = 640
c            energy low to high in eV
             do jk2=1,641
                 coveng(jk2) = engsnd(641-jk2+1)*1.e6
             enddo
          elseif ( icode .eq. 2) then
             icoveng = 770
c            energy low to high in eV
             do jk2=1,771
                 coveng(jk2) = engsnd(771-jk2+1)*1.e6
             enddo
          elseif ( icode .eq. 3) then
             icoveng = 89
c            energy low to high in eV
             do jk2=1,90
                 coveng(jk2) = engsnd(90-jk2+1)*1.e6
             enddo
          elseif ( icode .eq. 4) then
             icoveng = 48
c            energy low to high in eV
             do jk2=1,49
                 coveng(jk2) = engsnd(49-jk2+1)*1.e6
             enddo
          elseif ( icode .eq. 5) then
             icoveng = 175
c            energy low to high in eV
             do jk2=1,176
                 coveng(jk2) = engsnd(176-jk2+1)*1.e6
             enddo
          endif
c
c         over-ride for 89-group case
c
c          if ( icovreng == 89) then 
c            icoveng = 89
c           energy low to high in eV
c            do jk2=1,90
c                coveng(jk2) = engsnd(90-jk2+1)*1.e6
c            enddo
c          endif
c
c         expand stddev data to 640/770 grid
c
          xarray(1) = 0.0
          yarray(1) = 0.0
          do jk1=1,icovreng
             xarray(jk1+1) = covreng(jk1)
             yarray(jk1+1) = stdrdev1(jk1)
          enddo
          xarray(icovreng+2) = covreng(icovreng+1)
          yarray(icovreng+2) = 0.0
c          write (6,5612) (xarray(jk), yarray(jk), jk=1,icovreng+2)
c 5612     format (1x, 'covread interp array: ', 8g14.7) 
          do jk=1,icoveng
c            interpolate y histogram with x
             stddev1(jk) = fitmd(coveng(jk), icovreng+2, 
     1       xarray, yarray, 1)
c             write (6,2671) jk, coveng(jk), stddev1(jk)
c 2671        format (1x, 'covread interp: ', i5, 2g14.7)
          enddo
          xarray(1) = 0.0
          yarray(1) = 0.0
          do jk1=1,icovreng
             xarray(jk1+1) = covreng(jk1)
             yarray(jk1+1) = stdrdev2(jk1)
          enddo
          xarray(icovreng+2) = covreng(icovreng+1)
          yarray(icovreng+2) = 0.0
          do jk=1,icoveng
c            interpolate y histogram with x
             stddev2(jk) = fitmd(coveng(jk), icovreng+2, 
     1       xarray, yarray, 1)
          enddo
c
c         expand covrsp data to 640/770 grid
c
          xarray(1) = 0.0
          yarray(1) = 0.0
          do jk1=1,icovreng
             xarray(jk1+1) = covreng(jk1)
             if ( icov(1) .eq. 4) then 
c               convert to differential number
                yarray(jk1+1) = covrrsp2(jk1)/(covreng(jk1+1) - 
     1          covreng(jk1))
             else
                yarray(jk1+1) = covrrsp2(jk1)
             endif
          enddo
          xarray(icovreng+2) = covreng(icovreng+1)
          yarray(icovreng+2) = 0.0
          do jk=1,icoveng
c            interpolate y histogram with x
             covrsp2(jk) = fitmd(coveng(jk), icovreng+2, xarray,
     1       yarray, 1)
             if ( icov(1) .eq. 4 .and. icovfmt .ne. 7) then 
c            convert back to number fraction
                covrsp2(jk) = covrsp2(jk)*(coveng(jk+1) - coveng(jk))
             else
                continue
             endif
          enddo
          if ( icov(1) .eq. 4) then 
c
c            spectrum not response interpolation scheme 
c
             write (nt6, 8924) icovreng, icoveng
 8924        format (1x, 'Permit spectrum interpolation scheme ', 2i5)
c
c
          endif
          if ( icovreng .ne. icoveng) then 
             xarray(1) = 0.0
             yarray(1) = 0.0
             do jk1=1,icovreng
                xarray(jk1+1) = covreng(jk1)
                if ( icov(1) .eq. 4) then 
c                  convert to differential number
                   yarray(jk1+1) = covrrsp1(jk1)/(covreng(jk1+1) - 
     1             covreng(jk1))
                else
                   yarray(jk1+1) = covrrsp1(jk1)
                endif
             enddo
             xarray(icovreng+2) = covreng(icovreng+1)
             yarray(icovreng+2) = 0.0
             do jk=1,icoveng
c               interpolate y histogram with x
                covrsp1(jk) = fitmd(coveng(jk), icovreng+2, xarray,
     1          yarray, 1)
                if ( icov(1) .eq. 4) then 
c                  convert back to number fraction
                   covrsp1(jk) = covrsp1(jk)*(coveng(jk+1) - coveng(jk))
                else
                   continue
                endif
             enddo
c            enforce spectrum normalization in interpolated form
             xsum = 0.0
             do jk=1,icoveng
                xsum = xsum + covrsp1(jk)
             enddo
             if ( icon(9) .le. 0) then 
               write (6,2831) xsum
 2831          format (1x, 'covread sum/normalization = ', g14.7)
             endif
             if ( xsum .eq. 0.0) xsum = 1.0
             if (icon(6) .eq. 1) then 
c
c             inhibit normalization - because this is a response
c
               write (6, 3489) xsum, icon(6)
 3489          format (1x, 'Renormalization in covread inhibited due to',
     &              ' icon(6) setting ', g14.7, 2x, i5)
               xsum = 1.0
c
             endif
             do jk=1,icoveng
                covrsp1(jk) = covrsp1(jk)/xsum
             enddo
          else
             do jk = 1, icoveng
                covrsp1(jk) = covrrsp1(jk)
             enddo
          endif
          if ( icon(9) .le. 0) then 
            write (6,5614) (covrsp1(jr), jr=1,icoveng)
 5614       format (1x, 'covread convert back: ', 5g14.7)
          endif
c
c         Check for replacement of response with 640/770 fine group response
c
          if ( icov(2) .eq. 1) then 
c
c            Read in new array name and replace default response
c 
             if ( iself .eq. 0) then 
                read (nt5,*) fmt3, 
     1             imode, iselect_material
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
     1           file_name(1:iflast)
                if ( icov(1) .eq. 4) then 
c                  spectrum over-ride NOT response here
                   new_file_name = xoptical(1:leopt)//kdir(1:kblank2)//
     1              file_name(1:iflast)
                endif
                if ( icov(1) .eq. 4) then 
                   call filein(imode, new_file_name, ierr, ireverse)
                else
                   call filein_response(imode, new_file_name, 
     &              ierr, ireverse)
                endif
c                write (nt6, 6712) ilead,iflast, new_file_name
c 6712           format (1x, 'covread filename length check: ', 2i8,
c     &             2x, a)
                if ( ilead .gt. 0) then
                   ovr = file_name(ilead:iflast)
                else
                   ovr = file_name(ilead+1:iflast)
                endif
                write (nt6,5723) file_name(ilead+1:iflast)
 5723           format (1x, 'Override function:  ', a75)
c  
c               save response data
c
                isum = 2
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
c               check for self-shielding factor on over-ride response
c
                 if ( icov(1) .ne. 4) then
c                    do not apply for a spectrum, only for a resposne 
                    if ( icon(13) .eq. 2) then 
                       read(nt5,*) file_name
                       file_corr = file_name
                       kstatus = 1
                       if ( icon(9) < 0) then 
                       write (6,8712) kstatus, icoveng, icovreng, 
     &                   new_file_name
 8712                    format (1x, 'COVREAD correct call', 4x, 3i5,
     &                       1x, a)
                       endif
                       call correct(array_hld(1,isum), file_name)
                    endif
                 endif
                 zap = abs(coveng(icoveng)*1.e-6 - 
     1                 emid_hld(nenergy,isum))/
     1                (coveng(icoveng)*1.e-6)
               if ( nenergy .ne. icoveng .or. zap .gt. 0.1 ) then 
                   write (nt6, 3612) nenergy, icoveng, icovreng, 
     &              emid_hld(1,isum),
     1              emid_hld(nenergy,isum), coveng(1)*1.e-6,
     1              coveng(icoveng)*1.e-6, zap
 3612              format (1x, 'Covariance override energy mis-match ', 
     1             3i8, 6g14.7)
               endif
               do jc1=1,icoveng
                   covrsp1(jc1) = array_hld(jc1,isum)
                   covrsp2(jc1) = array_hld(jc1,isum)
               enddo
c
c               now check for non-zero response with zero stddev
c
               ilooperr = 0
c                find first non-zero stddev
               stddef = 0.0
               do jc1=icoveng,1,-1
                   if ( stddef .eq. 0.0 .and. 
     1                  stddev1(jc1) .ne. 0.0) then
                      stddef = stddev1(jc1)
                   endif
               enddo
               do jc1=icoveng,1,-1
                  if ( covrsp1(jc1) .gt. 1.e-32 .and. 
     1              stddev1(jc1) .eq. 0) then 
                     ilooperr = ilooperr + 1
c                     write (6,4519) jc1, ilooperr, stddev1(jc1), 
c     &                   covrsp1(jc1)
c 4519                format (1x, 'covread stddev check: ', 2i6,2g14.7)
c                     find nearest non-zero std. dev.
                      stddef2 = 0.0
                      do jk=1,25
                         iup = jc1 + jk
                         if ( stddef2 .eq. 0 .and. 
     1                        iup .lt. icoveng) then 
                            if ( stddev1(iup) .ne. 0.0) 
     1                          stddef2 = stddev1(iup)
                         endif
                         idown = jc1 - jk
                         if ( stddef2 .eq. 0 .and. idown .ge. 1) then 
                            if ( stddev1(idown) .ne. 0.0) 
     1                          stddef2 = stddev1(idown)
                         endif
                      enddo
                       if ( stddef2 .eq. 0.0) stddef2 = stddef
                       if ( icon(9) .lt. 2) 
     1                    write (nt6, 8729) ilooperr, jc1, coveng(jc1), 
     1                    covrsp1(jc1),stddev1(jc1), stddef2
 8729                  format (1x, 'Zero stddev with 7 real response ',
     1                 2i5, 4g14.7)
                       stddev1(jc1) = stddef2
                       stddev2(jc1) = stddef2
                   endif
               enddo
               if ( ilooperr .gt. 1) then 
                  write (6,8731) ilooperr
 8731             format (1x, 'Zero stddev with 31 real response ',
     &            i5, 
     1           ' times')
               endif
             else
c               non-sysmmetric covariance - do the work twice
c               do covrsp1a
                read (nt5,*) fmt3, 
     1             imode, iselect_material
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
     1           file_name(1:iflast)
                if ( icon(9) < 0) then 
                   ilocate = 1
                   write (6,5621) ilocate, new_file_name
 5621              format (1x, 'covread before filein call ',i5,2x,a)
                endif
                call filein(imode, new_file_name, ierr, ireverse)
c  
c               save response data
c
                isum = 2
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
                if ( icon(13) .eq. 2) then 
                   write (nt6, 3761)
                   stop 'Self-Shielding-3'
                endif
                 zap = abs(coveng(icoveng)*1.e-6 - 
     1                     emid_hld(nenergy,isum))/
     1                (coveng(icoveng)*1.e-6)
                if ( nenergy .ne. icoveng .or. zap .gt. 0.01 ) then 
                   write (nt6, 3612) nenergy, icoveng, emid_hld(1,isum),
     1              emid_hld(nenergy,isum), coveng(1)*1.e-6,
     1              coveng(icoveng)*1.e-6, zap
                endif
                do jc1=1,icoveng
                   covrsp1(jc1) = array_hld(jc1,isum)
                enddo
c
c               now check for non-zero response with zero stddev
c
                ilooperr = 0
c                find first non-zero stddev
                stddef = 0.0
               do jc1=icoveng,1,-1
                   if ( stddef .eq. 0.0 .and. 
     1                 stddev1(jc1) .ne. 0.0) then
                      stddef = stddev1(jc1)
                   endif
                enddo
                do jc1=icoveng,1,-1
                   if ( covrsp1(jc1) .gt. 1.e-32 .and. 
     1               stddev1(jc1) .eq. 0) then 
                      ilooperr = ilooperr + 1
c                      find nearest non-zero std. dev.
                       stddef2 = 0.0
                       do jk=1,25
                          iup = jc1 + jk
                          if ( stddef2 .eq. 0 .and. 
     1                         iup .lt. icoveng) then 
                             if ( stddev1(iup) .ne. 0.0) 
     1                           stddef2 = stddev1(iup)
                          endif
                          idown = jc1 - jk
                          if ( stddef2 .eq. 0 .and. idown .ge. 1) then 
                             if ( stddev1(idown) .ne. 0.0) 
     1                           stddef2 = stddev1(idown)
                          endif
                       enddo
                       if ( stddef2 .eq. 0.0) stddef2 = stddef
                       if ( icon(9) .lt. 2) 
     1                    write (nt6, 8229) ilooperr, jc1, coveng(jc1), 
     1                    covrsp1(jc1),stddev1(jc1), stddef2
 8229                  format (1x, 'Zero stddev with 2 real response ',
     1                 2i5, 4g14.7)
                       stddev1(jc1) = stddef2
                   endif
                enddo
                if ( ilooperr .gt. 1) then 
                   write (6,8731) ilooperr
                endif
c               do response covrsp2
                read (nt5,*) fmt3, 
     1             imode, iselect_material
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
     1           file_name(1:iflast)
                if ( icon(9) < 0) then 
                   ilocate = 2
                   write (6,5623) ilocate, new_file_name
                endif
                call filein(imode, new_file_name, ierr, ireverse)
c  
c               save response data
c
                isum = 2
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
c          Not implemented for self-shielding correction
c
                if ( icon(13) .eq. 2) then 
                   write (nt6, 3761)
 3761              format (1x, 'Self-shielding not permitted',
     1            ' in this case !')
                   stop 'Self-Shielding-2'
                endif
                 zap = abs(coveng(icoveng)*1.e-6 - 
     1                 emid_hld(nenergy,isum))/
     1                (coveng(icoveng)*1.e-6)
                if ( nenergy .ne. icoveng .or. zap .gt. 0.01 ) then 
                   write (nt6, 3612) nenergy, icoveng, emid_hld(1,isum),
     1              emid_hld(nenergy,isum), coveng(1)*1.e-6,
     1              coveng(icoveng)*1.e-6, zap
                endif
                do jc1=1,icoveng
                   covrsp2(jc1) = array_hld(jc1,isum)
                enddo
c
c               now check for non-zero response with zero stddev
c
                ilooperr = 0
c                find first non-zero stddev
                stddef = 0.0
               do jc1=icoveng,1,-1
                   if ( stddef .eq. 0.0 .and. 
     1                  stddev2(jc1) .ne. 0.0) then
                      stddef = stddev2(jc1)
                   endif
                enddo
                do jc1=icoveng,1,-1
                   if ( covrsp1(jc1) .gt. 1.e-32 .and. 
     1               stddev1(jc1) .eq. 0) then 
                      ilooperr = ilooperr + 1
c                      find nearest non-zero std. dev.
                       stddef2 = 0.0
                       do jk=1,25
                          iup = jc1 + jk
                          if ( stddef2 .eq. 0 .and. 
     1                         iup .lt. icoveng) then 
                             if ( stddev2(iup) .ne. 0.0) 
     1                           stddef2 = stddev2(iup)
                          endif
                          idown = jc1 - jk
                          if ( stddef2 .eq. 0 .and. idown .ge. 1) then 
                             if ( stddev1(idown) .ne. 0.0) 
     1                           stddef2 = stddev2(idown)
                          endif
                       enddo
                       if ( stddef2 .eq. 0.0) stddef2 = stddef
                       if ( icon(9) .lt. 2) 
     1                    write (nt6, 8329) ilooperr, jc1, coveng(jc1), 
     1                    covrsp2(jc1),stddev2(jc1), stddef2
 8329                  format (1x, 'Zero stddev with 3 real response ',
     1                 2i5, 4g14.7)
                       stddev2(jc1) = stddef2
                   endif
                enddo
                if ( ilooperr .gt. 1) then 
                   write (6,8731) ilooperr
                endif
             endif
          elseif (icov(2) .eq. 0) then 
c
c             No change
c
              continue
          else
              write (nt6, 3712) icov(2)
 3712         format (1x, 'Illegal ICOV(2) value = ', i5)
              stop 'ICOV-2'
          endif
c         expand correlation data to 640/770 grid
c
c
c         construct covariance mapping
c
          do j1=1,icoveng
             itran(j1) = 0
c            eng = midpoint of new energy grid
             eng = 0.5*(coveng(j1)+coveng(j1+1))
             do j2=1,icovreng
c               find interval of original energy grid that corresponds to eng
                if ( eng .ge. covreng(j2) .and. 
     1               eng .lt. covreng(j2+1) ) then
c                    write (6,1076) j1, j2, eng, covreng(j2), 
c     1                  covreng(j2+1)
c 1076               format (1x, 'Correlation match ', 2i6, 3g14.7)
                    itran(j1) = j2
                endif  
             enddo
          enddo
          do j1=1,icoveng
             do j2=j1,icoveng
                if ( icon(9) .lt. 0) then
                  write (6,1078) j1, j2,itran(j1), itran(j2), icoveng, 
     &                  corr(itran(j1), itran(j2))
 1078             format (1x, 'error check: ', 6i8, 1x, g14.6)
                endif
                cor(j1,j2) = corr(itran(j1), itran(j2))
                cor(j2,j1) = cor(j1,j2)
             enddo
          enddo
          do jk1=1,icoveng
              dif1 = stddev1(jk1)*0.01*covrsp1(jk1)
              do jk2=1,icoveng
                 dif2 = stddev2(jk2)*0.01*covrsp2(jk2)
                 dif = dif1*dif2
                 cov(jk1,jk2) = cor(jk1,jk2)*dif
              enddo
          enddo
c         ensure correct diagonal values - for iself=0 case
          if ( iself .eq. 0) then 
             do jk1=1,icoveng
                do jk2=1,icoveng
                  if ( jk1 .eq. jk2 .and. cor(jk1,jk2) 
     1                     .ne. 0.0) then 
                      zap1 = abs(cor(jk1,jk2) - 1.0)
                      zap2 = abs(cor(jk1,jk2) - 100.)
                      if ( zap1 .gt. 1.e-3 .and. zap2 .gt. 1.e-1) then 
                         write (nt6,8931) jk1,jk2, cor(jk1,jk2), 
     1                   cov(jk1,jk2), dif1, dif2, zap
 8931                    format (1x, 'COR .NE. 1.0 check value ',
     1                   2i5, 5g14.7)
c                         stop 'COR ERROR 1'
                         cor(jk1,jk2) = 1.0
                      endif 
                 endif
                enddo
             enddo
c            fix zero response entries
             do jk1=1,icoveng
                cor(jk1,jk1) = 1.0
                if ( abs(covrsp1(jk1)) .le. 1.0e-32 .or. 
     1               abs(covrsp2(jk1)) .le. 1.0e-32 ) then
                      cor(jk1,jk1) = 0.0
                      cov(jk1,jk1) = 0.0
                endif
             enddo
          endif
c
c        Check correlation matrix for eigenvalues
c
c         call eigen_out(cor, icoveng)
c
c        Check for over-ride of covariance matrix
c
         if ( icov(3) .eq. 0) then 
            continue
         elseif (icov(3) .eq. 2) then 
c
c           override generated covariance with an uncorrolated matrix
c 
               do jk1=1,icoveng
                  dif1 = stddev1(jk1)*0.01*covrsp1(jk1)
                  do jk2 = 1,icoveng
                      if ( jk1 .eq. jk2) then 
                          cor(jk1,jk1) = 1.0
                          dif2 = stddev2(jk2)*0.01*covrsp2(jk2)
                          cov(jk1,jk1) = dif1*dif2
                          if ( covrsp1(jk1) .eq. 0.0 .or. 
     1                         covrsp2(jk2) .eq. 0.0 ) then
                                       cor(jk1,jk1) = 0.0
                                       cov(jk1,jk1) = 0.0
                          endif
                      else
                          cov(jk1,jk2) = 0.0
                          cor(jk1,jk2) = 0.0
                      endif
                  enddo
               enddo
         elseif (icov(3) .eq. 3) then 
c
c           override generated covariance with an totally corrolated matrix
c 
               do jk1=1,icoveng
                  dif1 = stddev1(jk1)*0.01*covrsp1(jk1)
                  do jk2 = 1,icoveng
                      if ( jk1 .eq. jk2) then 
                          cor(jk1,jk1) = 1.0
                          dif2 = stddev2(jk2)*0.01*covrsp2(jk2)
                          cov(jk1,jk1) = dif1*dif2
                          if ( covrsp1(jk1) .eq. 0.0 .or. 
     1                         covrsp2(jk2) .eq. 0.0 ) then
                                       cor(jk1,jk1) = 0.0
                                       cov(jk1,jk1) = 0.0
                          endif
                      else
                          cov(jk1,jk2) = dif1*dif2
                          cor(jk1,jk2) = 1.0
                      endif
                  enddo
               enddo
          elseif (icov(3) .eq. 1) then
c
c              use default generation procedure
c 
               continue
          else 
               write (6,1452) icov(3)
 1452          format (1x, 'Illegal icov field ', i5)
               stop 'ICOV-3'
          endif
c
c         If data is projected onto new energy grid, save the projection
c         unless scaling is to be done, icov(1)=3
c
         if ( icon(9) < 0) then 
            write (6,4590) icov(1)
 4590       format (1x, 'COVREAD projection write option ', i5)
         endif
         if ( icov(1) .ne. 3) then 
            nfile = 38
            ifmt = lnblnk(fmt)
            jblank2 = lnblnk(jdir)
c            name = optical(1:mblank2)//jdir(1:jblank2)//
c     1             'covar/'//job(1:jblank3)//'.snlcov'
c            name = optical(1:mblank2)//jdir(1:jblank2)//
c     1             'covar/'//fmt(1:ifmt)//'.snlcov'
            name = optical(1:mblank2)//
     1             'covar/'//fmt(1:ifmt)//'.snlcov'
            lend = lnblnk(name)
            if ( icon(9) < 0) then 
              write (6,4189) name(1:lend), optical(1:mblank2), 
     &          jdir(1:jblank2), fmt(1:ifmt)
 4189         format (1x, 'snlcov file write: name = ', a,/,
     &                1x, '                optical = ', a,/,
     &                1x, '                   jdir = ', a,/,
     &                1x, '                    fmt = ', a,/)
            endif
             open(unit=nfile, form='formatted',
     1            file=name(1:lend),
     2            status='unknown', iostat=ilook, err=1909)
c           write snlcov type file
            write (nfile,8011) icoveng
 8011       format (i10)
            write (nfile,8012) (coveng(jk),jk=1,icoveng+1)
 8012       format (1x, 8(g14.7,1x))
            write (nfile,8012) (covrsp1(jk),jk=1,icoveng)
            write (nfile,8012) (stddev1(jk),jk=1,icoveng)
            write (nfile,8012) (covrsp2(jk),jk=1,icoveng)
            write (nfile,8012) (stddev2(jk),jk=1,icoveng)
            write (nfile,8012) ((cov(jk1,jk2), jk1=1,icoveng), 
     1         jk2=1,icoveng)
            write (nfile,8012) ((cor(jk1,jk2), jk1=1,icoveng), 
     1         jk2=1,icoveng)
            close (unit=nfile)
         endif
c
      elseif ( icode .eq. 0) then 
          if ( icon(9) < 0) then 
             write (6,561) icode, icovreng
 561         format (1x, 'icode zero default fill ', 2i15)
          endif
c
c         leave default energy grid
c
          icoveng = icovreng
          coveng(icoveng+1) = covreng(icoveng+1)
          do jk1=1,icovreng
            covrsp1(jk1) = covrrsp1(jk1)
            stddev1(jk1) = stdrdev1(jk1)
            covrsp2(jk1) = covrrsp2(jk1)
            stddev2(jk1) = stdrdev2(jk1)
            coveng(jk1) = covreng(jk1)
            do jk2=1,icovreng
                if ( icon(9) .lt. 0) then 
                  write (6,7843) jk1, jk2, corr(jk1,jk2)
 7843             format (1x, 'DEBUG corr ', 2i5, g14.7)
                endif
                cor(jk1,jk2) = corr(jk1,jk2)
                cov(jk1,jk2) = covr(jk1,jk2)
            enddo
          enddo
      else
           write (nt6, 7881) icode
7881      format (1x, 'ICODE format not defined = ', i5)
 7801      format (1x, 'COVREAD format not defined = ', i5)
           stop 'COVREAD-1'
      endif
      if ( icon(9) .le. 0) then 
         write (6,6614) (covrsp1(jr), jr=1,icoveng)
 6614    format (1x, 'covread exit: ', 5g14.7)
      endif
      if (icon(9) <= 0) then 
          write (6,6711)
 6711     format (1x, 'EXIT COVREAD ')
      endif
      return
1909  continue
      length = len(fmt)
      write (nt6,910) ilook, name(1:lend)
910   format (1x, '*** covread error ', i5,
     1    ' opening file ', a)
      stop 'covar read error'
      end
