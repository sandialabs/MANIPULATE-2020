      subroutine cov_combine
      common /io/ nt5, nt6, nfile, nplot, npun
      common /datain/ nenergy, energy(1001), array(1000,15),
     1        emid(1001)
      common /pltlab/ lab_lead(15), lab_trail(15), lab_nchar(15),
     1 lab_len(15), icurve(15), lab_file(15)
      dimension cov_part(1001,1001), coveng_part(1001),
     1   covrsp1_part(1000), covrsp2_part(1000), cor_part(1001,1001), 
     2   stddev1_part(1000), stddev2_part(1000)
      dimension cov_sum(1001, 1001), cor_sum(1001, 1001), 
     1          rsp_sum(1000), std_sum(1000), unc_sum(1000), 
     2          eng(1000), cross_section(1000)
      integer icovfmt_part, icode_part, icoveng_part
      character*250 fmt_part
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
      real icov_cov
      dimension covr(1001,1001), covreng(1001), corr(1001,1001),
     1    covrrsp1(1000), covrrsp2(1000),
     1    stdrdev1(1000),
     1    stdrdev2(1000), icov_cov(1001,1001)
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
      equivalence (icov_1, lcov_1), (icov_2, lcov_2), 
     1     (icov_3, lcov_3),   (icov_4, lcov_4)
c      logical ok
c     character*4 inull
      character*4 lcov_1, lcov_2, lcov_3, lcov_4
      character*80 label
      logical pos_def
      real, allocatable    :: input_replace(:,:)
      if (icon(9) < 0) then 
          write (6,6712) 
 6712     format (1x, 'COV_COMBINE entered ')
      endif
      ename = 'opt'
      lename = lnblnk(ename)
      call getenv(ename(1:lename), optical)
      mblank2 = lnblnk(optical)
      ename = 'job'
      lename = lnblnk(ename)
      call getenv(ename(1:lename), job)
      jblank3= lnblnk(job)
      read (nt5, *) number_of_files, scale, outfile
      if ( number_of_files .gt. 5) then 
          write (6,9845) number_of_files
 9845     format (1x, 'number of files in cov_combine exceeds 5 ')
          stop 'number_of_files'
      endif
c
c     zero out summation files
c
      do jk1=1,1000
        std_sum(jk1) = 0.0
        unc_sum(jk1) = 0.0
        rsp_sum(jk1) = 0.0
        do jk2=1,1000
          cov_sum(jk1,jk2) = 0.0
          cor_sum(jk1,jk2) = 0.0
        enddo
      enddo
c
c     Form combination
c
      do ifile = 1,number_of_files
         read (nt5,*) icovfmt_part, file_name,icode_part,weight_part
         if ( icon(9) < 0) then 
            write (6, 1691) icovfmt_part, icode_part, weight_part, 
     1         file_name
1691        format (1x, 'MANIPULATE call to cov_combine file ',
     &         'function ', 2i5, g14.6, 1x, a)
         endif
c
         fmt_part = file_name
c        only SNLRML *.lsl format is accepted here
         if ( icovfmt_part .ne. 2 .and. icovfmt_part .ne. 5) then 
           write (6,780) ifile, icovfmt_part
 780       format (1x, 'Illegal icovfmt_part in cov_combine ', 2i5)
           stop 'cov-comb-1'
         endif
c        only default energy grid is accepted here
         if ( icode_part .ne. 0) then 
           write (6,781) ifile, icode_part
 781       format (1x, 'Illegal icode_part in cov_combine ', 2i5)
           icode_part = 0
         endif
c
         call covread(icovfmt_part, fmt_part, icoveng_part, 
     1           coveng_part, covrsp1_part, stddev1_part, 
     2           covrsp2_part, stddev2_part, cov_part, cor_part, 
     3           icode_part)
c
c        Check input quantities
c
         if ( icon(9) < 0) then 
            write (6,75) (coveng_part(jk), jk=1,icoveng_part)
 75         format (1x, 'energy ', 5g14.7)
            write (6,73) (covrsp1_part(jk), jk=1,icoveng_part)
 73         format (1x, 'xsec   ', 5g14.7)
            write (6,91) (stddev1_part(jk)*0.01*covrsp1_part(jk), 
     1                  jk=1,icoveng_part)
 91         format (1x, 'unc_pt ', 5g14.7)
            write (6,71) (stddev1_part(jk), jk=1,icoveng_part)
 71         format (1x, 'stddev ', 5g14.7)
            do jk1=1,icoveng_part
             write (6,76) (cor_part(jk1,jk2), jk2=1,icoveng_part)
 76          format (1x, 'cor:   ', 5g14.7) 
            enddo
         endif
c
c           Form covariance matrix
c
         do jk1=1,icoveng_part
            do jk2=1,icoveng_part
              dummx = covrsp1_part(jk1)*covrsp1_part(jk2)*
     1               stddev1_part(jk1)*stddev1_part(jk2)*0.01*0.01
              cov_part(jk1,jk2) = cor_part(jk1,jk2)*dummx
              if ( icon(9) < -3) then 
                  write (6,81) jk1, jk2, dummx, cor_part(jk1,jk2),
     1                cov_part(jk1,jk2)
 81               format (1x, 'cov construct: ', 2i5, 3g14.7)
              endif
            enddo
         enddo
c
         if (icon(9) < 0) then
           do jk1=1,icoveng_part
             write (6,77) (cov_part(jk1,jk2), jk2=1,icoveng_part)
 77          format (1x, 'cov:   ', 5g14.7) 
           enddo
         endif
c
c        Now perform weighted combination if covariance matrices
c                - result stored in cov_sum
c                - we add the weighted uncertainties, so we convert
c                  the covariance back to the unc by taking a sqrt
c                  and sum.  Then we take square to get a covariance back.
c                  Assumption is that the components are uncorrelated.
c
         do jk1=1,icoveng_part
c
c          Force common response/corss-section, so do not weight and sum
c
c           rsp_sum(jk1) = rsp_sum(jk1) + weight_part*covrsp1_part(jk1)
           rsp_sum(jk1) = covrsp1_part(jk1)
c
c          No, the std-dev is NOT the weighted sum of the component std_dev
c          - it depends upon the correlation
c
c           unc_sum(jk1) = unc_sum(jk1) + 
c     1          weight_part*0.01*stddev1_part(jk1)*covrsp1_part(jk1)
           do jk2=1,icoveng_part
              cov_sum(jk1,jk2) = cov_sum(jk1,jk2) + 
     1                weight_part*weight_part*cov_part(jk1,jk2)
           enddo
c
c          set std-dev based on covariance
c
           if ( cov_sum(jk1,jk1) .gt. 0.0) then
              unc_sum(jk1) = sqrt(cov_sum(jk1,jk1))
           else
              unc_sum(jk1) = 0.0
           endif
c
         enddo
c
c        make sure the energy grids are the same for the components
c        also make sure that the baseline cross sections are the same
c
         if ( ifile == 1) then
            npoints = icoveng_part 
            do jk1=1,icoveng_part
              eng(jk1) = coveng_part(jk1)
              cross_section(jk1) = covrsp1_part(jk1)
            enddo
            eng(icoveng_part+1) = coveng_part(icoveng_part+1)
         else
           if ( icoveng_part /= npoints) then 
             write (6,249) ifile, npoints, icoveng_part
 249         format (1x, 'ERROR in common energy points ', 3i8)
             stop 'cov_cmb_pts'
           endif
           do jk1=1, icoveng_part
              if ( coveng_part(jk1) .ne. eng(jk1)) then
                 write (6,250) ifile, jk1, coveng_part(jk1), eng(jk1)
 250             format (1x, 'ERROR in common energy grid ', 
     1             2i8, 2g14.7)
                 stop 'cov_cmb_eng' 
              endif
           enddo 
           do jk1=1, icoveng_part
              if ( covrsp1_part(jk1) .ne. cross_section(jk1)) then
                 write (6,251) ifile, jk1, covrsp1_part(jk1), 
     &                  cross_section(jk1)
 251             format (1x, 'ERROR in common cross section ', 
     1             2i8, 2g14.7)
                 stop 'cov_cmb_rsp' 
              endif
           enddo 
         endif
c
      enddo
      do jk1=1,icoveng_part
        if ( rsp_sum(jk1) .ne. 0.0) then
           std_sum(jk1) = unc_sum(jk1)/rsp_sum(jk1)*100.
        else
           std_sum(jk1) = 0.0
        endif
      enddo
c      do jk1=1,icoveng_part
c        do jk2=1,icoveng_part
c           cov_sum(jk1,jk2) = cov_sum(jk1,jk2)*cov_sum(jk1,jk2)
c        enddo
c      enddo
      do jk1=1,icoveng_part
        do jk2=1,icoveng_part
c          dummx = rsp_sum(jk1)*rsp_sum(jk2)*
c     1               std_sum(jk1)*std_sum(jk2)*0.01*0.01
          dummx = unc_sum(jk1)*unc_sum(jk2)
          if ( dummx == 0.0) dummx = 1.0
          cor_sum(jk1,jk2) = cov_sum(jk1,jk2)/dummx
        enddo
      enddo
      if ( icon(9) < 0) then 
        write (6,871) number_of_files
 871    format (1x, 'cov_combination has started: ', i5)
        write (6,173) (rsp_sum(jk), jk=1,icoveng_part)
173     format (1x, 'sum:    ', 5g14.7)
        write (6,273) (unc_sum(jk), jk=1,icoveng_part)
273     format (1x, 'uncsum  ', 5g14.7)
        write (6,171) (std_sum(jk), jk=1,icoveng_part)
171     format (1x, 'stdsum: ', 5g14.7)
        do jk1=1,icoveng_part
         write (6,176) (cor_sum(jk1,jk2), jk2=1,icoveng_part)
176      format (1x,'cor_sum:', 5g14.7) 
        enddo
        do jk1=1,icoveng_part
         write (6,177) (cov_sum(jk1,jk2), jk2=1,icoveng_part)
177      format (1x,'cov_sum:', 5g14.7) 
        enddo
      endif
c
c     apply overall weighting factor
c
c      do jk1=1,npoints
c
c      do not scale response/corss-section
c
c        rsp_sum(jk1) = rsp_sum(jk1)*scale
c        do jk2=1,npoints
c          cov_sum(jk1,jk2) = cov_sum(jk1,jk2)*scale
c        enddo
c      enddo
c
c     Check correlation matrix for eigenvalues
c
      write (6,6734)
 6734 format (/,/,1x, 
     1   'Eigenvalues for combined covariance matrix: ')
      allocate( input_replace(1001, 1001),
     $      STAT=iStat)
      if (iStat/=0) then 
        write (6,25001) iStat
        stop 'Memory-allocation error in covread.'
      endif
25001 format ('Subroutine covread memory-(de)allocation error: ',
     &       2i5)
      ipass = 0
      input_replace = 0.0
c
c    Note: eigenvalue non-negativity check inhibited for combined covariacne option 
c
      if ( icon(1) .eq. 14) then
        go to 9871
      endif
c
      call eigen_out(cor_sum, npoints, input_replace, pos_def, ipass, 
     &    eigen_save)
c
c     Check for over-ride of covariance matrix
c
      deallocate(input_replace, STAT=iStat)
      if (iStat/=0) then 
        write (6,25000) iStat
25000   format ('Subroutine covread memory-(de)allocation error: ',
     &          2i5)
        stop 'Memory-(de)allocation error in eigen_out.'
      endif
 9871 continue
c
c     Write lsl-format combined file
c
      write (6, 7832) 
 7832 format (/,/,1x, 
     1 'Output combined covariance data in LSL - ',
     1  'format ',/,/)
c                                                                    
c Title cards                                                       
c                                                                  
      open (unit=78, file=outfile,
     &      status='unknown')
      write(78,8562)                                                
 8562 format(    '*COR    (LIBRARY)    (MAT.#)    (TEMP)K')         
c                                                                    
c Energy grid                                                         
c                                                                   
      write(78,8996)                                                 
 8996 format('*Number of Energies plus 1')                      
      write(78,8995) npoints+1                                          
      write(78,8993)                                                 
 8993 format('*Energy Grid ( eV )')          
 8995 format(i5)                             
      write(78,8990) (eng(jk),jk=1,npoints+1)  
c                                                  
c Cross section                                                       
c                                                                    
      write(78,8998)                                                 
 8998 format('*Cross Section (barn)')                           
      write(78,7454) (rsp_sum(jk), jk=1,npoints)                                                     
c                                                                   
c Standard deviation                                               
c                                                                   
      write(78,8991)                                            
 8991 format('*% Standard Deviation')                            
      write(78,8990) (std_sum(i1),i1=1,npoints)                           
c                                                                  
c Correlation coefficients                                        
c                                                                  
      write(78,8992)                                              
 8992 format('*Correlation Coefficient -- Upper Triangular')       
      do i1 = 1,npoints       
         write(78,8990) ((100.0*cor_sum(i1,i2)),i2=i1,npoints)              
      end do                                                       
 8990 format((1x,5(1pe14.7,1x)))                                    
 7454 format((1x,1p5e14.7))                                         
 7990 format((1x,15(i4,1x)))                                   
      close (unit=78)
c
c   plot std. dev. in plot format
c
      jblank4 = lnblnk(outfile)
      open (unit=78, file=outfile(1:jblank4)//'.std_pct_plt',
     &      status='unknown')
      write (78, 1205) eng(1)*1.e-6, 
     &        std_sum(1)*1.e-3
1205  format (2x, g14.7, 2x, g14.7)
      do jk=1, npoints
        write (78,1205) eng(jk)*1.e-6,std_sum(jk)
        write (78,1205) eng(jk+1)*1.e-6,std_sum(jk)
      enddo
      write (78, 1205) eng(npoints+1)*1.e-6, 
     &                 std_sum(npoints)*1.e-3
      close (unit=78)
c
c   plot correlation matrix in plot format
c
      open (unit=78, file=outfile(1:jblank4)//'.corplt',
     &      status='unknown')
      do jk1=1, npoints
        write (78,2295) eng(jk1)*1.e-6,eng(jk1)*1.e-6,
     &           (cor_sum(jk1,jk2), 
     &           jk2 = 1, npoints), cor_sum(jk1,npoints)
 2295   format (1x, g12.4, 2x, 191(',', g12.4, 2x))
        if ( icon(9) .lt. 0) then 
          if ( jk1 .lt. 5) then
             write (6,3491) (cor_sum(jk1,jk2), jk2 = 1, 5)
 3491        format (1x, 'corplt look: ', 5g12.4)
          endif
        endif
      enddo
      write (78,2295) eng(npoints+1)*1.e-6, 
     &     eng(npoints+1)*1.e-6,
     &     (cor_sum(npoints,jk2), jk2 = 1, npoints), 
     &      cor_sum(npoints,npoints)
      close (unit=78)
c
c    end of LSL rwrite logic                                
c
      return
      end
