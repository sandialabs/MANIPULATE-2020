      subroutine cov_norm(icovreng, covrrsp1, stdrdev1, input_replace,
     &                ipass)
      implicit none
      integer :: icovreng, ipass 
      integer :: iuc, jk1, jk2, il1, il2, iflag, jp, inc
      real :: flunorm, offset1, offset2, sum, error_max, dummx
      real, dimension(1001, 1001) :: input_replace 
      real, dimension(1000, 1000) :: abs_cov, abs2_cov
      real, dimension(1000) :: covrrsp1, stdrdev1, row_norm
c      real icon
      integer jcon
      common /guide/ jcon(40)
c
c     Enforce a fluence normalization on a spectrum covariance matrix 
c
c     form absolute covariance values
c
      if ( jcon(9) <= 0) then 
         write (6,6723) ipass, jcon(9)
 6723    format (1x, 'cov_norm entry: ', 2i5)
      endif
      inc = 0
      do jk1=1,icovreng
            do jk2=jk1, icovreng
               inc = inc + 1
c
c              Convert relative covariance to absolute correlation
c
               dummx = covrrsp1(jk1)*covrrsp1(jk2)*stdrdev1(jk1)*
     1                 stdrdev1(jk2)*0.01*0.01
               abs_cov(jk1,jk2) = input_replace(jk1,jk2)*dummx
               if ( jcon(9) <= -1) then 
                 write (6,445) jk1, jk2, input_replace(jk1,jk2), 
     &                 dummx,
     &                 abs_cov(jk1,jk2),  
     &                 covrrsp1(jk1), covrrsp1(jk2), stdrdev1(jk1), 
     &                 stdrdev1(jk2)
 445           format (1x, 'input_replace debug: ', 2i4, 7g14.7)
               endif
               abs_cov(jk2,jk1) = abs_cov(jk1,jk2)
               if ( abs(input_replace(jk1,jk2)) .gt. 1.0001) then 
                  write (6,7824) jk1, jk2, input_replace(jk1, jk2), 
     &                covrrsp1(jk1), covrrsp1(jk2), stdrdev1(jk1), 
     &                stdrdev1(jk2), abs_cov(jk1,jk2), dummx
 7824             format (1x, 'Warning: illegal init rel cov: ',
     &                    2i4, 8g14.7)
               endif
            enddo
      enddo
c
c     Check covariance normalizaton for the original absolute 
c        covarianced matrix
c
      iflag = 0
      error_max = 0.0
      do jk1 = 1,icovreng
            row_norm(jk1) = 0.0
            do jk2 = 1,icovreng
              row_norm(jk1) = row_norm(jk1) + abs_cov(jk1,jk2)
            enddo
            error_max = max(error_max, row_norm(jk1))
            if ( abs(row_norm(jk1)) .gt. 1.e-5) then 
                iflag = iflag + 1
                if ( jcon(9) < = -1) then
                   write (6,6827) jk1, row_norm(jk1)
 6827              format (1x, 'Initial renormalized row norm error ', 
     1             i5, 2x, g14.7)
                endif
            endif
      enddo   
      if ( iflag > 0) write (6,7390) ipass, iflag, error_max
 7390 format (1x, 'Initial absolute covariance renormalization ', 
     &         i5, ' had ', i4, ' normalization problems ',/, 1x,
     &         'Largest row normalization = ', g14.7)   
c
c     In case of a normalization error on covariance constraint,
c         e.g. sum over any row or column is not zero,
c         treat these as the covariance for an unnormalized
c         spectrum and apply the normalization condition
c         using derived formula from Smith, pg. 140
c
c     Only apply the renorm if there is a large error_max
      if ( error_max > 1.e-5) then
      flunorm = 1.0
        do jk1 = 1, icovreng
           do jk2 = 1, icovreng
              sum = 0.0
              do il1 = 1, icovreng
                 do il2 = 1, icovreng
                     offset1 = 0.0
                     offset2 = 0.0
                     if ( jk1 .eq. il1) offset1 = flunorm
                     if ( jk2 .eq. il2) offset2 = flunorm
                     sum = sum + (offset1-covrrsp1(jk1))
     &                          *abs_cov(il1,il2)
     &                          *(offset2-covrrsp1(jk2))
                 enddo
              enddo
              abs2_cov(jk1, jk2) = sum/flunorm**4
               if ( jcon(9) <= -1) then 
                 write (6,645) jk1, jk2, abs2_cov(jk1,jk2)
 645             format (1x, 'abs2_cov debug: ', 2i4, 6g14.7)
               endif
           enddo
        enddo 
      else
        abs2_cov = abs_cov
      endif
c
c      When the covariance matrix was modified, the standard deviations
c      may have been slightly modified.  Thus, update this data to 
c      be consistent.  Note this is the relative std. dev.
c
      do jp = 1,icovreng
            if ( jcon(9) <= -1) then 
              write (6,3982) jp,  abs_cov(jp,jp), 
     &          covrrsp1(jp)
 3982         format (1x, 'outflu sqrt: ', i4, 6g14.7)
           endif
           if ( covrrsp1(jp) /= 0.0) then 
              stdrdev1(jp) = sqrt(abs2_cov(jp,jp))*100./covrrsp1(jp)
           else
              stdrdev1(jp) = 100.
           endif
           if ( jcon(9) <= -1) then 
              write (6,2982) jp, stdrdev1(jp), abs2_cov(jp,jp), 
     &          covrrsp1(jp)
 2982         format (1x, 'renormed stddev: ', i4, 6g14.7)
           endif
      enddo
c
c     
c     check normalization condition for modified covariance values
c     e.g. sum of elements in any row/column is zero 
c     ENDF requires this to a precision of 1.E-5
c
      iflag = 0
      error_max = 0.0
      do jk1 = 1,icovreng
            row_norm(jk1) = 0.0
            do jk2 = 1,icovreng
              row_norm(jk1) = row_norm(jk1) + abs2_cov(jk1,jk2)
            enddo
            error_max = max(error_max, row_norm(jk1))
            if ( abs(row_norm(jk1)) .gt. 1.e-5) then 
                iflag = iflag + 1
                write (6,16827) jk1, row_norm(jk1)
16827           format (1x, 'renormalized row norm error ', 
     1         i5, 2x, g14.7)
            endif
      enddo   
      if ( iflag > 0) write (6,7290) iflag, error_max
 7290 format (1x, 'Absolute covariance renormalization had ',
     &         i4, ' normalization problems ',/, 1x,
     &         'Largest row normalization = ', g14.7)   
c
c    reform relative covariance elements with new normalization constraint
c
c
c     derive new relative covariance matrix 
c     (rel_corr) and output renormalized covariance 
c     files in LSL format
c
      inc = 0
      do jk1=1,icovreng
            do jk2=jk1, icovreng
               inc = inc + 1
c
c              Convert absolute covariance to relative correlation
c
               dummx = covrrsp1(jk1)*covrrsp1(jk2)*stdrdev1(jk1)*
     1                 stdrdev1(jk2)*0.01*0.01
               if ( dummx /= 0.0) then 
                  input_replace(jk1,jk2) = abs2_cov(jk1,jk2)/dummx
               else
                  if ( jk1 == jk2) then 
                      input_replace(jk1,jk2) = 1.0
                  else
                      input_replace(jk1,jk2) = 0.0
                  endif
               endif
               if ( abs(input_replace(jk1,jk2)) .gt. 1.0001) then 
                  write (6,7823) jk1, jk2, input_replace(jk1, jk2), 
     &                covrrsp1(jk1), covrrsp1(jk2), stdrdev1(jk1), 
     &                stdrdev1(jk2), abs2_cov(jk1,jk2), dummx
 7823             format (1x, 'Warning: illegal rel cov: ',
     &                    2i4, 8g14.7)
               endif
               if ( jcon(9) <= -1) then 
                  write (6,1445) jk1, jk2, input_replace(jk1,jk2), 
     &                 dummx,
     &                 abs_cov(jk1,jk2),  
     &                 covrrsp1(jk1), covrrsp1(jk2), stdrdev1(jk1), 
     &                 stdrdev1(jk2)
1445              format (1x, 'rel_cov cov_norm debug: ', 2i4, 7g14.7)
                endif
               input_replace(jk2,jk1) = input_replace(jk1,jk2)
            enddo
      enddo
      if ( jcon(9) <= 0) then 
         write (6,6724) ipass,iflag,jcon(9)
 6724    format (1x, 'cov_norm exit: ', 3i5)
      endif
c
      return
      end
