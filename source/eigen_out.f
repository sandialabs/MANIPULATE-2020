C********1*********2*********3*********4*********5*********6*********7*C
C<html>
C<body>
C<h3>                                                                      
C Module Name:   eigen_out
C</h3>
C<pre>
C Module Type:   subroutine
C
C Revision Date:  $Date: 2010/03/22 23:13:41 $
C CVS Version:    $Revision: 1.7 $
C
C Major Modification History:
C   03/17/2009   JEC   Added SLATEC support for solution of
C                      eigenvalues and eigenvectors.
C   07/28/2006   PJG   Original documented version
C                       
C Purpose:  This subroutine determines the eigenvalues and 
C           eigenvectors for the relative correlation matrix.
C           It is built upon the SLATEC library since SNL has
c           distribution rights for this mathematical library.
C
C Description:  The section contains documentation on methods,
C     references, etc...........................................end
C     .............................................................
C
C Invocation:  CALL eigen_out(cor,icovreng,input_replace,pos_def,ipass, 
C     &             eigen_save)
C               
C Arguments:   
C
C        icovreng     (IN) dimension of useful subsection of array
c        cor          (IN) matrix
C                         
C Subroutines Called:
C    TYPE            ENTRY                 
C    SUBROUTINE      TBD
C
C Common block variables set or changed:
C      Module "XXX.cmn"
C          variable = description
C      Module "YYY.cmn"
C          variable = description
C                    
C                                                                      
C</pre>
C</body>
C</html>
C*****7**1*********2*********3*********4*********5*********6*********7*C   
C=======================================================================
C
      subroutine eigen_out(cor,icovreng,input_replace,pos_def,ipass, 
     &             eigen_save)
C.. Implicits .. 
C 
C.. Implicit Declarations .. 
      implicit none
C 
C.. Include Statements ..
c      include "genpar.par"
c      include "mpif.h"
C.. Parameters
      integer, parameter :: corSize = 1001
C.. Formal Arguments .. 
      real, intent(inout), dimension(corSize, corSize) :: cor, 
     &                     input_replace
      real, intent(inout), dimension(1001) :: eigen_save
      integer, intent(in) :: icovreng
      integer, intent(inout) :: ipass
      logical, intent(inout) :: pos_def
      integer :: icon
      common /guide/ icon(40)
C
C  
C.. Local Arrays
c
      real, allocatable    :: TempVector(:), corWorkingArray(:,:), 
     &          evect(:,:), diag_new(:,:), diag(:,:), evect_t(:,:), 
     &          input_dup(:,:),  temp1(:,:), 
     &          temp2(:,:)
      complex, allocatable :: Eigenvalue(:), VectorOut(:,:)
      real :: max_eigen
      integer :: istage, j1, j2
c
C 
C.. Include Statements ..
C 
C.. External Calls .. 
C 
c      include "block.cmn"
c      include "uqsrc.cmn"
c      include "mpicom.cmn"
C 
C.. Local Scalars .. 
C
      character(LEN=80) ::  alf
      real sumx, correction
      integer :: jk1, jk2, jk, ndim, iStat, iError, ineg
      logical  AllZero
      integer :: mynum, i, iflag
      character*80 :: mpi_message
c
      data mynum /1/
C  
c    Allocate internal arrays.
c
c      write (6,45)
c 45   format (1x, 'Bypass eigen_out for now')
c      return
      if ( icon(9) .lt. 0) then 
        write (6,5619) ipass
 5619   format (1x, 'eigen_out entry ', 2i5)
      endif
      iflag = 0
      pos_def = .true.
      allocate(VectorOut(icovreng,icovreng),
     $         corWorkingArray(icovreng,icovreng),
     $         Eigenvalue(icovreng), 
     $         TempVector(2*icovreng), STAT=iStat)
      if (iStat/=0) then 
        write (6,25000) iStat
        stop 'Memory-(de)allocation error in eigen_out.'
      endif
      input_replace = 0.0
c
      if (icovreng > corSize) then
        write(6,104) corSize
 104    format('Covariance dimension exceeds maximum value of ',i4)
        stop 'Covariance dimension exceeds maximum value.'
      endif
      
      if ( ipass == 1) write (6,101) icovreng
 101  format (/,/,1x, 'Relative correlation matrix (cor)',
     &    ' has dimension = ', i5,/,/)
c     
      do i = 1,icovreng
        if (cor(i,i)==0.) then
          cor(i,i) = 1.
          if (i>1) write(6,105) i
        endif      
      enddo
  105 format ("Warning: Element ",i5," along diagonal of ",
     &  "relative correlation matrix has been changed ",
     &  "from 0.0 to 1.0.")
     
      corWorkingArray = cor(1:icovreng,1:icovreng)

      istage = 3
      iError = 0
      write (49,371) ipass,istage
 371  format (1x, 'eigen_out entrance with ipass/istage = ', 2i5)
      write (49, 372) size(corWorkingArray,1),size(corWorkingArray,1),
     &      size(VectorOut,1), iError, icovreng
 372  format (1x, 'SGEEV input ', 5i10)
      do j1 = 1, icovreng
         write (49, 373) j1, (corWorkingArray(j1,j2), 
     &         j2=1,icovreng)
 373     format (1x,'corWorking Array: ', i5, 1x, 5g14.6,
     &          (1x,'                  ', 5x, 1x, 5g14.6))
      enddo
      do j1 = 1,icovreng
         write (49, 473) j1, (cor(j1,j2), 
     &         j2=1,icovreng)
 473     format (1x,'cor Array:        ', i5, 1x, 5g14.6,
     &          (1x,'                  ', 5x, 1x, 5g14.6))
      enddo

      call SGEEV (corWorkingArray, size(corWorkingArray,1),
     $    size(corWorkingArray,1), Eigenvalue, VectorOut,
     $    size(VectorOut,1), TempVector, 1, iError)
      
      if (iError/=0) then
        write (6,*) "Eigenvalues for covariance matrix do not converge."
c        stop 'Eigenvalues for covariance matrix do not converge.'
      endif
c
      if ( ipass < 10 .or. ipass == 750) then
         if ( icon(9) < 2) then 
             write (6,102) ( real(Eigenvalue(jk)), jk=1,
     $               ubound(Eigenvalue,DIM=1) )
         endif
 102     format (/,1x, "Eigenvalues: ", 5g14.7,/,
     &       (  1x, "             ", 5g14.7)   )
c         do jk2=1,5 !ubound(VectorOut,DIM=2)
         do jk2=1,ubound(VectorOut,DIM=2)
            if ( icon(9) < 2) then 
              write (6,103) jk2, (real(VectorOut(jk1, jk2)), jk1 = 1,
     $        ubound(VectorOut,DIM=1))
            endif
 103       format (/,  1x, "Eigenvector #",i4, 2x, 5g14.7,/,
     &         (  1x, "             ",2x, 2x, 5g14.7)    )
         enddo
      endif
c
c    Find maximum eigenvalue
c
      max_eigen = 0.0
      do jk=1,ubound(Eigenvalue,DIM=1)
        max_eigen = max(max_eigen, real(Eigenvalue(jk)) )
      enddo
c
c     Check for positive definite matrix
c
      pos_def = .true.
      ineg = 0
      do jk=1,ubound(Eigenvalue,DIM=1)
c        if (real(Eigenvalue(jk)) <= 0.0 .or. 
        if (real(Eigenvalue(jk)) <= -1.e-6*max_eigen .or. 
     &      abs(imag(Eigenvalue(jk))) > 1.e-6*max_eigen) then
            pos_def = .false.
            ineg = ineg + 1
            if ( ipass < 10 .or. ipass == 750) then
              if ( icon(9) < 2) then 
                 write (6,6723) jk, real(Eigenvalue(jk)),
     &                imag(Eigenvalue(jk))
              endif
 6723         format (1x, 'Erroneous eignevalue: ', i5, 2g14.7)
            endif
        endif
        eigen_save(jk) = real(Eigenvalue(jk))
      enddo
c      
      if (.not. pos_def .and. icon(1) .ne. 14) then
        if ( ipass < 10 .or. ipass == 750) write (6,190) ineg, ipass
 190    format (1x, 'Covariance matrix FAILS positive definite test.'
     &     , 2i5)
c        stop 'Covariance matrix FAILS positive definite test.'
c
         if ( iflag <= -1) then
           do jk2=1, ubound(VectorOut,DIM=2)
             write (6,9103) jk2, (real(VectorOut(jk1, jk2)), jk1 = 1,
     $       ubound(VectorOut,DIM=1))
9103         format (/,  1x, "Eigenvector_a #",i4, 2x, 5g14.7,/,
     &         (  1x, "             ",2x, 2x, 5g14.7)    )
           enddo
         endif
c
         allocate(evect(icovreng,icovreng),
     &            diag_new(icovreng,icovreng), 
     &            diag(icovreng,icovreng), 
     &            evect_t(icovreng,icovreng), 
     &            input_dup(icovreng,icovreng),  
     &            temp1(icovreng, icovreng), 
     &            temp2(icovreng, icovreng),
     &            STAT=iStat)
         if (iStat/=0) then 
           write (6,25000) iStat
           stop 'Memory-(de)allocation error in diag eigen_out.'
         endif
c        Form matrix of eigenvector - evect
c           matrix stored row x column
c           place eigenvectors in columns
         do jk=1,ubound(Eigenvalue,DIM=1)
           evect(jk,:) = real(VectorOut(jk, :))
         enddo
c        Form transpose of eigenvector matrix = evect_t
         do jk=1,ubound(Eigenvalue,DIM=1)
           evect_t(:,jk) = real(VectorOut(jk, :))
         enddo
c        Form matrix of eigenvalues = diag
         diag = 0
         do jk=1,ubound(Eigenvalue,DIM=1)
            diag(jk,jk) = real(Eigenvalue(jk))
         enddo
c        Form identical input matrix = evect*diag*evect_t
         do jk1= 1, ubound(Eigenvalue,DIM=1) 
            do jk2 = 1, ubound(Eigenvalue,DIM=1)
               temp1(jk1, jk2) = SUM( diag(jk1,:)*evect_t(:,jk2)  )
            enddo
         enddo
         do jk1= 1, ubound(Eigenvalue,DIM=1) 
            do jk2 = 1, ubound(Eigenvalue,DIM=1)
               input_dup(jk1,jk2)=SUM(evect(jk1,:)*temp1(:,jk2))
            enddo
         enddo
         if( iflag <= -2) then 
           do jk1 = 1, ubound(Eigenvalue,DIM=1)
              write (6,982) jk1, (input_dup(jk1, jk2), jk2=1,
     &                      ubound(Eigenvalue,DIM=1))
 982          format (1x, 'Input dup: ', i5, /, 
     &            (1x, '               ', 5x, 5g14.7) )
           enddo
         endif
c
c        Form replacement eigenvector matrix = diag_new
         diag_new = 0
         do jk=1,ubound(Eigenvalue,DIM=1)
            diag_new(jk,jk) = max(1.e-35, real(Eigenvalue(jk)) )
         enddo
c        Form replacement input matrix = evect*diag_new*evect_t
         do jk1= 1, ubound(Eigenvalue,DIM=1) 
            do jk2 = 1, ubound(Eigenvalue,DIM=1)
               temp2(jk1, jk2) = SUM( diag_new(jk1,:)*evect_t(:,jk2))
            enddo
         enddo
         if ( icon(9) .lt. 0) then 
           write (6,2981) 
 2981      format (1x, 'eigen_out input_replace set ')
         endif
         do jk1= 1, ubound(Eigenvalue,DIM=1) 
            do jk2 = 1, ubound(Eigenvalue,DIM=1)
               input_replace(jk1,jk2)=SUM(evect(jk1,:)*temp2(:,jk2))
            enddo
         enddo
c        enforce proper range on the reformed correlation matrix
c           scale based on diagonal values
         do jk1 = 1, ubound(Eigenvalue,DIM=1)
            correction = 1.0/input_replace(jk1,jk1)
            do jk2 = 1, ubound(Eigenvalue,DIM=1)
              input_replace(jk1,jk2)=input_replace(jk1,jk2)*correction
            enddo
         enddo
c           ensure interval (-1, 1)
         do jk1 = 1, ubound(Eigenvalue,DIM=1)
           do jk2 = 1, ubound(Eigenvalue,DIM=1)
             if (input_replace(jk1,jk2)>1.0) 
     &             input_replace(jk1,jk2)=1.0
             if (input_replace(jk1,jk2)<-1.0) 
     &             input_replace(jk1,jk2)=-1.0
           enddo
         enddo
         if( iflag <= -1) then 
           do jk1 = 1, ubound(Eigenvalue,DIM=1)
              write (6,882) jk1, (input_replace(jk1, jk2), jk2=1,
     &                      ubound(Eigenvalue,DIM=1))
 882          format (1x, 'Input replace: ', i5, /, 
     &            (1x, '               ', 5x, 5g14.7) )
           enddo
         endif
c
c        Output replacement correlation matrix
c
         deallocate(evect,
     &            diag_new, diag, evect_t, 
     &            input_dup, STAT=iStat)
         if (iStat/=0) then 
           write (6,25000) iStat
           stop 'Memory-(de)allocation error in diag eigen_out.'
         endif                
c
      elseif (.not. pos_def) then
c
c       For icon(1)=14, set default input_replace since it is reset
c       and set ipass to 750 to stop looping
c
        ipass = 750
        input_replace = cor
c
      endif
c
      deallocate(corWorkingArray, Eigenvalue, VectorOut, TempVector,
     $          STAT=iStat)
      if (iStat/=0) then 
        write (6,25000) iStat
        stop 'Memory-(de)allocation error in eigen_out.'
      endif
c      write (6,782)
 782  format (1x, 'Force run stop at return from eigen_out')
c      stop 'development stop'
      
25000   format ('Subroutine eigen_out memory-(de)allocation error: ',
     &          2i5)
     
      return
      end subroutine eigen_out
