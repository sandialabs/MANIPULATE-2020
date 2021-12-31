      subroutine prune (material, ilead, itrail, length, nchar)
      character*(*) material
      character*1 blank 
      data blank /' '/
c
c        determine number of leading and trailing blanks
c
         nchar = len(material)
         ilead = 0
         do i=1,nchar
           if ( material(i:i) .eq. blank) then 
               ilead = ilead + 1
           else
               go to 6012
           endif
         enddo
 6012    continue
         itrail = 0 
         do i=1,nchar
           if ( material(nchar+1-i:nchar+1-i) .eq. blank) then 
               itrail = itrail + 1
           else
               go to 7012
           endif
         enddo
 7012    continue
         length =(nchar-itrail) - (ilead+1) + 1
      return
      end
