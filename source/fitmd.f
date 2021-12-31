      function fitmd( x, nx, xarray, yarray,icode)
      dimension xarray(nx), yarray(nx)
c
c     determine is xarray is increasing or deceasing
c
      ipick = 1
 101  continue
      if ( xarray(ipick+1) .gt. xarray(ipick)) then 
            idirection = 1
      elseif ( xarray(ipick+1) .lt. xarray(ipick)) then 
            idirection = -1
      else
            ipick = ipick + 1
            if ( ipick .ge. nx) then
c              all values are the same - return an answer
               if ( x .eq. xarray(1)) then 
                   fitmd = yarray(1)
                   return
               endif
               write (6,901) 
901            format (1x, 'error in fitmd - constant x ')
               stop 'fitmd-cnst'
            endif            
            go to 101
      endif
c
c     check array bounds
c
      if ( (x .lt. xarray(1) .and. idirection .eq.  1) .or. 
     1      (x .gt. xarray(1) .and. idirection .eq. -1)) then
            fitmd = 0.0
            return
      endif      
      if ( (x .lt. xarray(nx) .and. idirection .eq. -1) .or. 
     1      (x .gt. xarray(nx) .and. idirection .eq.  1)) then
            fitmd = 0.0
            return
      endif
c     bracket the array points
c
      ipick = 0
      if ( idirection .eq. 1) then 
         do i=1,nx
            place = x - xarray(i)
            if ( place .lt. 0.0) go to 909
            ipick = i
         enddo
      elseif ( idirection .eq. -1) then
         do i=1,nx
            place = xarray(i) - x
            if ( place .lt. 0.0) go to 909
            ipick = i
         enddo
      endif
      stop 'fitmd-err'
909   continue
c
c     perform interpolation
c
      ipick1 = ipick + 1
      call terp1(xarray(ipick), yarray(ipick), xarray(ipick1),
     1           yarray(ipick1), x, y, icode)
      fitmd = y
      return
      end

