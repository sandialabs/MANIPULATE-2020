      subroutine terp1 (x1,y1,x2,y2,x,y,i)
c    ******************************************************************
c     interpolate one point
c     (x1,y1) and (x2,y2) are the end points of the line
c     (x,y) is the interpolated point
c     i is the interpolation code
c     ******************************************************************
      common/terp6/thr6
c
      if (y1.eq.y2) go to 100
      if ( i .lt. 1 .or. i .gt. 5) then 
          write (6,901) i
 901      format (1x, '*** error in interpolation code *** ',
     1     i5)
          stop 'terp-err'
      endif
      go to (100,200,300,400,500),i
c
c     ***y is constant
  100 y=y1
      return
c
c     ***y is linear in x
  200 y=y1+(x-x1)*(y2-y1)/(x2-x1)
      return
c
c     ***y is linear in ln(x)
  300 y=y1+alog(x/x1)*(y2-y1)/alog(x2/x1)
      return
c
c     ***ln(y) is linear in x
  400 y=y1*exp((x-x1)*alog(y2/y1)/(x2-x1))
      return
c
c     ***ln(y) is linear in ln(x)

  500 if (y1.eq.0.0) go to 100
      y=y1*exp(alog(x/x1)*alog(y2/y1)/alog(x2/x1))
      return
      end

