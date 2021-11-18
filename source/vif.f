c      subroutine vif(koma,ok)
c      dimension alf(80),ilf(80),item(41),ytem(41),kar(48),
c     1          part(2),sign(2),blot(4),blan(3),add(47)
c      common/rv/alf,ytem,last
c      equivalence (alf(1),ilf(1)),(item(1),ytem(1))
c      logical dec,new,enew,ok
c      data kar/1h ,1h,,1h0,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1h.,1h+,
c     1 1h-,1he,1ha,1hb,1hc,1hd,1hf,1hg,1hh,1hi,1hj,1hk,1hl,1hm,
c     2 1hn,1ho,1hp,1hq,1hr,1hs,1ht,1hu,1hv,1hw,1hx,1hy,1hz,1h/,
c     3 1h(,1h),1h=,1h$,1h*,1h'/
c       data add/4h,,,,,4h0000,4h1111,4h2222,4h3333,4h4444,4h5555,4h6666,
c     1 4h7777,4h8888,4h9999,4h....,4h++++,4h----,4heeee,4haaaa,4hbbbb,
c     2 4hcccc,4hdddd,4hffff,4hgggg,4hhhhh,4hiiii,4hjjjj,4hkkkk,4hllll,
c     3 4hmmmm,4hnnnn,4hoooo,4hpppp,4hqqqq,4hrrrr,4hssss,4htttt,4huuuu,
c     4 4hvvvv,4hwwww,4hxxxx,4hyyyy,4hzzzz,4h1111,4h((((,4h)))),4h====,
c     5 4h$$$$,4h****,4h''''/
c       data blot/zff000000,z00ff0000,z0000ff00,z000000ff/
c       data blan/z00404040,z00004040,z00000040/
c      external and,or
c
c                         initialize all data items
c
c      do 1 jit=1,41
c    1 ytem(jit)=0.
c      l=1
c      jit=0
c
c                test whether item is numeric (first character
c                is digit, decimal, or sign), or alphabetic
c                (any other first character) or blank or comma
c
c    5 do 10 j=1,15
c      jch=j
c      if(ilf(l).eq.kar(jch))go to 200
c   10 continue
c
c                alphabetic item (word)
c
c      j1=16
c   12 do 30 kk=1,4
c      kb=kk-1
c
c               determine kkth character
c
c      do 15 j=j1,48
c      jch=j
c      if(ilf(l).eq.kar(jch))go to 202
c   15 continue
c
c                write error message for inscrutable character
c
c      write(3,20) l
c   20 format('0data card contains an inscrutable character in column',
c     1  i3,'.')
c      ok=.false.
c      return
c  202 if(jch-2)32,21,25
c
c           check whether comma is admissible character or spacer
c
c   21 if(koma-2)32,25,25
c
c                   add character to kkth position in word
c
c   25 jad=jch-1
c      if(kk.eq.1) jit=jit+1
c       temp=and( add(jad),blot(kk))
c       ytem(jit)=or(ytem(jit),temp)
c   29 l=l+1
c   30 j1=1
c      goto34
c   32 ytem(jit)=or(ytem(jit),blan(kb))
c   34 if(l.eq.81) go to 110
c      if(ilf(l).eq.kar(1)) go to 105
c      if(ilf(l).eq.kar(2))go to 204
c      go to 12
c  200 if(jch-2)105,105,35
c  204 if(koma-2)105,12,12
c
c                          numeric item (data value)
c
c   35 jit=jit+1
c
c        initialize parts of data item, signs, and decimal multiplier
c
c      dpl=0.
c      do 40 ip=1,2
c      part(ip)=0.
c   40 sign(ip)=1.
c
c                 first part of data item (significant value)
c
c      ip=1
c
c          initialize flags for decimal, new item, and new exponent
c
c      dec=.false.
c      new=.true.
c   45 enew=.false.
c
c         test whether character is digit, decimal, sign, or letter e
c
c      if(jch-13) 50,55,65
c
c              advance decimal multiplier, if decimal is flagged
c
c   50 if(dec)dpl=dpl-1.
c
c                            identify digit value
c
c      val=jch-3
c
c                 add digit to significant value or exponent
c
c      part(ip)=10.*part(ip)+val
c      go to 60
c
c                                flag decimal
c
c   55 dec=.true.
c   60 if(l-80)80,100,100
c
c         test whether sign is for significant value or for exponent
c
c   65 if(new)go to 70
c
c                     second part of data item (exponent)
c
c      ip=2
c
c                       stop decimal multiplier advance
c
c      dec=.false.
c
c                          flag entry into exponent
c
c      enew=.true.
c
c          test for negative significant value or negative exponent
c
c   70 if(jch.eq.15) sign(ip)=-1
c      if(l.lt.80) go to 80
c
c        write error message for item ending with a sign or the letter e
c
c      write(3,75) kar(jch)
c   75 format('0last data item on card ends in last card column with the 
c     1character ',a1,', which is not permitted.')
c      ok=.false.
c      return
c
c                            flag entry into item
c
c   80 new=.false.
c
c                          advance to next character
c
c      l=l+1
c      do 85 j=1,16
c      jch=j
c
c          determine whether character is or is not a blank or comma
c
c      if(ilf(l).eq.kar(jch))go to 206
c   85 continue
c
c                     write error message for numeric item
c                        containing unallowed character
c
c      write(3,90) l
c   90 format('0data card contains numeric data item which includes a cha
c     1racter other than a'/' numeric digit, decimal, sign, or the letter
c     2 e.  this character is in column',i3,'.')
c      ok=.false.
c      return
c  206 if(jch-2)95,100,45
c
c                if blank, test whether end of item or imbedded
c                 in exponent immediately following letter e
c
c   95 if(.not.enew) go to 100
c
c             continue, if blank immediately follows the letter e
c
c      enew=.false.
c      go to 60
c
c              calculate exponent, including decimal multiplier
c
c  100 expo=sign(2)*part(2)+dpl
c
c                    calculate value of numeric data item
c
c      ytem(jit)=sign(1)*part(1)*10.**expo
c
c                            test for end of card
c
c  105 if(l.eq.80) go to 110
c
c                          advance to next character
c
c      l=l+1
c      go to 5
c  110 last=jit
c      ok=.true.
c      return
c      end
      subroutine vif(koma,ok)
      implicit real*4(a-h,o-z), integer*4(i-n)
      integer*4 blot,blan,add
      logical*1 dec,new,enew,ok
      dimension alf(80),ilf(80),item(41),ytem(41),kar(48),
     1          part(2),sign(2),blot(4),blan(3),add(47)
      common/rv/alf,ytem,last
      character*1 cchar
      character*1 char_new
      integer char
      equivalence (alf,ilf),(item,ytem)
      data kar/1h ,1h,,1h0,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9,1h.,1h+,
     1 1h-,1he,1ha,1hb,1hc,1hd,1hf,1hg,1hh,1hi,1hj,1hk,1hl,1hm,
     2 1hn,1ho,1hp,1hq,1hr,1hs,1ht,1hu,1hv,1hw,1hx,1hy,1hz,1h/,
     3 1h(,1h),1h=,1h$,1h*,1h-/
       data add/4h,,,,,4h0000,4h1111,4h2222,4h3333,4h4444,4h5555,4h6666,
     1 4h7777,4h8888,4h9999,4h....,4h++++,4h----,4heeee,4haaaa,4hbbbb,
     2 4hcccc,4hdddd,4hffff,4hgggg,4hhhhh,4hiiii,4hjjjj,4hkkkk,4hllll,
     3 4hmmmm,4hnnnn,4hoooo,4hpppp,4hqqqq,4hrrrr,4hssss,4htttt,4huuuu,
     4 4hvvvv,4hwwww,4hxxxx,4hyyyy,4hzzzz,4h1111,4h((((,4h)))),4h====,
     5 4h$$$$,4h****,4h----/
      data blot/'377'o,'177400'o,'77600000'o,'37700000000'o/
      data blan/'04010020000'o,'04010000000'o,'04000000000'o/
c
c                         initialize all data items
c
      do 1 jit=1,41
    1 ytem(jit)=0.
      l=1
      jit=0
c
c
c                test whether item is numeric (first character
c                is digit, decimal, or sign), or alphabetic
c                (any other first character) or blank or comma
c
    5 do 10 j=1,15
      jch=j
      if(ilf(l).ne.kar(jch)) goto 10
      if(jch-2) 105,105,35
   10 continue
c
c                alphabetic item (word)
c
      j1=16
   12 do 30 kk=1,4
      kb=kk-1
c
c               determine kkth character
c
      do 15 j=j1,48
      jch=j
      if(ilf(l).eq.kar(jch))go to 202
   15 continue
c
c                write error message for inscrutable character
c
c
      write(6,20) l
   20 format('0data card contains an inscrutable character in column',
     1  i3,'.')
      ok=.false.
      return
  200 if(jch-2)105,105,35
  202 if(jch-2)32,21,25
c
c           check whether comma is admissible character or spacer
c
   21 if(koma-2)32,25,25
c
c                   add character to kkth position in word
c
   25 jad=jch-1
      if(kk.eq.1) jit=jit+1
      item(jit)=ior(item(jit),iand(add(jad),blot(kk)))
   29 l=l+1
   30 j1=1
      goto 34
   32 item(jit)=ior(item(jit),blan(kb))
   34 if(l.eq.81) go to 110
      if(ilf(l).eq.kar(1)) go to 105
      if(ilf(l).eq.kar(2))go to 204
      go to 12
  204 if(koma-2)105,12,12
c
c                          numeric item (data value)
c
   35 jit=jit+1
c
c        initialize parts of data item, signs, and decimal multiplier
c
      dpl=0.
      do 40 ip=1,2
      part(ip)=0.
   40 sign(ip)=1.
c
c                 first part of data item (significant value)
c
      ip=1
c
c          initialize flags for decimal, new item, and new exponent
c
      dec=.false.
      new=.true.
   45 enew=.false.
c
c         test whether character is digit, decimal, sign, or letter e
c
      if(jch-13) 50,55,65
c
c              advance decimal multiplier, if decimal is flagged
c
   50 if(dec)dpl=dpl-1.
c
c                            identify digit value
c
      val=jch-3
c
c                 add digit to significant value or exponent
c
      part(ip)=10.*part(ip)+val
      go to 60
c
c                                flag decimal
c
   55 dec=.true.
   60 if(l-80)80,100,100
c
c         test whether sign is for significant value or for exponent
c
   65 if(new)go to 70
c
c                     second part of data item (exponent)
c
      ip=2
c
c                       stop decimal multiplier advance
c
      dec=.false.
c
c                          flag entry into exponent
c
      enew=.true.
c
c          test for negative significant value or negative exponent
c
   70 if(jch.eq.15) sign(ip)=-1
      if(l.lt.80) go to 80
c
c        write error message for item ending with a sign or the letter e
c
      write(6,75) kar(jch)
   75 format('0last data item on card ends in last card column with the
     1character ',a1,', which is not permitted.')
      ok=.false.
      return
c
c                            flag entry into item
c
   80 new=.false.
c
c                          advance to next character
c
      l=l+1
      do 85 j=1,16
      jch=j
c
c          determine whether character is or is not a blank or comma
c
      if(ilf(l).eq.kar(jch))go to 206
   85 continue
c
c                     write error message for numeric item
c                        containing unallowed character
c
      write(6,90) l
   90 format('0data card contains numeric data item which includes a cha
     1racter other than a'/' numeric digit, decimal, sign, or the letter
     2 e.  this character is in column',i3,'.')
      ok=.false.
      return
  206 if(jch-2)95,100,45
c
c                if blank, test whether end of item or imbedded
c                 in exponent immediately following letter e
c
   95 if(.not.enew) go to 100
c
c             continue, if blank immediately follows the letter e
c
      enew=.false.
      go to 60
c
c              calculate exponent, including decimal multiplier
c
  100 expo=sign(2)*part(2)+dpl
c
c                    calculate value of numeric data item
c
      ytem(jit)=sign(1)*part(1)*10.**expo
c
c                            test for end of card
c
  105 if(l.eq.80) go to 110
c
c                          advance to next character
c
      l=l+1
      go to 5
  110 last=jit
      ok=.true.
      return
      end
