       program manipulate
c
c      ***********************************************
c      program to manipulate-2020 njoy output cross section 
c      data.  
c
c      ************************************************
cpjg
cpjg 6/27/2021 - Changes made to legacy manipulate code to rehost under Linux and GitHub
cpjg
cje
cje changes to this routine made 06 march 1992:
cje  
cje    1. elimination of calls to routine that were specifically
cje       designed around the matxs ( either ascii or binary )
cje 	  files.
cje    2. for the icon(1)=3 option only a 
cje       groupr tape will be read and tabular output will be
cje       be presented for the reaction rate of interest.
cje       the new subroutine is grprin and produces an output
cje       echoed to the standard output.
cje
cpjg
cpjg changes made to add icon(1)=-5 statistics option
cpjg  
c
c      control parameters:
c
c         icon(1)    type of action
c              = 9       end  
c              = 1       combine cross section elements - logic unclear, PJG 10/17/2012
c              = 2       prepare component plot file format output
c              = 3       display a matxsr file from njoy
c              =-3          same as 3 for xxx.plt_89 file
c              = 4       multiple plot file option
c              = 5       sum and plot components
c              =-5       perform statistical ananlysis on a series of functions  
c              = 6       fold source and response
c              = 7       weight and plot function
c              =-7       reciprocal weight and plot function
c              = 8       file difference
c              =-8       covariance manipulation options
c              =10       interpolate matxsr file
c              =11       interpolate tabulated file
c              =12       expand by interpolation a tabulated file
c              =13       expand an analytic function in 640 group structure
c              =14       combine covariance files in LSL format (common energy structure) - NEW
c              =15       extract PKA recoil spectrum for given reaction/energy - NEW
c              =16       extract integrated PKA recoil spectrum for as a function of energy - NEW
c
c         icon(2)    input energy grid/data format
c              = 0       njoy sand-ii 620/640/725/770 point neutron cross section
c                        (energies in ev, response data to be multiplied
c                         by 10**-3 to yield kerma units of mev-mb)
c              = 1       j. kelly sand-ii pun spectrum format
c                        (default quantity differential number spectrum)
c                          compute iselect_material = 1 number fraction
c                                                     2 energy fraction
c                                                     3 differential number fraction
c                                                     4 differential energy fraction
c                                                     5 integral number
c              = 2       same as 1
c              = 3       flash x-ray curr spectral input
c              = 4       list directed (plot output)
c              = 5       njoy sand-ii /770 point neutron pka recoil spectrum
c                        (energies in ev, response data to be multiplied
c                         by 10**-3 to kev)
c              = 6       njoy matxsr interface file - arbitrary structure
c              = 9       user defined - xxx.flux file 
c              =10       640 group response structure = xxx.res
c              =11       arbitrary response structure = xxx.res
c
c         icon(3)    display/output special option
c              = 0       no action
c              = 1       punch in damout format
c              = 2       punch in tdam format
c
c         icon(4)    print energy order
c              = 0       punch high to low
c              = 1        punch low to high
c  
c         icon(5)    sand structure
c              = 0        640 groups - SAND-II
c              = 1        770 groups - SAND-II extended
c              = 2         89 groups - NuGET neutron
c              = 3         48 groups - NuGET gamma
c              = 5        175 group - Vitamin-J
c              = 6        725 group - IAEA 
c
c         icon(6)    inhibit renormalization
c              = 0       normal
c              = 1       response, inhibit renormalization in filein
c
c         icon(7)    NJOY interface location
c              = 0       sync_Projects-Linux/NJOY2016 default
c              = 1       sync_sandialabs/NJOY2016 selection
c
c         icon(9)    print inhibit
c             = 0        default print
c             =>0        inhibit print by this level
c             < 0        enhance print (debug option)
c
c         icon(10)   label replacement
c             = 0        no action
c             = 1        read new legend information
c
c         icon(11)   special action
c             = 0        not used
c             = 1        on pka plot, divide entry 1 by entry 3
c                        to produce average recoil energy - 
c                        when icon(2) = 5
c             = 2        on icon(1)=7 combine, punch plt-format
c                        interface file on indicated dir:filename
c                        (tbd)
c          icon(12)   reverse order of data read - only partially implemented
c             = 0        no effect
c             = 1        reverse
c          icon(13)   attenuation function
c             = 1        if equal to 1 and icon(1)=6, then
c                        use new function for exponential attenuation
c                        with coefficient alpha 
c             = 2        if equal to 2 and icon(1)=6, then 
c                        use sensor self-shielding correction factor
c          icon(14)   new title
c             = 0        use first card in input stream
c             = 1        read-in new title
c
c          icon(15)   flip order of data for punch
c             = 0         do nothing
c             + 1        invert order
c
c         icon(16)    re-direct to look for pka files
c             = 0         do nothing
c             = 1         look in njoy/output-pka directory
c             = 2         look in NJOY9-2012 for groupr/pka files 
c
c        icon(17)     point interpolation for grprin option
c             = 0         no, not used
c             = 1         yes, interpolate
c
c        icon(18)     precision control for sum/difference
c             = 0         no action
c             = 1         force 0.0 response for differences < 1.e-4 of maximum element
c
c        icon(19)     inversion of output energy order
c             = 0         no change
c             = 1         invert order 
c
c      ************************************************
c
c      nenergy        number of energy bins
c      energy(i)      energy grid
c      array(i,j)     input data array -
c                              if source = number fraction
c                              if kerma =  mev-mb
c                              if other =  ???
c                        i = energy point
c                        j = material/reaction
c 
c
c      ************************************************
c
       external stat
       integer iresp_call
       character*80 title
       character*80  label,  mcnploc
       character*250 file_name, file_name2
       character*3 number_file
       character*80 job
       character*100 file_id
cje
cje addition character assignments.
cje
       character*80 outfile, putfile
       character*100 file_corr
       integer spc_tag_end, xsec_tag_end
       character*50 spc_tag, xsec_tag
       common /tag/ spc_tag, xsec_tag, spc_tag_end, xsec_tag_end
       common /sself/ file_corr
       common /whatever/ outfile
cje
       common /datain/ nenergy, energy(1001), array(1000,15)
     1         , emid(1001)
       dimension range(1000)
       dimension percent(1000)
       common /fpnew/ alpha
       common /datahld/ nenergy_hld(15), energy_hld(1001,15), 
     1          array_hld(1000,15)
     1         , emid_hld(1001,15)
       common /guide/ icon(40)
       common /io/ nt5, nt6, nfile, nplot, npun
       common /pltlab/ lab_lead(15), lab_trail(15), lab_nchar(15),
     1 lab_len(15), icurve(15), lab_file(15)
       character*80 lab_file
       character*20 lab, labtdam(3)
      character*145 idir, jdir, kdir, ldir
      character*250 new_file_name
      character*106 optical, xoptical
      character*15 ename
c      character*55 outfile
      common /location / optical, idir, jdir, kdir
      character*250 name
       common /misc/ lab(20), imaterial, noption, moption
      common /statsum/ arr_mean(1000), arr_sq(1000), 
     1  arr_std(1000),
     1  arr_var(1000)
      common /statinfo/ efactor, number_of_files, xoptical
      common /acov1/ icov(40),iself,ovr
      character*250 ovr
      character*256 epath
      common /new_io/ epath
       data (lab(jk),jk=1,3) /
     2      '    total kerma     ', '  kinematic limit   ',
     1           'displacement damage ' /
       data labtdam /'total kerma', 'kinematic limit',
     1           'displacement damage' /
       data icurve /1,13,4,2,3,5,6,7,8,9,10,11,12,14,15/
 1     format (a80)
 5434  continue
c       stop 'ok'
       imat2 = 1
       ename = 'job'
       lename = lnblnk(ename)
       call getenv(ename(1:lename), job)
       write (*,8623) job
 8623  format (1x, 'picked up exported job name = ', a)
       ename = 'opt'
       lename = lnblnk(ename)
       leopt = lnblnk(xoptical)
       call getenv(ename(1:lename), xoptical)
       leopt = lnblnk(xoptical)
       write (*,8698) leopt, xoptical(1:leopt)
 8698  format (1x, 'xoptical length ', i5, ' set = ', a) 
c       if (xoptical .ne. '/odsk1' .and. 
c     1     xoptical .ne. '/odsk2' .and. 
c     2     xoptical(1:leopt) .ne. '/esata' .and.
c     3     xoptical(1:leopt) .ne. '/esata/codes') then 
c           write (*,5013) xoptical
c           xoptical=''
c 5013      format (1x, 'blank optical default filled', 
c     1     1x, a6)
c       endif
       leopt = lnblnk(xoptical)
c       jdir = '/app/manipulate-2/response/'
       jdir = 'response/'
       jblank2 = lnblnk(jdir)
c       kdir = '/app/manipulate-2/spectrum/'
       kdir = 'spectrum/'
       kblank2 = lnblnk(kdir)
       nt6 = 6
       nfile = 20
       nt5 = 5
       npun = 17
       nplot = 77
cje
cje open pun and pun2 files ( 17, and 27 ) for output
cje
       open (unit=nt5,file='manipulate.inp',status='old')
       open (unit=nt6,file='manipulate.out',status='unknown')
       open (unit=39,file='manipulate.ext',status='unknown')
       open(unit=npun,file='pun',status='unknown')
       open(unit=27,file='pun2',status='unknown')
cpjg
c
       open (unit=30, file='fold.increment',status='unknown')
       open (unit=31, file='fold.sum', status='unknown')
c
cpjg
c
cje
       read(nt5,1,end=101) title
       spc_tag = "default"
       xsec_tag = "default"
       spc_tag_end = lnblnk(spc_tag)
       xsec_tag_end = lnblnk(xsec_tag)
 101   continue
       close (unit=20,err=1023)
 1023  continue
       write (nt6,2301) title
2301   format (1h1, /,/,5x, 'MANIPULATE Version 2020',/, 5x, 
     1         'title: ', a80,/,/)
       read (nt5,102,end=999) (icon(jk),jk=1,40)
       if ( icon(14) .eq. 1) then 
c            read (nt5,345,end=304) title
c
c            Modify, instead of new title, read tags
c
             read (nt5, *) spc_tag, xsec_tag
             spc_tag_end = lnblnk(spc_tag)
             xsec_tag_end = lnblnk(xsec_tag)
c
 345        format (a80)
       endif
 304   continue
       write (nt6,2302) (icon(jk),jk=1,40)
2302   format (1x, 'control parameters: ', 40i2)
       scale = 0.0
102    format (40i2)
       if ( icon(1) .eq. -8) then 
c
c       covariance manipulation modules
c
           if ( icon(5) .eq. 0) then
               nenergy = 640
           elseif ( icon(5) .eq. 1) then 
               nenergy = 770
           elseif ( icon(5) .eq. 2) then 
               nenergy = 89
           elseif ( icon(5) .eq. 5) then 
               nenergy = 175
           elseif ( icon(5) .eq. 6) then 
               nenergy = 725
           else
               write (6,5245) nenergy, icon(5)
 5245          format (1x, 'icon(5) error ', 2i6)
               stop 'icon(5)'
           endif
           read (nt5, *) (icov(jk), jk=1,3)
           if( icov(2) .ne. 0) then 
                write (nt6, 4502) icov(2)
 4502           format (1x, '*** Force fine definition of response',
     1          ' function: icov(2) =', i5)
           endif
           if( icov(3) .ne. 0) then 
                if ( icov(3) .eq. 1) then 
                    write (nt6, 4603) icov(3)
                elseif ( icov(3) .eq. 2) then 
                    write (nt6, 4503) icov(3)
                elseif ( icov(3) .eq. 3) then 
                    write (nt6, 4504) icov(3)
                else
                    write (nt6, 4703)
                    stop 'ICOV-3'
                endif
 4603           format (1x, '*** Force uncorrelated std. dev.',
     1          ' e.g. native energy grid diagonal covariance matrix:',
     1          ' icov(3) =', i5)
 4503           format (1x, '*** Force uncorrelated std. dev.',
     1          ' e.g. 640 bin diagonal covariance matrix:',
     2          ' icov(3) =', i5,/,
     3          1x,' This is only used for testing,the result has',
     4          ' no physical meaning')
 4504           format (1x, '*** Force totally correlated std. dev.',
     1          ' e.g. 640 bin full covariance matrix:',
     2          ' icov(3) =', i5,/,
     3          1x,' This is only used for testing,the result has',
     4          ' no physical meaning')
 4703           format (1x, '*** Forces uncorrelated std. dev.',
     1          'illegal icov(3) option: icov(3) = ', i5)
           endif
c
c          set flag is this is a covariance call for a response
c             else set to 0 for a source term
c
           iresp_call = 1
           if ( icon(9) < 0) then 
              write (6, 691) iresp_call
 691          format (1x, 'MANIPULATE call to covlate - for src/rsp ',
     &           'function ', i5)
           endif
           call covlate(iresp_call)
           go to 101
       endif
       if ( icon(1) .eq. 1) then 
c
c          combine cross sections
c
            go to 101
       elseif (icon(1) .eq. 14) then 
c
c           combine LSL-fomat covariance files
c
            read (nt5, *) (icov(jk), jk=1,3)
            if( icov(2) .ne. 0) then 
                write (nt6, 4502) icov(2)
            endif
            if( icov(3) .ne. 0) then 
                if ( icov(3) .eq. 5) then 
                    write (nt6, 4603) icov(3)
                else
                    write (nt6, 4703)
                    stop 'ICOV-3a'
                endif
            endif
            if (icon(9) < 0) then 
                write (6,7812) icon(1), (icov(jk), jk=1,3)
 7812           format (1x, 'MANIPULATE enter ICON 14 grab: ', 4i5)
            endif
            call cov_combine
            go to 101
c
       elseif ( icon(1) .eq. 15) then
c
c           extract PKA spectrum
c
            read(nt5,*) file_name,irxn,outfile, select_energy
            call pka_grprin(file_name,irxn,select_energy)
            go to 101
       elseif ( icon(1) .eq. 16) then
c
c           extract integated PKA spectrum
c
            read(nt5,*) file_name,irxn,outfile,scale
            call integral_pka_grprin(file_name,irxn)
            go to 101
c
       elseif (icon(1) .eq. 3 .or. icon(1) .eq. -3) then 
c           read matxsr file and report  - option 3
c           interface same logic with a plt_89 file - option -3
            noption = 1
cje         read(nt5,*) file_name
cje         call prune (file_name, ilead, itrail, length, nchar)
            if ( icon(1) .eq. 3) then 
cje
cje if icon(1) = 3, then read the xs data into array using the 
cje groupr ouput tape option from njoy.
cje
               read(nt5,*) file_name,irxn,outfile,scale
               call grprin(file_name,irxn)
cje            open(unit=nfile, form='unformatted',
cje  1              file=file_name, status='old', err=909)
cje            call pmatxs(nfile,nt6)
            elseif (icon(1) .eq. -3) then
cje
cje add read here so that other option follows old format.
c   hardwire for 640 group[ without header number
cje
               read(nt5,*) file_name,irxn,outfile,scale
cje            call prune (file_name, ilead, itrail, length, nchar)
cje
               nlen = lnblnk(file_name)
               name = xoptical(1:leopt)//
     2          jdir(1:jblank2)//file_name(1:nlen)
               open(unit=nfile, form='formatted',
     1              file=name, status='old', err=909)
               imaterial = 1
c               read (nfile,*) nenergy
               if ( icon(5) .eq. 0) then
                   nenergy = 640
               elseif ( icon(5) .eq. 1) then 
                   nenergy = 770
               elseif ( icon(5) .eq. 2) then 
                   nenergy = 89
               elseif ( icon(5) .eq. 5) then 
                   nenergy = 175
               elseif ( icon(5) .eq. 6) then 
                   nenergy = 725
               else
                   write (6,7245) nenergy, icon(5)
 7245              format (1x, 'icon(5) error ', 2i6)
                   stop 'icon(5)'
               endif
               do jkl=1,nenergy
                 do jkp = 1,10
                    array(jkl,jkp) = 0.0
                 enddo
               enddo
               if ( icon(12) .eq. 0) then
                   read (nfile,*) 
     1                (emid(jkl), array(jkl,1), jkl=1,nenergy) 
               elseif ( icon(12) .eq. 1) then
                   read (nfile,*) 
     1                (emid(jkl), array(jkl,1), jkl=nenergy,1,-1) 
               else
                   stop 'icon(12) err'
               endif
               close (unit=nfile)
            endif
            if ( icon(3) .eq. 1) then 
c
c              punch data in damout format
cje            close (unit=nfile)
cje            read (nt5,*) scale
cje            if ( icon(1) .eq. 3) then 
cje               open(unit=nfile, form='unformatted',
cje  1                file=file_name
cje  2                 , status='old', iostat=ilook, err=909)
cje            elseif ( icon(1) .eq. -3) then 
cje                continue
c                  open(unit=nfile, form='formatted',
c     1                 file=file_name, status='old', err=909)
cje            endif
cje            if ( icon(1) .eq. 3) then 
cje               do jkl=1,nenergy
cje                 do jkp = 1,10
cje                    array(jkl,jkp) = 0.0
cje                 enddo
cje               enddo
cje this will be replaced by the groupr input call equivalent
cje               call matxsin(nfile,nt6)
cje            elseif ( icon(1) .eq. -3) then 
cje               continue
cje            endif
c
c               write (npun,2341) ((array(jk,ik)*scale,ik=1,1),
c     1                             jk=nenergy,1,-1)
               write (npun, 5745) title(1:60)
               write (npun,2341) ((array(jk,ik)*scale,ik=1,1),
     1                             jk=1,nenergy)
 2341          format (5(e12.5,4x))
c              special punch format without zeros
c              modify to allow for zeros being turned into small numbers
c
               ilow = 1
               ihigh = 641
c              special look loop
                   jkp = 0
 3045              continue
                   jkp = jkp + 1
                   if ( array(jkp,1) .le. 1.e-32) then 
                       ilow = jkp
                       go to 3045
                   else
                       continue
                   endif
                   jkp = 641
 3046              continue
                   jkp = jkp - 1
                   if ( array(jkp,1) .le. 1.e-32) then 
                       ihigh = jkp
                       go to 3046
                   else
                       continue
                   endif
               write (27, 5745) title(1:60)
 5745          format (1x, a60)
               ilow_x = ilow
               ihigh_x = ihigh - 1
c               ihigh_x1 = 640 - ihigh_x + 1
c               ilow_x1  = 640 - ilow_x  + 1
               ihigh_x1 = ilow_x 
               ilow_x1  = ihigh_x
               write (27,6190)  outfile
 6190          format (1x, 'file = ', a)
               write (27,6923)  ihigh_x1, ilow_x1
 6923          format (1x, i5, 1x, i5)
c               write (27,2341) (( 
c     1          array(jk,ik)*scale,
c     1          ik=1,1),  jk=ihigh_x,ilow_x,-1)
c
               write (27,2341) (( 
     1          array(jk,ik)*scale,
     1          ik=1,1),  jk=ilow_x,ihigh_x)
c              add option 2 here also
cje            do jkl=1,nenergy
cje                emid(jkl) = 0.5*(energy(jkl) + energy(jkl+1))
cje            enddo
cje            do jkl=1,4
cje               npun2 = 31
cje               dummy = 0.0
cje               if ( jkl .eq. 2) go to 9558
cje               if ( jkl .eq. 1) then 
cje                  file_name2 = 'njoy89_jnk:xxxx_ker.dat'
cje               elseif ( jkl .eq. 3) then 
cje                  file_name2 = 'njoy89_jnk:xxxx_lmt.dat'
cje               elseif ( jkl .eq. 4) then 
cje                  file_name2 = 'njoy89_jnk:xxxx_dsp.dat'
cje              endif
cje               open(unit=npun2, form='formatted',
cje  1                file=file_name2
cje  2               , status='new', iostat=ilook, err=909)
cje               write (npun2, 2745) lab(jkl), title(1:60), nenergy, dummy
cje               write (npun2,2741) ( emid(jk)*1.e-6,  
cje  1             array(jk,jkl)*scale*1.e-3,
cje  1              jk=nenergy,1,-1)
cje               close (unit=npun2)
cje 9558             continue
cje            enddo
cje
cje print out the reaction on a separate tape for plotting 
cje with templegraph.  the tape will have the name outfile.###
cje with ### being the reaction of interest.
cje 
               i1 = int(irxn/100)
               i2 = int((irxn-i1*100)/10)
               i3 = int(irxn-i1*100-i2*10)
               iblank = lnblnk(outfile)
               if (i1 .ne. 0) then 
               outfile = outfile(1:iblank)//'.'//char(48+i1)//
     1         char(48+i2)//char(48+i3)
               elseif ( i2 .ne. 0) then 
               outfile = outfile(1:iblank)//'.'//
     1         char(48+i2)//char(48+i3)
               else
               outfile = outfile(1:iblank)//'.'//
     1         char(48+i3)
               endif
               ename = 'opt'
                lename = lnblnk(ename)
c               call getenv(ename(1:lename), optical)
c               if (optical .eq. ' ') then 
c                 optical='/odsk1'
c                 write (*,9013) optical
c 9013            format (1x, 'blank optical default filled', 
c     1           1x, a6)
c               endif
c               idir = '/app/manipulate-2/response/'
c               iblank2 = lnblnk(idir)
                efactor = 1.0
                if ( emid(1) .gt. 1.e6 .or.
     1               emid(nenergy) .gt. 1.e6) then 
                     efactor = 1.e-6
                endif
               if ( icon(9) .lt. 0) then 
                  write (6, 7412) leopt, jblank2, xoptical(1:leopt), 
     &                jdir(1:jblank2), outfile
 7412             format (1x, 'Main File 8: ', 2i5, /, (5x, a ))
               endif
               open(unit=8,file=xoptical(1:leopt)//
     2          jdir(1:jblank2)//outfile
     1         ,status='unknown')
               if ( icon(12) .eq. 1) then 
                  do i = nenergy,1,-1
                    if ( array(i,1) .le.0.0) then 
                                array(i,1) = 1.e-33
                    endif
                    write(8,2341)emid(i)*efactor,array(i,1)*scale
                  end do
               else
                  do i = 1,nenergy
                    if ( array(i,1) .le.0.0) then 
                                array(i,1) = 1.e-33
                    endif
                    write(8,2341) emid(i)*efactor,array(i,1)*scale
                  end do
               endif
               close(unit=8)
cje
c
c              plot files
c
               write (nplot,4312)
               write (nplot, 4313) 
               write (nplot, 4314)
               write (nplot, 8001) 
               lab(1) = '    total kerma     '
               lab(2) = 'dummy'
               lab(3) = '  kinematic limit   '
               lab(4) = 'displacement damage '
               lab(5) = 'dummy'
               lab(6) = 'dummy'
               lab(7) = 'dummy'
               if ( moption .le. 4) then 
                   moption = 4
               endif
               do jkl=1,moption
                  write (nplot,3001) lab(jkl)
                  write (nplot,2741) ( emid(jk)*1.e-6,  
     1             array(jk,jkl)*scale*1.e-3,
     1              jk=1,nenergy)
               enddo
               write (nplot, 8014)
            elseif ( icon(3) .eq. 2) then 
c
c              punch data in tdam format
               close (unit=nfile)
               read (nt5,*) scale
               open(unit=nfile, form='unformatted',
     1             file=file_name
     2              , status='old', iostat=ilook, err=909)
cje            call matxsin(nfile,nt6)
               do jkl=1,nenergy
                   emid(jkl) = 0.5*(energy(jkl) + energy(jkl+1))
               enddo
               do jkl=1,4
                  npun2 = 31
                  dummy = 0.0
                  if ( jkl .eq. 2) go to 9551
                  if ( jkl .eq. 1) then 
                     file_name2 = 'njoy89_jnk:xxxx_ker.dat'
                   elseif ( jkl .eq. 3) then 
                     file_name2 = 'njoy89_jnk:xxxx_lmt.dat'
                  elseif ( jkl .eq. 4) then 
                     file_name2 = 'njoy89_jnk:xxxx_dsp.dat'
                 endif
                  open(unit=npun2, form='formatted',
     1                file=file_name2
     2               , status='new', iostat=ilook, err=909)
                  write (npun2, 2745) lab(jkl), title(1:60), 
     1                nenergy, dummy
 2745             format (1x, a20, a60,/, 
     1            1x, ' 1 0 ', i3, ' 95.0 2 1 ', g14.7)
                  write (npun2,2741) ( emid(jk)*1.e-6,  
     1             array(jk,jkl)*scale*1.e-3,
     1              jk=nenergy,1,-1)
 2741             format (2(e12.5,4x))
                  close (unit=npun2)
 9551             continue
               enddo
c
c              plot files
c
               write (nplot,4312)
               write (nplot, 4313) 
               write (nplot, 4314)
               write (nplot, 8001) 
               lab(1) = '    total kerma     '
               lab(2) = 'dummy'
               lab(3) = '  kinematic limit   '
               lab(4) = 'displacement damage '
               lab(5) = 'dummy'
               lab(6) = 'dummy'
               lab(7) = 'dummy'
               if ( moption .le. 4) then 
                   moption = 4
               endif
               do jkl=1,moption
                  write (nplot,3001) lab(jkl)
                  write (nplot,2741) ( emid(jk)*1.e-6,  
     1             array(jk,jkl)*scale*1.e-3,
     1              jk=1,nenergy)
               enddo
               write (nplot, 8014)
            endif
            go to 101
       elseif ( icon(1) .eq. 2) then 
c
c           prepare plot file
c
c           input file name
c
            read(nt5,*) file_name
            call prune (file_name, ilead, itrail, length, nchar)
            ierr = 0
            call filein(icon(2), file_name, ierr, iselect_material)
c
c           replace file_name with new label
c
            if ( icon(10) .eq. 1) then 
               read(nt5,*) label
               file_name = label
               call prune (file_name, ilead, itrail, length, nchar)
            endif
c
c           output plot header info
c
c
c        write plot direction file header data -
c
            if ( icon(2) .eq. 5) then 
               nmat = 3
               write (nplot,4212)  file_name(ilead+1:nchar-itrail),
     1                      file_name(ilead+1:nchar-itrail)
               write (nplot, 4213) 
               write (nplot, 4214)
               lab(3) = 'total cross section -?'
               lab(2) = 'lower kin. kerma - ?'
               lab(1) = 'pka spectrum (kev) '
             elseif ( icon(2) .eq. 10 .or. icon(2) .eq. 11) then
               nmat = 1
               write (nplot,4212)  file_name(ilead+1:nchar-itrail),
     1                      file_name(ilead+1:nchar-itrail)
               write (nplot, 4213) 
               write (nplot, 4214)
               lab(1) = 'response'
               lab(2) = 'not used'
               lab(3) = 'not used '
             else
               nmat = 3
               write (nplot,4012)  file_name(ilead+1:nchar-itrail),
     1                      file_name(ilead+1:nchar-itrail)
               write (nplot, 4013) 
               write (nplot, 4014)
               lab(1) = '    total kerma     '
               lab(2) = '  kinematic limit   '
               lab(3) = 'displacement damage '
            endif
 4012       format ('**file**.',/,
     1           'generate a fancy plot.',/,
     2           'generate a fancy x log plot.',/,
     3           'x axis label text is "neutron energy (mev)".',/,
     4           'title text is "',a10,
     5           ' damage functions".',/,
     5           'y axis label text is "',a10,
     5                 ' damage (mev-mb)".',/,
     6           'x axis style is duplex.',/,
     7           'x axis origin is 2.00.',/,
     8           'x axis tick marks 5.',/,
     8           'x tick marks mode reversed .',/,
     8           'y tick marks mode reversed .',/,
     9           'y axis style is duplex.',/,
     a           'x log.')
 4013       format ('x grid.',/,
     1           'y log .',/,
     1           'x room 2.0 .',/,
     1           'y room 2.0 .',/,
     1           'y axis height 0.20   .',/,
     1           'x axis height 0.20   .',/,
     1           'y axis number mode commas .',/,
     2           'y grid.',/,
     3           'curve 1 symbol count 0.',/,
     3           'curve 2 symbol count 0.',/,
     3           'curve 3 symbol count 0.',/,
     3           'title style is triplex.',/,
     4           'frame on.',/,
     2           'curve 1 texture solid.',/,
     2           'curve 2 texture dashed.',/,
     2           'curve 3 texture 7.',/,
     2           'x integerize no. ',/,
     2           'y integerize no. ',/,
     3           'curve 1 thickness is 5.',/,
     3           'curve 2 thickness is 5.',/,
     3           'curve 3 thickness is 5.')
 4014       format ('legend frame on. ',/,
     1           'legend x origin 2.5  .',/,
     2           'legend y origin 5.2  .',/,
     3           'curve 1 label "total kerma", solid. ',/,
     3           'curve 2 label "kinematic kerma limit", dashed. ',/,
     3           'curve 3 label "displacement damage", 4 . ')
 4212       format ('**file**.',/,
     1           'generate a fancy plot.',/,
     2           'generate a fancy x log plot.',/,
     3           'x axis label text is "incident n',
     3           'eutron energy (mev)".',/,
     4           'title text is "',a10,
     5           ' pka spectrum".',/,
     5           'y axis label text is "',a10,
     5                 ' recoil energy (kev)".',/,
     6           'x axis style is duplex.',/,
     7           'x axis origin is 2.00.',/,
     8           'x axis tick marks 5.',/,
     8           'x tick marks mode reversed .',/,
     8           'y tick marks mode reversed .',/,
     9           'y axis style is duplex.',/,
     a           'x log.')
 4213       format ('x grid.',/,
     1           'y log .',/,
     1           'y axis height 0.20   .',/,
     1           'x axis height 0.20   .',/,
     1           'y axis number mode commas .',/,
     2           'y grid.',/,
     1           'x room 2.0 .',/,
     1           'y room 2.0 .',/,
     3           'curve 1 symbol count 0.',/,
     3           'curve 2 symbol count 0.',/,
     3           'curve 3 symbol count 0.',/,
     3           'title style is triplex.',/,
     4           'frame on.',/,
     2           'x integerize no. ',/,
     2           'y integerize no. ',/,
     2           'curve 1 texture solid.',/,
     2           'curve 2 texture dashed.',/,
     2           'curve 3 texture 7.',/,
     3           'curve 1 thickness is 5.',/,
     3           'curve 2 thickness is 5.',/,
     3           'curve 3 thickness is 5.')
 4214       format ('legend frame on. ',/,
     1           'legend x origin 2.5  .',/,
     2           'legend y origin 5.2  .',/,
     3           'curve 1 label "total cross section", 4 . ',/,
     3           'curve 2 label "lower kin. kerma", dashed. ',/,
     3           'curve 3 label "pka spectrum (kev)", solid . ')
 4312       format ('**file**.',/,
     1           'generate a fancy plot.',/,
     2           'generate a fancy x log plot.',/,
     3           'x axis label text is "n',
     3           'eutron energy (mev)".',/,
     4           'title text is "',
     5           'neutron damage functions".',/,
     5           'y axis label text is "',
     5                 ' damage kerma (mev-mb)".',/,
     6           'x axis style is duplex.',/,
     7           'x axis origin is 2.00.',/,
     8           'x axis tick marks 5.',/,
     8           'x tick marks mode reversed .',/,
     8           'y tick marks mode reversed .',/,
     9           'y axis style is duplex.',/,
     a           'x log.')
 4313       format ('x grid.',/,
     1           'y log .',/,
     1           'y axis height 0.20   .',/,
     1           'x axis height 0.20   .',/,
     1           'y axis number mode commas .',/,
     2           'y grid.',/,
     1           'x room 2.0 .',/,
     1           'y room 2.0 .',/,
     3           'curve 1 symbol count 0.',/,
     3           'curve 2 symbol count 0.',/,
     3           'curve 3 symbol count 0.',/,
     3           'title style is triplex.',/,
     4           'frame on.',/,
     2           'curve 1 texture solid.',/,
     2           'curve 2 texture dashed.',/,
     2           'curve 3 texture 7.',/,
     2           'x integerize no. ',/,
     2           'y integerize no. ',/,
     3           'curve 1 thickness is 5.',/,
     3           'curve 2 thickness is 5.',/,
     3           'curve 3 thickness is 5.')
 4314       format ('legend frame on. ',/,
     1           'legend x origin 2.5  .',/,
     2           'legend y origin 5.2  .',/,
     3           'curve 1 label "total kerma", 4 . ',/,
     3           'curve 2 label "lower kin. kerma", dashed. ',/,
     3           'curve 3 label "displacement kerma", solid . ')
c           output data arrays
c
            write (nplot, 8001) 
 8001       format (1x, 'input data.')
            if ( icon(2) .eq. 6) then 
                 itype = 1
            else
                 itype = 0
            endif
c
c           form recoil spectrum if pka option - icon(2)=5
c           and flagged by icon(11)=1
c
c           modified pka for new matxse format - but not lower kinematic limit
            if ( icon(2) .eq. 5 .and. icon(11) .eq. 1) then 
               do im = 1,nenergy
                   array(im,1) = array(im,6)/array(im,2)*1.e-6
               enddo
            endif
c
            do imat = 1,nmat
               write (nplot,3001) lab(imat)
3001           format (1x, '"', a20, ' "' )
               if ( itype .eq. 1 ) then
                  alow = array(1,imat)*0.1
                  ahigh = array(nenergy,imat)*0.1
                  write (nplot,3002) energy(1), alow, 
     1                (energy(ie), 
     1                array(ie,imat), energy(ie+1), 
     1                array(ie,imat), 
     1                ie=1,nenergy), energy(ie+1), 
     1                ahigh
               else
                  write (nplot,3002) (emid(ie), 
     1                array(ie,imat),
     1                ie=1,nenergy)
               endif
c               write (nplot,3002) (emid(ie), array(ie,imat),
c     1             ie=1,nenergy)
 3002           format (2x, g14.7, 2x, ' ', 2x, g14.7)
            enddo
c           output end-of-plot info
c
c
c        end plot description files
c
         write (nplot, 8014)
 8014    format ('end of data.',/,
     1           'go.')
            go to 101
       elseif ( icon(1) .eq. 4 .or. icon(1) .eq. 5 
     1     .or. icon(1) .eq. -5) then 
c
c           prepare multiple plot file
c
            read (nt5, *) number_of_files, scale, outfile
            itype = 1
            if ( number_of_files .gt. 15) then 
                write (6,9845) number_of_files
 9845           format (1x, 'number of files exceeds limit of 10 ',
     1           i10 )
                 stop 'number_of_files'
            endif
            do ifile = 1,number_of_files
c
c              input file name
c
               imat2 = 1
               if ( icon(1) .eq. 4) then
                  fraction = 1.0
                  read(nt5,*) file_name, imode, iselect_material
                  imat2 = iselect_material
               elseif ( icon(1) .eq. 5 .or. icon(1) .eq. -5) then
                  read(nt5,*) file_name, imode, iselect_material,
     1            fraction
               endif
               if ( iselect_material .le. 0) iselect_material = 1
               call prune (file_name, ilead, itrail, length, nchar)
               ierr = 0
               call filein(imode, file_name, ierr, iselect_material)
c
c              set flag for hisogram plot for number spectrum only
c
               if ( itype .eq. 1 .and. iselect_material .eq. 1 
     1            .and. imode .eq. 2) then 
                   itype = 1
               else
                   itype = 0
               endif
c
c              replace file_name with new label
c
               if ( icon(10) .eq. 1) then 
                  read(nt5,*) label
                  file_name = label
                  call prune (file_name, ilead, itrail, length, nchar)
               endif
               lab_file(ifile) = file_name
               lab_lead(ifile) = ilead
               lab_trail(ifile) = itrail
               lab_nchar(ifile) = nchar
               lab_len(ifile) = length
               nenergy_hld(ifile) = nenergy
c               write (6,18923) ifile, nenergy
c18923          format (1x, 'DEBUG: nenergy ', 2i5)
               if ( nenergy .le. 0) then 
                write (6, 3471) nenergy
3471            format (1x, 'ERROR: nenergy = 0, probable file exist',
     1                ' issue ', i5)
                stop 'nenergy'
               endif
               if ( ifile .eq. 1) then 
                   ematch_first = emid(1)
                   ematch_last = emid(nenergy)
                   ist = 1
               else
                   if ( ematch_first .eq. emid(1) .and. 
     1                  ematch_last .eq. emid(nenergy)) then 
                        ist = 1
c                   elseif ( ematch_first770 .eq. emid(1) .and. 
c     1                  ematch_last770 .eq. emid(nenergy)) then 
c                        ist = 1
                   elseif ( ematch_first .eq. emid(nenergy) .and. 
     1                  ematch_last .eq. emid(1)) then 
                        ist = 2
c                   elseif ( ematch_first770 .eq. emid(nenergy) .and. 
c     1                  ematch_last770 .eq. emid(1)) then 
c                        ist = 2
                   else
                         write (6,8923) nenergy, icon(5),
     2                    ematch_first, ematch_last, 
     1                    emid(1), emid(nenergy) 
c     2                    , ematch_first770, ematch_last770
 8923                    format (1x, '*** ERROR in energy grid ',
     1                   'match *** ',  2i5, 6g20.10)
                         stop 'egrid'
                   endif
               endif
               if ( icon(9) .lt. 0) then 
                 write (6,7853) ist, ifile, nenergy, imat2
 7853            format (1x, '*** MANIPULATE nenergy check ', 4i5)
               endif
               imat2_abs = iabs(imat2)
               if ( ist .eq. 1) then 
                  do ip = 1,nenergy
                     emid_hld(ip, ifile) = emid(ip)
                     array_hld(ip,ifile) = array(ip,imat2_abs)*
     1                  fraction
                  enddo
               elseif ( ist .eq. 2) then 
                  do ip = 1,nenergy
                     ipr = nenergy - ip + 1
                     emid_hld(ip, ifile) = emid(ipr)
                     array_hld(ip,ifile) = array(ipr,imat2_abs)*
     1                  fraction
                  enddo
               endif
c               write (6,7823) ist, ematch_first, ematch_last, emid(1), 
c     1         emid(nenergy), emid_hld(1,1), array_hld(1,1),
c     2         emid_hld(1,2), array_hld(1,2)
c 7823          format (1x, '*** debug *** ', i5, 14g14.7)
            enddo
            if (icon(1) .eq. -5) then
                 call stat
            endif
            if ( icon(1) .eq. 5) then 
c
c              make sum file - assume same energy grid
c
               isum = number_of_files + 1
               do ip=1,nenergy
                  emid_hld(ip,isum) = emid_hld(ip,1)
                  sum = 0.0
                  xmax_element = 0.0
                  do im=1,number_of_files
                    xmax_element = max(xmax_element, 
     &                   abs(array_hld(ip,im)))
                    sum = sum + array_hld(ip,im)
                  enddo
                  if ( icon(18) .eq. 1) then 
c 
c                   use precision control flag to eliminate trucated significance
c
                    if ( abs(sum)/xmax_element .lt. 1.e-4) sum = 0.0
c
                  endif
                  array_hld(ip,isum) = sum
                  if ( icon(9) .lt. 0) write (6,8911) ip, isum, sum, 
     &               (array_hld(ip,im), im=1,number_of_files), 
     &               array_hld(ip,isum)
 8911             format (1x, 'sum computation: ', 2i5, 18g14.7)
               enddo
c               write (6,8823) isum,  emid_hld(1,1), 
c     1          emid_hld(1,2), emid_hld(1,isum),
c     1          array_hld(1,1), 
c     1          array_hld(1,2), array_hld(1,isum)
c 8823          format (1x, '*** debxx *** ', i5, 14g14.7)
               number_of_files = number_of_files + 1
               file_name = 'cumulative sum'
               call prune (file_name, ilead, itrail, length, nchar)
c               lab_file(isum) = file_name
c               lab_lead(isum) = ilead
c               lab_trail(isum) = itrail
c               lab_nchar(isum) = nchar
c               lab_len(isum) = length
               nenergy_hld(isum) = nenergy            
c
cxxx
c
               if ( icon(3) .eq. 2) then 
c
c                 punch data in tdam format
                  scale = 1.0
                  jkl = 1
                  nenergy = nenergy_hld(isum)
                  do jf=1,nenergy
                      emid(jf) = emid_hld(jf,isum)*1.e+6
                      array(jf,jkl) = array_hld(jf,isum)*1.e+3
                  enddo
                  npun2 = 31
                  dummy = 0.0
                  if ( jkl .eq. 1) then 
                     file_name2 = 'ker.dat'
                   elseif ( jkl .eq. 2) then 
                     file_name2 = 'lmt.dat'
                  elseif ( jkl .eq. 3) then 
                     file_name2 = 'dsp.dat'
                  endif
                  open(unit=npun2, form='formatted',
     1                file=file_name2
     2               , status='new', iostat=ilook, err=909)
                  write (npun2, 2745) lab(jkl), title(1:60), 
     1              nenergy, dummy
                  if ( icon(19) .eq. 0) then 
                     write (npun2,2741) ( emid(jk)*1.e-6,  
     1                array(jk,jkl)*scale*1.e-3,
     1                 jk=nenergy,1,-1)
                  else
c                    reverse order
                     write (npun2,2741) ( emid(jk)*1.e-6,  
     1                array(jk,jkl)*scale*1.e-3,
     1                 jk=1,nenergy)
                  endif
                  close (unit=npun2)
c
c              plot files
c
c                write (nplot,4312)
c               write (nplot, 4313) 
c               write (nplot, 4314)
c               write (nplot, 8001) 
c               lab(1) = 'tdam total kerma     '
c               lab(2) = 'tdam kinematic limit   '
c               lab(3) = 'tdam displacement damage '
c                  write (nplot,3001) lab(jkl)
                open (unit=npun2, form='formatted',
     1             file='component_'//char(48+jkl), status='new', 
     2             iostat=ilook, err=909)
                  if ( icon(19) .eq. 0) then 
                     write (npun2,2741) ( emid(jk)*1.e-6,  
     1                array(jk,jkl)*scale,
     1                 jk=nenergy,1,-1)
                  else
                     write (npun2,2741) ( emid(jk)*1.e-6,  
     1                array(jk,jkl)*scale,
     1                 jk=1,nenergy)
                  endif
c               write (nplot, 8014)
               elseif ( icon(3) .eq. 1) then 
c
c              special punch format without zeros
c              modify to allow for zeros being turned into small numbers
c
               scale = 1.0
               ilow = 1
               ihigh = 641
               if ( icon(5) .eq. 2) ihigh = 90
c              special look loop
                   jkp = 0
 3945              continue
                   jkp = jkp + 1
                   write (6,4512) jkp,isum, array_hld(jkp,isum)
 4512              format (1x, 'test isum low: ', 2i5, g14.7)
                   if ( array_hld(jkp,isum) .le. 1.e-32) then 
                       ilow = jkp
                       if ( ilow < ihigh) go to 3945
                   else
                       continue
                   endif
                   jkp = 641
                   if ( icon(5) .eq. 2) jkp = 90
 3946              continue
                   jkp = jkp - 1
                   write (6,5512) jkp,isum, array_hld(jkp,isum)
 5512              format (1x, 'test isum high: ', 2i5, g14.7)
                   if ( array_hld(jkp,isum) .le. 1.e-32) then 
                       ihigh = jkp
                       if ( ihigh > 1) go to 3946
                   else
                       continue
                   endif
               write (27, 5745) title(1:60)
c 5745          format (1x, a60)
               ilow_x = ilow
               ihigh_x = ihigh - 1
               ihigh_x1 = ilow_x 
               ilow_x1  = ihigh_x
               write (27,6190)  outfile
c 6190          format (1x, 'file = ', a)
               write (27,6923)  ihigh_x1, ilow_x1
c 6923          format (1x, i5, 1x, i5)
               write (27,2341) (( 
     1          array_hld(jk,ik)*scale,
     1          ik=isum,isum),  jk=ilow_x,ihigh_x)
               endif
c
cxxx
c
            endif
c
c           output plot header info
c
c
c        write plot direction file header data -
c
c         write (nplot,5012) 
c         write (nplot, 5013) 
c 5012    format ('**file**.',/,
c     1           'generate a fancy plot.',/,
c     2           'generate a fancy x log plot.',/,
c     3           'x axis label text is "neutron energy (mev)".',/,
c     4           'title text is "',
c     5           'composite plot".',/,
c     5           'y axis label text is "',
c     5                 ' material damage (mev-mb)".',/,
c     6           'x axis style is duplex.',/,
c     7           'x axis origin is 2.00.',/,
c     8           'x axis tick marks 5.',/,
c     8           'x tick marks mode reversed .',/,
c     8           'y tick marks mode reversed .',/,
c     9           'y axis style is duplex.',/,
c     a           'x log.')
c 5013    format ('x grid.',/,
c     1           'y log .',/,
c     1           'y axis height 0.20   .',/,
c     1           'x axis height 0.20   .',/,
c     1           'y axis number mode commas .',/,
c     2           'y grid.',/,
c     4           'frame on.',/,
c     3           'title style is triplex.')
c         do ij=1,number_of_files
c            ibeg = lab_lead(ij) + 1
c            iend = lab_nchar(ij)-lab_trail(ij)
c            lenout = lab_len(ij)
c            write (nplot,7013) ij, ij, ij, icurve(ij),  
c     2      ij, lab_file(ij)(ibeg:iend),
c     3      icurve(ij) 
c         enddo
c 7013    format (
c     3           'curve ', i2, ' symbol count 0.',/,
c     3           'curve ', i2, ' thickness is 5.',/,
c     4           'curve ', i2, ' texture is ', i2, ' . ',/,
c     3           'curve ', i2, ' label "', a<lenout>,
c     3           '", ', i2, '  . ')
c         write (nplot, 5014)
c 5014    format ('legend frame on. ',/,
c     1           'legend x origin 7.0  .',/,
c     2           'legend y origin 5.2  .')
c
c           output data arrays
c
c            write (nplot, 8001) 
            do imat = 1,number_of_files
             if ( imat .lt. number_of_files) then 
                putfile='plot.'//char(48+imat)
                open(unit=nplot,file=putfile
     1         ,status='unknown')
             else
                number_file='sum'
                putfile='plot.'//number_file
c                idir = '/app/manipulate-2/response/'
                idir = 'response/'
                iblank2 = lnblnk(idir)
                if ( icon(9) < 0) then
                   write (*, 5634)  xoptical(1:leopt), idir(1:iblank2),
     &                             outfile
 5634              format (1x, 'sum file name: xoptical = ', a,/,
     &                     1x, '               idir     = ', a,/,
     &                     1x, '               outfile  = ', a)
                endif
                open(unit=nplot,
     1          file=xoptical(1:leopt)//idir(1:iblank2)//outfile
     1         ,status='unknown')
             endif
                efactor = 1.0
                if ( energy(1) .gt. 1.e6 .or.
     1               energy(nenergy) .gt. 1.e6) then 
                     efactor = 1.e-6
                endif
c               length = lab_len(imat)
c               ibeg = lab_lead(imat) + 1
c               iend = lab_nchar(imat)-lab_trail(imat)
c               write (nplot,3101) 
c     1            lab_file(imat)(ibeg:iend)
c3101           format (1x, '"', a<length>, ' "' )
               do jkl=1,nenergy_hld(imat)
                  if ( array_hld(jkl,imat)*scale .le. 0.0) then 
                        array_hld(jkl,imat) = 1.e-33
                   endif
               enddo
               if ( itype .eq. 1 ) then
                  alow = array_hld(1,imat)*0.1*scale
                  ahigh = array_hld(nenergy,imat)*0.1*scale
                  write (nplot,3002) energy(1)*efactor, alow, 
     1                (energy(ie)*efactor, 
     1                array_hld(ie,imat)*scale, energy(ie+1)*efactor, 
     1                array_hld(ie,imat)*scale, 
     1                ie=1,nenergy_hld(imat)), energy(ie+1)*efactor, 
     1                ahigh
               else
                   if ( icon(4) .eq. 1) then 
                   write (nplot,3002) (emid_hld(ie,imat)*efactor, 
     1                array_hld(ie,imat)*scale,
     1                ie=1,nenergy_hld(imat))
                   else
                   write (nplot,3002) (emid_hld(ie,imat)*efactor, 
     1                array_hld(ie,imat)*scale,
     1                ie=nenergy_hld(imat),1,-1)
                   endif
              endif
               close (unit=nplot)
            enddo
c           output end-of-plot info
c
c
c        end plot description files
c
c         write (nplot, 8014)
c
c       Special punch section to support lsl format output
c
         go to 101
       elseif ( icon(1) .eq. 6) then
c           
c           read in source files
c
            fraction = 1.0
            read(nt5,*) file_name, imode, iselect_material
            ireverse= iselect_material
            if ( iselect_material.le.0 ) iselect_material = 1
            call prune (file_name, ilead, itrail, length, nchar)
            ierr = 0
            iflast = lnblnk(file_name)
            kdir = 'spectrum/'
            kblank2 = lnblnk(kdir)
c            write (6,4598) ilead, itrail, length, nchar, 
c     &           file_name(1:iflast)
 4598       format (1x, 'Option 6 filein call: ', 4i5,/,
     &              1x, '                     =',a,'=')
            new_file_name = xoptical(1:leopt)//kdir(1:kblank2)//
     1       file_name(1:iflast)
            call filein(imode, new_file_name, ierr, ireverse)
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
c           read in response files
c
            fraction = 1.0
            read(nt5,*) file_name, imode, iselect_material
            ireverse = iselect_material
            if ( iselect_material .le. 0) iselect_material = 1
            call prune (file_name, ilead, itrail, length, nchar)
            ierr = 0
            iflast = lnblnk(file_name)
            new_file_name = xoptical(1:leopt)//jdir(1:jblank2)//
     1       file_name(1:iflast)
c           must inhibit response renormalization here
            call filein_response(imode, new_file_name, ierr, ireverse)
c
c           save response data
c
            isum = 2
            lab_file(isum) = file_name
            lab_lead(isum) = ilead
            lab_trail(isum) = itrail
            lab_nchar(isum) = nchar
            lab_len(isum) = length
            nenergy_hld(isum) = nenergy
            dif1 = abs( (emid(1) - emid_hld(1,1)) /emid(1))
            dif2 = abs( (emid(1) - emid_hld(nenergy_hld(1),1))
     1                  /emid(1))
            if ( dif1 .lt. 1.e-3) then
               do ip=1,nenergy
                  emid_hld(ip,isum) = emid(ip)
                  array_hld(ip,isum) = array(ip,iselect_material)
               enddo
            elseif ( dif2 .lt. 1.e-3) then
               do ip=1, nenergy
                  irev = nenergy - ip + 1
                  emid_hld(ip,isum) = emid(irev)
                  array_hld(ip,isum) = array(irev,iselect_material)
               enddo
            else
                write (nt6,8561) emid(1), emid_hld(1,1),
     1          emid(nenergy), emid_hld(nenergy_hld(1),1)
8561            format (1x, '*** first/last energy grid error before',
     &          ' fold *** ',/,
     1          5x, 4g14.7)
                stop 'energy-grid'
            endif
c
c           apply attenuation factor if requested
c
           if ( icon(13) .eq. 1) then 
              read(nt5,*) alpha
              fraction = 1.0
              read(nt5,*) file_name, imode, iselect_material
              if ( iselect_material .le. 0) iselect_material = 1
              call prune (file_name, ilead, itrail, length, nchar)
              ierr = 0
              call filein(imode, file_name, ierr, iselect_material)
c
c             save attenuation function
c
              isum = 3
              lab_file(isum) = file_name
              lab_lead(isum) = ilead
              lab_trail(isum) = itrail
              lab_nchar(isum) = nchar
              lab_len(isum) = length
              nenergy_hld(isum) = nenergy
              if ( emid(1) .eq. emid_hld(1,2)) then
                 do ip=1,nenergy
                    emid_hld(ip,isum) = emid(ip)
                    array_hld(ip,isum) = array(ip,iselect_material)
                 enddo
              elseif ( emid(1) .eq. emid_hld(nenergy_hld(1),2)) then
                 do ip=1, nenergy
                    irev = nenergy - ip + 1
                    emid_hld(ip,isum) = emid(irev)
                    array_hld(ip,isum) = array(irev,iselect_material)
                 enddo
              else
                  write (nt6,8261) emid(1), emid_hld(1,2),
     1            emid_hld(nenergy_hld(1),2)
8261              format (1x, '*** energy grid error in attenua',
     1            'tion *** ',/,
     1            5x, 3g14.7)
                  stop 'energy-grid'
              endif
c
c             apply factor to response
c
              do ip=1,nenergy
                  fact2 =-alpha*array_hld(ip,3)
c                  if ( fact2 .lt. 75) then
c                     act = 0.0
c                  else
                     fact = exp(fact2)
c                  endif
                  array_hld(ip,2) = array_hld(ip,2)*fact
              enddo
           endif
c
c           apply a self-shielding correction factor if requested
c
           if ( icon(13) .eq. 2) then 
               read(nt5,*) file_name
               file_corr = file_name
               if ( icon(9) < 0) then 
                   write (6,8712)
 8712              format (1x, 'MANIPULATE correct call')
               endif
               call correct(array_hld(1,2), file_name)
           endif
c              fraction = 1.0
c               iselect_material = 1
cc              read(nt5,*) file_name, imode, iselect_material
cc              if ( iselect_material .le. 0) iselect_material = 1
cc              call prune (file_name, ilead, itrail, length, nchar)
cc              ierr = 0
cc              call filein(imode, file_name, ierr, iselect_material)
cc
cc
cc             special self-shielding/cover MCNP correction factor
cc             read in low to high in 640 bin structure
cc
c              isum = 3
cc              label = self(i)
cc              write(name,310) reac(i,1),reac(i,2)
cc              call full_name(name,label,ireaction,lcov)
cc             ename = 'mcnploc'
cc             lename = lnblnk(ename)
cc             call getenv(ename(1:lename), mcnploc)
c              mcnploc = '/bx4/app/sandii/correction-interface'
c              lnloc = lnblnk(mcnploc)
c              iendx = lnblnk(file_name)
cc              leself = lnblnk(self(i))
cc              iendy = lnblnk(lcov_pjg)
c              call prune(file_name, ilead_1, itrail_1, 
c     1           length_1, nchar_1) 
cc              call prune(lcov_pjg, ilead_2, itrail_2, length_2, 
cc     1              nchar_2) 
c              file_id = mcnploc(1:lnloc)//'/'
c     1                 //file_name(ilead_1+1:nchar_1-itrail_1)//
c     1                 '.summary'
cc     1            //lcov_pjg(ilead_2+1:nchar_2-itrail_2)
cc     2            //'-'
cc     2            //self(i)(1:leself)//'.summary'
cc              write (6,8923)  file_id
c 8923         format (1x,'apply special self-shielding/cover correction',
c     1         /, 15x, ' file = ', a100)
                 jstatus = 1
c                open (unit=29, 
c     1              file=file_id, status='old', err=9099)
c                read (29, 2871) header
c 2871           format (a100)
c                read (29,*) (array_hld(ikl,isum), ikl=1,640)
cc                write (6,3321) (array_hld(ikl,isum),ikl=1,640)
c 3321           format ('ss-fct   ', 5(1pe13.6,1x),/, 
c     1          (   '         ', 5(1pe13.6,1x) )    )
c                close (unit=29)
cc
cc             save self-shielding correction function
cc
c              lab_file(isum) = file_name
c              lab_lead(isum) = ilead
c              lab_trail(isum) = itrail
c              lab_nchar(isum) = nchar
c              lab_len(isum) = length
c              nenergy_hld(isum) = nenergy
cc              if ( emid(1) .eq. emid_hld(1,2)) then
cc                 do ip=1,nenergy
cc                    emid_hld(ip,isum) = emid(ip)
cc                    array_hld(ip,isum) = array(ip,iselect_material)
cc                 enddo
cc              elseif ( emid(1) .eq. emid_hld(nenergy_hld(1),2)) then
cc                 do ip=1, nenergy
cc                    irev = nenergy - ip + 1
cc                    emid_hld(ip,isum) = emid(irev)
cc                    array_hld(ip,isum) = array(irev,iselect_material)
cc                 enddo
cc              else
cc                  write (nt6,5261) emid(1), emid_hld(1,2),
cc     1            emid_hld(nenergy_hld(1),2)
cc5261              format (1x, '*** energy grid error in attenua',
cc     1            'tion *** ',/,
cc     1            5x, 3g14.7)
cc                  stop 'energy-grid'
cc              endif
cc
cc             apply self-shielding correction factor to response
cc
c              do ip=1,nenergy
c                  fact2 =array_hld(ip,3)
c                  array_hld(ip,2) = array_hld(ip,2)*fact2
c              enddo
c           endif
cc
c           fold together data
c
            call fold
            go to 101
       elseif ( iabs(icon(1)) .eq. 7) then 
            number_of_files = 2
            itype = 1
            read (nt5,*) itytpe, scale_x, outfile
            do ifile = 1,number_of_files
c
c              input file name
c
c               fraction = 1.0
               read(nt5,*) file_name, imode, iselect_material, fraction
               if ( iselect_material .le. 0) iselect_material = 1
               call prune (file_name, ilead, itrail, length, nchar)
               ierr = 0
               iflast = lnblnk(file_name)
               new_file_name = xoptical(1:leopt)//jdir(1:jblank2)//
     1         file_name(1:iflast)
               call filein(imode, new_file_name, ierr, iselect_material)
c               call filein(imode, file_name, ierr, iselect_material)
c
c              set flag for histogram plot for number spectrum only
c
c               if ( itype .eq. 1 .and. iselect_material .eq. 1 
c     1            .and. imode .eq. 2) then 
c                   itype = 1
c               else
c                   itype = 0
c               endif
c
c              replace file_name with new label
c
c               if ( icon(10) .eq. 1) then 
c                  read(nt5,*) label
c                  file_name = label
c                  call prune (file_name, ilead, itrail, length, nchar)
c               endif
c               lab_file(ifile) = file_name
c               lab_lead(ifile) = ilead
c               lab_trail(ifile) = itrail
c               lab_nchar(ifile) = nchar
c               lab_len(ifile) = length
               nenergy_hld(ifile) = nenergy
               do ip = 1,nenergy
                  emid_hld(ip, ifile) = emid(ip)
                  array_hld(ip,ifile) = array(ip,1)*
     1               fraction
               enddo
            enddo
c
c           make weighted file - assume same energy grid
c
            isum = 3
            do ip=1,nenergy
               emid_hld(ip,isum) = emid(ip)
               if ( icon(1) .eq. 7) then 
                  sum = array_hld(ip,1)*array_hld(ip,2)
               elseif ( icon(1) .eq. -7) then
                  if ( array_hld(ip,2) .ne. 0.0) then 
                      sum = array_hld(ip,1)/array_hld(ip,2)
                  else
                       sum = -9.999e-30
                  endif
               endif
               array_hld(ip,isum) = sum
            enddo
            number_of_files = number_of_files + 1
c            file_name = 'weighted curve'
c            call prune (file_name, ilead, itrail, length, nchar)
c            lab_file(isum) = file_name
c            lab_lead(isum) = ilead
c            lab_trail(isum) = itrail
c            lab_nchar(isum) = nchar
c            lab_len(isum) = length
            nenergy_hld(isum) = nenergy            
c
c           output plot header info
c
c
c        write plot direction file header data -
c
c         write (nplot,5012) 
c         write (nplot, 5013) 
c         do ij=1,number_of_files
c            ibeg = lab_lead(ij) + 1
c            iend = lab_nchar(ij)-lab_trail(ij)
c            lenout = lab_len(ij)
c            write (nplot,7013) ij, ij, ij, icurve(ij),  
c     2      ij, lab_file(ij)(ibeg:iend),
c     3      icurve(ij) 
c         enddo
c         write (nplot, 5014)
c
c           output data arrays
c
c            write (nplot, 8001) 
c            do imat = 1,number_of_files
c               length = lab_len(imat)
c               ibeg = lab_lead(imat) + 1
c               iend = lab_nchar(imat)-lab_trail(imat)
c               write (nplot,3101) 
c     1            lab_file(imat)(ibeg:iend)
c               if ( itype .eq. 1 ) then
c                  alow = array_hld(1,imat)*0.1
c                  ahigh = array_hld(nenergy,imat)*0.1
c                  write (nplot,3002) energy(1), alow, 
c     1                (energy(ie), 
c     1                array_hld(ie,imat), energy(ie+1), 
c     1                array_hld(ie,imat), 
c     1                ie=1,nenergy_hld(imat)), energy(ie+1), 
c     1                ahigh
c               else
c                  write (nplot,3002) (emid_hld(ie,imat), 
c     1                array_hld(ie,imat),
c     1                ie=1,nenergy_hld(imat))
c               endif
c            enddo
c           output end-of-plot info
c
c
c        end plot description files
c
c         write (nplot, 8014)
c         if ( icon(11) .eq. 2) then 
c
c        punch interface plt-format file
c
c         read(nt5,*) file_name
         iblank = lnblnk(outfile)
         efactor = 1.0
         if ( emid(1) .gt. 1.e6 .or.
     1        emid(nenergy) .gt. 1.e6) then 
              efactor = 1.e-6
         endif
         open(unit=8,file=xoptical(1:leopt)//jdir(1:jblank2)//outfile
     1   ,status='unknown')
         if ( icon(12) .eq. 1) then 
            do i = nenergy,1,-1
              write(8,2341)emid(i)*efactor,array_hld(i,3)*scale_x
            end do
         else
            do i = 1,nenergy
              write(8,2341)emid(i)*efactor,array_hld(i,3)*scale_x
            end do
         endif
         close(unit=8)
c         open(unit=17, form='formatted',
c     1         file=file_name
c     2         , status='old', iostat= ilook, err=909)
c         write(17,8923) number_of_files 
c         if = number_of_file
c 8923    format (/, i5)
c         write (nplot,3002) (emid_hld(ie,if), 
c     1                array_hld(ie,if),
c     1                ie=1,nenergy_hld(if))
c         endif
         go to 101
       elseif ( icon(1) .eq. 8) then 
c
c           file difference
c
c           input file name
c
            ifile = 1
            fraction = 1.0
            read(nt5,*) file_name, imode, iselect_material
            call prune (file_name, ilead, itrail, length, nchar)
            ierr = 0
            call filein(imode, file_name, ierr, iselect_material)
c
c           save file
c
               nenergy_hld(ifile) = nenergy
               do ip = 1,nenergy
                  emid_hld(ip, ifile) = emid(ip)
                  array_hld(ip,ifile) = array(ip, iselect_material)*
     1               fraction
               enddo
            ifile = 2
            read(nt5,*) file_name2, imode, iselect_material
            call prune (file_name2, ilead, itrail, length, nchar)
            ierr = 0
            call filein(imode, file_name2, ierr, iselect_material)
c
c           save file
c
               nenergy_hld(ifile) = nenergy
               do ip = 1,nenergy
                  emid_hld(ip, ifile) = emid(ip)
                  array_hld(ip,ifile) = array(ip, iselect_material)*
     1               fraction
               enddo
c
c          difference
c
          if ( nenergy_hld(1) .ne. nenergy_hld(2)) then 
             write (nt6, 7601) nenergy_hld(1), 
     1                         nenergy_hld(2)
 7601        format (1x, '*** energy difference option requires ',
     1       'same number of points ', 2i8)
              stop 'points'
          endif
          write (nt6, 7602) file_name, file_name2
 7602     format (1h1, /,/,4x, 'file difference ',/,
     1     15x, 'file 1 = ', a50,/,
     2     15x, 'file 2 = ', a50,/,/,/)
          permax = 0.0
          do jk=1,nenergy_hld(1)
            res = array_hld(jk,1)
            if ( abs(res) .le. 1.e-33) then 
                 res = array_hld(jk,2)
            endif
            if ( abs(res) .le. 1.e-33) then 
                 res = 1.0
            endif
            res = res/100.
            percent(jk) = (array_hld(jk,1) - array_hld(jk,2))/res
            if ( abs(percent(jk)) .gt. permax) then 
                  permax = abs(percent(jk))
            endif
          enddo
          write (nt6, 7603) (emid_hld(jk,1), 
     1       array_hld(jk,1), array_hld(jk,2),
     3       percent(jk), 
     2       jk = 1,nenergy_hld(1)) 
 7603     format (5x, g14.7, 10x, g14.7, 5x, g14.7, 10x, g14.7)
          write (nt6, 7604) permax
 7604     format (/,/,/, 5x, 'max percent diff = ', g14.7)
       elseif (icon(1) .eq. 10 .or. icon(1) .eq. 11) then 
c           read matxsr file and report  - option 10
c           interface same logic with a plt_89 file - option -3
            noption = 1
            read(nt5,*) file_name
            call prune (file_name, ilead, itrail, length, nchar)
            if ( icon(1) .eq. 10) then 
               open(unit=nfile, form='unformatted',
     1              file=file_name, status='old', err=909)
cje            call pmatxs(nfile,nt6)
            elseif (icon(1) .eq. 11) then
               write (nt6,739) file_name
 739           format (1x, 'icon(1)=11 file_name check: ', a)
               open(unit=nfile, form='formatted',
     1              file=file_name, status='old', err=909)
               imaterial = 1
               read (nfile,*) nenergy
               do jkl=1,nenergy
                 do jkp = 1,10
                    array(jkl,jkp) = 0.0
                 enddo
               enddo
               if ( icon(12) .eq. 0) then
                   read (nfile,*) 
     1                (emid(jkl), array(jkl,1), jkl=1,nenergy) 
               elseif ( icon(12) .eq. 1) then
                   read (nfile,*) 
     1                (emid(jkl), array(jkl,1), jkl=nenergy,1,-1) 
               else
                   stop 'icon(12) err'
               endif
            endif
            if ( icon(3) .eq. 1) then 
c
c              punch data in damout format
               close (unit=nfile)
               read (nt5,*) scale
               if ( icon(1) .eq. 10) then 
                  open(unit=nfile, form='unformatted',
     1                file=file_name
     2                 , status='old', iostat=ilook, err=909)
               elseif ( icon(1) .eq. 11) then 
                   continue
c                  open(unit=nfile, form='formatted',
c     1                 file=file_name, status='old', err=909)
               endif
               if ( icon(1) .eq. 10) then 
                  do jkl=1,nenergy
                    do jkp = 1,10
                       array(jkl,jkp) = 0.0
                    enddo
                  enddo
cje               call matxsin(nfile,nt6)
                  do jkl=1,nenergy
                      emid(jkl) = 0.5*(energy(jkl) + energy(jkl+1))
                  enddo
               elseif ( icon(1) .eq. 11) then 
                  continue
               endif
            elseif ( icon(3) .eq. 2) then 
               close (unit=nfile)
               read (nt5,*) scale
               open(unit=nfile, form='unformatted',
     1             file=file_name
     2              , status='old', iostat=ilook, err=909)
cje            call matxsin(nfile,nt6)
               do jkl=1,nenergy
                   emid(jkl) = 0.5*(energy(jkl) + energy(jkl+1))
               enddo
            endif
c
c           interpolate on function
c              idirection = 1, interpolate on energy to get y
c              idirection = 2, interpolate on y to get energy
c
c              method = 1, y constant
c                       2, x linear in y
c                       3, x linear in ln(y)
c                       4, ln(x) linear in y
c                       5, ln(x) linear in ln(y)
c
            read (5,*) method, idirection, ipoint, value
            if ( ipoint .lt. 0) then 
                  ipoint = 1
            elseif ( ipoint .gt. 10) then 
                  ipoint = 1
            endif
            do i=1,nenergy
                range(i) = array(i,ipoint)
            enddo
            if ( idirection .eq. 1) then 
                 result = fitmd(value, nenergy, emid(1), 
     1                     array(1,ipoint), method)
            elseif ( idirection .eq. 2) then 
                 result = fitmd(value, nenergy, array(1,ipoint),
     1                     emid(1), method)
            endif
            write (6,9128) file_name, 
     1             method, idirection, value, result*scale
 9128       format (1x, 'file_name = ', a,/,
     1              1x, 'interpolation method    = ', i5,/,
     1              1x, '              direction = ', i5,/,
     2              1x, '              value     = ', 1pg14.7,/,
     3              1x, '              result    = ', 1pg14.7,/)
            go to 101
       elseif (icon(1) .eq. 12 .or. icon(1) .eq. 13 ) then 
            noption = 1
            if ( icon(5) .eq. 0) then
                nenergy = 640
            elseif ( icon(5) .eq. 1) then 
                nenergy = 770
            elseif ( icon(5) .eq. 1) then 
                   nenergy = 89
            elseif ( icon(5) .eq. 5) then 
                   nenergy = 175
            elseif ( icon(5) .eq. 6) then 
                   nenergy = 725
            else
                write (6,7245) nenergy, icon(5)
                stop 'icon(5)-2'
            endif
            idirection = 0
            if ( icon(1) .eq. 13) go to 3821
            read(nt5,*) file_name
            call prune (file_name, ilead, itrail, length, nchar)
               nlen = lnblnk(file_name)
               name = xoptical(1:leopt)//
     2          jdir(1:jblank2)//file_name(1:nlen)
               write (nt6,7391) xoptical(1:leopt), jdir(1:jblank2), 
     &            file_name(1:nlen), name
7391           format (1x, 'icon(1)=12 file_name check: ',/,
     &         4x, 'xoptical   ', a,/,
     &         4x, 'jdir       ', a,/,
     &         4x, 'file_name  ', a,/,
     &         4x, 'name       ', a,/)
               open(unit=nfile, form='formatted',
     1              file=name, status='old', err=909)
               imaterial = 1
               read (nfile,*) nenergy
 3821          continue
               do jkl=1,nenergy
                 do jkp = 1,10
                    array(jkl,jkp) = 0.0
                 enddo
               enddo
               if ( icon(1) .eq. 13) go to 3822
               if ( icon(12) .eq. 0) then
                   read (nfile,*) 
     1                (emid(jkl), array(jkl,1), jkl=1,nenergy) 
               elseif ( icon(12) .eq. 1) then
                   read (nfile,*) 
     1                (emid(jkl), array(jkl,1), jkl=nenergy,1,-1) 
               else
                   stop 'icon(12) err'
               endif
            if ( icon(3) .eq. 1) then 
c
c              punch data in damout format
               close (unit=nfile)
               read (nt5,*) scale
               if ( icon(1) .eq. 10) then 
                  open(unit=nfile, form='unformatted',
     1                file=file_name
     2                 , status='old', iostat=ilook, err=909)
               elseif ( icon(1) .eq. 11) then 
                   continue
c                  open(unit=nfile, form='formatted',
c     1                 file=file_name, status='old', err=909)
               endif
            elseif ( icon(3) .eq. 2) then 
               close (unit=nfile)
               read (nt5,*) scale
               open(unit=nfile, form='unformatted',
     1             file=file_name
     2              , status='old', iostat=ilook, err=909)
cje            call matxsin(nfile,nt6)
               do jkl=1,nenergy
                   emid(jkl) = 0.5*(energy(jkl) + energy(jkl+1))
               enddo
            endif
c
c           interpolate on function
c              idirection = 1, interpolate on energy to get y
c              idirection = 2, interpolate on y to get energy
c
c              method = 1, y constant
c                       2, x linear in y
c                       3, x linear in ln(y)
c                       4, ln(x) linear in y
c                       5, ln(x) linear in ln(y)
c
            read (5,*) method, idirection
 3822       continue
            write (6,8733)
c           fill default names for disk retrieves
c
            jdir = '\sync_sandialabs\MANIPULATE-2020\response\'
            jblank2 = lnblnk(jdir)
            kdir = '\sync_sandialabs\MANIPULATE-2020\spectrum\'
            kblank2 = lnblnk(kdir)
            ldir = '\sync_sandialabs\MANIPULATE-2020\'
            lblank2 = lnblnk(ldir)
cdebug            write (*,9012) xoptical(1:leopt)//jdir(1:jblank2)
9012        format (1x, 'opening directory location: ', a41)
 8733       format (1x, '640 grid interface function ',/,/)
            if (icon(5) .eq. 0) then
               write (nt6,7491) xoptical(1:leopt), jdir(1:jblank2)
7491           format (1x, 'grid interface file_name check: ',/,
     &         4x, 'xoptical   ', a,/,
     &         4x, 'jdir       ', a,/)
               open(unit=nfile, form='formatted',
     1           file=jdir(1:jblank2)//'sand641.nrg'
     2           , status='old', iostat=ilook, err=909)
               ilmt = 641
            elseif ( icon(5) .eq. 1) then
               open(unit=nfile, form='formatted',
     1           file=xoptical(1:leopt)//jdir(1:jblank2)//'sand771.nrg'
     2           , status='old', iostat=ilook, err=909)
               ilmt = 771
            endif
            read(nfile, 781) title
 781        format (a80)
            read (nfile,*) ienergy
c           if ( icon(12) .ne. 1) then 
c              read (nfile,*) (energy(jk), jk=ilmt,1,-1)
c           else
c              read (nfile,*) (energy(jk), jk=1,ilmt)
c           endif
            read (nfile,*) (energy(jk), jk=ilmt,1,-1)
            close (unit=nfile)
            do jkl=1,ilmt
               emid_hld(jkl,1) = 0.5*(energy(jkl) + energy(jkl+1))
            enddo
c
c
            open (unit=nfile, status='unknown', 
     1             form='formatted',
     1             file='file.interpolation' )
c            write (nfile,2128) file_name, 
c     1             method, idirection
             write (nfile, 3453)
 3453        format (5x, '640')
c
            do jkl=1,ilmt
               value = emid_hld(jkl,1)
               x = value
               if ( idirection .eq. 1) then 
                    result = fitmd(value, nenergy, emid(1), 
     1                        array(1,1), method)
               elseif ( idirection .eq. 2) then 
                    result = fitmd(value, nenergy, array(1,1),
     1                        emid(1), method)
               endif
               if ( icon(1) .eq. 13) then 
c
c               user supplied functional representation
c
c               option 1 - fit to gaas neutron damage efficiency
c
                 a0 = 0.872670
                 a1 = -0.187469
                 a2 = 1.237178e-7
                 a3 = -0.060753
c                energy in kev for equation
                 x = value*1.e3
                 result = a0 + a1*log10(x) + a2*(x**2)*log10(x) + 
     1                    a3*(log10(x)**2)
                 if ( x .le. 1.e-1) then 
                     result = 1.0
                 endif
                 if ( x .ge. 500.0) then 
                     result = 0.01
                 endif
               endif
               write (6,9629) x, result
               write (nfile,9669) x, result
 9629          format (2x, 1pg14.7, '   ,   ', 1pg14.7)
 9669          format (2x, 1pg14.7, '      ', 1pg14.7)
            enddo
            close (unit=nfile)
            write (6,2128) file_name, 
     1             method, idirection
 2128       format (1x, 'file_name = ', a50,/,
     1              1x, 'interpolation method    = ', i5,/,
     1              1x, '              direction = ', i5,/)
            go to 101
       elseif ( icon(1) .eq. 9) then
            close (unit=39)
            close (unit=nt6)
            stop 'normal'
       else
            write (nt6,9991) icon(1)
 9991       format (1x, 'icon(1)=',i2,
     1      ' option not implemented ')
            stop 'not-implemented'
       endif
 999   continue
       stop 'normal'
 909   continue
       length = len(file_name)
       write (nt6,910) file_name(1:length)
910    format (1x, '*** error opening file ', a50)
       go to 101
 9099  continue
       write (nt6, 6723) jstatus, file_id
 6723  format (1x, 'Error in MANIPULATE file open ', i5, 3x, a90)
       stop 'Open'
       end
