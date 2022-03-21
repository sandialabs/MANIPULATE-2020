       subroutine spectra_format

       common /datain/ nenergy, energy(1001), array(1000,15)

     1         , emid(1001)

       common /datahld/ nenergy_hld(15), energy_hld(1001,15),

     1          array_hld(1000,15)

     1         , emid_hld(1001,15)

       character*80 comment

       common /labx/ icomment, comment(20)

       common /guide/ icon(40)

       common /io/ nt5, nt6, nfile, nplot



       write (nt6, 6672)

 6672  format (1x, '*** spectra-format entererd ***')

c

c

c         convert differential number spectrum into number fraction

c

          sum = 0.0

          do jk=1,nenergy

             array(jk,1) = array(jk,1)*(energy(jk+1)-energy(jk))

     1                     *1.e+6

             sum = sum + array(jk,1)

             emid(jk) = 0.5*(energy(jk) + energy(jk+1))

          enddo

          if ( sum .le. 0.0) sum = 1.0

c         number fraction and integral number

          zsum = 0.0

          do jk=nenergy, 1, -1

             array(jk,1) = array(jk,1)/sum

             zsum = zsum + array(jk,1)

             array(jk,5) = zsum

          enddo

c         energy fraction

          esum = 0.0

          do jk=1,nenergy

             array(jk,2) = array(jk,1)*emid(jk)*1.e+6

             esum = esum + array(jk,2)

          enddo

          do jk=1,nenergy

             array(jk,2) = array(jk,2)/esum

          enddo

c         differential number fraction

          do jk=1,nenergy

             array(jk,3) = array(jk,1)/(energy(jk+1)-energy(jk))

     1                     /1.e+6

          enddo

c         differential energy fraction

          do jk=1,nenergy

             array(jk,4) = array(jk,2)/(energy(jk+1)-energy(jk))

     1                     /1.e+6

          enddo

c

c        output spectra

c

          if ( icon(9) .lt. 1) then

            write (nt6,8723) (jk, comment(jk),jk=1,icomment)

 8723       format (1x, 'comment ', i2, ': ', a80)

            avg = 0.0

            do jk=1,nenergy

                avg = avg + array(jk,1)*emid(jk)

            enddo

            write (nt6, 6623) avg

 6623       format (1x, 'average energy = ', g14.7)

            write (nt6, 8724)

 8724       format (/,/,/,/)

            write (nt6, 9011)

 9011        format (5x, 'lower', 10x, ' upper', 12x,

     1           'mid', 11x, 'number',

     1            11x, 'energy', 7x, 'differential',

     2            5x, 'differential',5x, 'integral', /,

     2            5x, 'energy', 9x, 'energy', 10x, 'energy',

     3            9x, 'fraction',

     3            9x, 'fraction', 9x, 'number', 10x, 'energy',

     4           10x, 'number', /,

     4            6x, 'mev', 12x, 'mev', 13x, 'mev',/ )

             do jk=1,nenergy

                 write (nt6, 9012) energy(jk), energy(jk+1), emid(jk),

     1           (array(jk,im), im=1,5)

 9012            format (1x, 8(g14.7, 2x))

             enddo

c

c            write out an interface file for the rebin code

c

             open(unit=47,file='rebin_energy',status='unknown')

             write (47, 8912) (energy(jk),jk=nenergy+1,1,-1)

             write (47, 8912) (array(jk,1),jk=nenergy,1,-1)

 8912        format (2x, 6(g14.7, 2x) )

             close (unit=47)



           endif

      return

      end

