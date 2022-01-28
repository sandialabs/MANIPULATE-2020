      subroutine stat

c      character assignments

       include "location.cmn"

       character*155 filename

       character*250 outfile

       common /whatever/ outfile

       common /datain/ nenergy, energy(1001), array(1000,15)

     1         , emid(1001)

       common /fpnew/ alpha

       common /datahld/ nenergy_hld(15), energy_hld(1001,15),

     1          array_hld(1000,15)

     1         ,emid_hld(1001,15)

       common /guide/ icon(40)

       common /io/ nt5, nt6, nfile, nplot, npun

       common /pltlab/ lab_lead(15), lab_trail(15), lab_nchar(15),

     1 lab_len(15), icurve(15), lab_file(15)

       character*80 lab_file

       character*20 lab

      character*106 xoptical

       common /misc/ lab(20), imaterial, noption, moption

      common /statsum/ arr_mean(1000), arr_sq(1000), arr_std(1000)

     1   ,arr_var(1000)

      common /statinfo/ efactor, number_of_files, xoptical

      do jk=1,1000

          arr_mean(jk) = 0.0

          arr_sq(jk) = 0.0

          arr_std(jk) = 0.0

          arr_var(jk) = 0.0

      enddo

      do ie = 1,nenergy

         sum = 0.0

         sum2 = 0.0

         do jk=1,number_of_files

            sum = sum + array_hld(ie,jk)

            sum2 = sum2 + array_hld(ie,jk)*array_hld(ie,jk)

         enddo

         arr_mean(ie) = sum/number_of_files

         arr_sq(ie) = sum2

      enddo

      do ie=1,nenergy

         zap = (arr_sq(ie) -

     1           number_of_files*arr_mean(ie)*arr_mean(ie))/

     1       (number_of_files - 1)

         if ( zap .lt. 0.0) zap = 0.0

         arr_var(ie) = sqrt(zap)*100.0

         if ( arr_mean(ie) .ne. 0.0) then

             arr_std(ie) = arr_var(ie)/arr_mean(ie)

         elseif (arr_mean(ie) .eq. 0.0 .and. arr_var(ie)

     1                        .eq. 0.0) then

             arr_std(ie) = 0.0

         else

             arr_std(ie) = 1.e+30

         endif

      enddo

c

c      punch output files

c

       leopt = lnblnk(xoptical)

       jblank2 = lnblnk(jdir)

       jout = lnblnk(outfile)

       filename = xoptical(1:leopt)//jdir(1:jblank2)//outfile(1:jout)

     1       //'.avg'

      write (6,9023) filename

 9023  format (1x, 'stat filename open = ', a150)

       open(unit=nplot,

     1 file=filename

     1 ,status='unknown')

      if ( icon(4) .eq. 1) then

           write (nplot,3002) (emid_hld(ie,1),

     1         arr_mean(ie),

     1         ie=1,nenergy)

           else

               write (nplot,3002) (emid_hld(ie,1),

     1         arr_mean(ie),

     1         ie=nenergy,1,-1)

 3002           format (2x, g14.7, 2x, ' ', 2x, g14.7)

      endif

      close (unit=nplot)

c

       filename = xoptical(1:leopt)//jdir(1:jblank2)//outfile(1:jout)

     1       //'.pctstd'

      write (6,9023) filename

       open(unit=nplot,

     1 file=filename

     1 ,status='unknown')

      if ( icon(4) .eq. 1) then

           write (nplot,3002) (emid_hld(ie,1),

     1         arr_std(ie),

     1         ie=1,nenergy)

           else

               write (nplot,3002) (emid_hld(ie,1),

     1         arr_std(ie),

     1         ie=nenergy,1,-1)

      endif

      close (unit=nplot)

c

      return

      end

