cmake_minimum_required(VERSION 2.8.12.2)
project(manipulate)
enable_language ( Fortran )

#Choose release or debug flags:
if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
  set( Release_FLAGS "-pipe -O3 -fno-math-errno -fno-signaling-nans -fno-trapv -fno-trapping-math -ffpe-summary='none' -fbranch-target-load-optimize2 -fdelete-dead-exceptions -fgcse-las -fgcse-sm -fipa-pta -fira-loop-pressure -fsched-pressure -fsched-spec-load -ftracer -w" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Release_FLAGS} -std=legacy" )
else ()
  set( Debug_FLAGS "-pipe -Og -fsignaling-nans -Wuninitialized -Wno-unused-parameter" )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Debug_FLAGS} -std=legacy -fbacktrace -fcheck=all -Wall -Wextra -Wno-align-commons -Wimplicit-procedure" )
endif ( CMAKE_BUILD_TYPE STREQUAL "Release" )
message( STATUS "Compiler used flags are: " ${CMAKE_Fortran_FLAGS} )

set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_Fortran_FLAGS}" )

set( SRC_FILES filein.f
               spectra_format.f
               grspin.f
               fold.f
               fitmd.f
               terp1.f
               manipulate.f
               prune.f
               grprin.f
               stat.f
               covlate.f
               covread.f
               force_upper.f
               correct.f
               spectra_out.f
               filein_response.f
               eigen_out.f
               balanc.f
               balbak.f
               cdiv.f
               fdump.f
               hqr.f
               hqr2.f
               i1mach.f
               j4save.f
               orthes.f
               ortran.f
               scopy.f
               scopym.f
               sgeev.f
               xercnt.f
               xerhlt.f
               xermsg.f
               xerprn.f
               xersve.f
               df.f
               xgetua.f
               cov_norm.f
               cov_combine.f
               pka_grprin.f
               integral_pka_grprin.f
               location.cmn )

# add the executable
add_executable( manipulate ${SRC_FILES} )
set_target_properties( manipulate PROPERTIES LINKER_LANGUAGE Fortran)
add_definitions( ${CMAKE_FORTRAN_FLAGS} )
install( TARGETS manipulate DESTINATION ${CMAKE_INSTALL_DIR} )
