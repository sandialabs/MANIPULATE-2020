@ECHO OFF
SET template=%1
SET master=%2
SET action=%3
if not defined template goto USAGE
perl manipulate_loop.prl %1 %2 %3
set saveerr=%ERRORLEVEL%
cd ..\input
if %saveerr% GTR 0 goto ERROR
ECHO Normal_
atch_file_termination
goto END
:USAGE
ECHO Usage: BATCH Manipulate_loop template master_list action_file
goto END
:ERROR
ECHO Batch file error occurred
:END
