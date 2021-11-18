@ECHO OFF
SET template=%1
SET master=%2
SET action=%3
if not defined template goto USAGE
SET LOOPROOT=\sync_Linux\Manipulate-2020
SET XSECROOT=\sync_Linux\xsec
SET LOOPPATH=\sync_Linux\Manipulate-2020\input
SET SYSTEMTYPE=dos
cd \sync_Linux\Manipulate-2010\input
perl \sync_Linux\Manipulate-2020\input\manipulate_loop.prl %1 %2 %3
set saveerr=%ERRORLEVEL%
cd \sync_Linux\Manipulate-2020\input
if %saveerr% GTR 0 goto ERROR
ECHO Normal_batch_file_termination-Bye
goto END
:USAGE
ECHO Usage: BATCH Manipulate_loop template master_list action_file
goto END
:ERROR
ECHO Batch file error occurred
:END
