@ECHO OFF
SET job=%1
if not defined job goto USAGE
SET GENROOT=\Projects-Git
SET GENPATH=\Projects-Git\Manipulate-2010
SET SYSTEMTYPE=dos
cd \sync_Projects-Git\Manipulate-2010\snl-work
perl \sync_Projects-Git\Manipulate-2010\input\manipulate-2010.prl %1 
set saveerr=%ERRORLEVEL%
cd \sync_Projects-Git\Manipulate-2010\input
if %saveerr% GTR 0 goto ERROR
ECHO Normal_Manipulate-2010_batch_file_termination-Bye!
goto END
:USAGE
ECHO Usage: Manipulate-2010 jobname
goto END
:ERROR
ECHO Manipulate-2010 Batch file error occurred
:END
