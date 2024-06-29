@echo off
setlocal

:: aliases for compilers
set "intel_fcomp=ifort"
set "gcc_fcomp=gfortran"

:: compiler flags
set "INTEL_FFLAGS=-O3"
set "GCC_FFLAGS=-O3 -fno-finite-math-only"


:: SHA256 for Fortran sources to be compiled (with \r terminating chars!)
set "EEASISS_IMPORTED_ROUTINES_SHA256=a3de67cdad38a5ed5980aa4be5f2d86a8362f2e5d1b68c12524eb26cd89926a5"
set "EEASISS_SHA256=919ebb06e72dff8f914b392d66fec6ffc2413b8f0f714db4b6682a7562784e3e"
set "CHECKSUM_VALIDATED=False"


IF /I "%1"=="intel" GOTO intel
IF /I "%1"=="gcc" GOTO gcc
IF /I "%1"=="clean" GOTO clean
GOTO error


:intel
	CALL :correct_sha256
    if %CHECKSUM_VALIDATED%==True (
        %intel_fcomp% eeasisss_src/imported_routines.f90 eeasisss_src/eeasisss.f90 -o eeasisss %INTEL_FFLAGS%
    )
    call :_remove_temp_files
    GOTO :EOF

:gcc
	CALL :correct_sha256
	%gcc_fcomp% eeasisss_src/imported_routines.f90 eeasisss_src/eeasisss.f90 -o eeasisss %GCC_FFLAGS%
	call :_remove_temp_files
    GOTO :EOF

:correct_sha256
    call :check_sha eeasisss_src/imported_routines.f90 , %EEASISS_IMPORTED_ROUTINES_SHA256%
    if %CHECKSUM_VALIDATED%==True (
        call :check_sha eeasisss_src/eeasisss.f90 , %EEASISS_SHA256%
    )
    del lf_fixed.txt 2>nul
	GOTO :EOF

:clean
    del eeasisss.exe 2>nul
    call :_remove_temp_files
	GOTO :EOF

:error
    IF "%1"=="" (
        ECHO make: *** No target specified. Use "make intel", "make gcc", or "make clean". Stop.
    ) ELSE (
        ECHO make: *** Unkown target '%1%'. Stop.
    )
    GOTO :EOF




:: Arguments are file and expected checksum
:check_sha
    :: The next one is for expansion of the for loop
    ::setlocal
    set fname=%~1
    set expected_sha=%~2
    echo Validating checksum for file '%fname%'

    :: First, remove carriage-return characters
    call :crlf_to_lf %fname% lf_fixed.txt

    setlocal EnableDelayedExpansion

    :: Collect the SHA256 with certUtil. To do so, parse the output
    :: of the certUtil call. The checksum is the first token
    set idx=0
    for /f %%F in ('certutil -hashfile lf_fixed.txt SHA256') do (
        set "out!idx!=%%F"
        set /a idx += 1
    )
    if /i %expected_sha%==%out1% (
        endlocal & set "CHECKSUM_VALIDATED=True"
        echo     OK.
    ) else (
        endlocal & set "CHECKSUM_VALIDATED=False"
        echo     *** Checksum mismatch!
        echo     Expected= %expected_sha%
        echo     Found   = %out1%
    )
    del lf_fixed.txt 2>nul
    GOTO :EOF

:: Turn carriage-return&line-feed pairs to line-feed only, leaving
:: all the line-feed characters to themselves.
:: Arguments are name of file with CRLF to be fixed and name of
:: file where the fixed version should be written to
:crlf_to_lf
    powershell.exe -noninteractive -NoProfile -ExecutionPolicy Bypass -Command "& {[IO.File]::WriteAllText('%~2', ([IO.File]::ReadAllText('%~1') -replace \"`r`n\", \"`n\"))};"
    GOTO :EOF

:_remove_temp_files
    :: Remove mod and obj files created during
    :: compilation. We only need the executable
    del "*.mod" 2>nul
    del "*.obj" 2>nul
    GOTO :EOF
