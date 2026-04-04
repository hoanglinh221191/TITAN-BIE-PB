@echo off
setlocal
set "REPO_ROOT=%~dp0"
if "%REPO_ROOT:~-1%"=="\" set "REPO_ROOT=%REPO_ROOT:~0,-1%"

set "PROFILE=native"
set "BUILD_TARGET=all"
set "ARG1=%~1"
set "ARG2=%~2"

if /I "%ARG1%"=="native" (
    set "PROFILE=native"
    if not "%ARG2%"=="" set "BUILD_TARGET=%ARG2%"
) else if /I "%ARG1%"=="portable" (
    set "PROFILE=portable"
    if not "%ARG2%"=="" set "BUILD_TARGET=%ARG2%"
) else if not "%ARG1%"=="" (
    set "BUILD_TARGET=%ARG1%"
)

if /I not "%PROFILE%"=="native" if /I not "%PROFILE%"=="portable" (
    echo Unsupported build profile: %PROFILE%
    echo Usage: compile.bat [native^|portable] [all^|release^|clean^|print-config]
    exit /b 2
)

echo Building Titan ^(profile=%PROFILE%, target=%BUILD_TARGET%^)...
if /I "%BUILD_TARGET%"=="all" (
    echo Cleaning previous objects for a profile-correct rebuild...
    C:\msys64\ucrt64\bin\mingw32-make.exe PROFILE=%PROFILE% clean
    if %errorlevel% neq 0 (
        echo.
        echo Cleanup failed!
        if /I not "%TITAN_NO_PAUSE%"=="1" pause
        exit /b %errorlevel%
    )
)

REM Invoke MSYS2 make with a pinned build profile.
C:\msys64\ucrt64\bin\mingw32-make.exe PROFILE=%PROFILE% %BUILD_TARGET% -j1

if %errorlevel% neq 0 (
    echo.
    echo Compilation failed!
    if /I not "%TITAN_NO_PAUSE%"=="1" pause
    exit /b %errorlevel%
)

if /I "%BUILD_TARGET%"=="all" (
    set "PROFILE_BIN_DIR=%REPO_ROOT%\test\bin\%PROFILE%"
    powershell -NoProfile -Command ^
        "$out='%REPO_ROOT%\test\bin\%PROFILE%';" ^
        "New-Item -ItemType Directory -Force -Path $out | Out-Null;" ^
        "Copy-Item -LiteralPath '%REPO_ROOT%\titan.exe' -Destination (Join-Path $out 'titan.exe') -Force"
    if %errorlevel% neq 0 (
        echo.
        echo Failed to snapshot titan.exe to %PROFILE_BIN_DIR%
        if /I not "%TITAN_NO_PAUSE%"=="1" pause
        exit /b %errorlevel%
    )
)

echo.
echo Compilation finished!
if /I "%BUILD_TARGET%"=="all" echo Profile snapshot: %PROFILE_BIN_DIR%\titan.exe
if /I not "%TITAN_NO_PAUSE%"=="1" pause
