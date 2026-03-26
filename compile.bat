@echo off
echo Building Titan...
:: Gọi công cụ make của MSYS2 với cờ -j (chạy song song đa luồng để tận dụng CPU)
C:\msys64\ucrt64\bin\mingw32-make.exe -j1

if %errorlevel% neq 0 (
    echo.
    echo Compilation failed!
    pause
    exit /b %errorlevel%
)

echo.
echo Compilation finished!
pause