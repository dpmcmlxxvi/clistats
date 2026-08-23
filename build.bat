@echo off
setlocal

set INSTALL_DIR=%~dp0install

if exist build rmdir /s /q build
mkdir build
cd build

cmake .. -DCMAKE_INSTALL_PREFIX="%INSTALL_DIR%"
if errorlevel 1 goto :error

cmake --build . --config Release
if errorlevel 1 goto :error

cmake --build . --config Release --target install
if errorlevel 1 goto :error

echo Build succeeded. Executable installed to %INSTALL_DIR%
cd ..
exit /b 0

:error
echo Build failed.
cd ..
exit /b 1