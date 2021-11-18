@echo off
echo.
echo Usage forAllCelfiles.bat [CELFILE_DIR] [RESULT_DIR] [Larpack parameters]
echo  example: forAllCelfiles.bat D:\data\cel D:\results -c d:\data\lib\HG133.bpmap
echo
echo Remember to change YOURDIR in this file to you local Larpack directory.
echo.
cd %1
for %%X in (*.cel) do (cd %2 & mkdir %%X & cd %%X && P:\YOURDIR\LarpackVc8.exe -i %1\%%X %3 %4 %5 %6 %7 %8 %9 && cd %1)
