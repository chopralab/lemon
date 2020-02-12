trap { Write-Error $_; Exit 1 }

git clone -q https://github.com/chemfiles/chemfiles.git C:\chemfiles

C:\Python35-x64\python.exe C:\lemon\scripts\appveyor\windows_build_wheels.py
