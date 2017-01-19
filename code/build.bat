@echo off

set CompilerFlags=/nologo /Zi /FC /GR- /EHa- /Za /WX /W4 /wd4001 /wd4100 /wd4996

IF NOT EXIST build mkdir build
pushd build
cl %CompilerFlags% /c ..\code\matrix.c
cl %CompilerFlags% /c ..\code\stdMatrix.c
cl %CompilerFlags% /c ..\code\main.c
cl %CompilerFlags% /c ..\code\mainTeste.c
cl %CompilerFlags% main.obj stdMatrix.obj matrix.obj
cl %CompilerFlags% mainTeste.obj matrix.obj
popd