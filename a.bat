del /q toweb\*.*
if "%2" == "end" goto compress
set compilerName=%userprofile%\emsdk\emsdk\upstream\emscripten\node_modules\google-closure-compiler-java\compiler.jar
set compilerOptions=--compilation_level ADVANCED_OPTIMIZATIONS --isolation_mode IIFE --externs=custom-externs.js
set commonFlags=-Os --no-entry -Wall -s DYNAMIC_EXECUTION=0 -s SUPPORT_ERRNO=0 -s ASSERTIONS=0 -s NO_FILESYSTEM=1 --js-library lib.js --pre-js pre.js -finline-functions -DNDEBUG
set jsCommon=%commonFlags% -s WASM=0 -s SINGLE_FILE=1 -s TEXTDECODER=0 -s INCOMING_MODULE_JS_API=['preRun','noInitialRun'] -s WASM_ASYNC_COMPILATION=0 -s ENVIRONMENT='worker' --closure 1 --memory-init-file 0
set wasmCommon=%commonFlags% -s WASM=1 %commonFlags% -D_USING64BITS_ 
if "%1" == "" goto end
del *.wasm
del *00*js

rem ==================== GENERATION OF ASM.JS ===============================
set fsquaresFiles=expression.c parseexpr.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c output.c bignbr.c showtime.c inputstr.c batch.c gcdrings.c fft.c
set fsquaresOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=33554432
cmd /c emcc %jsCommon% %fsquaresFiles% fsquares.c tsquares.c %fsquaresOptions% -o fsquaresW%1.js
if errorlevel 1 goto end

cmd /c emcc %jsCommon% %fsquaresFiles% fcubes.c %fsquaresOptions% -o fcubesW%1.js
if errorlevel 1 goto end

cmd /c emcc %jsCommon% %fsquaresFiles% tsqcubes.c tsquares.c %fsquaresOptions% -o tsqcubesW%1.js
if errorlevel 1 goto end

cmd /c emcc %jsCommon% %fsquaresFiles% contfrac.c %fsquaresOptions% -o contfracW%1.js
if errorlevel 1 goto end

set polfactFiles=expression.c parseexpr.c partition.c errors.c bigint.c linkedbignbr.c division.c baseconv.c karatsuba.c modmult.c sqroot.c rootseq.c quintics.c quinticsData.c bigrational.c output.c polynomial.c polyexpr.c multpoly.c divpoly.c fftpoly.c intpolfact.c polfact.c polfacte.c bignbr.c showtime.c inputstr.c fft.c
set polfactOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=254869504 -DPOLYEXPR=1
cmd /c emcc %jsCommon% %polfactFiles% %polfactOptions% -o polfactW%1.js
if errorlevel 1 goto end

set dilogFiles=expression.c parseexpr.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c dilog.c bignbr.c showtime.c inputstr.c fft.c
set dilogOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888
cmd /c emcc %jsCommon% %dilogFiles% %dilogOptions% -o dilogW%1.js
if errorlevel 1 goto end

set quadmodFiles=expression.c parseexpr.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c quadmod.c quadmodLL.c bignbr.c showtime.c inputstr.c fft.c
set quadmodOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888
cmd /c emcc %jsCommon% %quadmodFiles% %quadmodOptions% -o quadmodW%1.js
if errorlevel 1 goto end

set gaussianFiles=GaussExpr.c parseexpr.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c gaussian.c output.c bignbr.c showtime.c inputstr.c fft.c
set gaussianOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888
cmd /c emcc %jsCommon% %gaussianFiles% %gaussianOptions% -o gaussianW%1.js
if errorlevel 1 goto end

set ecmFiles=batch.c fft.c expression.c parseexpr.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c ecmfront.c sumSquares.c gcdrings.c bignbr.c showtime.c inputstr.c fromBlockly.c linkedbignbr.c
set ecmOptions=-DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -DUSING_BLOCKLY=1 -DENABLE_VERBOSE=1 -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr','_getFactorsAsciiPtr']" -s TOTAL_MEMORY=278396928
cmd /c emcc %jsCommon% %ecmFiles% %ecmOptions% -o ecmW%1.js
if errorlevel 1 goto end

set quadFiles=batch.c fft.c expression.c parseexpr.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c quad.c quadmodLL.c output.c bignbr.c showtime.c inputstr.c
set quadOptions=-DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456
cmd /c emcc %jsCommon% %quadFiles% %quadOptions% -o quadW%1.js
if errorlevel 1 goto end

rem ===================== GENERATION OF WASM ================================
cmd /c emcc %wasmCommon% -Dlang=0 %fsquaresFiles% fsquares.c tsquares.c %fsquaresOptions% -o fsquaresE.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %fsquaresFiles% fsquares.c tsquares.c %fsquaresOptions% -o fsquaresS.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %fsquaresFiles% fcubes.c %fsquaresOptions% -o fcubesE.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %fsquaresFiles% fcubes.c %fsquaresOptions% -o fcubesS.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %fsquaresFiles% tsqcubes.c tsquares.c %fsquaresOptions% -o tsqcubesE.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %fsquaresFiles% tsqcubes.c tsquares.c %fsquaresOptions% -o tsqcubesS.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %fsquaresFiles% contfrac.c %fsquaresOptions% -o contfracE.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %fsquaresFiles% contfrac.c %fsquaresOptions% -o contfracS.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %polfactFiles% %polfactOptions% -o polfactE.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %polfactFiles% %polfactOptions% -o polfactS.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %ecmFiles% %ecmOptions% -o ecmE.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %ecmFiles% %ecmOptions% -o ecmS.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %dilogFiles% %dilogOptions% -o dilogE.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %dilogFiles% %dilogOptions% -o dilogS.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %quadFiles% %quadOptions% -o quadE.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %quadFiles% %quadOptions% -o quadS.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %quadmodFiles% %quadmodOptions% -o quadmodE.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %quadmodFiles% %quadmodOptions% -o quadmodS.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %gaussianFiles% %gaussianOptions% -o gaussianE.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %gaussianFiles% %gaussianOptions% -o gaussianS.wasm
if errorlevel 1 goto end

java -jar %compilerName% %compilerOptions% --js intfwebw.js --js strings.js --js_output_file intWW.js
java -jar %compilerName% %compilerOptions% --js ecmfwebw.js --js strings.js --js_output_file ecmWW.js

java -jar %compilerName% -D app=0 %compilerOptions% --js interface.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js wizard.js --js worker.js --js_output_file fsquares.js
if errorlevel 1 goto end
copy FSQUARES.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\FSQUARES.HTM fsquares.js fsquaresE.wasm intWW.js
java -jar %compilerName% -D app=1 %compilerOptions% --js interface.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js wizard.js --js worker.js --js_output_file fsquares.js
if errorlevel 1 goto end
copy SUMCUAD.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\SUMCUAD.HTM fsquares.js fsquaresS.wasm intWW.js
java -jar %compilerName% -D app=2 %compilerOptions% --js interface.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js wizard.js --js worker.js --js_output_file fsquares.js
if errorlevel 1 goto end
copy FCUBES.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\FCUBES.HTM fsquares.js fcubesE.wasm intWW.js
java -jar %compilerName% -D app=3 %compilerOptions% --js interface.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js wizard.js --js worker.js --js_output_file fsquares.js
if errorlevel 1 goto end
copy SUMCUBOS.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\SUMCUBOS.HTM fsquares.js fcubesS.wasm intWW.js
java -jar %compilerName% -D app=4 %compilerOptions% --js interface.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js wizard.js --js worker.js --js_output_file fsquares.js
if errorlevel 1 goto end
copy CONTFRAC.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\CONTFRAC.HTM fsquares.js contfracE.wasm intWW.js
java -jar %compilerName% -D app=5 %compilerOptions% --js interface.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js wizard.js --js worker.js --js_output_file fsquares.js
if errorlevel 1 goto end
copy FRACCONT.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\FRACCONT.HTM fsquares.js contfracS.wasm intWW.js
java -jar %compilerName% -D app=6 %compilerOptions% --js interface.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js wizard.js --js worker.js --js_output_file fsquares.js
if errorlevel 1 goto end
copy TSQCUBES.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\TSQCUBES.HTM fsquares.js tsqcubesE.wasm intWW.js
java -jar %compilerName% -D app=7 %compilerOptions% --js interface.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js wizard.js --js worker.js --js_output_file fsquares.js
if errorlevel 1 goto end
copy TCUADCUB.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\TCUADCUB.HTM fsquares.js tsqcubesS.wasm intWW.js

java -jar %compilerName% -D lang=0 %compilerOptions% --js polyfact.js --js common.js --js cache.js --js calccode.js --js feedback.js --js worker.js --js_output_file polfactE.js
if errorlevel 1 goto end
copy POLFACT.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\POLFACT.HTM polfactE.js polfactE.wasm intWW.js
java -jar %compilerName% -D lang=1 %compilerOptions% --js polyfact.js --js common.js --js cache.js --js calccode.js --js feedback.js --js worker.js --js_output_file polfactS.js
if errorlevel 1 goto end
copy FACTPOL.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\FACTPOL.HTM polfactS.js polfactS.wasm intWW.js

java -jar %compilerName% -D lang=0 %compilerOptions% --js dislog.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js worker.js --js_output_file dilogE.js
if errorlevel 1 goto end
copy DILOG.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\DILOG.HTM dilogE.js dilogE.wasm ecmWW.js
java -jar %compilerName% -D lang=1 %compilerOptions% --js dislog.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js worker.js --js_output_file dilogS.js
if errorlevel 1 goto end
copy LOGDI.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\LOGDI.HTM dilogS.js dilogS.wasm ecmWW.js

java -jar %compilerName% -D lang=0 %compilerOptions% --js quadrmod.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js worker.js --js_output_file quadmodE.js
if errorlevel 1 goto end
copy QUADMOD.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\QUADMOD.HTM quadmodE.js quadmodE.wasm ecmWW.js
java -jar %compilerName% -D lang=1 %compilerOptions% --js quadrmod.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js worker.js --js_output_file quadmodS.js
if errorlevel 1 goto end
copy CUADMOD.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\CUADMOD.HTM quadmodS.js quadmodS.wasm ecmWW.js

java -jar %compilerName% -D lang=0 %compilerOptions% --js gauss.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js worker.js --js_output_file gaussianE.js
if errorlevel 1 goto end
copy GAUSSIAN.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\GAUSSIAN.HTM gaussianE.js gaussianE.wasm ecmWW.js
java -jar %compilerName% -D lang=1 %compilerOptions% --js gauss.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js worker.js --js_output_file gaussianS.js
if errorlevel 1 goto end
copy GAUSIANO.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\GAUSIANO.HTM gaussianS.js gaussianS.wasm ecmWW.js

java -jar %compilerName% -D lang=0 %compilerOptions% --js blocklyextern.js --js buttons.js --js ecmfront.js --js common.js --js calccode.js --js cache.js --js feedback.js --js wizard.js --js worker.js --js_output_file ecmE.js
if errorlevel 1 goto end
copy ECM.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\ECM.HTM ecmE.js ecmE.wasm ecmWW.js
java -jar %compilerName% -D lang=1 %compilerOptions% --js blocklyextern.js --js buttons.js --js ecmfront.js --js common.js --js calccode.js --js cache.js --js feedback.js --js wizard.js --js worker.js --js_output_file ecmS.js
if errorlevel 1 goto end
set emcc_DEBUG=
copy ECMC.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\ECMC.HTM ecmS.js ecmS.wasm ecmWW.js
copy calculatorSW.js calcSW.js

java -jar %compilerName% -D lang=0 %compilerOptions% --js quadr.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js worker.js --js_output_file quadE.js
if errorlevel 1 goto end
copy QUAD.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\QUAD.HTM quadE.js quadE.wasm ecmWW.js
java -jar %compilerName% -D lang=1 %compilerOptions% --js quadr.js --js common.js --js buttons.js --js cache.js --js calccode.js --js feedback.js --js worker.js --js_output_file quadS.js
if errorlevel 1 goto end
copy CUAD.HTM toweb
perl replaceEmbeddedJS.pl %1 toweb\CUAD.HTM quadS.js quadS.wasm ecmWW.js

java -jar %compilerName% %compilerOptions% --js dist.js --js common.js --js_output_file distE.js
copy distE.js distS.js
copy DIST.HTM toweb
perl replaceEmbeddedJS.pl 0000 toweb\DIST.HTM distS.js
copy DISTANCE.HTM toweb
perl replaceEmbeddedJS.pl 0000 toweb\DISTANCE.HTM distE.js

call w.bat
del *.wasm
:compress
call a1.bat %1
:end