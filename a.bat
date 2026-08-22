@if "%1" == "" goto :EOF
del /q toweb\*.*
@if "%2" == "end" goto compress
set optimization=--compilation_level ADVANCED_OPTIMIZATIONS
set compilerName=%userprofile%\emsdk\emsdk\upstream\emscripten\node_modules\google-closure-compiler-java\compiler.jar
set compilerOptions=%optimization% --language_out ECMASCRIPT_2015 --isolation_mode IIFE --externs=custom-externs.js --js androidextern.js --js commonNoAndroid.js --js supportsWasmSimd.js
set compilerOptionsAnd=%optimization% --language_out ECMASCRIPT_2015 --isolation_mode IIFE --externs=custom-externs.js --js androidextern.js --js commonAndroid.js
set compileFlags=-flto -Os -Wall -finline-functions -DNDEBUG
set commonLinkFlags=-flto -Os --no-entry -s ASSERTIONS=0 -s NO_FILESYSTEM=1 --js-library lib.js --pre-js pre.js
set jsCommon=%commonLinkFlags% -s WASM=0 -s DYNAMIC_EXECUTION=0 -s SINGLE_FILE=1 -s TEXTDECODER=2 -s INCOMING_MODULE_JS_API=['preRun','noInitialRun'] -s WASM_ASYNC_COMPILATION=0 -s ENVIRONMENT='worker' --closure 1
set wasmCommon=%commonLinkFlags% -s WASM=1 -D_USING64BITS_ -msimd128
del *.wasm
del *00*js

rem ===================== FILES TO USE DURING COMPILATION =====================
set fsquaresFiles=expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c MontgomeryMult.c sqroot.c output.c bignbr.c showtime.c inputstr.c batch.c gcdrings.c fft.c fftUtil.c
set fsquaresOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']"
set fsquaresJS=--js interface.js --js config.js --js common.js --js buttons.js --js feedback.js --js wizard.js

set polfactFiles=expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c linkedbignbr.c division.c baseconv.c karatsuba.c modmult.c MontgomeryMult.c sqroot.c rootseq.c lineareq.c quadraticeq.c cubiceq.c quartics.c quintics.c quinticsData.c bigrational.c output.c polynomial.c polyexpr.c multpoly.c divpoly.c fftpoly.c intpolfact.c modpolfact.c polfact.c polfacte.c bignbr.c showtime.c inputstr.c fft.c fftUtil.c
set polfactOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -DPOLYEXPR=1
set polfactJS=--js polyfact.js --js common.js --js feedback.js

set dilogFiles=expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c MontgomeryMult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c dilog.c bignbr.c showtime.c inputstr.c fft.c fftUtil.c
set dilogOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']"
set dilogJS=--js dislog.js --js config.js --js common.js --js buttons.js --js feedback.js

set quadmodFiles=expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c MontgomeryMult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c quadmod.c quadmodLL.c bignbr.c showtime.c inputstr.c fft.c fftUtil.c
set quadmodOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']"
set quadmodJS=--js quadrmod.js --js config.js --js common.js --js buttons.js --js feedback.js

set gaussianFiles=GaussExpr.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c MontgomeryMult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c gaussian.c output.c bignbr.c showtime.c inputstr.c gcdrings.c fft.c fftUtil.c
set gaussianOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']"
set gaussianJS=--js gauss.js --js config.js --js common.js --js buttons.js --js feedback.js

set ecmFiles=batch.c fft.c fftUtil.c expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c MontgomeryMult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c ecmfront.c sumSquares.c gcdrings.c bignbr.c showtime.c inputstr.c fromBlockly.c linkedbignbr.c
set ecmOptions=-s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr','_getFactorsAsciiPtr']"
set ecmJS=--js blocklyextern.js --js buttons.js --js ecmfront.js --js config.js --js common.js --js feedback.js --js wizard.js 

set quadFiles=batch.c fft.c fftUtil.c expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c MontgomeryMult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c quad.c quadmodLL.c output.c bignbr.c showtime.c inputstr.c
set quadOptions=-s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']"
set quadJS=--js quadr.js --js config.js --js common.js --js buttons.js --js feedback.js

rem ================================= COMPILE =================================
@call :compile %1 en
@call :compile %1 es
goto :generate_glue_code
@if errorlevel 1 exit /b 1
:compile
perl ClassicalMult_wasm.pl > ClassicalMult_wasm.S
perl internationalize.pl string_%2.txt string\strings.h
cmd /c emcc %compileFlags% %fsquaresFiles% fsquares.c ClassicalMult.c tsquares.c %jsCommon% %fsquaresOptions% -o toweb\fsquaresW%1%2.js
@if errorlevel 1 exit /b 1
cmd /c emcc %compileFlags% %fsquaresFiles% fsquares.c ClassicalMult_wasm.S tsquares.c %wasmCommon% %fsquaresOptions% -o fsquares_%2.wasm
@if errorlevel 1 exit /b 1

cmd /c emcc %compileFlags% %fsquaresFiles% fcubes.c ClassicalMult.c %jsCommon% %fsquaresOptions% -o toweb\fcubesW%1%2.js
@if errorlevel 1 exit /b 1
cmd /c emcc %compileFlags% %fsquaresFiles% fcubes.c ClassicalMult_wasm.S %wasmCommon% %fsquaresOptions% -o fcubes_%2.wasm
@if errorlevel 1 exit /b 1

cmd /c emcc %compileFlags% %fsquaresFiles% tsqcubes.c ClassicalMult.c tsquares.c %jsCommon% %fsquaresOptions% -o toweb\tsqcubesW%1%2.js
@if errorlevel 1 exit /b 1
cmd /c emcc %compileFlags% %fsquaresFiles% tsqcubes.c ClassicalMult_wasm.S tsquares.c %wasmCommon% %fsquaresOptions% -o tsqcubes_%2.wasm
@if errorlevel 1 exit /b 1

cmd /c emcc %compileFlags% %fsquaresFiles% contfrac.c ClassicalMult.c %jsCommon% %fsquaresOptions% -o toweb\contfracW%1%2.js
@if errorlevel 1 exit /b 1
cmd /c emcc %compileFlags% %fsquaresFiles% contfrac.c ClassicalMult_wasm.S %wasmCommon% %fsquaresOptions% -o contfrac_%2.wasm
@if errorlevel 1 exit /b 1

cmd /c emcc %compileFlags% %polfactFiles% ClassicalMult.c %jsCommon% %polfactOptions% -o toweb\polfactW%1%2.js
@if errorlevel 1 exit /b 1
cmd /c emcc %compileFlags% %polfactFiles% ClassicalMult_wasm.S %wasmCommon% %polfactOptions% -o polfact_%2.wasm
@if errorlevel 1 exit /b 1

cmd /c emcc %compileFlags% %dilogFiles% ClassicalMult.c %jsCommon% %dilogOptions% -o toweb\dilogW%1%2.js
@if errorlevel 1 exit /b 1
cmd /c emcc %compileFlags% %dilogFiles% ClassicalMult_wasm.S %wasmCommon% %dilogOptions% -o dilog_%2.wasm
@if errorlevel 1 exit /b 1

cmd /c emcc %compileFlags% %quadmodFiles% ClassicalMult.c %jsCommon% %quadmodOptions% -o toweb\quadmodW%1%2.js
@if errorlevel 1 exit /b 1
cmd /c emcc %compileFlags% %quadmodFiles% ClassicalMult_wasm.S %wasmCommon% %quadmodOptions% -o quadmod_%2.wasm
@if errorlevel 1 exit /b 1

cmd /c emcc %compileFlags% -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 %quadFiles% ClassicalMult.c %jsCommon% %quadOptions% -o toweb\quadW%1%2.js
@if errorlevel 1 exit /b 1
cmd /c emcc %compileFlags% -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 %quadFiles% ClassicalMult_wasm.S %wasmCommon% %quadOptions% -o quad_%2.wasm
@if errorlevel 1 exit /b 1

cmd /c emcc %compileFlags% %gaussianFiles% ClassicalMult.c %jsCommon% %gaussianOptions% -o toweb\gaussianW%1%2.js
@if errorlevel 1 exit /b 1
cmd /c emcc %compileFlags% %gaussianFiles% ClassicalMult_wasm.S %wasmCommon% %gaussianOptions% -o gaussian_%2.wasm
@if errorlevel 1 exit /b 1

cmd /c emcc %compileFlags% -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -DUSING_BLOCKLY=1 -DENABLE_VERBOSE=1 %ecmFiles% ClassicalMult.c %jsCommon% %ecmOptions% -o toweb\ecmW%1%2.js
@if errorlevel 1 exit /b 1
cmd /c emcc %compileFlags% -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -DUSING_BLOCKLY=1 -DENABLE_VERBOSE=1 %ecmFiles% ClassicalMult_wasm.S %wasmCommon% %ecmOptions% -o ecm_%2.wasm
@if errorlevel 1 exit /b 1

del ClassicalMult_wasm.S
del obj.o
@goto :EOF

:generate_glue_code
rem ========================= GENERATION OF GLUE CODE ========================
java -jar %compilerName% %compilerOptions% --js intfwebw.js --js commonwebw.js --js_output_file intWW.js
@if errorlevel 1 exit /b 1

java -jar %compilerName% -D app=0 %compilerOptions% --js cache.js --js calccode.js --js firstLineNoAndroid.js %fsquaresJS% --js worker.js --js_output_file WebGlue.js
@if errorlevel 1 exit /b 1
java -jar %compilerName% -D app=0 %compilerOptionsAnd% --js firstLineAndroid.js %fsquaresJS% --js workerAndroid.js --js_output_file AndroidGlue.js
@if errorlevel 1 exit /b 1
@call :generate_glue_code_subr FSQUARES.HTM fsquares.html SUMCUAD.HTM sumcuad.html fsquares %1

java -jar %compilerName% -D app=2 %compilerOptions% --js cache.js --js calccode.js --js firstLineNoAndroid.js %fsquaresJS% --js worker.js --js_output_file WebGlue.js
@if errorlevel 1 exit /b 1
java -jar %compilerName% -D app=2 %compilerOptionsAnd% --js firstLineAndroid.js %fsquaresJS% --js workerAndroid.js --js_output_file AndroidGlue.js
@if errorlevel 1 exit /b 1
@call :generate_glue_code_subr FCUBES.HTM fcubes.html SUMCUBOS.HTM sumcubos.html fcubes %1

java -jar %compilerName% -D app=4 %compilerOptions% --js cache.js --js calccode.js --js firstLineNoAndroid.js %fsquaresJS% --js worker.js --js_output_file WebGlue.js
@if errorlevel 1 exit /b 1
java -jar %compilerName% -D app=4 %compilerOptionsAnd% --js firstLineAndroid.js %fsquaresJS% --js workerAndroid.js --js_output_file AndroidGlue.js
@if errorlevel 1 exit /b 1
@call :generate_glue_code_subr CONTFRAC.HTM contfrac.html FRACCONT.HTM fraccont.html contfrac %1

java -jar %compilerName% -D app=6 %compilerOptions% --js cache.js --js calccode.js --js firstLineNoAndroid.js %fsquaresJS% --js worker.js --js_output_file WebGlue.js
@if errorlevel 1 exit /b 1
java -jar %compilerName% -D app=6 %compilerOptionsAnd% --js firstLineAndroid.js %fsquaresJS% --js workerAndroid.js --js_output_file AndroidGlue.js
@if errorlevel 1 exit /b 1
@call :generate_glue_code_subr TSQCUBES.HTM tsqcubes.html TCUADCUB.HTM tcuadcub.html tsqcubes %1

java -jar %compilerName% %compilerOptions% --js ecmfwebw.js --js commonwebw.js --js_output_file intWW.js
@if errorlevel 1 exit /b 1

java -jar %compilerName% -D android=0 %compilerOptions% --js cache.js --js calccode.js %polfactJS% --js worker.js --js_output_file WebGlue.js
@if errorlevel 1 exit /b 1
java -jar %compilerName% -D android=1 %compilerOptionsAnd% %polfactJS% --js workerAndroid.js --js_output_file AndroidGlue.js
@if errorlevel 1 exit /b 1
@call :generate_glue_code_subr POLFACT.HTM polfact.html FACTPOL.HTM factpol.html polfact %1

java -jar %compilerName% %compilerOptions% --js cache.js --js calccode.js %dilogJS% --js worker.js --js_output_file WebGlue.js
@if errorlevel 1 exit /b 1
java -jar %compilerName% %compilerOptionsAnd% %dilogJS% --js workerAndroid.js --js_output_file AndroidGlue.js
@if errorlevel 1 exit /b 1
@call :generate_glue_code_subr DILOG.HTM dilog.html LOGDI.HTM logdi.html dilog %1

java -jar %compilerName% %compilerOptions% --js cache.js --js calccode.js %quadmodJS% --js worker.js --js_output_file WebGlue.js
@if errorlevel 1 exit /b 1
java -jar %compilerName% %compilerOptionsAnd% %quadmodJS% --js workerAndroid.js --js_output_file AndroidGlue.js
@if errorlevel 1 exit /b 1
@call :generate_glue_code_subr QUADMOD.HTM quadmod.html CUADMOD.HTM cuadmod.html quadmod %1

java -jar %compilerName% %compilerOptions% --js cache.js --js calccode.js --js firstLineNoAndroid.js %gaussianJS% --js worker.js --js_output_file WebGlue.js
@if errorlevel 1 exit /b 1
java -jar %compilerName% %compilerOptionsAnd% --js firstLineAndroid.js %gaussianJS% --js workerAndroid.js --js_output_file AndroidGlue.js
@if errorlevel 1 exit /b 1
@call :generate_glue_code_subr GAUSSIAN.HTM gaussian.html GAUSIANO.HTM gausiano.html gaussian %1

java -jar %compilerName% %compilerOptions% --js cache.js --js calccode.js --js ecmNoAndroid.js --js firstLineNoAndroid.js %ecmJS% --js worker.js --js_output_file WebGlue.js
@if errorlevel 1 exit /b 1
java -jar %compilerName% %compilerOptionsAnd% --js ecmAndroid.js --js firstLineAndroid.js %ecmJS% --js workerAndroid.js --js_output_file AndroidGlue.js
@if errorlevel 1 exit /b 1
@call :generate_glue_code_subr ECM.HTM ecm.html ECMC.HTM ecmc.html ecm %1
copy calculatorSW.js toweb\calcSW.js

java -jar %compilerName% %compilerOptions% --js cache.js --js calccode.js %quadJS% --js worker.js --js_output_file WebGlue.js
@if errorlevel 1 exit /b 1
java -jar %compilerName% %compilerOptionsAnd% %quadJS% --js workerAndroid.js --js_output_file AndroidGlue.js
@if errorlevel 1 exit /b 1
@call :generate_glue_code_subr QUAD.HTM quad.html CUAD.HTM cuad.html quad %1

java -jar %compilerName% %compilerOptions% --js dist.js --js common.js --js_output_file distE.js
@if errorlevel 1 exit /b 1
copy distE.js distS.js
copy DIST.HTM toweb
perl replaceEmbeddedJS.pl 0000 toweb\DIST.HTM distS.js
copy DISTANCE.HTM toweb
perl replaceEmbeddedJS.pl 0000 toweb\DISTANCE.HTM distE.js
@goto :generate_graphic_apps

:generate_glue_code_subr
copy WebGlue.js WebGlueBak.js
copy AndroidGlue.js AndroidGlueBak.js
copy %1 toweb
copy %1 assets\%2
perl replaceEmbeddedJS.pl %6 toweb\%1 WebGlue.js %5_en.wasm intWW.js
perl replaceEmbeddedJSAnd.pl %6 assets\%2 AndroidGlue.js calc_privacy.html
copy %3 toweb
copy %3 assets\%4
perl replaceEmbeddedJS.pl %6 toweb\%3 WebGlueBak.js %5_es.wasm intWW.js
perl replaceEmbeddedJSAnd.pl %6 assets\%4 AndroidGlueBak.js privacidad_calc.html
@goto :EOF

:generate_graphic_apps
@call w.bat
del *.wasm
:compress
@call a1.bat %1
