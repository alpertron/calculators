del /q toweb\*.*
@if "%2" == "end" goto compress 
set compilerName=%userprofile%\emsdk\emsdk\upstream\emscripten\node_modules\google-closure-compiler-java\compiler.jar
set compilerOptions=--compilation_level ADVANCED_OPTIMIZATIONS --isolation_mode IIFE --externs=custom-externs.js --js androidextern.js --js commonNoAndroid.js
set compilerOptionsAnd=--compilation_level ADVANCED_OPTIMIZATIONS --isolation_mode IIFE --externs=custom-externs.js --js androidextern.js --js commonAndroid.js
set commonFlags=-Os --no-entry -Wall -s DYNAMIC_EXECUTION=0 -s SUPPORT_ERRNO=0 -s ASSERTIONS=0 -s NO_FILESYSTEM=1 --js-library lib.js --pre-js pre.js -finline-functions -DNDEBUG
set jsCommon=%commonFlags% -s WASM=0 -s SINGLE_FILE=1 -s TEXTDECODER=0 -s INCOMING_MODULE_JS_API=['preRun','noInitialRun'] -s WASM_ASYNC_COMPILATION=0 -s ENVIRONMENT='worker' --closure 1 --memory-init-file 0
set wasmCommon=%commonFlags% -s WASM=1 %commonFlags% -D_USING64BITS_ 
@if "%1" == "" goto end
del *.wasm
del *00*js

rem ==================== GENERATION OF ASM.JS ===============================
set fsquaresFiles=expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c output.c bignbr.c showtime.c inputstr.c batch.c gcdrings.c fft.c
set fsquaresOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=33554432
set fsquaresJS=--js interface.js --js config.js --js common.js --js buttons.js --js feedback.js --js wizard.js
cmd /c emcc %jsCommon% %fsquaresFiles% fsquares.c tsquares.c %fsquaresOptions% -o fsquaresW%1.js
@if errorlevel 1 goto end

cmd /c emcc %jsCommon% %fsquaresFiles% fcubes.c %fsquaresOptions% -o fcubesW%1.js
@if errorlevel 1 goto end

cmd /c emcc %jsCommon% %fsquaresFiles% tsqcubes.c tsquares.c %fsquaresOptions% -o tsqcubesW%1.js
@if errorlevel 1 goto end

cmd /c emcc %jsCommon% %fsquaresFiles% contfrac.c %fsquaresOptions% -o contfracW%1.js
@if errorlevel 1 goto end

set polfactFiles=expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c linkedbignbr.c division.c baseconv.c karatsuba.c modmult.c sqroot.c rootseq.c lineareq.c quadraticeq.c cubiceq.c quartics.c quintics.c quinticsData.c bigrational.c output.c polynomial.c polyexpr.c multpoly.c divpoly.c fftpoly.c intpolfact.c modpolfact.c polfact.c polfacte.c bignbr.c showtime.c inputstr.c fft.c
set polfactOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=255066112 -DPOLYEXPR=1
set polfactJS=--js polyfact.js --js common.js --js feedback.js
cmd /c emcc %jsCommon% %polfactFiles% %polfactOptions% -o polfactW%1.js
@if errorlevel 1 goto end

set dilogFiles=expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c dilog.c bignbr.c showtime.c inputstr.c fft.c
set dilogOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888
set dilogJS=--js dislog.js --js config.js --js common.js --js buttons.js --js feedback.js
cmd /c emcc %jsCommon% %dilogFiles% %dilogOptions% -o dilogW%1.js
@if errorlevel 1 goto end

set quadmodFiles=expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c quadmod.c quadmodLL.c bignbr.c showtime.c inputstr.c fft.c
set quadmodOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888
set quadmodJS=--js quadrmod.js --js config.js --js common.js --js buttons.js --js feedback.js
cmd /c emcc %jsCommon% %quadmodFiles% %quadmodOptions% -o quadmodW%1.js
@if errorlevel 1 goto end

set gaussianFiles=GaussExpr.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c gaussian.c output.c bignbr.c showtime.c inputstr.c gcdrings.c fft.c
set gaussianOptions=-s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888
set gaussianJS=--js gauss.js --js config.js --js common.js --js buttons.js --js feedback.js
cmd /c emcc %jsCommon% %gaussianFiles% %gaussianOptions% -o gaussianW%1.js
@if errorlevel 1 goto end

set ecmFiles=batch.c fft.c expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c ecmfront.c sumSquares.c gcdrings.c bignbr.c showtime.c inputstr.c fromBlockly.c linkedbignbr.c
set ecmOptions=-DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -DUSING_BLOCKLY=1 -DENABLE_VERBOSE=1 -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr','_getFactorsAsciiPtr']" -s TOTAL_MEMORY=280363008
set ecmJS=--js blocklyextern.js --js buttons.js --js ecmfront.js --js config.js --js common.js --js feedback.js --js wizard.js 
cmd /c emcc %jsCommon% %ecmFiles% %ecmOptions% -o ecmW%1.js
@if errorlevel 1 goto end

set quadFiles=batch.c fft.c expression.c parseexpr.c partition.c errors.c copyStr.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c ecm.c siqs.c siqsLA.c quad.c quadmodLL.c output.c bignbr.c showtime.c inputstr.c
set quadOptions=-DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=261095424
set quadJS=--js quadr.js --js config.js --js common.js --js buttons.js --js feedback.js
cmd /c emcc %jsCommon% %quadFiles% %quadOptions% -o quadW%1.js
@if errorlevel 1 goto end

rem ===================== GENERATION OF WASM ================================
cmd /c emcc %wasmCommon% -Dlang=0 %fsquaresFiles% fsquares.c tsquares.c %fsquaresOptions% -o fsquaresE.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %fsquaresFiles% fsquares.c tsquares.c %fsquaresOptions% -o fsquaresS.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %fsquaresFiles% fcubes.c %fsquaresOptions% -o fcubesE.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %fsquaresFiles% fcubes.c %fsquaresOptions% -o fcubesS.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %fsquaresFiles% tsqcubes.c tsquares.c %fsquaresOptions% -o tsqcubesE.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %fsquaresFiles% tsqcubes.c tsquares.c %fsquaresOptions% -o tsqcubesS.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %fsquaresFiles% contfrac.c %fsquaresOptions% -o contfracE.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %fsquaresFiles% contfrac.c %fsquaresOptions% -o contfracS.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %polfactFiles% %polfactOptions% -o polfactE.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %polfactFiles% %polfactOptions% -o polfactS.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %ecmFiles% %ecmOptions% -o ecmE.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %ecmFiles% %ecmOptions% -o ecmS.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %dilogFiles% %dilogOptions% -o dilogE.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %dilogFiles% %dilogOptions% -o dilogS.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %quadFiles% %quadOptions% -o quadE.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %quadFiles% %quadOptions% -o quadS.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %quadmodFiles% %quadmodOptions% -o quadmodE.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %quadmodFiles% %quadmodOptions% -o quadmodS.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=0 %gaussianFiles% %gaussianOptions% -o gaussianE.wasm
@if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 %gaussianFiles% %gaussianOptions% -o gaussianS.wasm
@if errorlevel 1 goto end

java -jar %compilerName% %compilerOptions% --js intfwebw.js --js strings.js --js_output_file intWW.js
@if errorlevel 1 goto end
java -jar %compilerName% %compilerOptions% --js ecmfwebw.js --js strings.js --js_output_file ecmWW.js
@if errorlevel 1 goto end

java -jar %compilerName% -D app=0 %compilerOptions% --js cache.js --js calccode.js %fsquaresJS% --js worker.js --js_output_file fsquaresE.js
@if errorlevel 1 goto end
java -jar %compilerName% -D app=0 %compilerOptionsAnd% %fsquaresJS% --js workerAndroid.js --js_output_file fsquaresEA.js
@if errorlevel 1 goto end
copy FSQUARES.HTM toweb
copy FSQUARES.HTM assets\fsquares.html
perl replaceEmbeddedJS.pl %1 toweb\FSQUARES.HTM fsquaresE.js fsquaresE.wasm intWW.js
perl replaceEmbeddedJSAnd.pl %1 assets\fsquares.html fsquaresEA.js calc_privacy.html
java -jar %compilerName% -D app=1 %compilerOptions% --js cache.js --js calccode.js %fsquaresJS% --js worker.js --js_output_file fsquaresS.js
@if errorlevel 1 goto end
java -jar %compilerName% -D app=1 %compilerOptionsAnd% %fsquaresJS% --js workerAndroid.js --js_output_file fsquaresSA.js
@if errorlevel 1 goto end
copy SUMCUAD.HTM toweb
copy SUMCUAD.HTM "assets\sumcuad.html"
perl replaceEmbeddedJS.pl %1 toweb\SUMCUAD.HTM fsquaresS.js fsquaresS.wasm intWW.js
perl replaceEmbeddedJSAnd.pl %1 "assets\sumcuad.html" fsquaresSA.js privacidad_calc.html

java -jar %compilerName% -D app=2 %compilerOptions% --js cache.js --js calccode.js %fsquaresJS% --js worker.js --js_output_file fsquaresE.js
@if errorlevel 1 goto end
java -jar %compilerName% -D app=2 %compilerOptionsAnd% %fsquaresJS% --js workerAndroid.js --js_output_file fsquaresEA.js
@if errorlevel 1 goto end
copy FCUBES.HTM toweb
copy FCUBES.HTM assets\fcubes.html
perl replaceEmbeddedJS.pl %1 toweb\FCUBES.HTM fsquaresE.js fcubesE.wasm intWW.js
perl replaceEmbeddedJSAnd.pl %1 assets\fcubes.html fsquaresEA.js calc_privacy.html
java -jar %compilerName% -D app=3 %compilerOptions% --js cache.js --js calccode.js %fsquaresJS% --js worker.js --js_output_file fsquaresS.js
@if errorlevel 1 goto end
java -jar %compilerName% -D app=3 %compilerOptionsAnd% %fsquaresJS% --js workerAndroid.js --js_output_file fsquaresSA.js
@if errorlevel 1 goto end
copy SUMCUBOS.HTM toweb
copy SUMCUBOS.HTM "assets\sumcubos.html"
perl replaceEmbeddedJS.pl %1 toweb\SUMCUBOS.HTM fsquaresS.js fcubesS.wasm intWW.js
perl replaceEmbeddedJSAnd.pl %1 "assets\sumcubos.html" fsquaresSA.js privacidad_calc.html

java -jar %compilerName% -D app=4 %compilerOptions% --js cache.js --js calccode.js %fsquaresJS% --js worker.js --js_output_file fsquaresE.js
@if errorlevel 1 goto end
java -jar %compilerName% -D app=4 %compilerOptionsAnd% %fsquaresJS% --js workerAndroid.js --js_output_file fsquaresEA.js
@if errorlevel 1 goto end
copy CONTFRAC.HTM toweb
copy CONTFRAC.HTM assets\contfrac.html
perl replaceEmbeddedJS.pl %1 toweb\CONTFRAC.HTM fsquaresE.js contfracE.wasm intWW.js
perl replaceEmbeddedJSAnd.pl %1 assets\contfrac.html fsquaresEA.js calc_privacy.html
java -jar %compilerName% -D app=5 %compilerOptions% --js cache.js --js calccode.js %fsquaresJS% --js worker.js --js_output_file fsquaresS.js
@if errorlevel 1 goto end
java -jar %compilerName% -D app=5 %compilerOptionsAnd% %fsquaresJS% --js workerAndroid.js --js_output_file fsquaresSA.js
@if errorlevel 1 goto end
copy FRACCONT.HTM toweb
copy FRACCONT.HTM "assets\fraccont.html"
perl replaceEmbeddedJS.pl %1 toweb\FRACCONT.HTM fsquaresS.js contfracS.wasm intWW.js
perl replaceEmbeddedJSAnd.pl %1 "assets\fraccont.html" fsquaresSA.js privacidad_calc.html

java -jar %compilerName% -D app=6 %compilerOptions% --js cache.js --js calccode.js %fsquaresJS% --js worker.js --js_output_file fsquaresE.js
@if errorlevel 1 goto end
java -jar %compilerName% -D app=6 %compilerOptionsAnd% %fsquaresJS% --js workerAndroid.js --js_output_file fsquaresEA.js
@if errorlevel 1 goto end
copy TSQCUBES.HTM toweb
copy TSQCUBES.HTM assets\tsqcubes.html
perl replaceEmbeddedJS.pl %1 toweb\TSQCUBES.HTM fsquaresE.js tsqcubesE.wasm intWW.js
perl replaceEmbeddedJSAnd.pl %1 assets\tsqcubes.html fsquaresEA.js calc_privacy.html
java -jar %compilerName% -D app=7 %compilerOptions% --js cache.js --js calccode.js %fsquaresJS% --js worker.js --js_output_file fsquaresS.js
java -jar %compilerName% -D app=7 %compilerOptionsAnd% %fsquaresJS% --js workerAndroid.js --js_output_file fsquaresSA.js
@if errorlevel 1 goto end
copy TCUADCUB.HTM toweb
copy TCUADCUB.HTM "assets\tcuadcub.html"
perl replaceEmbeddedJS.pl %1 toweb\TCUADCUB.HTM fsquaresS.js tsqcubesS.wasm intWW.js
perl replaceEmbeddedJSAnd.pl %1 "assets\tcuadcub.html" fsquaresSA.js privacidad_calc.html

java -jar %compilerName% -D lang=0 %compilerOptions% --js cache.js --js calccode.js %polfactJS% --js worker.js --js_output_file polfactE.js
@if errorlevel 1 goto end
java -jar %compilerName% -D lang=0 %compilerOptionsAnd% %polfactJS% --js workerAndroid.js --js_output_file polfactEA.js
@if errorlevel 1 goto end
copy POLFACT.HTM toweb
copy POLFACT.HTM assets\polfact.html
perl replaceEmbeddedJS.pl %1 toweb\POLFACT.HTM polfactE.js polfactE.wasm ecmWW.js
perl replaceEmbeddedJSAnd.pl %1 assets\polfact.html polfactEA.js calc_privacy.html
java -jar %compilerName% -D lang=1 %compilerOptions% --js cache.js --js calccode.js %polfactJS% --js worker.js --js_output_file polfactS.js
@if errorlevel 1 goto end
java -jar %compilerName% -D lang=1 %compilerOptionsAnd% %polfactJS% --js workerAndroid.js --js_output_file polfactSA.js
@if errorlevel 1 goto end
copy FACTPOL.HTM toweb
copy FACTPOL.HTM "assets\factpol.html"
perl replaceEmbeddedJS.pl %1 toweb\FACTPOL.HTM polfactS.js polfactS.wasm ecmWW.js
perl replaceEmbeddedJSAnd.pl %1 "assets\factpol.html" polfactSA.js privacidad_calc.html

java -jar %compilerName% -D lang=0 %compilerOptions% --js cache.js --js calccode.js %dilogJS% --js worker.js --js_output_file dilogE.js
@if errorlevel 1 goto end
java -jar %compilerName% -D lang=0 %compilerOptionsAnd% %dilogJS% --js workerAndroid.js --js_output_file dilogEA.js
@if errorlevel 1 goto end
copy DILOG.HTM toweb
copy DILOG.HTM assets\dilog.html
perl replaceEmbeddedJS.pl %1 toweb\DILOG.HTM dilogE.js dilogE.wasm ecmWW.js
perl replaceEmbeddedJSAnd.pl %1 assets\dilog.html dilogEA.js calc_privacy.html
java -jar %compilerName% -D lang=1 %compilerOptions% --js cache.js --js calccode.js %dilogJS% --js worker.js --js_output_file dilogS.js
@if errorlevel 1 goto end
java -jar %compilerName% -D lang=1 %compilerOptionsAnd% %dilogJS% --js workerAndroid.js --js_output_file dilogSA.js
@if errorlevel 1 goto end
copy LOGDI.HTM toweb
copy LOGDI.HTM "assets\logdi.html"
perl replaceEmbeddedJS.pl %1 toweb\LOGDI.HTM dilogS.js dilogS.wasm ecmWW.js
perl replaceEmbeddedJSAnd.pl %1 "assets\logdi.html" dilogSA.js privacidad_calc.html

java -jar %compilerName% -D lang=0 %compilerOptions% --js cache.js --js calccode.js %quadmodJS% --js worker.js --js_output_file quadmodE.js
@if errorlevel 1 goto end
java -jar %compilerName% -D lang=0 %compilerOptionsAnd% %quadmodJS% --js workerAndroid.js --js_output_file quadmodEA.js
@if errorlevel 1 goto end
copy QUADMOD.HTM toweb
copy QUADMOD.HTM assets\quadmod.html
perl replaceEmbeddedJS.pl %1 toweb\QUADMOD.HTM quadmodE.js quadmodE.wasm ecmWW.js
perl replaceEmbeddedJSAnd.pl %1 assets\quadmod.html quadmodEA.js calc_privacy.html
java -jar %compilerName% -D lang=1 %compilerOptions% --js cache.js --js calccode.js %quadmodJS% --js worker.js --js_output_file quadmodS.js
@if errorlevel 1 goto end
java -jar %compilerName% -D lang=1 %compilerOptionsAnd% %quadmodJS% --js workerAndroid.js --js_output_file quadmodSA.js
@if errorlevel 1 goto end
copy CUADMOD.HTM toweb
copy CUADMOD.HTM "assets\cuadmod.html"
perl replaceEmbeddedJS.pl %1 toweb\CUADMOD.HTM quadmodS.js quadmodS.wasm ecmWW.js
perl replaceEmbeddedJSAnd.pl %1 "assets\cuadmod.html" quadmodSA.js privacidad_calc.html

java -jar %compilerName% -D lang=0 %compilerOptions% --js cache.js --js calccode.js %gaussianJS% --js worker.js --js_output_file gaussianE.js
@if errorlevel 1 goto end
java -jar %compilerName% -D lang=0 %compilerOptionsAnd% %gaussianJS% --js workerAndroid.js --js_output_file gaussianEA.js
@if errorlevel 1 goto end
copy GAUSSIAN.HTM toweb
copy GAUSSIAN.HTM assets\gaussian.html
perl replaceEmbeddedJS.pl %1 toweb\GAUSSIAN.HTM gaussianE.js gaussianE.wasm ecmWW.js
perl replaceEmbeddedJSAnd.pl %1 assets\gaussian.html gaussianEA.js calc_privacy.html
java -jar %compilerName% -D lang=1 %compilerOptions% --js cache.js --js calccode.js %gaussianJS% --js worker.js --js_output_file gaussianS.js
@if errorlevel 1 goto end
java -jar %compilerName% -D lang=1 %compilerOptionsAnd% %gaussianJS% --js workerAndroid.js --js_output_file gaussianSA.js
@if errorlevel 1 goto end
copy GAUSIANO.HTM toweb
copy GAUSIANO.HTM "assets\gausiano.html"
perl replaceEmbeddedJS.pl %1 toweb\GAUSIANO.HTM gaussianS.js gaussianS.wasm ecmWW.js
perl replaceEmbeddedJSAnd.pl %1 "assets\gausiano.html" gaussianSA.js privacidad_calc.html

java -jar %compilerName% -D lang=0 %compilerOptions% --js cache.js --js calccode.js --js ecmNoAndroid.js %ecmJS% --js worker.js --js_output_file ecmE.js
@if errorlevel 1 goto end
java -jar %compilerName% -D lang=0 %compilerOptionsAnd% --js ecmAndroid.js %ecmJS% --js workerAndroid.js --js_output_file ecmEA.js
@if errorlevel 1 goto end
copy ECM.HTM toweb
copy ECM.HTM assets\ecm.html
perl replaceEmbeddedJS.pl %1 toweb\ECM.HTM ecmE.js ecmE.wasm ecmWW.js
perl replaceEmbeddedJSAnd.pl %1 assets\ecm.html ecmEA.js calc_privacy.html
java -jar %compilerName% -D lang=1 %compilerOptions% --js cache.js --js calccode.js --js ecmNoAndroid.js %ecmJS% --js worker.js --js_output_file ecmS.js
@if errorlevel 1 goto end
java -jar %compilerName% -D lang=1 %compilerOptionsAnd% --js ecmAndroid.js %ecmJS% --js workerAndroid.js --js_output_file ecmSA.js
@if errorlevel 1 goto end
java -jar %compilerName% %compilerOptions% --js ecmfwebw.js --js commonwebw.js --js_output_file ecmWW.js
copy ECMC.HTM toweb
copy ECMC.HTM "assets\ecmc.html"
perl replaceEmbeddedJS.pl %1 toweb\ECMC.HTM ecmS.js ecmS.wasm ecmWW.js
perl replaceEmbeddedJSAnd.pl %1 "assets\ecmc.html" ecmSA.js privacidad_calc.html
copy calculatorSW.js toweb\calcSW.js

java -jar %compilerName% -D lang=0 %compilerOptions% --js cache.js --js calccode.js %quadJS% --js worker.js --js_output_file quadE.js
@if errorlevel 1 goto end
java -jar %compilerName% -D lang=0 %compilerOptionsAnd% %quadJS% --js workerAndroid.js --js_output_file quadEA.js
@if errorlevel 1 goto end
copy QUAD.HTM toweb
copy QUAD.HTM assets\quad.html
perl replaceEmbeddedJS.pl %1 toweb\QUAD.HTM quadE.js quadE.wasm ecmWW.js
perl replaceEmbeddedJSAnd.pl %1 assets\quad.html quadEA.js calc_privacy.html
java -jar %compilerName% -D lang=1 %compilerOptions% --js cache.js --js calccode.js %quadJS% --js worker.js --js_output_file quadS.js
@if errorlevel 1 goto end
java -jar %compilerName% -D lang=1 %compilerOptionsAnd% %quadJS% --js workerAndroid.js --js_output_file quadSA.js
@if errorlevel 1 goto end
copy CUAD.HTM toweb
copy CUAD.HTM "assets\cuad.html"
perl replaceEmbeddedJS.pl %1 toweb\CUAD.HTM quadS.js quadS.wasm ecmWW.js
perl replaceEmbeddedJSAnd.pl %1 "assets\cuad.html" quadSA.js privacidad_calc.html

java -jar %compilerName% %compilerOptions% --js dist.js --js common.js --js_output_file distE.js
@if errorlevel 1 goto end
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