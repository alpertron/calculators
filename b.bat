set compilerName=%userprofile%\emsdk\emsdk\upstream\emscripten\node_modules\google-closure-compiler-java\compiler.jar
set compilerOptions=--compilation_level ADVANCED_OPTIMIZATIONS --language_in ECMASCRIPT5_STRICT --language_out ECMASCRIPT5_STRICT --externs=custom-externs.js
set commonFlags=--no-entry -Wall -s DYNAMIC_EXECUTION=0 -s SUPPORT_ERRNO=0 -s ASSERTIONS=0 -s NO_FILESYSTEM=1 --js-library lib.js --pre-js pre.js
set jsCommon=-Os -s WASM=0 %commonFlags% -s SINGLE_FILE=1 -s TEXTDECODER=0 -s MIN_IE_VERSION=11 -s INCOMING_MODULE_JS_API=['preRun','noInitialRun'] -s LEGACY_VM_SUPPORT=1 -s WASM_ASYNC_COMPILATION=0 -s ENVIRONMENT='worker' --closure 1 --memory-init-file 0
set wasmCommon=-Os -s WASM=1 %commonFlags% -D_USING64BITS_ 
if "%1" == "" goto :end
del *.wasm
del *00*js

rem ==================== GENERATION OF ASM.JS ===============================
set ecmFiles=batch.c fft.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c ecmfront.c gcdrings.c bignbr.c showtime.c inputstr.c
set ecmOptions=-DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -DENABLE_VERBOSE=1 -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s DOUBLE_MODE=1
cmd /c emcc %jsCommon% %ecmFiles% %ecmOptions% -o ecmW%1.js
if errorlevel 1 goto end

rem ===================== GENERATION OF WASM ================================

cmd /c emcc %wasmCommon% %ecmFiles% %ecmOptions% -Dlang=0 -o ecmE.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% %ecmFiles% %ecmOptions% -Dlang=1 -o ecmS.wasm
if errorlevel 1 goto end

set emcc_DEBUG=
java -jar %compilerName% %compilerOptions% --js ecmfwebw.js --js_output_file ecmWW.js
java -jar %compilerName% -D lang=0 %compilerOptions% --js ecmfront.js --js_output_file ecmE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 ECM.HTM ecmE.js ecmE.wasm
java -jar %compilerName% -D lang=1 %compilerOptions% --js ecmfront.js --js_output_file ecmS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 ECMC.HTM ecmS.js ecmS.wasm
copy calculatorSW.js calcSW.js

perl csp.pl
:end