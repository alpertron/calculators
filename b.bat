set compilerName=%userprofile%\emsdk\emsdk\upstream\emscripten\node_modules\google-closure-compiler-java\compiler.jar
set compilerOptions=--compilation_level ADVANCED_OPTIMIZATIONS --language_in ECMASCRIPT5_STRICT --language_out ECMASCRIPT5_STRICT --externs=custom-externs.js
set commonOptions=-Os -Wall -s WASM=0 -s TEXTDECODER=0 -s MIN_IE_VERSION=11 -s ASSERTIONS=0 -s NO_FILESYSTEM=1 --closure 1 --llvm-lto 1 --js-library lib.js --pre-js pre.js --memory-init-file 0
set wasmCommon=-Os -Wall -s WASM=1 -D_USING64BITS_ -flto -s ASSERTIONS=0 -s NO_FILESYSTEM=1
if "%1" == "" goto :end
del *.wasm
del *00*js

rem ==================== GENERATION OF ASM.JS ===============================
call emsdk activate latest-fastcomp > NUL

cmd /c emcc -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -DENABLE_VERBOSE=1 batch.c fft.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c ecmfront.c gcdrings.c bignbr.c showtime.c inputstr.c -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s DOUBLE_MODE=1 %commonOptions% -o ecmW%1.js
if errorlevel 1 goto end

rem ===================== GENERATION OF WASM ================================
call emsdk activate latest > NUL

cmd /c emcc %wasmCommon% -Dlang=0 -flto -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -DENABLE_VERBOSE=1 batch.c fft.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c ecmfront.c gcdrings.c bignbr.c showtime.c inputstr.c -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -o ecmE.wasm
if errorlevel 1 goto end

cmd /c emcc %wasmCommon% -Dlang=1 -flto -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -DENABLE_VERBOSE=1 batch.c fft.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c ecmfront.c gcdrings.c bignbr.c showtime.c inputstr.c -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -o ecmS.wasm
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