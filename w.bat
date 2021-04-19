set compilerName=%userprofile%\emsdk\emsdk\upstream\emscripten\node_modules\google-closure-compiler-java\compiler.jar
rem set compilerOptions2=--compilation_level WHITESPACE_ONLY --language_in ECMASCRIPT5_STRICT --language_out ECMASCRIPT5_STRICT --externs=custom-externs.js
set compilerOptions2=--compilation_level SIMPLE_OPTIMIZATIONS --language_in ECMASCRIPT5_STRICT --language_out ECMASCRIPT5_STRICT --externs=custom-externs.js
set commonOptions=--no-entry -Os -Wall -s WASM=0 -s TEXTDECODER=0 -s MIN_IE_VERSION=11 -s ASSERTIONS=0 -s NO_FILESYSTEM=1 -s WASM_ASYNC_COMPILATION=0 --pre-js preGraphics.js --closure 1 --memory-init-file 0
set wasmCommon=--no-entry -Os -Wall -s WASM=1 -D_USING64BITS_ -s ASSERTIONS=0 -s NO_FILESYSTEM=1 --js-library lib.js --pre-js pre.js
rem del *.wasm
rem del *00*js

rem ==================== GENERATION OF ASM.JS ===============================
cmd /c emcc ulam.c -s EXPORTED_FUNCTIONS="['_moveSpiral', '_drawPartialUlamSpiral', '_nbrChanged', '_getInformation', '_getPixels']" -s TOTAL_MEMORY=33554432 %commonOptions% -o ulamW.js
if errorlevel 1 goto end
cmd /c emcc gausspr.c -s EXPORTED_FUNCTIONS="['_moveGraphic', '_drawPartialGraphic', '_nbrChanged', '_getInformation', '_getPixels']" -s TOTAL_MEMORY=33554432 %commonOptions% -o gaussprW.js
if errorlevel 1 goto end

rem ===================== GENERATION OF WASM ================================
cmd /c emcc %wasmCommon% ulam.c -s EXPORTED_FUNCTIONS="['_moveSpiral', '_drawPartialUlamSpiral', '_nbrChanged', '_getInformation', '_getPixels']" -s TOTAL_MEMORY=33554432 -o ulam.wasm
if errorlevel 1 goto end
cmd /c emcc %wasmCommon% gausspr.c -s EXPORTED_FUNCTIONS="['_moveGraphic', '_drawPartialGraphic', '_nbrChanged', '_getInformation', '_getPixels']" -s TOTAL_MEMORY=33554432 -o gausspr.wasm
if errorlevel 1 goto end

perl generateTempJS.pl ulam.js ulamW.js ulamT.js moveSpiral drawPartialUlamSpiral nbrChanged getInformation getPixels
java -jar %compilerName% %compilerOptions2% --js ulamT.js --js_output_file ulamU.js
copy ulamU.js ulamV.js
perl replaceEmbeddedJS.pl 0000 ULAM.HTM ulamV.js ulam.wasm
perl replaceEmbeddedJS.pl 0000 EULAM.HTM ulamU.js ulam.wasm

perl generateTempJS.pl gausspr.js gaussprW.js gaussprT.js moveGraphic drawPartialGraphic nbrChanged getInformation getPixels
java -jar %compilerName% %compilerOptions2% --js gaussprT.js --js_output_file gaussprU.js
copy gaussprU.js gaussprV.js 
perl replaceEmbeddedJS.pl 0000 GAUSSPR.HTM gaussprV.js gausspr.wasm
perl replaceEmbeddedJS.pl 0000 PRGAUSS.HTM gaussprU.js gausspr.wasm

:end