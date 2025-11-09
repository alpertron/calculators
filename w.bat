set compilerName=%userprofile%\emsdk\emsdk\upstream\emscripten\node_modules\google-closure-compiler-java\compiler.jar
set compilerOptions2=--compilation_level ADVANCED_OPTIMIZATIONS --isolation_mode IIFE --externs=custom-externs.js
set commonOptions=--no-entry -Os -Wall -s WASM=0 -s TEXTDECODER=0 -s ASSERTIONS=0 -s NO_FILESYSTEM=1 -s WASM_ASYNC_COMPILATION=0 --pre-js preGraphics.js --closure 1 --memory-init-file 0
set wasmCommon=--no-entry -Os -Wall -s WASM=1 -D_USING64BITS_ -s ASSERTIONS=0 -s NO_FILESYSTEM=1 --js-library lib.js --pre-js pre.js
rem del *.wasm
rem del *00*js
rem --compilation_level WHITESPACE_ONLY
rem --compilation_level ADVANCED_OPTIMIZATIONS
rem ===================== GENERATION OF WASM ================================
cmd /c emcc %wasmCommon% ulam.c isprime.c MontMultGraphic.c graphics.c copyStr.c -s EXPORTED_FUNCTIONS="['_initUlam', '_moveGraphic', '_drawPartialGraphic', '_nbrChanged', '_getInformation', '_getPixels']" -s TOTAL_MEMORY=33554432 -o ulam.wasm
if errorlevel 1 goto end
cmd /c emcc %wasmCommon% gausspr.c isprime.c MontMultGraphic.c graphics.c -s EXPORTED_FUNCTIONS="['_initGaussPr', '_moveGraphic', '_drawPartialGraphic', '_nbrChanged', '_getInformation', '_getPixels']" -s TOTAL_MEMORY=33554432 -o gausspr.wasm
if errorlevel 1 goto end

copy /b ulam.js + common.js + strings.js + commonGraphics.js ulamT.js
java -jar %compilerName% %compilerOptions2% --js ulamT.js --js initGraphicNoAndroid.js --js_output_file ulamU.js
copy ulamU.js ulamV.js
copy ULAM.HTM toweb
copy ULAM.HTM assets\ulam.html
perl replaceEmbeddedJS.pl 0000 toweb\ULAM.HTM ulamV.js ulam.wasm
@if errorlevel 1 goto end
copy /b ulam.js + common.js + commonAndroid.js + stringsAndroid.js + commonGraphics.js ulamT.js
java -jar %compilerName% %compilerOptions2% --js ulamT.js --js initGraphicAndroid.js --js androidextern.js --js_output_file ulamEA.js
@if errorlevel 1 goto end
perl replaceEmbeddedJSAnd.pl 0000 "assets\ulam.html" ulamEA.js privacidad_calc.html
@if errorlevel 1 goto end
copy EULAM.HTM toweb
copy EULAM.HTM assets\eulam.html
perl replaceEmbeddedJS.pl 0000 toweb\EULAM.HTM ulamU.js ulam.wasm
@if errorlevel 1 goto end
java -jar %compilerName% %compilerOptions2% --js ulamT.js --js initGraphicAndroid.js --js androidextern.js --js_output_file ulamSA.js
@if errorlevel 1 goto end
perl replaceEmbeddedJSAnd.pl 0000 "assets\eulam.html" ulamSA.js privacidad_calc.html
@if errorlevel 1 goto end

copy /b gausspr.js + common.js + strings.js + commonGraphics.js gaussprT.js
java -jar %compilerName% %compilerOptions2% --js gaussprT.js --js initGraphicNoAndroid.js --js_output_file gaussprU.js
copy gaussprU.js gaussprV.js
copy GAUSSPR.HTM toweb
copy GAUSSPR.HTM assets\gausspr.html
perl replaceEmbeddedJS.pl 0000 toweb\GAUSSPR.HTM gaussprV.js gausspr.wasm
@if errorlevel 1 goto end
copy /b gausspr.js + common.js + commonAndroid.js + stringsAndroid.js + commonGraphics.js gaussprT.js
java -jar %compilerName% %compilerOptions2% --js gaussprT.js --js initGraphicAndroid.js --js androidextern.js --js_output_file gaussprEA.js
@if errorlevel 1 goto end
perl replaceEmbeddedJSAnd.pl 0000 "assets\gausspr.html" gaussprEA.js privacidad_calc.html
@if errorlevel 1 goto end
copy PRGAUSS.HTM toweb
copy PRGAUSS.HTM assets\prgauss.html
perl replaceEmbeddedJS.pl 0000 toweb\PRGAUSS.HTM gaussprU.js gausspr.wasm
java -jar %compilerName% %compilerOptions2% --js gaussprT.js --js initGraphicAndroid.js --js androidextern.js --js_output_file gaussprSA.js
@if errorlevel 1 goto end
perl replaceEmbeddedJSAnd.pl 0000 "assets\prgauss.html" gaussprSA.js privacidad_calc.html
del ulamT.js
del gaussprT.js

:end