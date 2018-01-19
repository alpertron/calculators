set compilerName="C:\Program Files\Emscripten\emscripten\tag-1.37.3\third_party\closure-compiler\compiler.jar"
if exist "C:\Program Files\Emscripten\emscripten\tag-1.37.3" goto compile
set compilerName="C:\Program Files\Emscripten\emscripten\1.37.21\third_party\closure-compiler\compiler.jar"
:compile
del fsquares????.js
del fsquaresW????.js
copy fsquares.wasm fsquares%1.wasm
java -jar %compilerName% --compilation_level WHITESPACE_ONLY --js interface.js --js_output_file fsquares%1.js
if errorlevel 1 goto end
cmd /c emcc -Os -DFSQUARES_APP=1 expression.c partition.c errors.c bigint.c division.c fsquares.c fcubes.c baseconv.c karatsuba.c modmult.c sqroot.c contfrac.c bignbr.c showtime.c inputstr.c batch.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=33554432 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o fsquaresW%1.js
if errorlevel 1 goto end

java -jar %compilerName% --compilation_level WHITESPACE_ONLY --js polyfact.js --js_output_file polfact%1.js
if errorlevel 1 goto end
cmd /c emcc -Os expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c polynomial.c polfact.c bignbr.c showtime.c inputstr.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=67108864 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o polfactW%1.js
if errorlevel 1 goto end

java -jar %compilerName% --compilation_level WHITESPACE_ONLY --js dislog.js --js_output_file dilog%1.js
if errorlevel 1 goto end
cmd /c emcc -Os expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c dilog.c bignbr.c showtime.c inputstr.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o dilogW%1.js
if errorlevel 1 goto end

java -jar %compilerName% --compilation_level WHITESPACE_ONLY --js quadrmod.js --js_output_file quadmod%1.js
if errorlevel 1 goto end
cmd /c emcc -Os expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c quadmod.c bignbr.c showtime.c inputstr.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o quadmodW%1.js
if errorlevel 1 goto end

java -jar %compilerName% --compilation_level WHITESPACE_ONLY --js gauss.js --js_output_file gaussian%1.js
if errorlevel 1 goto end
cmd /c emcc -Os GaussExpr.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c gaussian.c bignbr.c showtime.c inputstr.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o gaussianW%1.js
if errorlevel 1 goto end

del ecm????.js
del ecmW????.js
java -jar %compilerName% --compilation_level SIMPLE_OPTIMIZATIONS --js ecmfront.js --js_output_file ecm%1.js
if errorlevel 1 goto end
copy calculatorSW.js calcSW.js
copy ecm.wasm ecm%1.wasm
set EMCC_DEBUG=
cmd /c emcc -Os -Wall -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 batch.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c ecmfront.c bignbr.c showtime.c inputstr.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s NO_FILESYSTEM=1 -s DOUBLE_MODE=1 --memory-init-file 0 -o ecmW%1.js

java -jar %compilerName% --compilation_level WHITESPACE_ONLY --js dist.js --js_output_file dist%1.js

:end