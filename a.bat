set compilerName="C:\Program Files\Emscripten\emscripten\tag-1.37.3\third_party\closure-compiler\compiler.jar"
if exist "C:\Program Files\Emscripten\emscripten\tag-1.37.3" goto compile
set compilerName="C:\Program Files\Emscripten\emscripten\1.37.21\third_party\closure-compiler\compiler.jar"
:compile
del fsquares????.js
del fsquaresW????.js
copy fsquares.wasm fsquares%1.wasm
java -jar %compilerName% -D app=0 --compilation_level ADVANCED_OPTIMIZATIONS --js interface.js --js_output_file fsquares.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 FSQUARES.HTM fsquares.js
java -jar %compilerName% -D app=1 --compilation_level ADVANCED_OPTIMIZATIONS --js interface.js --js_output_file fsquares.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 SUMCUAD.HTM fsquares.js
java -jar %compilerName% -D app=2 --compilation_level ADVANCED_OPTIMIZATIONS --js interface.js --js_output_file fsquares.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 FCUBES.HTM fsquares.js
java -jar %compilerName% -D app=3 --compilation_level ADVANCED_OPTIMIZATIONS --js interface.js --js_output_file fsquares.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 SUMCUBOS.HTM fsquares.js
java -jar %compilerName% -D app=4 --compilation_level ADVANCED_OPTIMIZATIONS --js interface.js --js_output_file fsquares.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 CONTFRAC.HTM fsquares.js
java -jar %compilerName% -D app=5 --compilation_level ADVANCED_OPTIMIZATIONS --js interface.js --js_output_file fsquares.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 FRACCONT.HTM fsquares.js

cmd /c emcc -Os -DFSQUARES_APP=1 expression.c partition.c errors.c bigint.c division.c fsquares.c fcubes.c baseconv.c karatsuba.c modmult.c sqroot.c contfrac.c bignbr.c showtime.c inputstr.c batch.c gcdrings.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=33554432 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o fsquaresW%1.js
if errorlevel 1 goto end
java -jar %compilerName% --compilation_level ADVANCED_OPTIMIZATIONS --js intfwebw.js --js_output_file intWW.js

java -jar %compilerName% -D lang=0 --compilation_level ADVANCED_OPTIMIZATIONS --js polyfact.js --js_output_file polfactE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 POLFACT.HTM polfactE.js
java -jar %compilerName% -D lang=1 --compilation_level ADVANCED_OPTIMIZATIONS --js polyfact.js --js_output_file polfactS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 FACTPOL.HTM polfactS.js
cmd /c emcc -Os expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c polynomial.c polfact.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=67108864 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o polfactW%1.js
if errorlevel 1 goto end
copy polfact.wasm polfact%1.wasm

java -jar %compilerName% -D lang=0 --compilation_level ADVANCED_OPTIMIZATIONS --js dislog.js --js_output_file DILOGE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 DILOG.HTM DILOGE.js
java -jar %compilerName% -D lang=1 --compilation_level ADVANCED_OPTIMIZATIONS --js dislog.js --js_output_file DILOGS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 LOGDI.HTM DILOGS.js
cmd /c emcc -Os expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c dilog.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o dilogW%1.js
if errorlevel 1 goto end

java -jar %compilerName% -D lang=0 --compilation_level ADVANCED_OPTIMIZATIONS --js quadrmod.js --js_output_file QUADE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 QUADMOD.HTM QUADE.js
java -jar %compilerName% -D lang=1 --compilation_level ADVANCED_OPTIMIZATIONS --js quadrmod.js --js_output_file QUADS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 CUADMOD.HTM QUADS.js

cmd /c emcc -Os expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c quadmod.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o quadmodW%1.js
if errorlevel 1 goto end

java -jar %compilerName% -D lang=0 --compilation_level ADVANCED_OPTIMIZATIONS --js gauss.js --js_output_file gaussianE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 GAUSSIAN.HTM gaussianE.js
java -jar %compilerName% -D lang=1 --compilation_level ADVANCED_OPTIMIZATIONS --js gauss.js --js_output_file gaussianS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 GAUSIANO.HTM gaussianS.js
cmd /c emcc -Os GaussExpr.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c gaussian.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o gaussianW%1.js
if errorlevel 1 goto end

del ecm????.js
del ecmW????.js
java -jar %compilerName% -D lang=0 --compilation_level ADVANCED_OPTIMIZATIONS --js ecmfront.js --js_output_file ecmE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 ECM.HTM ECME.js
java -jar %compilerName% -D lang=1 --compilation_level ADVANCED_OPTIMIZATIONS --js ecmfront.js --js_output_file ecmS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 ECMC.HTM ECMS.js
perl generateServiceWorker.pl %1 calculatorSW.js calcSW.js
copy ecm.wasm ecm%1.wasm
set EMCC_DEBUG=
cmd /c emcc -Os -Wall -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 batch.c fft.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c ecmfront.c gcdrings.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s NO_FILESYSTEM=1 -s DOUBLE_MODE=1 --memory-init-file 0 -o ecmW%1.js
java -jar %compilerName% --compilation_level ADVANCED_OPTIMIZATIONS --js ecmfwebw.js --js_output_file ecmWW.js

java -jar %compilerName% -D lang=0 --compilation_level ADVANCED_OPTIMIZATIONS --js quadr.js --js_output_file quadE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 QUAD.HTM quadE.js
java -jar %compilerName% -D lang=1 --compilation_level ADVANCED_OPTIMIZATIONS --js quadr.js --js_output_file quadS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 CUAD.HTM quadS.js
copy quad.wasm quad%1.wasm
cmd /c emcc -Os -Wall -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 batch.c fft.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c quad.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s NO_FILESYSTEM=1 -s DOUBLE_MODE=1 --memory-init-file 0 -o quadW%1.js

java -jar %compilerName% --compilation_level WHITESPACE_ONLY --js dist.js --js_output_file dist%1.js

:end