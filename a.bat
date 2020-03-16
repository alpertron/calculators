set compilerName=%userprofile%\emsdk\emsdk\upstream\emscripten\node_modules\google-closure-compiler-java\compiler.jar
set compilerOptions=--compilation_level ADVANCED_OPTIMIZATIONS --language_in ECMASCRIPT5_STRICT --language_out ECMASCRIPT5_STRICT --externs=custom-externs.js
if "%1" == "" goto :end
del *.wasm
del *00*js

rem ==================== GENERATION OF ASM.JS ===============================
cmd /c emsdk activate latest-fastcomp > NUL
cmd /c emcc -Os -s WASM=0 -DFSQUARES_APP=1 expression.c partition.c errors.c bigint.c division.c fsquares.c fcubes.c baseconv.c karatsuba.c modmult.c sqroot.c contfrac.c output.c bignbr.c showtime.c inputstr.c batch.c gcdrings.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=33554432 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o fsquaresW%1.js
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=0 expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c rootseq.c quintics.c bigrational.c output.c polynomial.c intpolfact.c polfact.c polfacte.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o polfactW%1.js
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=0 expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c dilog.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o dilogW%1.js
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=0 expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c quadmod.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o quadmodW%1.js
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=0 GaussExpr.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c gaussian.c output.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o gaussianW%1.js
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=0 -Wall -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -DENABLE_VERBOSE=1 batch.c fft.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c ecmfront.c gcdrings.c bignbr.c showtime.c inputstr.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s NO_FILESYSTEM=1 -s DOUBLE_MODE=1 --memory-init-file 0 -o ecmW%1.js
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=0 -Wall -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 batch.c fft.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c quad.c output.c bignbr.c showtime.c inputstr.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s NO_FILESYSTEM=1 -s DOUBLE_MODE=1 --memory-init-file 0 -o quadW%1.js
if errorlevel 1 goto end

rem ===================== GENERATION OF WASM ================================
cmd /c emsdk activate latest > NUL
cmd /c emcc -Os -s WASM=1 -Dlang=0 -flto --llvm-lto 1 -Dapplic=0 -DFSQUARES_APP=1 -D_USING64BITS_ expression.c partition.c errors.c bigint.c division.c fsquares.c baseconv.c karatsuba.c modmult.c sqroot.c output.c bignbr.c showtime.c inputstr.c batch.c gcdrings.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=33554432 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o fsquaresE.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=1 -flto --llvm-lto 1 -Dapplic=0 -DFSQUARES_APP=1 -D_USING64BITS_ expression.c partition.c errors.c bigint.c division.c fsquares.c baseconv.c karatsuba.c modmult.c sqroot.c output.c bignbr.c showtime.c inputstr.c batch.c gcdrings.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=33554432 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o fsquaresS.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=0 -flto --llvm-lto 1 -Dapplic=1 -DFSQUARES_APP=1 -D_USING64BITS_ expression.c partition.c errors.c bigint.c division.c fsquares.c fcubes.c baseconv.c karatsuba.c modmult.c sqroot.c output.c bignbr.c showtime.c inputstr.c batch.c gcdrings.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=33554432 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o fcubesE.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=1 -flto --llvm-lto 1 -Dapplic=1 -DFSQUARES_APP=1 -D_USING64BITS_ expression.c partition.c errors.c bigint.c division.c fsquares.c fcubes.c baseconv.c karatsuba.c modmult.c sqroot.c output.c bignbr.c showtime.c inputstr.c batch.c gcdrings.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=33554432 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o fcubesS.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=0 -flto --llvm-lto 1 -Dapplic=2 -DFSQUARES_APP=1 -D_USING64BITS_ expression.c partition.c errors.c bigint.c division.c fsquares.c baseconv.c karatsuba.c modmult.c sqroot.c contfrac.c output.c bignbr.c showtime.c inputstr.c batch.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=33554432 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o contfracE.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=1 -flto --llvm-lto 1 -Dapplic=2 -DFSQUARES_APP=1 -D_USING64BITS_ expression.c partition.c errors.c bigint.c division.c fsquares.c baseconv.c karatsuba.c modmult.c sqroot.c contfrac.c output.c bignbr.c showtime.c inputstr.c batch.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=33554432 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o contfracS.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=0 -flto --llvm-lto 1 -D_USING64BITS_ expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c rootseq.c quintics.c bigrational.c output.c polynomial.c intpolfact.c polfact.c polfacte.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o polfactE.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=1 -flto --llvm-lto 1 -D_USING64BITS_ expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c rootseq.c quintics.c bigrational.c output.c polynomial.c intpolfact.c polfact.c polfacte.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o polfactS.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=0 -flto --llvm-lto 1 -D_USING64BITS_ -Wall -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -DENABLE_VERBOSE=1 batch.c fft.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c ecmfront.c gcdrings.c bignbr.c showtime.c inputstr.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s NO_FILESYSTEM=1 -s DOUBLE_MODE=1 --memory-init-file 0 -o ecmE.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=1 -flto --llvm-lto 1 -D_USING64BITS_ -Wall -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 -DENABLE_VERBOSE=1 batch.c fft.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c ecmfront.c gcdrings.c bignbr.c showtime.c inputstr.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s NO_FILESYSTEM=1 -s DOUBLE_MODE=1 --memory-init-file 0 -o ecmS.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=0 -flto --llvm-lto 1 -D_USING64BITS_ expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c dilog.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o dilogE.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=1 -flto --llvm-lto 1 -D_USING64BITS_ expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c dilog.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o dilogS.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=0 -flto --llvm-lto 1 -D_USING64BITS_ -Wall -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 batch.c fft.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c quad.c output.c bignbr.c showtime.c inputstr.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s NO_FILESYSTEM=1 -s DOUBLE_MODE=1 --memory-init-file 0 -o quadE.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=1 -flto --llvm-lto 1 -D_USING64BITS_ -Wall -DFACTORIZATION_FUNCTIONS=1 -DFACTORIZATION_APP=1 batch.c fft.c expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c quad.c output.c bignbr.c showtime.c inputstr.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_copyString','_getInputStringPtr']" -s TOTAL_MEMORY=268435456 -s NO_FILESYSTEM=1 -s DOUBLE_MODE=1 --memory-init-file 0 -o quadS.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=0 -flto --llvm-lto 1 -D_USING64BITS_  expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c quadmod.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o quadmodE.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=1 -flto --llvm-lto 1 -D_USING64BITS_  expression.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c quadmod.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o quadmodS.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=0 -flto --llvm-lto 1 -D_USING64BITS_  GaussExpr.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c gaussian.c output.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o gaussianE.wasm
if errorlevel 1 goto end

cmd /c emcc -Os -s WASM=1 -Dlang=1 -flto --llvm-lto 1 -D_USING64BITS_  GaussExpr.c partition.c errors.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c siqs.c gaussian.c output.c bignbr.c showtime.c inputstr.c fft.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork','_getInputStringPtr']" -s TOTAL_MEMORY=301989888 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o gaussianS.wasm
if errorlevel 1 goto end

java -jar %compilerName% -D app=0 %compilerOptions% --js interface.js --js_output_file fsquares.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 FSQUARES.HTM fsquares.js fsquaresE.wasm
java -jar %compilerName% -D app=1 %compilerOptions% --js interface.js --js_output_file fsquares.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 SUMCUAD.HTM fsquares.js fsquaresS.wasm
java -jar %compilerName% -D app=2 %compilerOptions% --js interface.js --js_output_file fsquares.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 FCUBES.HTM fsquares.js fcubesE.wasm
java -jar %compilerName% -D app=3 %compilerOptions% --js interface.js --js_output_file fsquares.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 SUMCUBOS.HTM fsquares.js fcubesS.wasm
java -jar %compilerName% -D app=4 %compilerOptions% --js interface.js --js_output_file fsquares.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 CONTFRAC.HTM fsquares.js contfracE.wasm
java -jar %compilerName% -D app=5 %compilerOptions% --js interface.js --js_output_file fsquares.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 FRACCONT.HTM fsquares.js contfracS.wasm

java -jar %compilerName% %compilerOptions% --js intfwebw.js --js_output_file intWW.js

java -jar %compilerName% -D lang=0 %compilerOptions% --js polyfact.js --js_output_file polfactE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 POLFACT.HTM polfactE.js polfactE.wasm
java -jar %compilerName% -D lang=1 %compilerOptions% --js polyfact.js --js_output_file polfactS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 FACTPOL.HTM polfactS.js polfactS.wasm

java -jar %compilerName% -D lang=0 %compilerOptions% --js dislog.js --js_output_file dilogE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 DILOG.HTM dilogE.js dilogE.wasm
java -jar %compilerName% -D lang=1 %compilerOptions% --js dislog.js --js_output_file dilogS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 LOGDI.HTM dilogS.js dilogS.wasm

java -jar %compilerName% -D lang=0 %compilerOptions% --js quadrmod.js --js_output_file quadmodE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 QUADMOD.HTM quadmodE.js quadmodE.wasm
java -jar %compilerName% -D lang=1 %compilerOptions% --js quadrmod.js --js_output_file quadmodS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 CUADMOD.HTM quadmodS.js quadmodS.wasm

java -jar %compilerName% -D lang=0 %compilerOptions% --js gauss.js --js_output_file gaussianE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 GAUSSIAN.HTM gaussianE.js gaussianE.wasm
java -jar %compilerName% -D lang=1 %compilerOptions% --js gauss.js --js_output_file gaussianS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 GAUSIANO.HTM gaussianS.js gaussianS.wasm

java -jar %compilerName% -D lang=0 %compilerOptions% --js ecmfront.js --js_output_file ecmE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 ECM.HTM ecmE.js ecmE.wasm
java -jar %compilerName% -D lang=1 %compilerOptions% --js ecmfront.js --js_output_file ecmS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 ECMC.HTM ecmS.js ecmS.wasm
copy calculatorSW.js calcSW.js

set emcc_DEBUG=
java -jar %compilerName% %compilerOptions% --js ecmfwebw.js --js_output_file ecmWW.js

java -jar %compilerName% -D lang=0 %compilerOptions% --js quadr.js --js_output_file quadE.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 QUAD.HTM quadE.js quadE.wasm
java -jar %compilerName% -D lang=1 %compilerOptions% --js quadr.js --js_output_file quadS.js
if errorlevel 1 goto end
perl replaceEmbeddedJS.pl %1 CUAD.HTM quadS.js quadS.wasm

java -jar %compilerName% --compilation_level WHITESPACE_ONLY --js dist.js --js_output_file dist%1.js

perl csp.pl
:end