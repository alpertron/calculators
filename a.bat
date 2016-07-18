java -jar "C:\Program Files\Emscripten\emscripten\1.35.0\third_party\closure-compiler\compiler.jar" --compilation_level ADVANCED_OPTIMIZATIONS --js interface.js --js_output_file fsquares.js
cmd /c emcc -Os expression.c bigint.c division.c fsquares.c fcubes.c baseconv.c karatsuba.c modmult.c sqroot.c contfrac.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork']" -s TOTAL_MEMORY=33554432 -s NO_BROWSER=1 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o fsquaresW.js

java -jar "C:\Program Files\Emscripten\emscripten\1.35.0\third_party\closure-compiler\compiler.jar" --compilation_level ADVANCED_OPTIMIZATIONS --js polyfact.js --js_output_file polfact.js
cmd /c emcc -Os expression.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c polynomial.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork']" -s TOTAL_MEMORY=33554432 -s NO_BROWSER=1 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o polfactW.js

java -jar "C:\Program Files\Emscripten\emscripten\1.35.0\third_party\closure-compiler\compiler.jar" --compilation_level ADVANCED_OPTIMIZATIONS --js dislog.js --js_output_file dilog.js
cmd /c emcc -Os expression.c bigint.c division.c baseconv.c karatsuba.c modmult.c sqroot.c factor.c dilog.c --llvm-lto 1 --js-library lib.js --pre-js pre.js -s EXPORTED_FUNCTIONS="['_doWork']" -s TOTAL_MEMORY=33554432 -s NO_BROWSER=1 -s NO_FILESYSTEM=1 --closure 1 --memory-init-file 0 -o dilogW.js
