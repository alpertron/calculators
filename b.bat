java -jar closure-compiler.jar ^
    --js='useBlockly.js' ^
    --js='.\Blockly\blocks\**.js' ^
    --js='.\Blockly\core\**.js' ^
    --js='.\Blockly\externs\svg-externs.js' ^
    --compilation_level ADVANCED_OPTIMIZATIONS ^
    --dependency_mode=PRUNE --entry_point=useBlockly ^
    --generate_exports ^
    --isolation_mode IIFE ^
    --jscomp_off visibility ^
    --js_output_file blockly.js

rem    --compilation_level SIMPLE_OPTIMIZATIONS ^
