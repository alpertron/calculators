java -jar closure-compiler.jar --js='useBlockly.js' ^
    --js='.\Blockly\blocks\**.js' ^
    --js='.\Blockly\core\**.js' ^
    --js='.\Blockly\generators\**.js' ^
    --compilation_level SIMPLE_OPTIMIZATIONS ^
    --dependency_mode=PRUNE --entry_point=useBlockly ^
    --generate_exports ^
    --isolation_mode IIFE ^
    --jscomp_off visibility ^
    --js_output_file blockly.js
