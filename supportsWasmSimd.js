/*
    This file is part of Alpertron Calculators.

    Copyright 2026 Dario Alejandro Alpern

    Alpertron Calculators is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Alpertron Calculators is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
*/
function supportsWasmSimd()
{
  // Minimal WASM module containing a v128 SIMD instruction.
  const simdModule = new Uint8Array(
  [
    0x00, 0x61, 0x73, 0x6d, // \0asm
    0x01, 0x00, 0x00, 0x00, // WASM version 1
  
    // Type section
    0x01, 0x05,
    0x01,
    0x60, 0x00, 0x01, 0x7b, // () -> v128
  
    // Function section
    0x03, 0x02,
    0x01, 0x00,
  
    // Code section
    0x0a, 0x0a,
    0x01,
    0x08,
    0x00,
    0x41, 0x00,             // i32.const 0
    0xfd, 0x0f,             // i8x16.splat
    0xfd, 0x62,             // i8x16.popcnt
    0x0b                    // end
  ]);

  return typeof(WebAssembly) === "object" &&
         WebAssembly["validate"](simdModule);
}
