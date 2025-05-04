/*
    This file is part of Alpertron Calculators.

    Copyright 2023 Dario Alejandro Alpern

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
/* global Android */
/* global asmGetInformation */
/* global asmDrawPartialGraphic */
/* global asmMoveGraphic */
/* global asmNbrChanged */
/* global bitsCanvas */
/* global buffer */
/* global getPixels */
/* global imgData */
/* global pixels */
/* global updateGraphic */

let startOffset;
function startLowLevelCode()
{
  let length, bytes;
  let info;
  let getPixels;
  let asmJSbuffer;
  let bufSize = 33554432;
  startOffset = 11000000|0;
  imgData = canvas.getContext("2d").createImageData(2048, 4096);  // 32 MB;
  buffer = imgData.data.buffer;
  bitsCanvas = new Uint8Array(buffer);
  if (typeof(WebAssembly) === "undefined")
  {                                      // Asm.js initialization.
    asmJSbuffer = new ArrayBuffer(bufSize);
    HEAPU8 = new Uint8Array(asmJSbuffer);    // Reserve 32 MB for asm.js variables and buffers.
    env = {"a": {"buffer": asmJSbuffer},
      "abort": function(_q)
               {
                 /* Not used*/
               },
    };
    asm = instantiate(env);  // Link asm.js module.
    getPixels = startLowLevelCodeCallback(asm);
    pixels = HEAPU8.subarray(getPixels());
    updateGraphic(center, 1);
  }
  else
  {                                      // WebAssembly initialization.
    wasm = get("wasmb64").text;
    while (wasm.charCodeAt(0) < 32)
    {
      wasm = wasm.substring(1);
    }    
    while (wasm.charCodeAt(wasm.length-1) < 32)
    {
      wasm = wasm.substring(0, wasm.length-1);
    }
    length = wasm.length * 3 / 4;
    if (wasm.charCodeAt(wasm.length - 1) === 61)
    {                                    // Base64 ending equal sign found.
      length--;
    }
    if (wasm.charCodeAt(wasm.length - 2) === 61)
    {                                    // Another base64 ending equal sign found.
      length--;
    }
    bytes = new Int8Array(length);
    // Decode Base64.
    b64decode(wasm, bytes);
    info = {};
    WebAssembly["instantiate"](bytes, info).then(function(results)
    {
      asm = results["instance"]["exports"];
      asmGetInformation = asm["getInformation"];
      asmDrawPartialGraphic = asm["drawPartialGraphic"];
      asmMoveGraphic = asm["moveGraphic"];
      asmNbrChanged = asm["nbrChanged"];
      
      HEAPU8 = new Uint8Array(asm["memory"]["buffer"]);
      pixels = HEAPU8.subarray(asm["getPixels"]());
      updateGraphic(center, 1);
    });
  }
}

function updateImageData(top, left, width, height)
{
  let rowNbr;
  let leftPixelSrc = top*8192+left*4;
  let leftPixelDest = (top+32)*8192+left*4;

  for (rowNbr=0; rowNbr<height; rowNbr++)
  {
    bitsCanvas.set(pixels.subarray(leftPixelSrc, leftPixelSrc + 4*width), leftPixelDest);
    leftPixelSrc += 8192;
    leftPixelDest += 8192;
  }
}

function copyStr(startOffset, str)
{
  let idx;
  for (idx = 0; idx < str.length; idx++)
  {
    HEAPU8[(startOffset+idx) >> 0] = str.charCodeAt(idx);
  }
  HEAPU8[(startOffset+idx) >> 0] = 0;
}

function drawPartialGraphic(xminDisp, xmaxDisp, yminDisp, ymaxDisp)
{
  asmDrawPartialGraphic(xminDisp, xmaxDisp, yminDisp, ymaxDisp, 0);
}

function setNewDimensionsForCanvas()
{
  let newDomRect = canvas.getBoundingClientRect();
  canvas.width = newDomRect.width;
  canvas.height = newDomRect.height;
}
