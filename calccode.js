/*
    This file is part of Alpertron Calculators.

    Copyright 2015 Dario Alejandro Alpern

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
/* global calcURLs */
/* global fileContents */
let workPar;
function getCalculatorCode(jsFileName, workerParameter)
{
  workPar = workerParameter;
  if (asmjs)
  {
    let req = new XMLHttpRequest();
    req.open("GET", jsFileName, true);
    req.responseType = "arraybuffer";
    req.onreadystatechange = function (_aEvt)
    {
      if (req.readyState === 4 && req.status === 200)
      {
        fileContents = /** @type {ArrayBuffer} */ (req.response);
        if (workPar)
        {
          callWorker(workPar);
        }
      }
    };
    req.send(null);
  }
  else
  {
    let wasm = get("wasmb64").text;
    while (wasm.charCodeAt(0) < 32)
    {
      wasm = wasm.substring(1);
    }    
    while (wasm.charCodeAt(wasm.length-1) < 32)
    {
      wasm = wasm.substring(0, wasm.length-1);
    }    
    let length = wasm.length*3/4;
    if (wasm.charCodeAt(wasm.length-1) === 61)
    {
      length--;
    }
    if (wasm.charCodeAt(wasm.length-2) === 61)
    {
      length--;
    }
    fileContents=new Int8Array(length);
    b64decode(wasm, fileContents);
    calcURLs.shift();  // Do not fetch Javascript file that will not be used.
  }
}
