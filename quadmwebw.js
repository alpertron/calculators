/*
    This file is part of Alpertron Calculators.

    Copyright 2018 Dario Alejandro Alpern

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
/* global convertToString */
/* global ptrToString */
let exports, HEAPU8, wasmLoaded;

let env =
{
  "databack": function(data)
  {
    self.postMessage(ptrToString(data));
  },
  "tenths": function()
  {
    return Math.floor(new Date().getTime() / 100);
  }
};

let info =
{
  "env": env
};  

self.onmessage = function (e)
{
  if (wasmLoaded)
  {
    convertToString(exports["getInputStringPtr"](), e.data[0]);
    exports["doWork"]();
  }
  WebAssembly["instantiate"](e.data[1], info).then(function(results)
  {
    wasmLoaded = 1;
    exports = results["instance"]["exports"];
    HEAPU8 = new Uint8Array(exports["memory"]["buffer"]);
    convertToString(exports["getInputStringPtr"](), e.data[0]);
    exports["doWork"]();
  }, function() {});
};
