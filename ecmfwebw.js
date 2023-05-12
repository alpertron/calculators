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
/* global convertToString */
/* global ptrToString */
let exports, HEAPU8, wasmLoaded;

let info =
{
  "env":
{
  "databack": function(data)
  {
    self.postMessage(ptrToString(data));
  },
  "startSkipTest": function()
  {
    self.postMessage("51");   // Show Skip Test button on screen
  },
  "endSkipTest": function()
  {
    self.postMessage("52");   // Hide Skip Test button from screen
  },
  "tenths": function()
  {
    return Math.floor(new Date().getTime() / 100);
  },
  "getCunn": function(data)
  {
    try
    {
      let req = new XMLHttpRequest();
      // Web worker protocol is blob:, so we need to change that to https: as appropriate.
      req.open("GET", "https://www.alpertron.com.ar/"+ptrToString(data), false);
      req.send(null);
      if (req.status === 200)
      {
        convertToString(exports["getFactorsAsciiPtr"](), req.responseText);
      }
    }
    catch (ex)
    {    // Nothing to do if the factors could not be retrieved.
    }
  }
}
};  

self.onmessage = function (e)
{
  if (typeof(WebAssembly) === "undefined")
  {
    return;
  }
  if (wasmLoaded)
  {
    convertToString(exports["getInputStringPtr"](), e.data[0]);
    exports["doWork"]();
    return;  
  }
  WebAssembly["instantiate"](e.data[1], info).then(function(results)
  {
    wasmLoaded = 1;
    exports = results["instance"]["exports"];
    HEAPU8 = new Uint8Array(exports["memory"]["buffer"]);
    convertToString(exports["getInputStringPtr"](), e.data[0]);
    exports["doWork"]();
  },
  function()
  {
  });
};
