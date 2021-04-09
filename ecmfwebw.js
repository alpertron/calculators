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
var exports, HEAPU8, wasmLoaded;

function ptrToString(ptr)
{
  var t=-1;
  var i = 0;
  var str="", outString="";
  do
  {
    for (i=0; i<1024; i++)
    {
      t = HEAPU8[((ptr++)>>0)];
      if (t === 0)
      {
        break;
      }
      if (t>=128)
      {
        t = ((t-192)<<6) + HEAPU8[((ptr++)>>0)] - 128;
      }
      str += String.fromCharCode(t);
    }
    outString += str;
    str = "";
  } while (t !== 0);
  outString += str;
  return outString;
}

function convertToString(ptr, str)
{
  var dest = ptr;
  var length = str.length;
  var i, t;
  for (i=0; i<length; i++)
  {
    t = str.charCodeAt(i);
    if (t<128)
    {
      HEAPU8[(dest++) >> 0] = t;
    }
    else
    {
      HEAPU8[(dest++) >> 0] = (t >> 6) + 192;
      HEAPU8[(dest++) >> 0] = (t & 63) + 128;
    }
  }
  HEAPU8[dest >> 0] = 0;
}

var info =
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
    var req = new XMLHttpRequest();
    // Web worker protocol is blob:, so we need to change that to https: as appropriate.
    req.open('GET', "https://www.alpertron.com.ar/"+ptrToString(data), false);
    req.send(null);
    if (req.status === 200)
    {
      convertToString(exports["getFactorsAsciiPtr"](), req.responseText);
    }
  }
}
};  

self.onmessage = function (e)
{
  var request;
  var afterKey;
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
    return;
  });
};
