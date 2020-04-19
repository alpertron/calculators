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
var exports, HEAPU8, wasmLoaded;
var env =
{
  "databack": function(data)
  {
    self.postMessage(PtrToString(data));
  },
  "tenths": function()
  {
    return Math.floor(new Date().getTime() / 100);
  }
};

var info =
{
  "env": env
};  

self.onmessage = function (e)
{
  var request;
  if (wasmLoaded)
  {
    ConvertToString(exports["getInputStringPtr"](), e.data[0]);
    exports["doWork"]();
    return;  
  }
  WebAssembly["instantiate"](e.data[1], info).then(function(results)
  {
    wasmLoaded = 1;
    exports = results["instance"]["exports"];
    HEAPU8 = new Uint8Array(exports["memory"]["buffer"]);
    ConvertToString(exports["getInputStringPtr"](), e.data[0]);
    exports["doWork"]();
    return;
  });
}

function PtrToString(ptr)
{
  var t=-1;
  var i = 0;
  var str="", outString="";
  do
  {
    for (i=0; i<1024; i++)
    {
      t = HEAPU8[((ptr++)>>0)];
      if (t==0)
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
  } while (t!=0);
  outString += str;
  return outString;
}

function ConvertToString(ptr, str)
{
  var dest = ptr;
  var length = str.length;
  var i, t;
  for (i=0; i<length; i++)
  {
    t = str.charCodeAt(i);
    if (t<128)
    {
      HEAPU8[dest++] = t;
    }
    else
    {
      HEAPU8[dest++] = (t >> 6) + 192;
      HEAPU8[dest++] = (t & 63) + 128;
    }
  }
  HEAPU8[dest] = 0;
}
