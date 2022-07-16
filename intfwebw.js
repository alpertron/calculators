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
let exports, HEAPU8, wasmLoaded;

function ptrToString(ptr)
{
  let t=-1;
  let i = 0;
  let str="", outString="";
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
  let dest = ptr;
  let length = str.length;
  let i, t;
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
  console.log("onmessage ecmfwebw 111: " + e.origin);
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
    try
    {
      exports["doWork"]();
    } catch (err)
    {
      self.postMessage("2<p>"+err.message+"</p>");
    }
  });
};
