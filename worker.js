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
/* global comingFromWorker */
/* global fileContents */
/* global get */
let blob;
let worker = null;
function endWorker()
{
  if (worker != null)
  {
    worker.terminate();
    worker = null;
  }
}

function callWorker(param)
{
  if (worker == null)
  {
    if (!blob)
    {
      if (typeof(WebAssembly) === "undefined")
      {    // Asm.js
        blob = new Blob([fileContents],{type: "text/javascript"});
      }
      else
      {    // WebAssembly
        blob = new Blob([get("worker").textContent],{type: "text/javascript"});
      }
    }   
    worker = new Worker(window.URL.createObjectURL(blob));
    worker.onmessage = comingFromWorker;
  }
  if (typeof(WebAssembly) === "undefined")
  {      // Asm.js
    worker.postMessage(param);
  }
  else
  {      // WebAssembly.
    worker.postMessage([param, fileContents]);
  }
}

function getVersionText()
{
  let langName = (typeof(WebAssembly) === "undefined")? "asm.js": "WebAssembly";
  return (lang ? " Esta es la versión "+langName+".</p>":
                 " This is the "+langName+" version.</p>");
}