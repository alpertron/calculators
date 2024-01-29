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
/* global dowork */
/* global fileName */
/* global getStorage */
/* global getVersionText */
/* global lang */
/* global hide */
/* global newCurveOrFactor */
/* global resultDirty */
/* global resultText */
/* global tofile */
function loadPolyCalc(firstChar, value)
{
  window.sessionStorage.setItem((firstChar === "M"? "F": "E"), value);
  window.location.replace(lang? "FACTPOL.HTM": "POLFACT.HTM");
}

function comingFromPolfact(value)
{
  let fromPolfact = window.sessionStorage.getItem("F");
  if (fromPolfact != null)
  {    // Number to factor coming from polynomial factorization calculator.
    window.sessionStorage.removeItem("F");
    value.value = fromPolfact;
    dowork(-2);    // Perform factorization.
  }
  fromPolfact = window.sessionStorage.getItem("E");
  if (fromPolfact != null)
  {    // Number to factor coming from polynomial factorization calculator.
    window.sessionStorage.removeItem("E");
    value.value = fromPolfact;
    dowork(0);     // Perform evaluation.
  }
  else
  {
    let search = window.location.search;
    if (search.startsWith("?q="))
    {
      value.value = decodeURIComponent(search.substring(3)).replace(/\{/g, "(").replace(/\}/g, ")");
      dowork(-2);
    }
    else
    {
      let ecmFactor = getStorage("ecmFactors");
      if (ecmFactor)
      {          // Continue factoring.
        value.value = ecmFactor.slice(0,ecmFactor.indexOf("="));
        newCurveOrFactor.value = getStorage("ecmCurve");
        dowork(-2);
        newCurveOrFactor.value = "";
      }
    }
  }
}

function showVersion(lang)  
{
  if ((typeof(Worker) === "undefined"))
  {    // Web workers not supported on this browser.
    resultDirty = true;
    resultText = (lang ? "<p>Esta calculadora necesita Web Workers. Por favor use otro navegador Web.</p>" :
                         "<p>This calculator requires Web Workers. Please use another Web browser.</p>");
    return null;
  }
  return getVersionText();
}

function downloadResult()
{
  hide("savefile");
  let fileBlob = new Blob([tofile], { type: "text/plain" });
  let fileUrl = URL.createObjectURL(fileBlob);
  let a = document.createElement("a");
  a.href = fileUrl;
  a.download = fileName;
  let clickHandler = function()
  {
    setTimeout(function()
    {
      URL.revokeObjectURL(fileUrl);
      this.removeEventListener("click", clickHandler);
    },
    150);
  };
  a.addEventListener("click", clickHandler, false);
  a.click();
}