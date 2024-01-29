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
/* global hide */
/* global fileName */
/* global tofile */
function loadPolyCalc(firstChar, value)
{
}

function comingFromPolfact(value)
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

function showVersion(lang)  
{
  return "";
}

function downloadResult()
{
  hide("savefile");
  Android.download(fileName, tofile);
}