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
/* global fileContents */
/* global get */
function setStorage(name, data)
{
  Android.setStorage(name, data);
}

function getStorage(name)
{
  return Android.getStorage(name);
}

function getCalculatorCode(fileName, workerParam)
{
  fileContents = 1;
}

function registerServiceWorker()
{
  get("appver").value = Android.getAppVer();
}

function endCalculation()
{
  Android.endCalculation();
}
