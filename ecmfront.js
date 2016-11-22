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
var worker = 0;
var app;
var oldWhiteSpace;
var oldOverflowX;
function document_getElementById(x)
{
  return document.getElementById(x);
}
function callWorker(param)
{
  if (!worker)
  {
  	worker = new Worker('ecmW.js?1911');
    worker.onmessage = function(e)
    { // First character of e.data is '1' for intermediate text
      // and it is '2' for end of calculation.
      if (e.data.substring(0, 1) == '4')
      {
        document_getElementById('status').innerHTML = e.data.substring(1);
      }
      else
      {
        document_getElementById('result').innerHTML = e.data.substring(1);
        if (e.data.substring(0, 1) == '2')
        {   // First character passed from web worker is '2'.
          document_getElementById('status').innerHTML = "";
          document_getElementById('eval').disabled = false;
          document_getElementById('factor').disabled = false;
          document_getElementById('stop').disabled = true;
        }
      }
    }
  }
  worker.postMessage(param);
}

function dowork(n)
{
  var app = parseInt(document_getElementById('app').value) + n;
  var res = document_getElementById('result');
  var valueText = document_getElementById('value').value;
  var digitGroup = document_getElementById('digits').value;
  document_getElementById('help').style.display = "none";
  res.style.display = "block";
  if (valueText == "")
  {
    res.innerHTML = (app & 1 ? "Por favor ingrese una expresión." :
                               "Please type an expression.");
    return;
  }
  document_getElementById('eval').disabled = true;
  document_getElementById('factor').disabled = true;
  document_getElementById('stop').disabled = false;
  res.innerHTML = (app & 1 ? "Factorizando la expresión..." :
                             "Factoring expression...");
  param = digitGroup + ',' + app + ',' + valueText + String.fromCharCode(0)
  callWorker(param);
}

window.onload = function ()
{
  var param;
  document_getElementById('stop').disabled = true;
  document_getElementById('eval').onclick = function ()
  {
    dowork(0);
  }
  document_getElementById('factor').onclick = function ()
  {
    dowork(2);
  }
  document_getElementById('stop').onclick = function ()
  {
    worker.terminate();
    worker = 0;
    document_getElementById('eval').disabled = false;
    document_getElementById('factor').disabled = false;
    document_getElementById('stop').disabled = true;
    document_getElementById('result').innerHTML = 
      (app & 1 ? "<p>Cálculo detenido por el usuario.</p>" :
                 "<p>Calculation stopped by user</p>");
    document_getElementById('status').innerHTML = "";
  }
  document_getElementById('helpbtn').onclick = function ()
  {
    document_getElementById('help').style.display = "block";
    document_getElementById('result').style.display = "none";
  }
  document_getElementById('batch').onchange = function ()
  {
    var entry = document_getElementById('entry');
    if (document_getElementById('batch').checked)
    {
      value = document_getElementById("value");
      oldWhiteSpace = value.style.whiteSpace;
      oldOverflowX = value.style.overflowX;
      entry.innerHTML = '<textarea id="value" rows="5" class="input" placeholder="One numerical expression or loop per line"></textarea>';
      value = document_getElementById("value");
      value.style.whiteSpace = "nowrap";
      value.style.overflowX = "scroll";
    }
    else
    {
      entry.innerHTML = '<input type="text" id="value" value="" placeholder="Number or numerical expression" class="input"/>';
      value = document_getElementById("value");
      value.style.whiteSpace = oldWhiteSpace;
      value.style.overflowX = oldOverflowX;
    }
  }

}

