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

function localStorage_setItem(name, data)
{
  localStorage.setItem(name, data);
}

function localStorage_getItem(name)
{
  return localStorage.getItem(name);
}

function styleButtons(style1, style2)
{
  document_getElementById('eval').style.display = style1;
  document_getElementById('factor').style.display = style1;
  document_getElementById('stop').style.display = style2;
  document_getElementById('more').style.display = style2;
}

function restartFactorization(type)
{
  document_getElementById('modal').style.display = "none";
  worker.terminate();
  worker = 0;
  dowork(type);
}

function callWorker(param)
{
  if (!worker)
  {
  	worker = new Worker('ecmW.js?0412');
    worker.onmessage = function(e)
    { // First character of e.data is '1' for intermediate text
      // and it is '2' for end of calculation.
	  // It is '9' for saving expression to factor into Web Storage.
      var firstChar = e.data.substring(0, 1);
      if (firstChar == '8')
      {
        localStorage_setItem("ecmFactors", e.data.substring(1));
        localStorage_setItem("ecmCurve", "");
      }
      else if (firstChar == '7')
      {
        localStorage_setItem("ecmCurve", e.data.substring(1));
      }
      else if (firstChar == '4')
      {
        document_getElementById('status').innerHTML = e.data.substring(1);
      }
      else
      {
        document_getElementById('result').innerHTML = e.data.substring(1);
        if (e.data.substring(0, 1) == '2')
        {   // First character passed from web worker is '2'.
          document_getElementById('status').innerHTML = "";
          styleButtons("inline", "none");  // Enable eval and factor
          document_getElementById('modal').style.display = "none";
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
  var charNull = String.fromCharCode(0);
  document_getElementById('help').style.display = "none";
  res.style.display = "block";
  if (valueText == "")
  {    // Nothing in input box.
    res.innerHTML = (app & 1 ? "<p>Por favor ingrese una expresión.</p>" :
                               "<p>Please type an expression.</p>");
    return;
  }
  if (typeof(Worker) === "undefined")
  {    // Web workers not supported on this browser.
    res.innerHTML = (app & 1 ? "<p>Esta calculadora necesita Web Workers. Por favor use otro navegador Web.</p>" :
                               "<p>This calculator requires Web Workers. Please use another Web browser.</p>");
    return;							   
  }
  styleButtons("none", "inline");  // Enable "more" and "stop" buttons
  res.innerHTML = (app & 1 ? "<p>Factorizando la expresión...</p>" :
                             "<p>Factoring expression...</p>");
  if (n < -2)
  {
	app += 6;   // Convert to factorization.
  }
  param = digitGroup + ',' + app + ',' + valueText + charNull +
          localStorage_getItem("ecmFactors");
  if (n == -1 || n == -2)
  {
    param += charNull + document_getElementById('curve').value;   // Append new factor or curve number.
  }
  if (n == -3 || n == -4)
  {
    param += '*' + document_getElementById('curve').value + '^1(2)';   // Append new factor or curve number.
  }
  callWorker(param + charNull);
}

window.onload = function ()
{
  var param;
  document_getElementById('eval').onclick = function ()
  {
    localStorage.setItem('ecmFactors','');
    dowork(document_getElementById('batch').checked? 4: 0);
  }
  document_getElementById('factor').onclick = function ()
  {
    localStorage.setItem('ecmFactors','');
    dowork(document_getElementById('batch').checked? 6: 2);
  }
  document_getElementById('more').onclick = function ()
  {
    document_getElementById('modal').style.display = "block";
  }
  document_getElementById('close').onclick = function ()
  {
    document_getElementById('modal').style.display = "none";
  }
  document_getElementById('ncurve').onclick = function ()
  {
	restartFactorization(-2);
  }
  document_getElementById('nfactor').onclick = function ()
  {
	restartFactorization(-4);
  }
  document_getElementById('curve').onkeypress = function(event)
  {
    return (event.charCode == 8 || event.charCode == 0) ? null : event.charCode >= 48 && event.charCode <= 57;
  }
  document_getElementById('stop').onclick = function ()
  {
    worker.terminate();
    worker = 0;
	styleButtons("inline", "none");  // Enable eval and factor
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
  window.onclick = function(event)
  {
	var modal = document_getElementById('modal');
    if (event.target == modal)
    {
      modal.style.display = "none";
    }
  }
  ecmFactor = localStorage_getItem("ecmFactors");
  if (ecmFactor)
  {          // Continue factoring.
    document_getElementById('value').value = ecmFactor.slice(0,ecmFactor.indexOf('='));
	document_getElementById('curve').value = localStorage.getItem("ecmCurve");
    dowork(-2);
	document_getElementById('curve').value = "";
  }
}
