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
(function(global)
{   // This method separates the name space from the Google Analytics code.
var worker = 0;
var app;
function document_getElementById(x)
{
  return document.getElementById(x);
}
function callWorker(param)
{
  if (!worker)
  {
  	worker = new Worker('quadmodW0020.js');
	worker.onmessage = function(e)
	{ // First character of e.data is '1' for intermediate text
      // and it is '2' for end of calculation.
	  document_getElementById('result').innerHTML = e.data.substring(1);
	  if (e.data.substring(0, 1) == '2')
	  {   // First character passed from web worker is '2'.
	    document_getElementById('solve').disabled = false;
	    document_getElementById('stop').disabled = true;
      }
	}
  }
  worker.postMessage(param);
}

function dowork(n)
{
  var param;
  var app = parseInt(document_getElementById('app').value) + n;
  var res = document_getElementById('result');
  var quadrText = document_getElementById('quad').value.trim();
  var linText = document_getElementById('lin').value.trim();
  var constText = document_getElementById('const').value.trim();
  var modText = document_getElementById('mod').value.trim();
  var digitGroup = document_getElementById('digits').value.trim();
  document_getElementById('help').style.display = "none";
  res.style.display = "block";
  var missing = "";
  if (quadrText == "")
  {
    missing = (app & 1 ? "coeficiente cuadrático." : "quadratic coefficient.");
  }
  if (linText == "")
  {
    missing = (app & 1 ? "coeficiente lineal." : "linear coefficient.");
  }
  if (constText == "")
  {
    missing = (app & 1 ? "término independiente." : "constant coefficient.");
  }
  if (modText == "")
  {
    missing = (app & 1 ? "módulo." : "modulus.");
  }
  if (missing != "")
  {
    res.innerHTML = (app & 1 ? "Por favor ingrese un número o expresión para el "+missing :
                               "Please type a number or expression for the "+missing);
    return;
  }
  document_getElementById('solve').disabled = true;
  document_getElementById('stop').disabled = false;
  res.innerHTML = (app & 1 ? "Resolviendo la ecuación cuadrática..." :
                             "Solving the quadratic equation...");
  param = digitGroup + ',' + app + ',' + quadrText + String.fromCharCode(0) + linText +
    String.fromCharCode(0) + constText +String.fromCharCode(0) + modText + String.fromCharCode(0);
  callWorker(param);
}

window.onload = function ()
{
  var param;
  document_getElementById('stop').disabled = true;
  document_getElementById('solve').onclick = function ()
  {
    dowork(0);
  }
  document_getElementById('stop').onclick = function ()
  {
    worker.terminate();
    worker = 0;
    document_getElementById('solve').disabled = false;
    document_getElementById('stop').disabled = true;
    document_getElementById('result').innerHTML = 
      (app & 1 ? "<p>Cálculo detenido por el usuario.</p>" :
                 "<p>Calculation stopped by user</p>");
  }
  document_getElementById('helpbtn').onclick = function ()
  {
    document_getElementById('help').style.display = "block";
    document_getElementById('result').style.display = "none";
  }
  if ('serviceWorker' in navigator)
  { // Attempt to register service worker.
    // There is no need to do anything on registration success or failure in this JavaScript module.
    navigator.serviceWorker.register('calcSW.js').then(function() {}, function() {});
  }
}
})(this);

if (typeof(window) !== "undefined")
{   // In main thread: register Google Analytics.
  addEventListener("load", function ()
  {
    (function(i,s,o,g,r,a,m){i["GoogleAnalyticsObject"]=r;i[r]=i[r]||function(){
      (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
      m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
      })(window,document,"script","https://www.google-analytics.com/analytics.js","ga");
  
    ga("create", "UA-4438475-1", "auto");
    ga("send", "pageview");
  });
}

