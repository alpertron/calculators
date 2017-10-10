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
    worker = new Worker('polfactW0014.js');
    worker.onmessage = function(e)
    { // First character of e.data is '1' for intermediate text
      // and it is '2' for end of calculation.
      document_getElementById('result').innerHTML = e.data.substring(1);
      if (e.data.substring(0, 1) == '2')
      {   // First character passed from web worker is '2'.
        document_getElementById('eval').disabled = false;
        document_getElementById('factor').disabled = false;
        document_getElementById('stop').disabled = true;
      }
    }
  }
  worker.postMessage(param);
}

function dowork(n)
{
  var app = parseInt(document_getElementById('app').value) + n + 
    (document_getElementById('sup').checked ? 4 : 0);
  var res = document_getElementById('result');
  var polyText = document_getElementById('poly').value;
  var modText = document_getElementById('mod').value;
  var digitGroup = document_getElementById('digits').value;
  document_getElementById('help').style.display = "none";
  res.style.display = "block";
  if (polyText == "")
  {
    res.innerHTML = (app & 1 ? "Por favor ingrese una expresión para el polinomio a evaluar." :
                               "Please type an expression for the polynomial to evaluate.");
    return;
  }
  if (modText == "")
  {
    res.innerHTML = (app & 1 ? "Por favor ingrese un número o expresión para el módulo." :
                               "Please type a number or expression for the modulus.");
    return;
  }
  document_getElementById('eval').disabled = true;
  document_getElementById('factor').disabled = true;
  document_getElementById('stop').disabled = false;
  res.innerHTML = (app & 1 ? "Factorizando el polinomio..." :
                             "Factoring polynomial...");
  param = digitGroup + ',' + app + ',' + modText + String.fromCharCode(0) + polyText +
  String.fromCharCode(0);
  callWorker(param);
}

window.onload = function ()
{
  var param;
  document_getElementById('stop').disabled = true;
  document_getElementById('eval').onclick = function ()
  {
    dowork(2);
  }
  document_getElementById('factor').onclick = function ()
  {
    dowork(0);
  }
  document_getElementById('stop').onclick = function ()
  {
    worker.terminate();
    worker = 0;
    document_getElementById('eval').disabled = false;
    document_getElementById('factor').disabled = false;
    document_getElementById('stop').disabled = true;
    document_getElementById('result').innerHTML = 
      (app & 1 ? "<p>Factorización detenida por el usuario.</p>" :
                 "<p>Factorization stopped by user</p>");
  }
  document_getElementById('helpbtn').onclick = function ()
  {
    document_getElementById('help').style.display = "block";
    document_getElementById('result').style.display = "none";
  }
  document_getElementById('poly').oninput = function ()
  {
    var input = document_getElementById('poly');
    var loc = input.value.length - input.selectionStart;
    input.value = input.value.replace(".", "x^");
	setTimeout(function()
	{
      loc = input.value.length - loc;
      input.selectionStart = loc;
      input.selectionEnd = loc;
	}, 30);		
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


