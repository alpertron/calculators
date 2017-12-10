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
function callWorker(param)
{
  if (!worker)
  {
    worker = new Worker('fsquaresW0026.js');
    worker.onmessage = function(e)
    {
      document.getElementById('result').innerHTML = e.data;
    }
  }
  worker.postMessage(param);
}

window.onload = function()
{
  var param;
  if (document.getElementById('input'))
  {
    document.getElementById('input').onkeypress = function(e)
    {
      var output, res;
	  var digitGroup = document.getElementById('digits').value;
      app = document.getElementById('app').value;
      res = document.getElementById('result');
      res.style.display = "block";
      var input = document.getElementById('input').value;
      if (!e) e = window.event;
      var keyCode = e.keyCode || e.which;
      if (keyCode == 13)
	  {  // Used pressed Enter key
	    output = document.getElementById('output')
	    if (input == "")
	    {
	      res.innerHTML = (app & 1 ? "Por favor ingrese un número o expresión." : "Please type a number or expression.");
          return;
	    }
	    if (app==0)
		{
	      res.innerHTML = "Computing sum of squares...";
		}
		else if (app==1)
		{
		  res.innerHTML = "Calculando suma de cuadrados...";
		}
		else if (app==2)
		{
		  res.innerHTML = "Computing sum of cubes...";
		}
		else
		{
		  res.innerHTML = "Calculando suma de cubos...";
		}
	    param = digitGroup + ',' + app + ',' + input + String.fromCharCode(0);
		callWorker(param);
	  }
    }
  }
  document.getElementById('calc').onclick = function()
  {
    var app, res, valueA, valueB, valueC, digitGroup;
    app = parseInt(document.getElementById('app').value);
    res = document.getElementById('result');
	res.style.display = "block";
	valueA = document.getElementById('num').value;
	if (valueA == "")
	{
	  if (app >= 4)
	  {
	    res.innerHTML = (app & 1 ? "Por favor ingrese un número o expresión para el numerador." :
	                            "Please type a number or expression for numerator.");
	  }
	  else
	  {
	    res.innerHTML = (app & 1 ? "Por favor ingrese un número o expresión." :
	                            "Please type a number or expression.");
	  }
	  return;
	}
	if (app >= 4)
	{
	  valueB = document.getElementById('delta').value;
	  if (valueB == "")
	  {
	    res.innerHTML = (app & 1 ? "Por favor ingrese un número o expresión para el argumento de la raíz cuadrada." :
                              "Please type a number or expression for square root argument.");
	    return;
	  }
	  valueC = document.getElementById('den').value;
	  if (valueC == "")
	  {
	    res.innerHTML = (app & 1 ? "Por favor ingrese un número o expresión para el denominador." :
                              "Please type a number or expression for denominator.");
	    return;
	  }
    }
	digitGroup = document.getElementById('digits').value;
	document.getElementById('help').style.display = "none";
	switch (app)
	{
    case 0:
      res.innerHTML = "Computing sum of squares...";
      break;
    case 1:
      res.innerHTML = "Calculando suma de cuadrados...";
      break;
    case 2:
      res.innerHTML = "Computing sum of cubes...";
      break;
    case 3:
      res.innerHTML = "Calculando suma de cubos...";
      break;
    case 4:
      res.innerHTML = "Computing continued fraction expansion...";
      break;
    default:
      res.innerHTML = "Calculando desarrollo en fracciones continuas...";
      break;
	}
	param = digitGroup + ',' + app + ',' + valueA + String.fromCharCode(0);
	if (app >= 4)
	{
	  param += valueB + String.fromCharCode(0) + valueC + String.fromCharCode(0);
	}
    callWorker(param);
  }
  document.getElementById('helpbtn').onclick = function()
  {
    document.getElementById('help').style.display = "block";
    document.getElementById('result').style.display = "none";
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

