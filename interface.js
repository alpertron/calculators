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
function callWorker(param)
{
  if (!worker)
  {
    worker = new Worker('fsquaresW.js?0108');
    worker.onmessage = function(e)
    {
      if (app < 4)
      {
        document.getElementById('output').value = e.data;
      }
      else
      {
	      document.getElementById('result').innerHTML = e.data;
      }
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
      var output;
	    var digitGroup = document.getElementById('digits').value;
      app = document.getElementById('app').value;
      var input = document.getElementById('input').value;
      if (!e) e = window.event;
      var keyCode = e.keyCode || e.which;
      if (keyCode == 13)
	    {  // Used pressed Enter key
	      output = document.getElementById('output')
	      if (input == "")
	      {
          output.value = (app&1?"Por favor ingrese un número o expresión.": "Please type a number or expression.");
          return;
	      }
	      if (app==0)
		    {
		      output.value = "Computing sum of squares...";
		    }
		    else if (app==1)
		    {
          output.value = "Calculando suma de cuadrados...";
		    }
		    else if (app==2)
		    {
		      output.value = "Computing sum of cubes...";
		    }
		    else
		    {
		      output.value = "Calculando suma de cubos...";
		    }
		    param = digitGroup + ',' + app + ',' + input;
		    callWorker(param);
	    }
    }
  }
  if (document.getElementById('calc'))
  {
    document.getElementById('calc').onclick = function()
    {
	    var app = document.getElementById('app').value;
  	  var res = document.getElementById('result');
	    var valueA = document.getElementById('num').value;
      var valueB = document.getElementById('delta').value;
      var valueC = document.getElementById('den').value;
	    var digitGroup = document.getElementById('digits').value;
	    document.getElementById('help').style.display = "none";
	    res.style.display = "block";
	    if (valueA == "")
      {
	      res.innerHTML = (app&1? "Por favor ingrese un número o expresión para el numerador." :
	                              "Please type a number or expression for numerator.");
        return;
	    }
	    if (valueB == "")
	    {
  	    res.innerHTML = (app&1? "Por favor ingrese un número o expresión para el argumento de la raíz cuadrada." :
	                              "Please type a number or expression for square root argument.");
	      return;
	    }
	    if (valueC == "")
	    {
  	    res.innerHTML = (app&1? "Por favor ingrese un número o expresión para el denominador." :
	                              "Please type a number or expression for denominator.");
	      return;
	    }
      res.innerHTML = (app&1? "Calculando desarrollo en fracciones continuas..." :
	                          "Computing continued fraction expansion...");
      param = digitGroup + ',' + app + ',' + valueA + String.fromCharCode(0) + valueB +
			String.fromCharCode(0) + valueC + String.fromCharCode(0);
      callWorker(param);
	  }
    document.getElementById('helpbtn').onclick = function()
	  {
      document.getElementById('help').style.display = "block";
      document.getElementById('result').style.display = "none";
	  }
  }
}

