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
var digits;
var config;

function get(x)
{
  return document.getElementById(x);
}

function setStorage(name, data)
{
  localStorage.setItem(name, data);
}

function getStorage(name)
{
  return localStorage.getItem(name);
}

function styleButtons(style1, style2)
{
  get("eval").style.display = style1;
  get("factor").style.display = style1;
  get("config").style.display = style1;
  get("stop").style.display = style2;
  get("more").style.display = style2;
}

function restartFactorization(type)
{
  get("modal-more").style.display = "none";
  worker.terminate();
  worker = 0;
  dowork(type);
}

function callWorker(param)
{
  if (!worker)
  {
    worker = new Worker("ecmW0023.js");
    worker.onmessage = function(e)
    { // First character of e.data is "1" for intermediate text
      // and it is "2" for end of calculation.
      // It is "9" for saving expression to factor into Web Storage.
      var firstChar = e.data.substring(0, 1);
      if (firstChar == "8")
      {
        setStorage("ecmFactors", e.data.substring(1));
        setStorage("ecmCurve", "");
      }
      else if (firstChar == "7")
      {
        setStorage("ecmCurve", e.data.substring(1));
      }
      else if (firstChar == "4")
      {
        get("status").innerHTML = e.data.substring(1);
      }
      else
      {
        get("result").innerHTML = e.data.substring(1);
        if (e.data.substring(0, 1) == "2")
        {   // First character passed from web worker is "2".
          get("status").innerHTML = "";
          styleButtons("inline", "none");  // Enable eval and factor
          get("modal-more").style.display = "none";
        }
      }
    };
  }
  worker.postMessage(param);
}

function dowork(n)
{
  app = parseInt(get("app").value) + n;
  var res = get("result");
  var valueText = get(config.substr(0,1)=="1"?"textarea":"value").value;
  var charNull = String.fromCharCode(0);
  var helphelp = get("helphelp");
  get("help").style.display = "none";
  helphelp.style.display = "block";
  helphelp.innerHTML = (app & 1 ? "<p>Aprieta el botón <strong>Ayuda</strong> para obtener ayuda para esta aplicación. Apriétalo de nuevo para retornar a la factorización.</p>":
                                  "<p>Press the <strong>Help</strong> button to get help about this application. Press it again to return to the factorization.</p>");
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
  param = digits + "," + app + "," + config + valueText + charNull +
          getStorage("ecmFactors");
  if (n == -1 || n == -2)
  {
    param += charNull + get("curve").value;   // Append new factor or curve number.
  }
  if (n == -3 || n == -4)
  {
    param += "*" + get("curve").value + "^1(2)";   // Append new factor or curve number.
  }
  callWorker(param + charNull);
}

function isBatch()
{   
  if (config.substr(0,1)=="1")
  {
    get("value").style.display = "none";
    get("lvalue").style.display = "none";
    get("textarea").style.display = "inline";
    get("ltextarea").style.display = "inline";
  }
  else
  {
    get("value").style.display = "inline";
    get("lvalue").style.display = "inline";
    get("textarea").style.display = "none";
    get("ltextarea").style.display = "none";
  }
}
window.onload = function ()
{
  var param, index;
  get("eval").onclick = function ()
  {
    setStorage("ecmFactors","");
    dowork(0);
  };
  get("factor").onclick = function ()
  {
    setStorage("ecmFactors","");
    dowork(2);
  };
  get("more").onclick = function ()
  {
    get("modal-more").style.display = "block";
  };
  get("config").onclick = function ()
  {
	get("digits").value = digits;
    get("batch").checked = (config.substr(0,1)=="1");
    get("verbose").checked = (config.substr(1,1)=="1");
    get("pretty").checked = (config.substr(2,1)=="1");
    get("cunnin").checked = (config.substr(3,1)=="1");  
    get("modal-config").style.display = "block";
  };
  get("close-config").onclick = function ()
  {
    get("modal-config").style.display = "none";
  };
  get("cancel-config").onclick = function ()
  {
    get("modal-config").style.display = "none";
  };
  get("save-config").onclick = function ()
  {
    oldconfig = config;
    config = (get("batch").checked? "1" :"0") +
             (get("verbose").checked? "1" :"0") +
             (get("pretty").checked? "1" :"0") +
             (get("cunnin").checked? "1" :"0");
    digits = get("digits").value;
    setStorage("ecmConfig", digits+","+config);
    if (config.substr(0,1) != oldconfig.substr(0,1))
    {
      isBatch();
    }
    get("modal-config").style.display = "none";
  };
  get("close-more").onclick = function ()
  {
    get("modal-more").style.display = "none";
  };
  get("ncurve").onclick = function ()
  {
    restartFactorization(-2);
  };
  get("nfactor").onclick = function ()
  {
    restartFactorization(-4);
  };
  get("curve").onkeypress = function(event)
  {
    return (event.charCode == 8 || event.charCode == 0) ? null : event.charCode >= 48 && event.charCode <= 57;
  };
  get("stop").onclick = function ()
  {
    worker.terminate();
    worker = 0;
    styleButtons("inline", "none");  // Enable eval and factor
    get("result").innerHTML =
      (app & 1 ? "<p>Cálculo detenido por el usuario.</p>" :
                 "<p>Calculation stopped by user</p>");
    get("status").innerHTML = "";
  };
  get("value").onkeydown = function (event)
  {
	if (event.keyCode == 13)
	{
	  event.preventDefault();        // Do not propagate Enter key.
      setStorage("ecmFactors","");   // Perform factorization.
      dowork(2);
	}
	return true;
  }
  get("helpbtn").onclick = function ()
  {
    var help = get("help");
	var helpStyle = help.style;
	var helphelpStyle = get("helphelp").style;
	var result = get("result");
	var resultStyle = result.style;
	if (helpStyle.display == "block" && result.innerHTML != "")		
	{
	  helpStyle.display = "none";
	  helphelpStyle.display = resultStyle.display = "block";
	}
	else
	{
      helpStyle.display = "block";
	  helphelpStyle.display = resultStyle.display = "none";
	}
  };
  window.onclick = function(event)
  {
    var modal = get("modal");
    if (event.target == modal)
    {
      modal.style.display = "none";
    }
  };
  digits = getStorage("ecmConfig");
  if (digits == null || digits == "")
  {
    digits = 6;
    config = "0010";
    setStorage("ecmConfig", digits+","+config);
  }
  else
  {
    index = digits.indexOf(",");
    if (index<0)
    {
      digits = 6;
      config = "0010";
      setStorage("ecmConfig", digits+","+config);
    }
    else
    {
      config = digits.substr(index+1);
      digits = digits.substr(0,index);
    }
  }
  isBatch();
  ecmFactor = getStorage("ecmFactors");
  if (ecmFactor)
  {          // Continue factoring.
    get("value").value = ecmFactor.slice(0,ecmFactor.indexOf("="));
    get("curve").value = getStorage("ecmCurve");
    dowork(-2);
    get("curve").value = "";
  }
  (function(i,s,o,g,r,a,m){i["GoogleAnalyticsObject"]=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,"script","https://www.google-analytics.com/analytics.js","ga");

  ga("create", "UA-4438475-1", "auto");
  ga("send", "pageview");
}
