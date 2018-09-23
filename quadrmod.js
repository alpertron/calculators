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
/** @define {number} */ var lang = 0;   // Use with Closure compiler.
(function(global)
{   // This method separates the name space from the Google Analytics code.
  var worker = 0;
  var blob;
  var workerParam;
  var fileContents = 0;
  function get(x)
  {
    return document.getElementById(x);
  }
  function callWorker(param)
  {
    if (!worker)
    {
      if (!blob)
      {
        blob = new Blob([new Uint8Array(fileContents)]);
      }
      worker = new Worker(window.URL.createObjectURL(blob));
      worker.onmessage = function(e)
      { // First character of e.data is '1' for intermediate text
        // and it is '2' for end of calculation.
        get('result').innerHTML = e.data.substring(1);
        if (e.data.substring(0, 1) == '2')
        {   // First character passed from web worker is '2'.
          get('solve').disabled = false;
          get('stop').disabled = true;
        }
      }
    }
    worker.postMessage(param);
  }

  function dowork(n)
  {
    var param;
    var res = get('result');
    var quadrText = get('quad').value.trim();
    var linText = get('lin').value.trim();
    var constText = get('const').value.trim();
    var modText = get('mod').value.trim();
    var digitGroup = get('digits').value.trim();
    get('help').style.display = "none";
    res.style.display = "block";
    var missing = "";
    if (quadrText == "")
    {
      missing = (lang? "coeficiente cuadrático." : "quadratic coefficient.");
    }
    if (linText == "")
    {
      missing = (lang? "coeficiente lineal." : "linear coefficient.");
    }
    if (constText == "")
    {
      missing = (lang? "término independiente." : "constant coefficient.");
    }
    if (modText == "")
    {
      missing = (lang? "módulo." : "modulus.");
    }
    if (missing != "")
    {
      res.innerHTML = (lang? "Por favor ingrese un número o expresión para el "+missing :
                                 "Please type a number or expression for the "+missing);
      return;
    }
    get('solve').disabled = true;
    get('stop').disabled = false;
    res.innerHTML = (lang? "Resolviendo la ecuación cuadrática..." :
                               "Solving the quadratic equation...");
    param = digitGroup + ',' + lang + ',' + quadrText + String.fromCharCode(0) + linText +
      String.fromCharCode(0) + constText +String.fromCharCode(0) + modText + String.fromCharCode(0);
    callWorker(param);
  }

  function endFeedback()
  {
    get("main").style.display = "block";
    get("feedback").style.display = "none";
    get("quad").focus();   
  }

  window.onload = function ()
  {
    var param;
    get('stop').disabled = true;
    get('solve').onclick = function ()
    {
      dowork(0);
    }
    get('stop').onclick = function ()
    {
      worker.terminate();
      worker = 0;
      get('solve').disabled = false;
      get('stop').disabled = true;
      get('result').innerHTML = 
        (lang? "<p>Cálculo detenido por el usuario.</p>" :
                   "<p>Calculation stopped by user</p>");
    }
    get('helpbtn').onclick = function ()
    {
      get('help').style.display = "block";
      get('result').style.display = "none";
    }
    get("formlink").onclick = function ()
    {
      get("main").style.display = "none";
      get("feedback").style.display = "block";
      get("formfeedback").reset();
      get("name").focus();
      return false;   // Do not follow the link.
    }
    get("formcancel").onclick = function ()
    {
      endFeedback();
    }
    get("formsend").onclick = function()
    {
      var userdata = get("userdata");
      if (get("adduserdata").checked)
      {
        userdata.value = "ax^2 + bx + c = 0 (mod n)" + 
                         "\na = " + get("quad").value + "\nb = " + get("lin").value +
                         "\nc = " + get("const").value + "\nn = " + get("mod").value;
      }
      else
      {
        userdata.value = "";      
      }
      var xhr = new XMLHttpRequest();
      xhr.onreadystatechange = function (event)
      {
        if (xhr.readyState == 4) 
        {             // XHR finished.
          if (xhr.status == 200)
          {           // PHP page loaded.
            alert(lang?"Comentarios enviados satisfactoriamente.": "Feedback sent successfully.");
          }
          else
          {           // PHP page not loaded.
            alert(lang?"No se pudieron enviar los comentarios.": "Feedback could not be sent.");
          }
          endFeedback();
        }
      };
      xhr.open("POST", (lang? "/enviomail.php": "/sendmail.php"), true);
      xhr.setRequestHeader('Content-type', 'application/x-www-form-urlencoded');
      var elements = get("formfeedback").elements;
      var contents = "";
      var useAmp = 0;
      for (var i = 0; i < elements.length; i++)
      {
        var element = elements[i];
        if (element.type == "radio" && element.checked == false)
        {
          continue;
        }
        if (element.name)
        {
          if (useAmp)
          {
            contents += '&';
          }
          contents += element.name + "=" + encodeURIComponent(element.value);
          useAmp++;
        }
      }
      xhr.send(contents);
      return false;   // Send form only through JavaScript.
    }
    if ('serviceWorker' in navigator)
    { // Attempt to register service worker.
      // There is no need to do anything on registration success or failure in this JavaScript module.
      navigator["serviceWorker"].register('calcSW.js').then(function() {}, function() {});
    }
  }
  var req = new XMLHttpRequest();
  req.open('GET', "quadmodW0000.js", true);
  req.responseType = "arraybuffer";
  req.onreadystatechange = function (aEvt)
  {
    if (req.readyState == 4 && req.status == 200)
    {
      fileContents = req.response;
      if (workerParam)
      {
        callWorker(workerParam);
      }
    }
  };
  req.send(null);
})(this);

// Register Google Analytics.
addEventListener("load", function ()
{
  (function(i,s,o,g,r,a,m){i["GoogleAnalyticsObject"]=r;i[r]=i[r]||function(){
    (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
    m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
    })(window,document,"script","https://www.google-analytics.com/analytics.js","ga");
  
  ga("create", "UA-4438475-1", "auto");
  ga("send", "pageview");
});
