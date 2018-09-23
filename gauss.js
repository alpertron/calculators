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
/** @define {number} */ var lang = 1;   // Use with Closure compiler.
(function(global)
{   // This method separates the name space from the Google Analytics code.
var worker = 0;
var app;
var blob;
var workerParam;
var fileContents = 0;
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
    { // First character of e.data is "1" for intermediate text
      // and it is "2" for end of calculation.
      // It is "9" for saving expression to factor into Web Storage.
      var firstChar = e.data.substring(0, 1);
      if (firstChar == "8")
      {
//        setStorage("ecmFactors", e.data.substring(1));
//        setStorage("ecmCurve", "");
      }
      else if (firstChar == "7")
      {
//        setStorage("ecmCurve", e.data.substring(1));
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
  app = lang + n;
  var res = get("result");
  var valueText = get("value").value;
  var helphelp = get("helphelp");
  get("help").style.display = "none";
  helphelp.style.display = "block";
  helphelp.innerHTML = (lang? "<p>Aprieta el botón <strong>Ayuda</strong> para obtener ayuda para esta aplicación. Apriétalo de nuevo para retornar a la factorización.</p>":
                              "<p>Press the <strong>Help</strong> button to get help about this application. Press it again to return to the factorization.</p>");
  res.style.display = "block";
  if (valueText == "")
  {
    res.innerHTML = (lang? "Por favor ingrese una expresión." :
                           "Please type an expression.");
    return;
  }
  styleButtons("none", "inline");  // Enable "more" and "stop" buttons
  res.innerHTML = (lang? "Factorizando la expresión..." :
                         "Factoring expression...");
  param = digits + "," + app + "," + valueText + String.fromCharCode(0);
  callWorker(param);
}

function endFeedback()
{
  get("main").style.display = "block";
  get("feedback").style.display = "none";
  get("value").focus();   
}

window.onload = function ()
{
  var param;
  get("eval").onclick = function ()
  {
    dowork(0);
  };
  get("factor").onclick = function ()
  {
    dowork(2);
  };
  get("stop").onclick = function ()
  {
    worker.terminate();
    worker = 0;
    styleButtons("inline", "none");  // Enable eval and factor
    get("result").innerHTML =
      (lang? "<p>Cálculo detenido por el usuario.</p>" :
             "<p>Calculation stopped by user</p>");
    get("status").innerHTML = "";
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
      userdata.value = "\n" + get("value").value + "\n" + get("result").innerHTML + "\n" + get("status").innerHTML;
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
  if ('serviceWorker' in navigator)
  { // Attempt to register service worker.
    // There is no need to do anything on registration success or failure in this JavaScript module.
    navigator["serviceWorker"].register('calcSW.js').then(function() {}, function() {});
  }
}
var req = new XMLHttpRequest();
req.open('GET', "gaussianW0000.js", true);
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
