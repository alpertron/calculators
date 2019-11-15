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
var asmjs = typeof(WebAssembly) === "undefined";
function get(x)
{
  return document.getElementById(x);
}
function callWorker(param)
{
  if (!worker)
  {
    if (asmjs)
    {    // Asm.js
      if (!blob)
      {
        blob = new Blob([new Uint8Array(fileContents)]);
      }
    }
    else
    {    // WebAssembly
      if (!blob)
      {
        blob = new Blob(Array.prototype.map.call(document.querySelectorAll('script[type=\'text\/js-worker\']'), function (oScript) { return oScript.textContent; }),{type: 'text/javascript'});
      }
    }    
    worker = new Worker(window.URL.createObjectURL(blob));
    worker.onmessage = function(e)
    { // First character of e.data is '1' for intermediate text
      // and it is '2' for end of calculation.
      get('result').innerHTML = e.data.substring(1);
      if (e.data.substring(0, 1) == '2')
      {   // First character passed from web worker is '2'.
        get('eval').disabled = false;
        get('factor').disabled = false;
        get('stop').disabled = true;
      }
    }
  }
  if (asmjs)
  {      // Asm.js
    worker.postMessage(param);
  }
  else
  {      // WebAssembly.
    worker.postMessage([param, fileContents]);
  }
}

function dowork(n)
{
  var app = lang + n + (get('sup').checked ? 4 : 0);
  var res = get('result');
  var polyText = get('poly').value;
  var modText = get('mod').value;
  var digitGroup = get('digits').value;
  get('help').style.display = "none";
  res.style.display = "block";
  if (polyText == "")
  {
    res.innerHTML = (lang? "Por favor ingrese una expresión para el polinomio a evaluar." :
                           "Please type an expression for the polynomial to evaluate.");
    return;
  }
  if (modText == "")
  {
    res.innerHTML = (lang? "Por favor ingrese un número o expresión para el módulo." :
                           "Please type a number or expression for the modulus.");
    return;
  }
  get('eval').disabled = true;
  get('factor').disabled = true;
  get('stop').disabled = false;
  res.innerHTML = (lang? "Factorizando el polinomio..." :
                         "Factoring polynomial...");
  param = digitGroup + ',' + app + ',' + modText + String.fromCharCode(0) + polyText +
  String.fromCharCode(0);
  if (!fileContents)
  {
    workerParam = param;
  }
  else
  {
    callWorker(param);
  }
}

function endFeedback()
{
  get("main").style.display = "block";
  get("feedback").style.display = "none";
  get("poly").focus();   
}

window.onload = function ()
{
  var param;
  get('stop').disabled = true;
  get('eval').onclick = function ()
  {
    dowork(2);
  }
  get('factor').onclick = function ()
  {
    dowork(0);
  }
  get('stop').onclick = function ()
  {
    worker.terminate();
    worker = 0;
    get('eval').disabled = false;
    get('factor').disabled = false;
    get('stop').disabled = true;
    get('result').innerHTML = 
      (lang? "<p>Factorización detenida por el usuario.</p>" :
             "<p>Factorization stopped by user</p>");
  }
  get('helpbtn').onclick = function ()
  {
    get('help').style.display = "block";
    get('result').style.display = "none";
  }
  get('poly').oninput = function ()
  {
    var input = get('poly');
    var loc = input.value.length - input.selectionStart;
    input.value = input.value.replace(".", "x^");
    setTimeout(function()
    {
      loc = input.value.length - loc;
      input.selectionStart = loc;
      input.selectionEnd = loc;
    }, 30);   
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
      userdata.value = get("poly").value + " (mod " + get("mod").value + ")";
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

  // Generate accordion.
  var acc = document.querySelectorAll("h2");
  var idx, x, y;

  for (idx = 0; idx < acc.length; idx++)
  {
    acc[idx].addEventListener("click", function()
    {
    // "active" means that panel is being displayed.
      this.children[0].classList.toggle("active");
      var panel = this.nextElementSibling;
      if (panel.style.display === "block")
      {
        panel.style.display = "none";
      }
      else
      {
        panel.style.display = "block";
      }
    });
  }
  if ('serviceWorker' in navigator)
  { // Attempt to register service worker.
    // There is no need to do anything on registration success or failure in this JavaScript module.
    navigator["serviceWorker"].register('calcSW.js').then(function() {}, function() {});
  }
}
var req = new XMLHttpRequest();
req.open('GET', (asmjs? "polfactW0000.js": "polfact0000.wasm"), true);
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
