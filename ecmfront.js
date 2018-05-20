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

// In order to reduce the number of files to read from Web server, this 
// Javascript file includes both the Javascript in the main thread and the 
// Javascript that drives WebAssembly on its own Web Worker.
(function(global)
{   // This method separates the name space from the Google Analytics code.
var worker = 0;
var app;
var digits;
var config;
var asmjsFileName = "ecm0050.js";
var wasmFileName = "ecm0050.wasm";
var asmjs = typeof(WebAssembly) === "undefined";

function msgRecvByWorker(e)
{
  var request;
  if (wasmLoaded)
  {
    ConvertToString(exports["getInputStringPtr"](), e.data);
    exports["doWork"]();
    return;  
  }
  request = new XMLHttpRequest();
  request.open('GET', wasmFileName);
  request.responseType = 'arraybuffer';
  request.send();

  request.onload = function()
  {
    if (request.status != 200)
    {
      return;
    }
    var bytes = request.response;
    WebAssembly["instantiate"](bytes, info).then(function(results)
    {
      wasmLoaded = 1;
      exports = results["instance"]["exports"];
      HEAPU8 = new Uint8Array(exports["memory"]["buffer"]);
      ConvertToString(exports["getInputStringPtr"](), e.data);
      exports["doWork"]();
      return;
    });
  };
}

function oneexpr()
{
  get("next").value = (app & 1? "Hecho": "Done");
  get("wzddesc").innerHTML = (app & 1? "Paso 1 de 1: Expresión a factorizar": "Step 1 of 1: Expression to factor");
  get("wzdexam").innerHTML = "&nbsp;";
  wizardTextInput = "";
  wizardStep = 9;
}

function PtrToString(ptr)
{
  var t=-1;
  var i = 0;
  var str="", outString="";
  do
  {
    for (i=0; i<1024; i++)
    {
      t = HEAPU8[((ptr++)>>0)];
      if (t==0)
      {
        break;
      }
      if (t>=128)
      {
        t = ((t-192)<<6) + HEAPU8[((ptr++)>>0)] - 128;
      }
      str += String.fromCharCode(t);
    }
    outString += str;
    str = "";
  } while (t!=0);
  outString += str;
  return outString;
}

function ConvertToString(ptr, str)
{
  var dest = ptr;
  var length = str.length;
  var i, t;
  for (i=0; i<length; i++)
  {
    t = str.charCodeAt(i);
    if (t<128)
    {
      HEAPU8[dest++] = t;
    }
    else
    {
      HEAPU8[dest++] = (t >> 6) + 192;
      HEAPU8[dest++] = (t & 63) + 128;
    }
  }
  HEAPU8[dest] = 0;
}

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
  get("openwizard").style.display = style1;
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
    worker = new Worker(asmjs? "ecmW0050.js": asmjsFileName);
    worker.onmessage = function(e)
    { // First character of e.data is:
      // "1" for intermediate output
      // "2" for end calculation
      // "4" for sending intermediate data
      // "6" for pausing calculation and showing the Continue button
      // "7" for saving curve number into local storage
      // "8" for saving input expression into local storage
      // "9" for sending data to console.
      var firstChar = e.data.substring(0, 1);
      if (firstChar == "9")
      {
        console.log(e.data.substring(1));
      }
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
        if (firstChar == "2" || firstChar == "6")
        {   // First character passed from web worker is "2".
          get("status").innerHTML = "";
          styleButtons("inline", "none");  // Enable eval and factor
          get("modal-more").style.display = "none";
          if (firstChar == "6")
          {
            get("cont").style.display = "block";
          }
        }
      }
    };
  }
  worker.postMessage(param);
}

function dowork(n)
{
  var param;
  app = parseInt(get("app").value) + n;
  var res = get("result");
  var valueText = get("value").value.replace(/\u2011/g, "-");
  var charNull = String.fromCharCode(0);
  var helphelp = get("helphelp");
  get("cont").style.display = "none";
  get("help").style.display = "none";
  helphelp.style.display = "block";
  helphelp.innerHTML = (app & 1 ? '<p class="pad">Aprieta el botón <strong>Ayuda</strong> para obtener ayuda para esta aplicación. Apriétalo de nuevo para retornar a la factorización. Los usuarios con teclado pueden presionar CTRL+ENTER para comenzar la factorización. Esta es la versión '+(asmjs? "asm.js": "WebAssembly")+".</p>":
                                  '<p class="pad">Press the <strong>Help</strong> button to get help about this application. Press it again to return to the factorization. Keyboard users can press CTRL+ENTER to start factorization. This is the '+(asmjs? "asm.js": "WebAssembly")+" version.</p>");
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

function selectLoop()
{   
  get("next").value = (app & 1 ? "Siguiente": "Next");
  get("wzddesc").innerHTML = (app & 1 ? "Paso 1 de 5: Valor inicial de x": "Step 1 of 5: Initial value of x");
  get("wzdexam").innerHTML = (app & 1? "No usar variables <var>x</var> o <var>c</var>. Ejemplo para números de Smith menores que 10000: <code>1</code>": 
                                       "Do not use variables <var>x</var> or <var>c</var>. Example for Smith numbers less than 10000: <code>1</code>");
  wizardStep = 1;
}
  
function wizardNext()
{
  var nextBtn = get("next");
  var wzdDescText = get("wzddesc");
  var wzdExamText = get("wzdexam");
  var wzdInput = get("wzdinput");
  var valueInput = get("value");
  nextBtn.disabled = true;
  switch (++wizardStep)
  {
    case 2:
      wizardTextInput += "x="+wzdInput.value;
      get("mode").style.display = "none";
      wzdDescText.innerHTML = (app & 1? "Paso 2 de 5: Valor de x para la nueva iteración": "Step 2 of 5: Value of x for new iteration");
      wzdExamText.innerHTML = (app & 1? "Variables <var>x</var> y/o <var>c</var> requeridas. Ejemplo para números de Smith menores que 10000: <code>x+1</code>":
                                        "Variables <var>x</var> and/or <var>c</var> required. Example for Smith numbers less than 10000: <code>x+1</code>");
      break;
    case 3:
      wizardTextInput += ";x="+wzdInput.value;
      wzdDescText.innerHTML = (app & 1? "Paso 3 de 5: Condición para finalizar el ciclo": "Step 3 of 5: End loop condition");
      wzdExamText.innerHTML = (app & 1? "Variables <var>x</var> y/o <var>c</var> requeridas. Ejemplo para números de Smith menores que 10000: <code>x&lt;10000</code>":
                                           "Variables <var>x</var> and/or <var>c</var> required. Example for Smith numbers less than 10000: <code>x&lt;10000</code>");
      break;
    case 4:
      wizardTextInput += ";"+wzdInput.value;
      wzdDescText.innerHTML = (app & 1? "Paso 4 de 5: Expresión a factorizar": "Step 4 of 5: Expression to factor");
      wzdExamText.innerHTML = (app & 1? "Variables <var>x</var> y/o <var>c</var> requeridas. Ejemplo para números de Smith menores que 10000: <code>x</code>":
                                           "Variables <var>x</var> and/or <var>c</var> required. Example for Smith numbers less than 10000: <code>x</code>");
      break;
    case 5:
      wizardTextInput += ";"+wzdInput.value;
      nextBtn.value = (app & 1? "Hecho": "Done");
      nextBtn.disabled = false;
      wzdDescText.innerHTML = (app & 1? "Paso 5 de 5: Condición para procesar la expresión": "Step 5 of 5: Process expression condition");
      wzdExamText.innerHTML = (app & 1? "Variables <var>x</var> y/o <var>c</var> requeridas. Ejemplo para números de Smith menores que 10000: <code>sumdigits(x,10) == sumdigits(concatfact(2,x),10) and not isprime(x)</code>":
                                        "Variables <var>x</var> and/or <var>c</var> required. Example for Smith numbers less than 10000: <code>sumdigits(x,10) == sumdigits(concatfact(2,x),10) and not isprime(x)</code>");
      break;
    case 6:
      if (wzdInput.value != "")
      {
        wizardTextInput += ";"+wzdInput.value;
      }
      valueInput.value = wizardTextInput;
      wizardStep = 0;
      get("hex").checked = get("hexW").checked;
      saveConfig();
      get("main").style.display = "block";
      get("wizard").style.display = "none";
      valueInput.focus();
      break;
    default:
      wizardStep = 0;
      valueInput.value = wzdInput.value;
      get("hex").checked = get("hexW").checked;
      saveConfig();
      get("main").style.display = "block";
      get("wizard").style.display = "none";
      valueInput.focus();
      break;
  } 
  if (wizardStep)
  {
    wzdInput.value = "";
    wzdInput.focus();
  }
}

function saveConfig()
{    
  config = "1" +   // Batch mode
           (get("verbose").checked? "1" :"0") +
           (get("pretty").checked? "1" :"0") +
           (get("cunnin").checked? "1" :"0") +
           (get("hex").checked? "1" :"0");
  digits = get("digits").value;
  setStorage("ecmConfig", digits+","+config);
}
    
function startUp()
{
  var param, index, ecmFactor;
  app = parseInt(get("app").value);
  get("value").wrap="off";
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
    get("verbose").checked = (config.substr(1,1)=="1");
    get("pretty").checked = (config.substr(2,1)=="1");
    get("cunnin").checked = (config.substr(3,1)=="1");  
    get("hex").checked = (config.substr(4,1)=="1");
    get("modal-config").style.display = "block";
  };
  get("openwizard").onclick = function ()
  {
    get("main").style.display = "none";
    get("wizard").style.display = "block";
    get("mode").style.display = "block";
    get("oneexpr").checked = true;
    get("next").disabled = true;
    get("wzdinput").value = "";
    get("wzdinput").focus();
    get("hexW").checked = (config.substr(4,1)=="1");
    get("decW").checked = (config.substr(4,1)!="1");
    oneexpr();
  };
  get("wzdinput").onkeydown = function (event)
  {
    if (event.keyCode == 10 || event.keyCode == 13)
    {
      event.preventDefault();          // Do not propagate Enter key.
      if (get("next").disabled == false)
      {                                // Next button is not disabled.
        wizardNext();                  // Perform same operation as if the user had pressed Next button.
      }
    }
    if (event.altKey)
    {                                  // User pressed ALT key.
      if (event.keyCode == 80)
      {                                // User pressed ALT-P.
        event.preventDefault();        // Do not propagate key.
        if (get("oneexpr").checked)
        {
          get("oneexpr").checked = false;
          get("loop").checked = true;
          selectLoop();
        }
        else
        {
          get("oneexpr").checked = true;
          get("loop").checked = false;
          oneexpr();
        }
      }
      else if (event.keyCode == 68)
      {                                // User pressed ALT-D.
        event.preventDefault();        // Do not propagate key.
        get("decW").checked = true;
        get("hexW").checked = false;
      }
      else if (event.keyCode == 72)
      {                                // User pressed ALT-H.
        event.preventDefault();        // Do not propagate key.
        get("decW").checked = false;
        get("hexW").checked = true;
      }
    }
    return true;
  };
  get("oneexpr").onclick = function ()
  {
    oneexpr();
  };
  get("loop").onclick = function ()
  {
  selectLoop();
  };
  get("next").onclick = function ()
  {
    wizardNext();
  };
  get("wzdinput").oninput = function ()
  {
    var inputValue = get("wzdinput").value;
    var nextBtn = get("next");
    if (inputValue != "")
    {         // User typed something on input box.
      if (wizardStep == 1 || wizardStep == 9 || (inputValue.lastIndexOf("x") >= 0 || inputValue.lastIndexOf("c") >= 0 ||
          inputValue.lastIndexOf("X") >= 0 || inputValue.lastIndexOf("C") >= 0))
      {       // At least one x or c. Indicate valid.
        nextBtn.disabled = false;
      }
      else
      {
        nextBtn.disabled = true;
      }
    }
    else if (wizardStep == 5)
    {         // Last step is optional, so empty input is valid.
      nextBtn.disabled = false;
    }
    else
    {         // For required input, empty input is invalid.
      nextBtn.disabled = true;
    }
  };
  get("cancel").onclick = function ()
  {
    get("main").style.display = "block";
    get("wizard").style.display = "none";
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
    saveConfig();
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
    if ((event.keyCode == 10 || event.keyCode == 13) && event.ctrlKey)
    {
      event.preventDefault();          // Do not propagate Enter key.
      setStorage("ecmFactors","");     // Perform factorization.
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
  get("continue").onclick = function ()
  {
    get("cont").style.display = "none";
    callWorker("C");  // Indicate worker that user pressed Continue button.
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
    config = "00100";
    setStorage("ecmConfig", digits+","+config);
  }
  else
  {
    index = digits.indexOf(",");
    if (index<0)
    {
      digits = 6;
      config = "00100";
      setStorage("ecmConfig", digits+","+config);
    }
    else
    {
      config = digits.substr(index+1);
      while (config.length < 5)
      {  // Convert legacy configuration.
        config += "0";
      }
      digits = digits.substr(0,index);
    }
  }
  ecmFactor = getStorage("ecmFactors");
  if (ecmFactor)
  {          // Continue factoring.
    get("value").value = ecmFactor.slice(0,ecmFactor.indexOf("="));
    get("curve").value = getStorage("ecmCurve");
    dowork(-2);
    get("curve").value = "";
  }
  
  if ('serviceWorker' in navigator)
  { // Attempt to register service worker.
    // There is no need to do anything on registration success or failure in this JavaScript module.
    navigator.serviceWorker.register('calcSW.js').then(function() {}, function() {});
  }
}
if (typeof(window) === "undefined")
{    // Inside Web Worker
  var wizardStep = 0;
  var wizardTextInput;
  var exports, HEAPU8, wasmLoaded;
  var env =
  {
    "databack": function(data)
    {
      self.postMessage(PtrToString(data));
    },
    "tenths": function()
    {
      return Math.floor(new Date().getTime() / 100);
    },
    "getCunn": function(data)
    {
      var req = new XMLHttpRequest();
      req.open('GET', PtrToString(data), false);
      req.send(null);
      if (req.status == 200)
      {
        ConvertToString(exports["getFactorsAsciiPtr"](), req.responseText);
      }
    }
  };

  var info =
  {
    "env": env
  };  

  global.addEventListener('message', msgRecvByWorker);
}
else
{    // Outside Web Worker
  if (!getStorage("ecmFactors"))
  {          // No factorization. Read factorization asm.js or wasm file in idle time.
    fetch(asmjs? asmFileName: wasmFileName).then(function(response) {return;}).catch(function(err) {});
  }
  addEventListener("load", startUp);
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
