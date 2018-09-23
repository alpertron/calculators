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
/** @define {number} */ var app = 0;   // Use with Closure compiler.
(function(global)
{   // This method separates the name space from the Google Analytics code.
var wizardStep = 0;
var wizardTextInput;
var worker = 0;
var fileContents = 0;
var hex = 0;
var blob;
var lang = app % 2;
var workerParam;
var asmjs = typeof(WebAssembly) === "undefined";

function get(x)
{
  return document.getElementById(x);
}

function styleButtons(style1, style2)
{
  get("calc").style.display = style1;
  if (app < 4)
  {    // Continued fraction applet does not use wizard.
    get("openwizard").style.display = style1;
  }
  if (get("stop") != null)
  {
    get("stop").style.display = style2;
  }
}

function callWorker(param)
{
  var helphelp = get("helphelp");
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
    { // First character of e.data is:
      // "1" for intermediate output
      // "2" for end calculation
      // "4" for sending data to status line
      // "6" for pausing calculation and showing the Continue button
      var firstChar = e.data.substring(0, 1);
      if (firstChar == "4")
      {
        get("status").innerHTML = e.data.substring(1);
      }
      else
      {
        get("result").innerHTML = e.data.substring(1);
        if (firstChar == "2" || firstChar == "6")
        {   // First character passed from web worker is "2".
          get("status").innerHTML = "";
          styleButtons("inline", "none");  // Enable buttons that must be enabled when applet is not running
          if (firstChar == "6")
          {
            get("cont").style.display = "block";
          }
        }
      }
    };
  }
  helphelp.style.display = "block";
  helphelp.innerHTML = (lang ? '<p class="pad">Aprieta el botón <strong>Ayuda</strong> para obtener ayuda para esta aplicación. Apriétalo de nuevo para retornar a esta pantalla. Los usuarios con teclado pueden presionar CTRL+ENTER para comenzar el cálculo. Esta es la versión '+(asmjs? "asm.js": "WebAssembly")+".</p>":
                               '<p class="pad">Press the <strong>Help</strong> button to get help about this application. Press it again to return to this screen. Keyboard users can press CTRL+ENTER to start calculation. This is the '+(asmjs? "asm.js": "WebAssembly")+" version.</p>");
  if (asmjs)
  {      // Asm.js
    worker.postMessage(param);
  }
  else
  {      // WebAssembly.
    worker.postMessage([param, fileContents]);
  }
}

function performCalc()
{
  var res, valueA, valueB, valueC, digitGroup;
  res = get('result');
  res.style.display = "block";
  valueA = get('num').value;
  if (valueA == "")
  {
    if (app >= 4)
    {
      res.innerHTML = (lang ? "Por favor ingrese un número o expresión para el numerador." :
                              "Please type a number or expression for numerator.");
    }
    else
    {
      res.innerHTML = (lang ? "Por favor ingrese un número o expresión." :
                              "Please type a number or expression.");
    }
    return;
  }
  if (app >= 4)
  {
    valueB = get('delta').value;
    if (valueB == "")
    {
      res.innerHTML = (lang ? "Por favor ingrese un número o expresión para el argumento de la raíz cuadrada." :
                              "Please type a number or expression for square root argument.");
      return;
    }
    valueC = get('den').value;
    if (valueC == "")
    {
      res.innerHTML = (lang ? "Por favor ingrese un número o expresión para el denominador." :
                              "Please type a number or expression for denominator.");
      return;
    }
  }
  digitGroup = get('digits').value;
  get('help').style.display = "none";
  if (app == 0)    // Closure compiler cannot optimize switch, so a series of "if" instructions is used.
  {
    res.innerHTML = "Computing sum of squares...";
  }
  else if (app == 1)
  {
    res.innerHTML = "Calculando suma de cuadrados...";
  }
  else if (app == 2)
  {
    res.innerHTML = "Computing sum of cubes...";
  }
  else if (app == 3)
  {
    res.innerHTML = "Calculando suma de cubos...";
  }
  else if (app == 4)
  {
    res.innerHTML = "Computing continued fraction expansion...";
  }
  else
  {
    res.innerHTML = "Calculando desarrollo en fracciones continuas...";
  }
  if (app >= 4)
  {
    hex = (get("converg").checked? 1: 0);
  }
  param = digitGroup + ',' + (app+hex*64) + ',' + valueA + String.fromCharCode(0);
  if (app >= 4)
  {
    param += valueB + String.fromCharCode(0) + valueC + String.fromCharCode(0);
  }
  else
  {
    styleButtons("none", "inline");  // Enable "stop" button
  }
  callWorker(param);
}

function oneexpr()
{
  get("next").value = (lang? "Hecho": "Done");
  get("wzddesc").innerHTML = (lang? "Paso 1 de 1: Expresión a factorizar": "Step 1 of 1: Expression to factor");
  get("wzdexam").innerHTML = "&nbsp;";
  wizardTextInput = "";
  wizardStep = 9;
}

function selectLoop()
{   
  get("next").value = (lang ? "Siguiente": "Next");
  get("wzddesc").innerHTML = (lang ? "Paso 1 de 5: Valor inicial de x": "Step 1 of 5: Initial value of x");
  get("wzdexam").innerHTML = (lang? "No usar variables <var>x</var> o <var>c</var>. Ejemplo para números de Smith menores que 10000: <code>1</code>": 
                                       "Do not use variables <var>x</var> or <var>c</var>. Example for Smith numbers less than 10000: <code>1</code>");
  wizardStep = 1;
}
  
function wizardNext()
{
  var nextBtn = get("next");
  var wzdDescText = get("wzddesc");
  var wzdExamText = get("wzdexam");
  var wzdInput = get("wzdinput");
  var valueInput = get("num");
  var textExample = (lang? "Variables <var>x</var> y/o <var>c</var> requeridas. Ejemplo para números de Smith menores que 10000: <code>":
                           "Variables <var>x</var> and/or <var>c</var> required. Example for Smith numbers less than 10000: <code>");
  nextBtn.disabled = true;
  switch (++wizardStep)
  {
    case 2:
      wizardTextInput += "x="+wzdInput.value;
      get("mode").style.display = "none";
      wzdDescText.innerHTML = (lang? "Paso 2 de 5: Valor de x para la nueva iteración": "Step 2 of 5: Value of x for new iteration");
      wzdExamText.innerHTML = textExample + "x+1</code>";
      break;
    case 3:
      wizardTextInput += ";x="+wzdInput.value;
      wzdDescText.innerHTML = (lang? "Paso 3 de 5: Condición para finalizar el ciclo": "Step 3 of 5: End loop condition");
      wzdExamText.innerHTML = textExample + "x&lt;10000</code>";
      break;
    case 4:
      wizardTextInput += ";"+wzdInput.value;
      wzdDescText.innerHTML = (lang? "Paso 4 de 5: Expresión a factorizar": "Step 4 of 5: Expression to factor");
      wzdExamText.innerHTML = textExample + "x</code>";
      break;
    case 5:
      wizardTextInput += ";"+wzdInput.value;
      nextBtn.value = (lang? "Hecho": "Done");
      nextBtn.disabled = false;
      wzdDescText.innerHTML = (lang? "Paso 5 de 5: Condición para procesar la expresión": "Step 5 of 5: Process expression condition");
      wzdExamText.innerHTML = textExample + "sumdigits(x,10) == sumdigits(concatfact(2,x),10) and not isprime(x)</code>";
      break;
    case 6:
      if (wzdInput.value != "")
      {
        wizardTextInput += ";"+wzdInput.value;
      }
      valueInput.value = wizardTextInput;
      wizardStep = 0;
      hex = (get("hexW").checked? 1: 0);
      get("main").style.display = "block";
      get("wizard").style.display = "none";
      valueInput.focus();
      break;
    default:
      wizardStep = 0;
      valueInput.value = wzdInput.value;
      hex = (get("hexW").checked? 1: 0);
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

function endFeedback()
{
  get("main").style.display = "block";
  get("feedback").style.display = "none";
  get(app<4?"value": "num").focus();
}

function startUp()
{
  var param;
  if (app<4)
  {
    get('num').onkeypress = function(e)
    {
      var output, res;
      var digitGroup = get('digits').value;
      res = get('result');
      res.style.display = "block";
      var input = get('num').value;
      if (!e) e = window.event;
      var keyCode = e.keyCode || e.which;
      if (keyCode == 13)
      {  // Used pressed Enter key
        output = get('result')
        if (input == "")
        {
          res.innerHTML = (lang ? "Por favor ingrese un número o expresión." : "Please type a number or expression.");
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
        styleButtons("none", "inline");  // Enable "stop" button
        callWorker(param);
      }
    }
  }
  get('calc').onclick = function()
  {
    performCalc();
  };
  if (app < 4)
  {    // Continued fraction applet does not use wizard.
    get("openwizard").onclick = function ()
    {
      get("main").style.display = "none";
      get("wizard").style.display = "block";
      get("mode").style.display = "block";
      get("oneexpr").checked = true;
      get("next").disabled = true;
      get("hexW").checked = (hex? true: false);
      get("decW").checked = (hex? false: true);
      get("wzdinput").value = "";
      get("wzdinput").focus();
      oneexpr();
    }
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
    }
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
  }
  if (get("stop") != null)
  {
    get("stop").onclick = function ()
    {
      worker.terminate();
      worker = 0;
      styleButtons("inline", "none");  // Enable buttons that have to be enabled when applet is not running.
      get("result").innerHTML =
        (lang ? "<p>Cálculo detenido por el usuario.</p>" :
                "<p>Calculation stopped by user</p>");
      get("status").innerHTML = "";
    }
  };
  if (get("continue") != null)
  {
    get("continue").onclick = function ()
    {
      get("cont").style.display = "none";
      callWorker("C");  // Indicate worker that user pressed Continue button.
    }
  }
  get("num").onkeydown = function (event)
  {
    if ((event.keyCode == 10 || event.keyCode == 13) && event.ctrlKey)
    {
      event.preventDefault();          // Do not propagate Enter key.
      performCalc();                   // Perform calculation.
    }
    return true;
  }
  get('helpbtn').onclick = function()
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
      if(app < 4)
      {   
        userdata.value = "\n" + get("num").value + "\n" + get("result").innerHTML + "\n" + get("status").innerHTML;
      }
      else
      {
        userdata.value = "\nnum = " + get("num").value + "\ndelta = " + get("delta").value + "\nden = " + get("den").value;         
      }
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
req.open('GET', (asmjs? "fsquaresW0000.js": "fsquares0000.wasm"), true);
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
addEventListener("load", startUp);
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
