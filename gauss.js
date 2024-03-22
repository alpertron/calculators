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
/* global callWorker */
/* global changeInputmode */
/* global clickFormLink */
/* global endCalculation */
/* global endWorker */
/* global formSend */
/* global generateFuncButtons */
/* global get */
/* global getConfig */
/* global getStorage */
/* global getCalculatorCode */
/* global hide */
/* global registerServiceWorker */
/* global setStorage */
/* global show */
/** @define {number} */ const lang = 1;   // Use with Closure compiler.
const debugEcm = false;
let app;
let digits;
let config;
let fileContents = 0;
let funcnames;
let parens;
let currentInputBox;
let verboseValue;
let prettyValue;
let CunninghamValue;
if (lang)
{
  funcnames =
  [
    "Suma,+,Resta,-,Multiplicación,*,División,/,Resto,%,Potencia,^,Resultado anterior,ans,Parte real,Re(,Parte imaginaria,Im(,Norma\n\nRe(z)^2 + Im(z)^2,Norm(,Número aleatorio\n\nPrimer argumento: mínimo valor del número aleatorio\nSegundo argumento: máximo valor del número aleatorio,Random(",
    "Máximo común divisor\n\nSe pueden usar uno o más argumentos,GCD(,Mínimo común múltiplo\n\nSe pueden usar uno o más argumentos,LCM(,¿El valor es primo?,IsPrime(",
    "Primo siguiente,N(,Primo anterior,B(",
    "Inverso modular\n\nPrimer argumento: valor\nSegundo argumento: módulo,ModInv(,Exponenciación modular\n\nPrimer argumento: base\nSegundo argumento: exponente\nTercer argumento: módulo,ModPow(",
    "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partición,P("
  ];
  parens = "Paréntesis izquierdo,(,Paréntesis derecho,),Unidad imaginaria,i,";
}
else
{
  funcnames =
  [
    "Sum,+,Subtraction,-,Multiplication,*,Division,/,Remainder,%,Power,^,Last answer,ans,Real part,Re(,Imaginary part,Im(,Norm\n\nRe(z)^2 + Im(z)^2,Norm(,Random number\n\nFirst argument: minimum value for random number\nSecond argument: maximum value for random number,Random(",
    "Greatest Common Divisor\n\nOne or more arguments can be used,GCD(,Least Common Multiple\n\nOne or more arguments can be used,LCM(,The value is prime?,IsPrime(",
    "Next prime after,N(,Last prime before,B(",
    "Modular inverse\n\nFirst argument: value\nSecond argument: modulus,ModInv(,Modular power\n\nFirst argument: base\nSecond argument: exponent\nThird argument: modulus,ModPow(",
    "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partition,P("
  ];
  parens = "Left parenthesis,(,Right parenthesis,),Imaginary unit,i,";
}

function getFuncNames()
{
  return funcnames;
}

function getParens()
{
  return parens;
}

function styleButtons(style1, style2)
{
  get("eval").style.display = style1;
  get("factor").style.display = style1;
  get("config").style.display = style1;
  get("functions").style.display = style1;
  get("stop").style.display = style2;
  get("more").style.display = style2;
}

function fromWorker(e)
{
  // First character of e is "1" for intermediate text
  // and it is "2" for end of calculation.
  // It is "9" for saving expression to factor into Web Storage.
  let firstChar = e.substring(0, 1);
  if (firstChar === "8" && debugEcm)
  {
    setStorage("ecmFactors", e.substring(1));
    setStorage("ecmCurve", "");
  }
  else if (firstChar === "7" && debugEcm)
  {
    setStorage("ecmCurve", e.substring(1));
  }
  else if (firstChar === "4")
  {
    get("status").innerHTML = e.substring(1);
  }
  else
  {
    get("result").innerHTML = e.substring(1);
    if (e.substring(0, 1) === "2")
    {   // First character passed from web worker is "2".
      get("status").innerHTML = "";
      styleButtons("inline", "none");  // Enable eval and factor
      hide("modal-more");
      endCalculation();
    }
  }
}

function saveConfig(fromWizard)
{
  config = "1" +   // Batch mode
           verboseValue +
           prettyValue +
           CunninghamValue +
           (get("hex").checked? "1" : "0") +
           (get("kbd")[1].selected? "1" : "0");
  digits = get("digits").value.trim();
  setStorage("ecmConfig", digits+","+config);
  changeInputmode(get("kbd")[1].selected);
}

function comingFromWorker(e)
{
  fromWorker(e.data);
}

function dowork(n)
{
  app = lang + n + (config.charAt(4) === "1"? 16: 0);
  let res = get("result");
  let valueText = get("value").value;
  let helphelp = get("helphelp");
  hide("help");
  show("helphelp");
  helphelp.innerHTML = (lang? "<p>Aprieta el botón <strong>Ayuda</strong> para obtener ayuda para esta aplicación. Apriétalo de nuevo para retornar a la factorización.</p>":
                              "<p>Press the <strong>Help</strong> button to get help about this application. Press it again to return to the factorization.</p>");
  show("result");
  if (valueText === "")
  {
    res.innerHTML = (lang? "Por favor ingrese una expresión." :
                           "Please type an expression.");
    return;
  }
  styleButtons("none", "inline");  // Enable "more" and "stop" buttons
  res.innerHTML = (lang? "Factorizando la expresión..." :
                         "Factoring expression...");
  let param = digits + "," + app + "," + valueText + String.fromCharCode(0);
  callWorker(param);
}

function getCalcURLs()
{
  return ["gaussianW0000.js",
          "gaussian.webmanifest", "gausiano.webmanifest", "gaussian-icon-1x.png", "gaussian-icon-2x.png", "gaussian-icon-4x.png", "gaussian-icon-180px.png", "gaussian-icon-512px.png", "favicon.ico"];
}

function getFormSendValue()
{
  get("userdata").value = "\n" + get("value").value + "\n" + get("result").innerHTML + "\n" + get("status").innerHTML;
}

function popstate(event)
{
  if (get("feedback").style.display == "block" ||
      get("sentOK").style.display == "block" ||
      get("notSent").style.display == "block")
  {     // End feedback.
    show("main");
    hide("feedback");
    hide("sentOK");
    hide("notSent");
    get("value").focus();   
  }
  else if (get("modal-config").style.display == "block")
  {     // End configuration mode.
    hide("modal-config");
  }
}

function startUp()
{
  get("btnSentOK").onclick = function()
  {
    history.back();
  }
  get("btnNotSent").onclick = function()
  {
    history.back();
  }
  get("eval").onclick = function()
  {
    dowork(0);
  };
  get("factor").onclick = function()
  {
    dowork(2);
  };
  get("stop").onclick = function()
  {
    endWorker();
    styleButtons("inline", "none");  // Enable eval and factor
    get("result").innerHTML =
      (lang? "<p>Cálculo detenido por el usuario.</p>" :
             "<p>Calculation stopped by user</p>");
    get("status").innerHTML = "";
  };
  get("more").onclick = function()
  {
    show("modal-more");
  };
  get("config").onclick = function()
  {
    get("digits").value = digits;
    get("pretty").checked = (config.charAt(2) === "1");
    get("cunnin").checked = (config.charAt(3) === "1");  
    get("hex").checked = (config.charAt(4) === "1");  
    get("kbd")[+config.charAt(5)].selected = "selected";
    show("modal-config");
  };
  get("close-config").onclick = function()
  {
    history.back();   // Close configuration mode.
  };
  get("cancel-config").onclick = function()
  {
    history.back();   // Close configuration mode.
  };
  get("save-config").onclick = function()
  {
    saveConfig(false);
    history.back();   // Close configuration mode.
  };
  get("close-more").onclick = function()
  {
    history.back();   // Close configuration mode.
  };
  get("helpbtn").onclick = function()
  {
    if (get("help").style.display === "block" && get("result").innerHTML !== "")
    {
      hide("help");
      show("helphelp");
      show("result");
    }
    else
    {
      show("help");
      hide("helphelp");
      hide("result");
    }
  };
  get("value").onkeydown = function(evt)
  {
    if (evt.ctrlKey && evt.key === "Enter")
    {
      dowork(2);
      evt.stopPropagation();
    }
  };
  get("funccat").onchange = function()
  {
    generateFuncButtons("funccat", "funcbtns");
  };
  get("formlink").onclick = clickFormLink;
  get("formcancel").onclick = function()
  {
    history.back();
  };
  get("comments").oninput = function(_event)
  {
    get("formsend").disabled = (get("comments").value === "");
  };
  get("formsend").onclick = formSend;
  window.onclick = function(event)
  {
    if (event.target === get("modal"))
    {
      hide("modal");
    }
  };
  get("config").onclick = function()
  {
    verboseValue = config.charAt(1);
    prettyValue = config.charAt(2);
    CunninghamValue = config.charAt(3);
    get("digits").value = digits;
    get("hex").checked = (config.charAt(4) === "1");
    history.pushState({id: 5}, "", location.href);
    get("kbd")[+config.charAt(5)].selected = "selected";
    show("modal-config");
  };
  getConfig();
  generateFuncButtons("funccat", "funcbtns");
  registerServiceWorker();
  currentInputBox = get("value");
};
getCalculatorCode("gaussianW0000.js", false);
window.addEventListener("load", startUp);
window.addEventListener("popstate", popstate);
window["fromWorker"] = fromWorker;
