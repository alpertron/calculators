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
/* global fillCache */
/* global formSend */
/* global initMenubarEvents */
/** @define {number} */ const lang = 1;   // Use with Closure compiler.
const debugEcm = false;
const asmjs = typeof(WebAssembly) === "undefined";
let worker = 0;
let app;
let blob;
let digits;
let config;
let fileContents = 0;
let funcnames;
let parens;
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
  get("functions").style.display = style1;
  get("stop").style.display = style2;
  get("more").style.display = style2;
}

function callWorker(param)
{
  if (!worker)
  {
    if (!blob)
    {
      if (asmjs)
      {    // Asm.js
        blob = new Blob([fileContents],{type: "text/javascript"});
      }
      else
      {    // WebAssembly
        blob = new Blob([get("worker").textContent],{type: "text/javascript"});
      }
    }   
    worker = new Worker(window.URL.createObjectURL(blob));
    worker.onmessage = function(e)
    { // First character of e.data is "1" for intermediate text
      // and it is "2" for end of calculation.
      // It is "9" for saving expression to factor into Web Storage.
      let firstChar = e.data.substring(0, 1);
      if (firstChar === "8" && debugEcm)
      {
        setStorage("ecmFactors", e.data.substring(1));
        setStorage("ecmCurve", "");
      }
      else if (firstChar === "7" && debugEcm)
      {
        setStorage("ecmCurve", e.data.substring(1));
      }
      else if (firstChar === "4")
      {
        get("status").innerHTML = e.data.substring(1);
      }
      else
      {
        get("result").innerHTML = e.data.substring(1);
        if (e.data.substring(0, 1) === "2")
        {   // First character passed from web worker is "2".
          get("status").innerHTML = "";
          styleButtons("inline", "none");  // Enable eval and factor
          get("modal-more").style.display = "none";
        }
      }
    };
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
  app = lang + n;
  let res = get("result");
  let valueText = get("value").value;
  let helphelp = get("helphelp");
  get("help").style.display = "none";
  helphelp.style.display = "block";
  helphelp.innerHTML = (lang? "<p>Aprieta el botón <strong>Ayuda</strong> para obtener ayuda para esta aplicación. Apriétalo de nuevo para retornar a la factorización.</p>":
                              "<p>Press the <strong>Help</strong> button to get help about this application. Press it again to return to the factorization.</p>");
  res.style.display = "block";
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

function endFeedback()
{
  get("main").style.display = "block";
  get("feedback").style.display = "none";
  get("value").focus();   
}

let calcURLs = ["gaussianW0000.js",
               "gaussian.webmanifest", "gausiano.webmanifest", "gaussian-icon-1x.png", "gaussian-icon-2x.png", "gaussian-icon-4x.png", "gaussian-icon-180px.png", "gaussian-icon-512px.png", "favicon.ico"];

function generateFuncButtons(optionCategory, funcButtons, inputId)
{
  let button;
  let catIndex;
  let funcbtns = get(funcButtons);
  let catnbr = get(optionCategory).selectedIndex;
  let funcname = (parens + funcnames[+catnbr]).split(",");
  // Append all buttons to document fragment instead of funcbtns
  // and finally append the fragment to funcbtns to minimize redraws.
  let fragment = document.createDocumentFragment();
  for (catIndex = 0; catIndex < funcname.length/2; catIndex++)
  {
    button = document.createElement("button");
    button.setAttribute("type", "button");        // Indicate this is a button, not submit.
    button.setAttribute("title", funcname[catIndex*2]);  // Text of tooltip.
    button.innerHTML = funcname[catIndex*2 + 1];         // Text of button.
    button.classList.add("funcbtn");
    button.onclick = function()
    {
      let input = get(inputId);
      input.focus();
      let start = input.selectionStart;
      input.value = input.value.substring(0, start) +
                    this.innerText +
                    input.value.substring(input.selectionEnd);
        // Place the caret at the end of the appended text.
      input.selectionStart = start + this.innerText.length;
      input.selectionEnd = input.selectionStart;
    };
    fragment.appendChild(button);
  }
  funcbtns.innerHTML = "";
  funcbtns.appendChild(fragment);
}

function getFormSendValue()
{
  get("userdata").value = "\n" + get("value").value + "\n" + get("result").innerHTML + "\n" + get("status").innerHTML;
}

window.onload = function()
{
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
    worker.terminate();
    worker = 0;
    styleButtons("inline", "none");  // Enable eval and factor
    get("result").innerHTML =
      (lang? "<p>Cálculo detenido por el usuario.</p>" :
             "<p>Calculation stopped by user</p>");
    get("status").innerHTML = "";
  };
  get("more").onclick = function()
  {
    get("modal-more").style.display = "block";
  };
  get("config").onclick = function()
  {
    get("digits").value = digits;
    get("batch").checked = (config.substring(1,2) === "1");
    get("verbose").checked = (config.substring(1,2) === "1");
    get("pretty").checked = (config.substring(2,3) === "1");
    get("cunnin").checked = (config.substring(3,4) === "1");  
    get("modal-config").style.display = "block";
  };
  get("close-config").onclick = function()
  {
    get("modal-config").style.display = "none";
  };
  get("cancel-config").onclick = function()
  {
    get("modal-config").style.display = "none";
  };
  get("save-config").onclick = function()
  {
    config = (get("batch").checked? "1" :"0") +
             (get("verbose").checked? "1" :"0") +
             (get("pretty").checked? "1" :"0") +
             (get("cunnin").checked? "1" :"0");
    digits = get("digits").value;
    setStorage("ecmConfig", digits+","+config);
    get("modal-config").style.display = "none";
  };
  get("close-more").onclick = function()
  {
    get("modal-more").style.display = "none";
  };
  get("helpbtn").onclick = function()
  {
    let help = get("help");
    let helpStyle = help.style;
    let helphelpStyle = get("helphelp").style;
    let result = get("result");
    let resultStyle = result.style;
    if (helpStyle.display === "block" && result.innerHTML !== "")
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
    generateFuncButtons("funccat", "funcbtns", "value");
  };
  get("formlink").onclick = function()
  {
    get("main").style.display = "none";
    get("feedback").style.display = "block";
    get("formfeedback").reset();
    get("name").focus();
    return false;   // Do not follow the link.
  };
  get("formcancel").onclick = function()
  {
    endFeedback();
  };
  get("formsend").onclick = formSend;
  window.onclick = function(event)
  {
    let modal = get("modal");
    if (event.target === modal)
    {
      modal.style.display = "none";
    }
  };
  digits = getStorage("ecmConfig");
  if (digits === null || digits === "")
  {
    digits = 6;
    config = "0010";
    setStorage("ecmConfig", digits+","+config);
  }
  else
  {
    let index = digits.indexOf(",");
    if (index<0)
    {
      digits = 6;
      config = "0010";
      setStorage("ecmConfig", digits+","+config);
    }
    else
    {
      config = digits.substring(index+1);
      digits = digits.substring(0,index);
    }
  }
  initMenubarEvents();
  generateFuncButtons("funccat", "funcbtns", "value");
  if ("serviceWorker" in navigator)
  { // Attempt to register service worker.
    // There is no need to do anything on registration success or failure in this JavaScript module.
    navigator["serviceWorker"]["register"]("calcSW.js").then(
              function()
              {
                /* Nothing to do */
              },
              function()
              {
                /* Nothing to do */
              });
    fillCache("cacheGAUSS");
  }
};
getCalculatorCode("gaussianW0000.js", false);
