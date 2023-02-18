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
/* global getCalculatorCode */
/* global initMenubarEvents */
/** @define {number} */ const lang = 1;   // Use with Closure compiler.
const exprTextEs = "Por favor ingrese un número o expresión para ";
const exprTextEn = "Please type a number or expression for the ";
const asmjs = typeof(WebAssembly) === "undefined";
let worker = 0;
let blob;
let fileContents = 0;
let result, dlog, stop, base, pow, mod, digits, main, help, helpbtn, formlink;
let feedback, formfeedback, name, formcancel, formsend, userdata;
let currentInputBox;
let funcnames;
let parens;
if (lang)
{
  funcnames =
  [
    "Suma,+,Resta,-,Multiplicación,*,División,/,Resto,%,Potencia,^,Resultado anterior,ans,Raíz cuadrada entera,sqrt(,Raíz entera\n\nPrimer argumento: radicando\nSegundo argumento: orden de la raíz,iroot(,Número aleatorio\n\nPrimer argumento: mínimo valor del número aleatorio\nSegundo argumento: máximo valor del número aleatorio,Random(,Valor absoluto,Abs(,Signo,Sign(",
    "Igual,=,Distinto,!=,Mayor,>,Menor o igual,<=,Menor,<,Mayor o igual,>=",
    "Y lógica, AND ,O lógica, OR ,O exclusiva, XOR ,Negación lógica, NOT ,Desplazamiento a la izquierda\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHL ,Desplazamiento a la derecha\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHR ",
    "Máximo común divisor\n\nSe pueden usar uno o más argumentos,GCD(,Mínimo común múltiplo\n\nSe pueden usar uno o más argumentos,LCM(,¿El valor es primo?,IsPrime(",
    "Primo siguiente,N(,Primo anterior,B(,Cantidad de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,NumDigits(,Suma de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,SumDigits(,Invertir dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,RevDigits(",
    "Inverso modular\n\nPrimer argumento: valor\nSegundo argumento: módulo,ModInv(,División modular\n\nPrimer argumento: dividendo\nSegundo argumento: divisor\nTercer argumento: módulo,ModDiv(,Exponenciación modular\n\nPrimer argumento: base\nSegundo argumento: exponente\nTercer argumento: módulo,ModPow(,Indicador de Euler,Totient(,Símbolo de Jacobi\n\nPrimer argumento: valor superior\nSegundo argumento: valor inferior,Jacobi(",
    "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partición,P("
  ];
  parens = "Paréntesis izquierdo,(,Paréntesis derecho,),";
}
else
{
  funcnames =
  [
    "Sum,+,Subtraction,-,Multiplication,*,Division,/,Remainder,%,Power,^,Last answer,ans,Integer square root,sqrt(,Integer root\n\nFirst argument: radicand\nSecond argument: root order,iroot(,Random number\n\nFirst argument: minimum value for random number\nSecond argument: maximum value for random number,Random(,Absolute value,Abs(,Sign,Sign(",
    "Equal,=,Not equal,!=,Greater,>,Not greater,<=,Less,<,Not less,>=",
    "Logic AND, AND ,Logic OR, OR ,Exclusive OR, XOR ,Logic NOT, NOT ,Shift left\n\nLeft operand: value to shift\nRight operand: number of bits, SHL ,Shift right\n\nLeft operand: value to shift\nRight operand: number of bits, SHR ",
    "Greatest Common Divisor\n\nOne or more arguments can be used,GCD(,Least Common Multiple\n\nOne or more arguments can be used,LCM(,The value is prime?,IsPrime(",
    "Next prime after,N(,Last prime before,B(,Number of digits\n\nFirst argument: value\nSecond argument: base,NumDigits(,Sum of digits\n\nFirst argument: value\nSecond argument: base,SumDigits(,Reverse digits\n\nFirst argument: value\nSecond argument: base,RevDigits(",
    "Modular inverse\n\nFirst argument: value\nSecond argument: modulus,ModInv(,Modular division\n\nFirst argument: dividend\nSecond argument: divisor\nThird argument: modulus,ModDiv(,Modular power\n\nFirst argument: base\nSecond argument: exponent\nThird argument: modulus,ModPow(,Totient,Totient(,Jacobi symbol\n\nFirst argument: upper value\nSecond argument: lower value,Jacobi(",
    "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partition,P("
  ];
  parens = "Left parenthesis,(,Right parenthesis,),";
}

function get(x)
{
  return document.getElementById(x);
}
function exprText(es, en)
{
  return lang? exprTextEs + es: exprTextEn + en;
}
function callWorker(param)
{
  if (!worker)
  {
    if (!blob)
    {
      if (asmjs)
      {    // Asm.js
        blob = new Blob([fileContents]);
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
      result.innerHTML = e.data.substring(1);
      if (e.data.substring(0, 1) === "2")
      {   // First character passed from web worker is "2".
        dlog.disabled = false;
        stop.disabled = true;
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
  const app = lang + n;
  const baseText = base.value;
  const powText = pow.value;
  const modText = mod.value;
  const digitGroup = digits.value;
  help.style.display = "none";
  result.style.display = "block";
  if (baseText === "")
  {
    result.innerHTML = exprText("la base.", "base.");
    return;
  }
  if (powText === "")
  {
    result.innerHTML = exprText("la potencia.", "power.");
    return;
  }
  if (modText === "")
  {
    result.innerHTML = exprText("el módulo.", "modulus.");
    return;
  }
  dlog.disabled = true;
  stop.disabled = false;
  result.innerHTML = (lang? "Calculando el logaritmo discreto..." :
                            "Computing discrete logarithm...");
  const param = digitGroup + "," + app + "," + baseText + String.fromCharCode(0) + powText +
  String.fromCharCode(0) + modText + String.fromCharCode(0);
  callWorker(param);
}

function endFeedback()
{
  main.style.display = "block";
  feedback.style.display = "none";
  base.focus();   
}

const calcURLs = ["dilogW0000.js",
                "dilog.webmanifest", "logdi.webmanifest", "dilog-icon-1x.png", "dilog-icon-2x.png", "dilog-icon-4x.png", "dilog-icon-180px.png", "dilog-icon-512px.png", "favicon.ico"];

function buttonClick()
{
  const input = currentInputBox;
  input.focus();
  const start = input.selectionStart;
  input.value = input.value.substring(0, start) +
                this.innerText +
                input.value.substring(input.selectionEnd);
    // Place the caret at the end of the appended text.
  input.selectionStart = start + this.innerText.length;
  input.selectionEnd = input.selectionStart;
}
    
function generateFuncButtons(optionCategory, funcButtons)
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
    button.onclick = buttonClick;
    fragment.appendChild(button);
  }
  funcbtns.innerHTML = "";
  funcbtns.appendChild(fragment);
}

function getFormSendValue()
{
  userdata.value = "\nBase = " + base.value + 
        (lang? "\nPotencia = ":"\npower = ") + pow.value +
        (lang? "\nMódulo = ": "\nModulus = ") + mod.value;
}

window.onload = function()
{
  result = get("result");
  dlog = get("dlog");
  stop = get("stop");
  base = get("base");
  pow = get("pow");
  mod = get("mod");
  digits = get("digits");
  help = get("help");
  main = get("main");
  helpbtn = get("helpbtn");
  formlink = get("formlink");
  feedback = get("feedback");
  formfeedback = get("formfeedback");
  name = get("name");
  formcancel = get("formcancel");
  formsend = get("formsend");
  userdata = get("userdata");
  stop.disabled = true;
  dlog.onclick = function()
  {
    dowork(0);
  };
  stop.onclick = function()
  {
    worker.terminate();
    worker = 0;
    dlog.disabled = false;
    stop.disabled = true;
    result.innerHTML = 
      (lang? "<p>Cálculo detenido por el usuario.</p>" :
             "<p>Calculation stopped by user</p>");
  };
  helpbtn.onclick = function()
  {
    help.style.display = "block";
    result.style.display = "none";
  };
  base.onfocus = function()
  {
    currentInputBox = base;
  };
  pow.onfocus = function()
  {
    currentInputBox = pow;
  };
  mod.onfocus = function()
  {
    currentInputBox = mod;
  };
  get("funccat").onchange = function()
  {
    generateFuncButtons("funccat", "funcbtns");
  };
  formlink.onclick = function()
  {
    main.style.display = "none";
    feedback.style.display = "block";
    formfeedback.reset();
    name.focus();
    return false;   // Do not follow the link.
  };
  formcancel.onclick = function()
  {
    endFeedback();
  };
  formsend.onclick = formSend;
  currentInputBox = base;
  generateFuncButtons("funccat", "funcbtns");
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
    fillCache("cacheDILOG");
  }
  initMenubarEvents();
};
getCalculatorCode("dilogW0000.js", false);
