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
/* global clickFormLink */
/* global endWorker */
/* global formSend */
/* global generateFuncButtons */
/* global get */
/* global getCalculatorCode */
/* global hide */
/* global registerServiceWorker */
/* global show */
/** @define {number} */ const lang = 1;   // Use with Closure compiler.
const exprTextEs = "Por favor ingrese un número o expresión para ";
const exprTextEn = "Please type a number or expression for the ";
let fileContents = 0;
let result, dlog, stop, base, pow, mod, digits, helpbtn, formlink;
let formfeedback, formcancel, formsend, userdata;
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

function getFuncNames()
{
  return funcnames;
}

function getParens()
{
  return parens;
}

function exprText(es, en)
{
  return lang? exprTextEs + es: exprTextEn + en;
}

function fromWorker(e)
{
  // First character of e is "1" for intermediate text
  // and it is "2" for end of calculation.
  result.innerHTML = e.substring(1);
  if (e.substring(0, 1) === "2")
  {   // First character passed from web worker is "2".
    dlog.disabled = false;
    stop.disabled = true;
  }
}

function comingFromWorker(e)
{
  fromWorker(e.data);
}

function dowork(n)
{
  const app = lang + n;
  const baseText = base.value;
  const powText = pow.value;
  const modText = mod.value;
  const digitGroup = digits.value;
  hide("help");
  show("result");
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

function getCalcURLs()
{
  return ["dilogW0000.js",
          "dilog.webmanifest", "logdi.webmanifest", "dilog-icon-1x.png", "dilog-icon-2x.png", "dilog-icon-4x.png", "dilog-icon-180px.png", "dilog-icon-512px.png", "favicon.ico"];
}

function getFormSendValue()
{
  userdata.value = "\nBase = " + base.value + 
        (lang? "\nPotencia = ":"\npower = ") + pow.value +
        (lang? "\nMódulo = ": "\nModulus = ") + mod.value;
}

function popstate(event)
{
  if (get("feedback").style.display == "block" ||
      get("sentOK").style.display == "block" ||
      get("notSent").style.display == "block")
  {         // End feedback.
    show("main");
    hide("feedback");
    hide("sentOK");
    hide("notSent");
    base.focus();   
  }
}

function startUp()
{
  result = get("result");
  dlog = get("dlog");
  stop = get("stop");
  base = get("base");
  pow = get("pow");
  mod = get("mod");
  digits = get("digits");
  helpbtn = get("helpbtn");
  formlink = get("formlink");
  formfeedback = get("formfeedback");
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
    endWorker();
    dlog.disabled = false;
    stop.disabled = true;
    result.innerHTML = 
      (lang? "<p>Cálculo detenido por el usuario.</p>" :
             "<p>Calculation stopped by user</p>");
  };
  get("btnSentOK").onclick = function()
  {
    history.back();
  }
  get("btnNotSent").onclick = function()
  {
    history.back();
  }
  helpbtn.onclick = function()
  {
    show("help");
    hide("result");
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
  formlink.onclick = clickFormLink;
  get("comments").oninput = function(_event)
  {
    get("formsend").disabled = (get("comments").value === "");
  };
  formcancel.onclick = function()
  {
    history.back();
  };
  formsend.onclick = formSend;
  currentInputBox = base;
  generateFuncButtons("funccat", "funcbtns");
  registerServiceWorker();
};
getCalculatorCode("dilogW0000.js", false);
window.addEventListener("load", startUp);
window.addEventListener("popstate", popstate);
window["fromWorker"] = fromWorker;
