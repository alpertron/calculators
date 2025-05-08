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
/* global completeFuncButtons */
/* global endCalculation */
/* global endWorker */
/* global formSend */
/* global generateFuncButtons */
/* global get */
/* global getCalculatorCode */
/* global getConfig */
/* global hide */
/* global registerServiceWorker */
/* global setStorage */
/* global show */
/** @define {number} */ const lang = 0;   // Use with Closure compiler.
let fileContents = 0;
let currentInputBox;
let funcnames;
let parens;
let verboseValue;
let prettyValue;
let CunninghamValue;
let hexValue;
let digits;
let config;
if (lang)
{
  funcnames =
  [
    "Suma,+,Resta,-,Multiplicación,*,División,/,Resto,%,Potencia,^,Resultado anterior,ans,Raíz cuadrada entera,sqrt1,Raíz entera\n\nPrimer argumento: radicando\nSegundo argumento: orden de la raíz,iroot2,Número aleatorio\n\nPrimer argumento: mínimo valor del número aleatorio\nSegundo argumento: máximo valor del número aleatorio,Random2,Valor absoluto,Abs1,Signo,Sign1",
    "Igual,=,Distinto,!=,Mayor,>,Menor o igual,<=,Menor,<,Mayor o igual,>=",
    "Y lógica, AND ,O lógica, OR ,O exclusiva, XOR ,Negación lógica, NOT ,Desplazamiento a la izquierda\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHL ,Desplazamiento a la derecha\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHR ",
    "Máximo común divisor\n\nSe pueden usar uno o más argumentos,GCD2,Mínimo común múltiplo\n\nSe pueden usar uno o más argumentos,LCM2,¿El valor es primo?,IsPrime1",
    "Primo siguiente,N1,Primo anterior,B1,Cantidad de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,NumDigits2,Suma de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,SumDigits2,Invertir dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,RevDigits2",
    "Parte entera del cociente\n\nPrimer argumento: dividendo\nSegundo argumento: divisor,FloorDiv2,Módulo\n\nPrimer argumento: valor\nSegundo argumento: módulo,Mod2,Inverso modular\n\nPrimer argumento: valor\nSegundo argumento: módulo,ModInv2,División modular\n\nPrimer argumento: dividendo\nSegundo argumento: divisor\nTercer argumento: módulo,ModDiv3,Exponenciación modular\n\nPrimer argumento: base\nSegundo argumento: exponente\nTercer argumento: módulo,ModPow3,Indicador de Euler,Totient1,Símbolo de Jacobi\n\nPrimer argumento: valor superior\nSegundo argumento: valor inferior,Jacobi2",
    "Factorial,!,Primorial,#,Fibonacci,F1,Lucas,L1,Partición,P1",
    "Suma,+,Resta,-,Multiplicación,*,División,/,Prefijo hex,0x,10,A,11,B,12,C,13,D,14,E,15,F"
  ];
  parens = "Paréntesis izquierdo,(,Paréntesis derecho,),";
}
else
{
  funcnames =
  [
    "Sum,+,Subtraction,-,Multiplication,*,Division,/,Remainder,%,Power,^,Last answer,ans,Integer square root,sqrt1,Integer root\n\nFirst argument: radicand\nSecond argument: root order,iroot2,Random number\n\nFirst argument: minimum value for random number\nSecond argument: maximum value for random number,Random2,Absolute value,Abs1,Sign,Sign1",
    "Equal,=,Not equal,!=,Greater,>,Not greater,<=,Less,<,Not less,>=",
    "Logic AND, AND ,Logic OR, OR ,Exclusive OR, XOR ,Logic NOT, NOT ,Shift left\n\nLeft operand: value to shift\nRight operand: number of bits, SHL ,Shift right\n\nLeft operand: value to shift\nRight operand: number of bits, SHR ",
    "Greatest Common Divisor\n\nOne or more arguments can be used,GCD2,Least Common Multiple\n\nOne or more arguments can be used,LCM2,The value is prime?,IsPrime1",
    "Next prime after,N1,Last prime before,B1,Number of digits\n\nFirst argument: value\nSecond argument: base,NumDigits2,Sum of digits\n\nFirst argument: value\nSecond argument: base,SumDigits2,Reverse digits\n\nFirst argument: value\nSecond argument: base,RevDigits2",
    "Integer part of quotient\n\nFirst argument: dividend\nSecond argument: divisor,FloorDiv2,Modulo\n\nFirst argument: value\nSecond argument: modulo,Mod2,Modular inverse\n\nFirst argument: value\nSecond argument: modulus,ModInv2,Modular division\n\nFirst argument: dividend\nSecond argument: divisor\nThird argument: modulus,ModDiv3,Modular power\n\nFirst argument: base\nSecond argument: exponent\nThird argument: modulus,ModPow3,Totient,Totient1,Jacobi symbol\n\nFirst argument: upper value\nSecond argument: lower value,Jacobi2",
    "Factorial,!,Primorial,#,Fibonacci,F1,Lucas,L1,Partition,P1",
    "Sum,+,Subtraction,-,Multiplication,*,Division,/,Hex prefix,0x,10,A,11,B,12,C,13,D,14,E,15,F"
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

function fromWorker(e)
{
  // First character of e is "1" for intermediate text
  // and it is "2" for end of calculation.
  get("result").innerHTML = e.substring(1);
  if (e.substring(0, 1) === "2")
  {   // First character passed from web worker is "2".
    get("solve").disabled = false;
    get("stop").disabled = true;
    endCalculation();
  }
}

function saveConfig()
{
  config = "1" +   // Batch mode
           verboseValue +
           prettyValue +
           CunninghamValue +
           hexValue +
           (get("kbd")[1].selected? "1" : "0");
  digits = get("digits").value.trim();
  setStorage("ecmConfig", digits+","+config);
  changeInputmode(get("kbd")[1].selected);
}

function comingFromWorker(e)
{
  fromWorker(e.data);
}

function dowork()
{
  let param;
  let res = get("result");
  let quadrText = get("quad").value.trim();
  let linText = get("lin").value.trim();
  let constText = get("const").value.trim();
  let modText = get("mod").value.trim();
  let digitGroup = digits;
  hide("help");
  show("result");
  let missing = "";
  if (quadrText === "")
  {
    missing = (lang? "coeficiente cuadrático." : "quadratic coefficient.");
  }
  if (linText === "")
  {
    missing = (lang? "coeficiente lineal." : "linear coefficient.");
  }
  if (constText === "")
  {
    missing = (lang? "término independiente." : "constant coefficient.");
  }
  if (modText === "")
  {
    missing = (lang? "módulo." : "modulus.");
  }
  if (missing !== "")
  {
    res.innerHTML = (lang? "Por favor ingrese un número o expresión para el "+missing :
                           "Please type a number or expression for the "+missing);
    return;
  }
  get("solve").disabled = true;
  get("stop").disabled = false;
  res.innerHTML = (lang? "Resolviendo la ecuación cuadrática..." :
                         "Solving the quadratic equation...");
  param = digitGroup + "," + lang + "," + quadrText + String.fromCharCode(0) + linText +
    String.fromCharCode(0) + constText +String.fromCharCode(0) + modText + String.fromCharCode(0);
  callWorker(param);
}

function getCalcURLs()
{
  return ["quadmodW0000.js",
          "quadmod.webmanifest", "cuadmod.webmanifest", "quadmod-icon-1x.png", "quadmod-icon-2x.png", "quadmod-icon-4x.png", "quadmod-icon-180px.png", "quadmod-icon-512px.png", "favicon.ico"];
}

function getFormSendValue()
{
  get("userdata").value = "ax^2 + bx + c = 0 (mod n)" + 
                     "\na = " + get("quad").value + "\nb = " + get("lin").value +
                     "\nc = " + get("const").value + "\nn = " + get("mod").value;
}

function popstate(event)
{
  if (get("feedback").style.display === "block" ||
      get("sentOK").style.display === "block" ||
      get("notSent").style.display === "block")
  {
    show("main");
    hide("feedback");
    hide("sentOK");
    hide("notSent");
    get("quad").focus();   
  }
  else if (get("modal-config").style.display === "block")
  {     // End configuration mode.
    hide("modal-config");
  }
}

function startUp()
{
  get("btnSentOK").onclick = function()
  {
    history.back();
  };
  get("btnNotSent").onclick = function()
  {
    history.back();
  };
  get("stop").disabled = true;
  get("solve").onclick = dowork;
  get("stop").onclick = function()
  {
    endWorker();
    get("solve").disabled = false;
    get("stop").disabled = true;
    get("result").innerHTML = 
      (lang? "<p>Cálculo detenido por el usuario.</p>" :
             "<p>Calculation stopped by user</p>");
  };
  get("helpbtn").onclick = function()
  {
    show("help");
    hide("result");
  };
  get("quad").onfocus = function()
  {
    currentInputBox = get("quad");
  };
  get("lin").onfocus = function()
  {
    currentInputBox = get("lin");
  };
  get("const").onfocus = function()
  {
    currentInputBox = get("const");
  };
  get("mod").onfocus = function()
  {
    currentInputBox = get("mod");
  };
  get("funccat").onchange = function()
  {
    generateFuncButtons("funccat", "funcbtns");
  };
  get("clrinput").onclick = function()
  {
    get("quad").value = "";
    get("lin").value = "";
    get("const").value = "";
    get("mod").value = "";
    get("quad").focus();
  };
  get("config").onclick = function()
  {
    verboseValue = config.charAt(1);
    prettyValue = config.charAt(2);
    CunninghamValue = config.charAt(3);
    hexValue = config.charAt(4);
    get("digits").value = digits;
    history.pushState({id: 5}, "", location.href);
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
    saveConfig();
    history.back();   // Close configuration mode.
  };
  getConfig();
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
  currentInputBox = get("quad");
  registerServiceWorker();
  completeFuncButtons("funcbtns");
}

window.addEventListener("load", startUp);
window.addEventListener("popstate", popstate);
getCalculatorCode("quadmodW0000.js", false);
window["fromWorker"] = fromWorker;
