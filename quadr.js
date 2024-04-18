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
/* global getCalculatorCode */
/* global getConfig */
/* global hide */
/* global registerServiceWorker */
/* global show */
/** @define {number} */ const lang = 1;   // Use with Closure compiler.
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
    "Factorial,!,Primorial,#,Fibonacci,F1,Lucas,L1,Partición,P1"
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
    "Factorial,!,Primorial,#,Fibonacci,F1,Lucas,L1,Partition,P1"
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
    get("steps").disabled = false;
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

function dowork(n)
{
  let param;
  let app = lang + n;
  let res = get("result");
  let coefAText = get("coefA").value.trim();
  let coefBText = get("coefB").value.trim();
  let coefCText = get("coefC").value.trim();
  let coefDText = get("coefD").value.trim();
  let coefEText = get("coefE").value.trim();
  let coefFText = get("coefF").value.trim();
  let digitGroup = get("digits").value.trim();
  hide("help");
  show("result");
  let missing = "";
  let zero = String.fromCharCode(0);
  if (coefAText === "")
  {
    missing = (lang? "coeficiente <var>a</var>." : "coefficient <var>a</var>.");
  }
  if (coefBText === "")
  {
    missing = (lang? "coeficiente <var>b</var>." : "coefficient <var>b</var>.");
  }
  if (coefCText === "")
  {
    missing = (lang? "coeficiente <var>c</var>." : "coefficient <var>c</var>.");
  }
  if (coefDText === "")
  {
    missing = (lang? "coeficiente <var>d</var>." : "coefficient <var>d</var>.");
  }
  if (coefEText === "")
  {
    missing = (lang? "coeficiente <var>e</var>." : "coefficient <var>e</var>.");
  }
  if (coefFText === "")
  {
    missing = (lang? "coeficiente <var>f</var>." : "coefficient <var>f</var>.");
  }
  if (missing !== "")
  {
    res.innerHTML = (lang? "Por favor ingrese un número o expresión para el "+missing :
                               "Please type a number or expression for the "+missing);
    return;
  }
  get("solve").disabled = true;
  get("steps").disabled = true;
  get("stop").disabled = false;
  res.innerHTML = (lang? "Resolviendo la ecuación cuadrática..." :
                             "Solving the quadratic equation...");
  param = digitGroup + "," + app + "," + coefAText + zero + coefBText + zero + coefCText + zero +
                                         coefDText + zero + coefEText + zero + coefFText + zero;
  callWorker(param);
}

function moveNext(e, curr, next)
{    
  let nextInput = get(next);
  if (e.key === "Enter" && curr.value.trim().length > 0)
  {
    e.preventDefault();
    nextInput.focus();
    nextInput.setSelectionRange(0, nextInput.value.length);
  }
}

function getCalcURLs()
{
  return ["quadW0000.js",
          "quad.webmanifest", "cuad.webmanifest", "quad-icon-1x.png", "quad-icon-2x.png", "quad-icon-4x.png", "quad-icon-180px.png", "quad-icon-512px.png", "favicon.ico"];
}

function getFormSendValue()
{
  get("userdata").value = "ax^2 + bxy + cy^2 + dx + ey + f = 0" +
                     "\na = " + get("coefA").value.trim() + "\nb = " + get("coefB").value.trim() +
                     "\nc = " + get("coefC").value.trim() + "\nd = " + get("coefD").value.trim() +
                     "\ne = " + get("coefE").value.trim() + "\nf = " + get("coefF").value.trim();  
}

function popstate(event)
{
  if (get("feedback").style.display == "block" ||
      get("sentOK").style.display == "block" ||
      get("notSent").style.display == "block")
  {           // End feedback.
    show("main");
    hide("feedback");
    hide("sentOK");
    hide("notSent");
    get("coefA").focus();   
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
  get("stop").disabled = true;
  get("solve").onclick = function()
  {
    dowork(0);
  };
  get("steps").onclick = function()
  {
    dowork(2);
  };
  get("stop").onclick = function()
  {
    endWorker();
    get("solve").disabled = false;
    get("steps").disabled = false;
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
  get("coefA").onkeydown = function(e)
  {
    moveNext(e, this, "coefB");
  };
  get("coefB").onkeydown = function(e)
  {
    moveNext(e, this, "coefC");
  };
  get("coefC").onkeydown = function(e)
  {
    moveNext(e, this, "coefD");
  };
  get("coefD").onkeydown = function(e)
  {
    moveNext(e, this, "coefE");
  };
  get("coefE").onkeydown = function(e)
  {
    moveNext(e, this, "coefF");
  };
  get("coefF").onkeydown = function(e)
  {
    if (e.key === "Enter" && e.target.value.trim().length > 0)
    {
      e.preventDefault();
      get("coefA").focus();
      dowork(0);
    }
  };
  get("coefA").onfocus = function()
  {
    currentInputBox = get("coefA");
  };
  get("coefB").onfocus = function()
  {
    currentInputBox = get("coefB");
  };
  get("coefC").onfocus = function()
  {
    currentInputBox = get("coefC");
  };
  get("coefD").onfocus = function()
  {
    currentInputBox = get("coefD");
  };
  get("coefE").onfocus = function()
  {
    currentInputBox = get("coefE");
  };
  get("coefF").onfocus = function()
  {
    currentInputBox = get("coefF");
  };
  get("funccat").onchange = function()
  {
    generateFuncButtons("funccat", "funcbtns");
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
  currentInputBox = get("coefA");
  generateFuncButtons("funccat", "funcbtns");
  registerServiceWorker();
};
window.addEventListener("load", startUp);
window.addEventListener("popstate", popstate);
getCalculatorCode("quadW0000.js", false);
window["fromWorker"] = fromWorker;
