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
/* global endCalculation */
/* global endWorker */
/* global formSend */
/* global generateFuncButtons */
/* global get */
/* global getCalculatorCode */
/* global hide */
/* global registerServiceWorker */
/* global show */
/** @define {number} */ const lang = 1;   // Use with Closure compiler.
let fileContents = 0;
let currentInputBox;
let funcnames;
let parens;
if (lang)
{
  funcnames =
  [
    "Suma,+,Resta,-,Multiplicación,*,División,/,Resto,%,Potencia,^,Raíz cuadrada entera,sqrt(,Raíz entera\n\nPrimer argumento: radicando\nSegundo argumento: orden de la raíz,iroot(,Número aleatorio\n\nPrimer argumento: mínimo valor del número aleatorio\nSegundo argumento: máximo valor del número aleatorio,Random(,Valor absoluto,Abs(,Signo,Sign(",
    "Igual,=,Distinto,!=,Mayor,>,Menor o igual,<=,Menor,<,Mayor o igual,>=",
    "Y lógica, AND ,O lógica, OR ,O exclusiva, XOR ,Negación lógica, NOT ,Desplazamiento a la izquierda\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHL ,Desplazamiento a la derecha\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHR ",
    "Máximo común divisor\n\nSe pueden usar uno o más argumentos,GCD(,Mínimo común múltiplo\n\nSe pueden usar uno o más argumentos,LCM(,¿El valor es primo?,IsPrime(,Cantidad de factores primos,NumFact(,menor divisor primo,MinFact(,mayor divisor primo,MaxFact(,Cantidad de divisores,NumDivs(,Suma de divisores,SumDivs(",
    "Primo siguiente,N(,Primo anterior,B(,Cantidad de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,NumDigits(,Suma de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,SumDigits(,Invertir dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,RevDigits(,Concatenar factores primos\n\nPrimer argumento: modo\n0: Primos no repetidos en forma ascendente\n1: Primos no repetidos en forma descendente\n2: Primos repetidos en forma ascendente\n3: Primos repetidos en forma descendente\nSegundo argumento: valor a factorizar,ConcatFact(",
    "Inverso modular\n\nPrimer argumento: valor\nSegundo argumento: módulo,ModInv(,División modular\n\nPrimer argumento: dividendo\nSegundo argumento: divisor\nTercer argumento: módulo,ModDiv(,Exponenciación modular\n\nPrimer argumento: base\nSegundo argumento: exponente\nTercer argumento: módulo,ModPow(,Indicador de Euler,Totient(,Símbolo de Jacobi\n\nPrimer argumento: valor superior\nSegundo argumento: valor inferior,Jacobi(",
    "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partición,P("
  ];
  parens = "Paréntesis izquierdo,(,Paréntesis derecho,),";
}
else
{
  funcnames =
  [
    "Sum,+,Subtraction,-,Multiplication,*,Division,/,Remainder,%,Power,^,Integer square root,sqrt(,Integer root\n\nFirst argument: radicand\nSecond argument: root order,iroot(,Random number\n\nFirst argument: minimum value for random number\nSecond argument: maximum value for random number,Random(,Absolute value,Abs(,Sign,Sign(",
    "Equal,=,Not equal,!=,Greater,>,Not greater,<=,Less,<,Not less,>=",
    "Logic AND, AND ,Logic OR, OR ,Exclusive OR, XOR ,Logic NOT, NOT ,Shift left\n\nLeft operand: value to shift\nRight operand: number of bits, SHL ,Shift right\n\nLeft operand: value to shift\nRight operand: number of bits, SHR ",
    "Greatest Common Divisor\n\nOne or more arguments can be used,GCD(,Least Common Multiple\n\nOne or more arguments can be used,LCM(,The value is prime?,IsPrime(,Number of prime factors,NumFact(,smallest prime divisor,MinFact(,greatest prime divisor,MaxFact(,Number of divisors,NumDivs(,Sum of divisors,SumDivs(",
    "Next prime after,N(,Last prime before,B(,Number of digits\n\nFirst argument: value\nSecond argument: base,NumDigits(,Sum of digits\n\nFirst argument: value\nSecond argument: base,SumDigits(,Reverse digits\n\nFirst argument: value\nSecond argument: base,RevDigits(,Concatenate prime factors\n\nFirst argument: Mode\n0: No repeated primes in ascending order\n1: No repeated primes in descending order\n2: Repeated primes in ascending order\n3: Repeated primes in descending order\nSecond argument: Value to factor,ConcatFact(",
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
