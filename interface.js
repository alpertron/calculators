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
/* global clearWizardTextInput */
/* global clickFormLink */
/* global endCalculation */
/* global endWorker */
/* global get */
/* global getConfig */
/* global getVersionText */
/* global formSend */
/* global generateFuncButtons */
/* global getCalculatorCode */
/* global hide */
/* global keyDownOnWizard */
/* global registerServiceWorker */
/* global selectLoop */
/* global setStorage */
/* global setWizardStep */
/* global show */
/* global typedOnWizard */
/* global wizardNext */
/** @define {number} */ const app = 0;   // Use with Closure compiler.
const lang = app % 2;
let fileContents = 0;
let currentInputBox;
let funcnames;
let parens;
let value;
let verboseValue;
let prettyValue;
let CunninghamValue;
let digits;
let config;

if (lang)
{
  if ((app === 4) || (app === 5))
  {  
    funcnames =
    [
      "Suma,+,Resta,-,Multiplicación,*,División,/,Resto,%,Potencia,^,Resultado anterior,ans,Raíz cuadrada entera,sqrt1,Raíz entera\n\nPrimer argumento: radicando\nSegundo argumento: orden de la raíz,iroot2,Número aleatorio\n\nPrimer argumento: mínimo valor del número aleatorio\nSegundo argumento: máximo valor del número aleatorio,Random2,Valor absoluto,Abs1,Signo,Sign1",
      "Igual,=,Distinto,!=,,,Mayor,>,Menor o igual,<=,Menor,<,Mayor o igual,>=",
      "Desplazamiento a la izquierda\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHL ,Desplazamiento a la derecha\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHR ,,,Y lógica, AND ,O lógica, OR ,O exclusiva, XOR ,Negación lógica, NOT ",
      ",,Máximo común divisor\n\nSe pueden usar uno o más argumentos,GCD2,Mínimo común múltiplo\n\nSe pueden usar uno o más argumentos,LCM2,¿El valor es primo?,IsPrime1",
      "Primo siguiente,N1,Primo anterior,B1,Concatenar factores primos\n\nPrimer argumento: modo\n0: Primos no repetidos en forma ascendente\n1: Primos no repetidos en forma descendente\n2: Primos repetidos en forma ascendente\n3: Primos repetidos en forma descendente\nSegundo argumento: valor a factorizar,ConcatFact2,,,Cantidad de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,NumDigits2,Suma de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,SumDigits2,Invertir dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,RevDigits2",
      "Parte entera del cociente\n\nPrimer argumento: dividendo\nSegundo argumento: divisor,FloorDiv2,Módulo\n\nPrimer argumento: valor\nSegundo argumento: módulo,Mod2,Inverso modular\n\nPrimer argumento: valor\nSegundo argumento: módulo,ModInv2,División modular\n\nPrimer argumento: dividendo\nSegundo argumento: divisor\nTercer argumento: módulo,ModDiv3,Exponenciación modular\n\nPrimer argumento: base\nSegundo argumento: exponente\nTercer argumento: módulo,ModPow3,Indicador de Euler,Totient1,Símbolo de Jacobi\n\nPrimer argumento: valor superior\nSegundo argumento: valor inferior,Jacobi2",
      "Factorial,!,,,Primorial,#,Fibonacci,F1,Lucas,L1,Partición,P1",
      "Suma,+,Resta,-,Multiplicación,*,División,/,,,Prefijo hex,0x,10,A,11,B,12,C,13,D,14,E,15,F"
    ];
  }
  else
  {
    funcnames =
    [
      "Suma,+,Resta,-,Multiplicación,*,División,/,Resto,%,Potencia,^,Resultado anterior,ans,Raíz cuadrada entera,sqrt1,Raíz entera\n\nPrimer argumento: radicando\nSegundo argumento: orden de la raíz,iroot2,Número aleatorio\n\nPrimer argumento: mínimo valor del número aleatorio\nSegundo argumento: máximo valor del número aleatorio,Random2,Valor absoluto,Abs1,Signo,Sign1",
      "Variable,x,Inicializar variable,x=,Variable es menor que,x<,,,Variable es menor o igual que,x<=,Contador,c,Contador es menor que,c<,Contador es menor o igual que,c<=,Separador,;",
      "Igual,=,Distinto,!=,,,Mayor,>,Menor o igual,<=,Menor,<,Mayor o igual,>=",
      "Desplazamiento a la izquierda\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHL ,Desplazamiento a la derecha\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHR ,,,Y lógica, AND ,O lógica, OR ,O exclusiva, XOR ,Negación lógica, NOT ",
      ",,Máximo común divisor\n\nSe pueden usar uno o más argumentos,GCD2,Mínimo común múltiplo\n\nSe pueden usar uno o más argumentos,LCM2,¿El valor es primo?,IsPrime1",
      "Primo siguiente,N1,Primo anterior,B1,Concatenar factores primos\n\nPrimer argumento: modo\n0: Primos no repetidos en forma ascendente\n1: Primos no repetidos en forma descendente\n2: Primos repetidos en forma ascendente\n3: Primos repetidos en forma descendente\nSegundo argumento: valor a factorizar,ConcatFact2,,,Cantidad de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,NumDigits2,Suma de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,SumDigits2,Invertir dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,RevDigits2",
      "Parte entera del cociente\n\nPrimer argumento: dividendo\nSegundo argumento: divisor,FloorDiv2,Módulo\n\nPrimer argumento: valor\nSegundo argumento: módulo,Mod2,Inverso modular\n\nPrimer argumento: valor\nSegundo argumento: módulo,ModInv2,División modular\n\nPrimer argumento: dividendo\nSegundo argumento: divisor\nTercer argumento: módulo,ModDiv3,Exponenciación modular\n\nPrimer argumento: base\nSegundo argumento: exponente\nTercer argumento: módulo,ModPow3,Indicador de Euler,Totient1,Símbolo de Jacobi\n\nPrimer argumento: valor superior\nSegundo argumento: valor inferior,Jacobi2",
      "Factorial,!,,,Primorial,#,Fibonacci,F1,Lucas,L1,Partición,P1",
      "Suma,+,Resta,-,Multiplicación,*,División,/,,,Prefijo hex,0x,10,A,11,B,12,C,13,D,14,E,15,F"
    ];
  }
  parens = "Paréntesis izquierdo,(,Paréntesis derecho,),Nueva línea,\u23CE,";
}
else
{
  if ((app === 4) || (app === 5))
  {  
    funcnames =
    [
      "Sum,+,Subtraction,-,Multiplication,*,Division,/,Remainder,%,Power,^,Last answer,ans,Integer square root,sqrt1,Integer root\n\nFirst argument: radicand\nSecond argument: root order,iroot2,Random number\n\nFirst argument: minimum value for random number\nSecond argument: maximum value for random number,Random2,Absolute value,Abs1,Sign,Sign1",
      "Equal,=,Not equal,!=,,,Greater,>,Not greater,<=,Less,<,Not less,>=",
      "Shift left\n\nLeft operand: value to shift\nRight operand: number of bits, SHL ,Shift right\n\nLeft operand: value to shift\nRight operand: number of bits, SHR ,,,Logic AND, AND ,Logic OR, OR ,Exclusive OR, XOR ,Logic NOT, NOT ",
      ",,Greatest Common Divisor\n\nOne or more arguments can be used,GCD2,Least Common Multiple\n\nOne or more arguments can be used,LCM2,The value is prime?,IsPrime1",
      "Next prime after,N1,Last prime before,B1,,,Number of digits\n\nFirst argument: value\nSecond argument: base,NumDigits2,Sum of digits\n\nFirst argument: value\nSecond argument: base,SumDigits2,Reverse digits\n\nFirst argument: value\nSecond argument: base,RevDigits2",
      "Integer part of quotient\n\nFirst argument: dividend\nSecond argument: divisor,FloorDiv2,Modulo\n\nFirst argument: value\nSecond argument: modulo,Mod2,Modular inverse\n\nFirst argument: value\nSecond argument: modulus,ModInv2,Modular division\n\nFirst argument: dividend\nSecond argument: divisor\nThird argument: modulus,ModDiv3,Modular power\n\nFirst argument: base\nSecond argument: exponent\nThird argument: modulus,ModPow3,Totient,Totient1,Jacobi symbol\n\nFirst argument: upper value\nSecond argument: lower value,Jacobi2",
      "Factorial,!,,,Primorial,#,Fibonacci,F1,Lucas,L1,Partition,P1",
      "Sum,+,Subtraction,-,Multiplication,*,Division,/,,,Hex prefix,0x,10,A,11,B,12,C,13,D,14,E,15,F"
    ];
  }
  else
  {
    funcnames =
    [
      "Sum,+,Subtraction,-,Multiplication,*,Division,/,Remainder,%,Power,^,Last answer,ans,Integer square root,sqrt1,Integer root\n\nFirst argument: radicand\nSecond argument: root order,iroot2,Random number\n\nFirst argument: minimum value for random number\nSecond argument: maximum value for random number,Random2,Absolute value,Abs1,Sign,Sign1",
      "Variable,x,Initialize variable,x=,Variable is less than,x<,,,Variable is less or equal than,x<=,Counter,c,Counter is less than,c<,Counter is less or equal than,c<=,Separator,;",
      "Equal,=,Not equal,!=,,,Greater,>,Not greater,<=,Less,<,Not less,>=",
      "Shift left\n\nLeft operand: value to shift\nRight operand: number of bits, SHL ,Shift right\n\nLeft operand: value to shift\nRight operand: number of bits, SHR ,,,Logic AND, AND ,Logic OR, OR ,Exclusive OR, XOR ,Logic NOT, NOT ",
      ",,Greatest Common Divisor\n\nOne or more arguments can be used,GCD2,Least Common Multiple\n\nOne or more arguments can be used,LCM2,The value is prime?,IsPrime1",
      "Next prime after,N1,Last prime before,B1,,,Number of digits\n\nFirst argument: value\nSecond argument: base,NumDigits2,Sum of digits\n\nFirst argument: value\nSecond argument: base,SumDigits2,Reverse digits\n\nFirst argument: value\nSecond argument: base,RevDigits2",
      "Integer part of quotient\n\nFirst argument: dividend\nSecond argument: divisor,FloorDiv2,Modulo\n\nFirst argument: value\nSecond argument: modulo,Mod2,Modular inverse\n\nFirst argument: value\nSecond argument: modulus,ModInv2,Modular division\n\nFirst argument: dividend\nSecond argument: divisor\nThird argument: modulus,ModDiv3,Modular power\n\nFirst argument: base\nSecond argument: exponent\nThird argument: modulus,ModPow3,Totient,Totient1,Jacobi symbol\n\nFirst argument: upper value\nSecond argument: lower value,Jacobi2",
      "Factorial,!,,,Primorial,#,Fibonacci,F1,Lucas,L1,Partition,P1",
      "Sum,+,Subtraction,-,Multiplication,*,Division,/,,,Hex prefix,0x,10,A,11,B,12,C,13,D,14,E,15,F"
    ];
  }
  parens = "Left parenthesis,(,Right parenthesis,),New line,\u23CE,";
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
  get("calc").style.display = style1;
  if (get("calc5") !== null)
  {
    get("calc5").style.display = style1;
  }
  if (get("calcN") !== null)
  {
    get("calcN").style.display = style1;
  }
  if ((app !== 4) && (app !== 5))
  {    // Continued fraction applet does not use wizard.
    get("openwizard").style.display = style1;
  }
  if (get("stop") !== null)
  {
    get("stop").style.display = style2;
  }
}

let calcURLs;

if (app < 2)
{
  calcURLs = ["fsquaresW0000.js",
               "fsquares.webmanifest", "sumcuad.webmanifest", "fsquares-icon-1x.png", "fsquares-icon-2x.png", "fsquares-icon-4x.png", "fsquares-icon-180px.png", "fsquares-icon-512px.png", "favicon.ico"];
}
else if (app < 4)
{
  calcURLs = ["fsquaresW0000.js",
               "fcubes.webmanifest", "sumcubos.webmanifest", "fcubes-icon-1x.png", "fcubes-icon-2x.png", "fcubes-icon-4x.png", "fcubes-icon-180px.png", "fcubes-icon-512px.png", "favicon.ico"];
}
else if (app < 6)
{
  calcURLs = ["fsquaresW0000.js",
               "contfrac.webmanifest", "fraccont.webmanifest", "contfrac-icon-1x.png", "contfrac-icon-2x.png", "contfrac-icon-4x.png", "contfrac-icon-180px.png", "contfrac-icon-512px.png", "favicon.ico"];
}
else
{
  calcURLs = ["fsquaresW0000.js",
               "tsqcubes.webmanifest", "tcuadcub.webmanifest", "tsqcubes-icon-1x.png", "tsqcubes-icon-2x.png", "tsqcubes-icon-4x.png", "tsqcubes-icon-180px.png", "tsqcubes-icon-512px.png", "favicon.ico"];
}

function getCalcURLs()
{
  return calcURLs;
}

function fromWorker(e)
{
  // First character of e is:
  // "1" for intermediate output
  // "2" for end calculation
  // "4" for sending data to status line
  // "6" for pausing calculation and showing the Continue button
  let firstChar = e.substring(0, 1);
  if (firstChar === "4")
  {
    get("status").innerHTML = e.substring(1);
  }
  else if (firstChar === "9")
  {
    console.log(e.substring(1)+ "\n");
  }
  else
  {
    get("result").innerHTML = e.substring(1);
    if (firstChar === "2" || firstChar === "6")
    {   // First character passed from web worker is "2" or "6".
      get("status").innerHTML = "";
      styleButtons("inline", "none");  // Enable buttons that must be enabled when applet is not running
      if (firstChar === "6")
      {
        show("cont");
      }
      else
      {
        endCalculation();
      }
    }
  }
}

function comingFromWorker(e)
{
  fromWorker(e.data);
}

function saveConfig(fromWizard)
{
  if (fromWizard)
  {
    get("hex").checked = get("hexW").checked;
  }
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

function performCalc(from)
{
  let valueA, valueB, valueC, digitGroup;
  let res = get("result");
  let id;
  show("result");
  valueA = get("num").value;
  if (valueA === "")
  {
    if ((app === 4) || (app === 5))
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
  if ((app === 4) || (app === 5))
  {
    valueB = get("delta").value;
    if (valueB === "")
    {
      res.innerHTML = (lang ? "Por favor ingrese un número o expresión para el argumento de la raíz cuadrada." :
                              "Please type a number or expression for square root argument.");
      return;
    }
    valueC = get("den").value;
    if (valueC === "")
    {
      res.innerHTML = (lang ? "Por favor ingrese un número o expresión para el denominador." :
                              "Please type a number or expression for denominator.");
      return;
    }
  }
  digitGroup = get("digits").value;
  hide("help");
  if (app === 0)    // Closure compiler cannot optimize switch, so a series of "if" instructions is used.
  {
    res.innerHTML = "Computing sum of squares...";
  }
  else if (app === 1)
  {
    res.innerHTML = "Calculando suma de cuadrados...";
  }
  else if (app === 4)
  {
    res.innerHTML = "Computing continued fraction expansion...";
  }
  else if (app === 5)
  {
    res.innerHTML = "Calculando desarrollo en fracciones continuas...";
  }
  res.scrollIntoView({behavior: "smooth"});
  let param = "";
  if ((app === 6) || (app === 7))
  {         // Sum of two squares and a power.
    param = from + ",";
  }
  let options = app + ((config.charAt(4) === "1")? 64: 0);
  if ((app === 4) || (app === 5))
  {         // Continued fractions.
    options += get("converg").checked? 32: 0;
  }
  param += digitGroup + "," + options + "," + valueA + String.fromCharCode(0);
  if ((app === 4) || (app === 5))
  {         // Continued fractions.
    param += valueB + String.fromCharCode(0) + valueC + String.fromCharCode(0);
  }
  else
  {
    styleButtons("none", "inline");  // Enable "stop" button
  }
  hide("cont");
  callWorker(param);
  let helphelp = get("helphelp");
  let langName = (typeof(WebAssembly) === "undefined")? "asm.js": "WebAssembly";
  show("helphelp");
  let versionText = getVersionText();
  helphelp.innerHTML = (lang ? "<p>Aprieta el botón <strong>Ayuda</strong> para obtener ayuda para esta aplicación. Apriétalo de nuevo para retornar a esta pantalla. Los usuarios con teclado pueden presionar CTRL+ENTER para comenzar el cálculo. "+versionText+"</p>":
                               "<p>Press the <strong>Help</strong> button to get help about this application. Press it again to return to this screen. Keyboard users can press CTRL+ENTER to start calculation. "+versionText+"</p>");
}

function oneexpr()
{
  get("next").value = (lang? "Hecho": "Done");
  get("wzddesc").innerHTML = (lang? "Paso 1 de 1: Expresión a factorizar": "Step 1 of 1: Expression to factor");
  get("wzdexam").innerHTML = "&nbsp;";
  clearWizardTextInput();
  setWizardStep(9);
}

function getFormSendValue()
{
  let userdata = get("userdata");
  if ((app !== 4) && (app !== 5))
  {     // Not continued fraction application.
    userdata.value = "\n" + get("num").value + "\n" + get("result").innerHTML + "\n" + get("status").innerHTML;
  }
  else
  {     // Continued fraction application.
    userdata.value = "\nnum = " + get("num").value + "\ndelta = " + get("delta").value + "\nden = " + get("den").value;         
  }
}

function popstate(event)
{
  if (get("feedback").style.display === "block" ||
      get("sentOK").style.display === "block" ||
      get("notSent").style.display === "block")
  {         // End feedback.
    show("main");
    hide("feedback");
    hide("sentOK");
    hide("notSent");
    get("num").focus();
  }
  else if (get("wizard") != null &&
           get("wizard").style.display === "block")
  {         // End wizard.
    show("main");
    hide("wizard");
    get("num").focus();
  }
  else if (get("modal-config").style.display === "block")
  {     // End configuration mode.
    hide("modal-config");
  }
}

function startUp()
{
  let param;
  value = get("num");
  if ((app !== 4) && (app !== 5))
  {    // Not continued fraction.
    get("num").onkeydown = function(e)
    {
      let digitGroup = get("digits").value;
      let res = get("result");
      show("result");
      let input = get("num").value;
      let keyCode = e.key;
      if (keyCode === "Enter")
      {  // Used pressed Enter key
        if (input === "")
        {
          res.innerHTML = (lang ? "Por favor ingrese un número o expresión." : "Please type a number or expression.");
          return;
        }
        if (app === 0)
        {
          res.innerHTML = "Computing sum of squares...";
        }
        else if (app === 1)
        {
          res.innerHTML = "Calculando suma de cuadrados...";
        }
        param = digitGroup + "," + app + "," + input + String.fromCharCode(0);
        styleButtons("none", "inline");  // Enable "stop" button
        callWorker(param);
      }
    };
  }
  get("btnSentOK").onclick = function()
  {
    history.back();
  };
  get("btnNotSent").onclick = function()
  {
    history.back();
  };
  get("calc").onclick = function()
  {
    performCalc(3);
  };
  if (get("calc5") !== null)
  {
    get("calc5").onclick = function()
    {
      performCalc(5);
    };
  }
  if (get("calcN") !== null)
  {
    get("calcN").onclick = function()
    {
      let expon = get("expon").value;
      if (!/^\d+$/.test(expon) || parseInt(expon, 10)%2 === 0)
      {
        get("result").innerHTML = lang?"<p>El exponente <var>n</var> debe ser impar.</p>":
                                       "<p>Exponent <var>n</var> must be odd.</p>";
      }
      else
      {
        performCalc(expon);
      }
    };
  }
  if ((app !== 4) && (app !== 5))
  {    // Continued fraction applet does not use wizard.
    get("openwizard").onclick = function()
    {
      get("exprwiz").innerHTML = get("expr").innerHTML;
      hide("main");
      show("wizard");
      get("wzdupper").style.display = "grid";
      get("oneexpr").checked = true;
      get("next").disabled = true;
      get("decW").checked = !get("hexW").checked;
      get("wzdinput").value = "";
      get("wzdinput").focus();
      history.pushState({id: 2}, "", location.href);
      oneexpr();
    };
    get("wzdinput").onkeydown = keyDownOnWizard;
    get("oneexpr").onclick = function()
    {
      oneexpr();
    };
    get("loop").onclick = function()
    {
      selectLoop();
    };
    get("next").onclick = wizardNext;
    get("wzdinput").oninput = typedOnWizard;
    get("cancel").onclick = function()
    {
      history.back();
    };
  }
  if (get("stop") !== null)
  {
    get("stop").onclick = function()
    {
      endWorker();
      styleButtons("inline", "none");  // Enable buttons that have to be enabled when applet is not running.
      get("result").innerHTML =
        (lang ? "<p>Cálculo detenido por el usuario.</p>" :
                "<p>Calculation stopped by user</p>");
      get("status").innerHTML = "";
    };
  }
  get("continue").onclick = function()
  {
    hide("cont");
    callWorker("C");  // Indicate worker that user pressed Continue button.
  };
  get("num").onkeydown = function(event)
  {
    let keyCode = event.key;
    if (keyCode === "Enter" && event.ctrlKey)
    {
      event.preventDefault();          // Do not propagate Enter key.
      performCalc(0);                  // Perform calculation.
    }
    return true;
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
  get("formlink").onclick = clickFormLink;
  get("formcancel").onclick = function()
  {
    history.back();
  };
  get("funccat").onchange = function()
  {
    generateFuncButtons("funccat", "funcbtns");
  };
  if (get("wzdfunccat"))
  {
    get("wzdfunccat").onchange = function()
    {
      generateFuncButtons("wzdfunccat", "wzdfuncbtns");
    };
  }
  get("clrinput").onclick = function()
  {
    if ((app === 4) || (app === 5))
    {
      get("delta").value = "";
      get("den").value = "";
    }
    get("num").value = "";
    get("num").focus();
  };
  get("num").onfocus = function()
  {
    currentInputBox = get("num");
  };
  if (get("wzdinput"))
  {
    get("wzdinput").onfocus = function()
    {
      currentInputBox = get("wzdinput");
    };
  }
  else
  {
    get("delta").onfocus = function()
    {
      currentInputBox = get("delta");
    };
    get("den").onfocus = function()
    {
      currentInputBox = get("den");
    };
  }
  get("comments").oninput = function(_event)
  {
    get("formsend").disabled = (get("comments").value === "");
  };
  get("formsend").onclick = formSend;
  currentInputBox = get("num");
  generateFuncButtons("funccat", "funcbtns");
  if (get("wzdfunccat"))
  {
    generateFuncButtons("wzdfunccat", "wzdfuncbtns");
  }
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
  getConfig();
  registerServiceWorker();
}
getCalculatorCode("fsquaresW0000.js", false);
window.addEventListener("load", startUp);
window.addEventListener("popstate", popstate);
window["fromWorker"] = fromWorker;
