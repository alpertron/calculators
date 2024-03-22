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
      "Suma,+,Resta,-,Multiplicación,*,División,/,Resto,%,Potencia,^,Resultado anterior,ans,Raíz cuadrada entera,sqrt(,Raíz entera\n\nPrimer argumento: radicando\nSegundo argumento: orden de la raíz,iroot(,Número aleatorio\n\nPrimer argumento: mínimo valor del número aleatorio\nSegundo argumento: máximo valor del número aleatorio,Random(,Valor absoluto,Abs(,Signo,Sign(",
      "Igual,=,Distinto,!=,Mayor,>,Menor o igual,<=,Menor,<,Mayor o igual,>=",
      "Y lógica, AND ,O lógica, OR ,O exclusiva, XOR ,Negación lógica, NOT ,Desplazamiento a la izquierda\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHL ,Desplazamiento a la derecha\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHR ",
      "Máximo común divisor\n\nSe pueden usar uno o más argumentos,GCD(,Mínimo común múltiplo\n\nSe pueden usar uno o más argumentos,LCM(,¿El valor es primo?,IsPrime(",
      "Primo siguiente,N(,Primo anterior,B(,Cantidad de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,NumDigits(,Suma de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,SumDigits(,Invertir dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,RevDigits(",
      "Inverso modular\n\nPrimer argumento: valor\nSegundo argumento: módulo,ModInv(,División modular\n\nPrimer argumento: dividendo\nSegundo argumento: divisor\nTercer argumento: módulo,ModDiv(,Exponenciación modular\n\nPrimer argumento: base\nSegundo argumento: exponente\nTercer argumento: módulo,ModPow(,Indicador de Euler,Totient(,Símbolo de Jacobi\n\nPrimer argumento: valor superior\nSegundo argumento: valor inferior,Jacobi(",
      "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partición,P("
    ];
  }
  else
  {
    funcnames =
    [
      "Suma,+,Resta,-,Multiplicación,*,División,/,Resto,%,Potencia,^,Resultado anterior,ans,Raíz cuadrada entera,sqrt(,Raíz entera\n\nPrimer argumento: radicando\nSegundo argumento: orden de la raíz,iroot(,Número aleatorio\n\nPrimer argumento: mínimo valor del número aleatorio\nSegundo argumento: máximo valor del número aleatorio,Random(,Valor absoluto,Abs(,Signo,Sign(",
      "Variable,x,Inicializar variable,x=,Variable es menor que,x<,Variable es menor o igual que,x<=,Contador,c,Contador es menor que,c<,Contador es menor o igual que,c<=,Separador,;",
      "Igual,=,Distinto,!=,Mayor,>,Menor o igual,<=,Menor,<,Mayor o igual,>=",
      "Y lógica, AND ,O lógica, OR ,O exclusiva, XOR ,Negación lógica, NOT ,Desplazamiento a la izquierda\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHL ,Desplazamiento a la derecha\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHR ",
      "Máximo común divisor\n\nSe pueden usar uno o más argumentos,GCD(,Mínimo común múltiplo\n\nSe pueden usar uno o más argumentos,LCM(,¿El valor es primo?,IsPrime(",
      "Primo siguiente,N(,Primo anterior,B(,Cantidad de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,NumDigits(,Suma de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,SumDigits(,Invertir dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,RevDigits(",
      "Inverso modular\n\nPrimer argumento: valor\nSegundo argumento: módulo,ModInv(,División modular\n\nPrimer argumento: dividendo\nSegundo argumento: divisor\nTercer argumento: módulo,ModDiv(,Exponenciación modular\n\nPrimer argumento: base\nSegundo argumento: exponente\nTercer argumento: módulo,ModPow(,Indicador de Euler,Totient(,Símbolo de Jacobi\n\nPrimer argumento: valor superior\nSegundo argumento: valor inferior,Jacobi(",
      "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partición,P("
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
      "Sum,+,Subtraction,-,Multiplication,*,Division,/,Remainder,%,Power,^,Last answer,ans,Integer square root,sqrt(,Integer root\n\nFirst argument: radicand\nSecond argument: root order,iroot(,Random number\n\nFirst argument: minimum value for random number\nSecond argument: maximum value for random number,Random(,Absolute value,Abs(,Sign,Sign(",
      "Equal,=,Not equal,!=,Greater,>,Not greater,<=,Less,<,Not less,>=",
      "Logic AND, AND ,Logic OR, OR ,Exclusive OR, XOR ,Logic NOT, NOT ,Shift left\n\nLeft operand: value to shift\nRight operand: number of bits, SHL ,Shift right\n\nLeft operand: value to shift\nRight operand: number of bits, SHR ",
      "Greatest Common Divisor\n\nOne or more arguments can be used,GCD(,Least Common Multiple\n\nOne or more arguments can be used,LCM(,The value is prime?,IsPrime(",
      "Next prime after,N(,Last prime before,B(,Number of digits\n\nFirst argument: value\nSecond argument: base,NumDigits(,Sum of digits\n\nFirst argument: value\nSecond argument: base,SumDigits(,Reverse digits\n\nFirst argument: value\nSecond argument: base,RevDigits(",
      "Modular inverse\n\nFirst argument: value\nSecond argument: modulus,ModInv(,Modular division\n\nFirst argument: dividend\nSecond argument: divisor\nThird argument: modulus,ModDiv(,Modular power\n\nFirst argument: base\nSecond argument: exponent\nThird argument: modulus,ModPow(,Totient,Totient(,Jacobi symbol\n\nFirst argument: upper value\nSecond argument: lower value,Jacobi(",
      "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partition,P("
    ];
  }
  else
  {
    funcnames =
    [
      "Sum,+,Subtraction,-,Multiplication,*,Division,/,Remainder,%,Power,^,Last answer,ans,Integer square root,sqrt(,Integer root\n\nFirst argument: radicand\nSecond argument: root order,iroot(,Random number\n\nFirst argument: minimum value for random number\nSecond argument: maximum value for random number,Random(,Absolute value,Abs(,Sign,Sign(",
      "Variable,x,Initialize variable,x=,Variable is less than,x<,Variable is less or equal than,x<=,Counter,c,Counter is less than,c<,Counter is less or equal than,c<=,Separator,;",
      "Equal,=,Not equal,!=,Greater,>,Not greater,<=,Less,<,Not less,>=",
      "Logic AND, AND ,Logic OR, OR ,Exclusive OR, XOR ,Logic NOT, NOT ,Shift left\n\nLeft operand: value to shift\nRight operand: number of bits, SHL ,Shift right\n\nLeft operand: value to shift\nRight operand: number of bits, SHR ",
      "Greatest Common Divisor\n\nOne or more arguments can be used,GCD(,Least Common Multiple\n\nOne or more arguments can be used,LCM(,The value is prime?,IsPrime(",
      "Next prime after,N(,Last prime before,B(,Number of digits\n\nFirst argument: value\nSecond argument: base,NumDigits(,Sum of digits\n\nFirst argument: value\nSecond argument: base,SumDigits(,Reverse digits\n\nFirst argument: value\nSecond argument: base,RevDigits(",
      "Modular inverse\n\nFirst argument: value\nSecond argument: modulus,ModInv(,Modular division\n\nFirst argument: dividend\nSecond argument: divisor\nThird argument: modulus,ModDiv(,Modular power\n\nFirst argument: base\nSecond argument: exponent\nThird argument: modulus,ModPow(,Totient,Totient(,Jacobi symbol\n\nFirst argument: upper value\nSecond argument: lower value,Jacobi(",
      "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partition,P("
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
  let param = "";
  if ((app === 6) || (app === 7))
  {         // Sum of two squares and a power.
    param = from + ",";
  }
  let options = app + ((config.charAt(4) === '1')? 64: 0);
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
  if (get("feedback").style.display == "block" ||
      get("sentOK").style.display == "block" ||
      get("notSent").style.display == "block")
  {         // End feedback.
    show("main");
    hide("feedback");
    hide("sentOK");
    hide("notSent");
    get("num").focus();
  }
  else if (get("wizard") != null &&
           get("wizard").style.display == "block")
  {         // End wizard.
    show("main");
    hide("wizard");
    get("num").focus();
  }
  else if (get("modal-config").style.display == "block")
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
  }
  get("btnNotSent").onclick = function()
  {
    history.back();
  }
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
      show("mode");
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
