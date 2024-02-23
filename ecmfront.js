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
/* global clearWizardTextInput */
/* global clickFormLink */
/* global comingFromPolfact */
/* global completeFuncButtons */
/* global endWorker */
/* global formSend */
/* global generateFuncButtons */
/* global get */
/* global getCalculatorCode */
/* global getStorage */
/* global hide */
/* global keyDownOnWizard */
/* global loadPolyCalc */
/* global registerServiceWorker */
/* global selectLoop */
/* global setStorage */
/* global setWizardStep */
/* global show */
/* global showVersion */
/* global typedOnWizard */
/* global useBlockly */
/* global wizardNext */
/** @define {number} */ const lang = 1;        // Use with Closure compiler.
const points=[0,6, 2,9, 4,0, 5,6, 7,1, 8,0, 13,9, 14,9, 15,7, 16,7, 17,0, 18,13, 20,5, 22,10, 23,12, 24,6, 27,7];
let fileContents = null;
let app;
let digits;
let config;
let fromFile;
let tofile;
let fileName;
let workerParam;
let bmodeLoaded = 0;
let statusText = "";
let resultText = "";
let divisorsDirty = false;
let statusDirty = false;
let resultDirty = false;
let blocklyLoaded = 0;
let scriptsLoaded = 0;
let script1;
let script2;
let funcnames;
let parens;
let currentInputBox;

// DOM resources
let value;
let btnNext;
let btnEval;
let btnPrime;
let btnFactor;
let btnConfig;
let btnFromFile;
let btnBlocklyMode;
let btnOpenWizard;
let btnMore;
let btnToFile;
let btnStop;
let chkCunningham;
let chkDecW;
let chkHex;
let chkHexW;
let chkPretty;
let chkVerbose;
let divResult;
let getFile;
let newCurveOrFactor;
let wzdDescText;
let wzdExamText;
let wzdInput;

if (lang)
{
  funcnames =
  [
    "Suma,+,Resta,-,Multiplicación,*,División,/,Resto,%,Potencia,^,Resultado anterior,ans,Raíz cuadrada entera,sqrt(,Número aleatorio\n\nPrimer argumento: mínimo valor del número aleatorio\nSegundo argumento: máximo valor del número aleatorio,Random(,Valor absoluto,Abs(,Signo,Sign(,Variable,x,Contador,c",
    "Igual,=,Distinto,!=,Mayor,>,Menor o igual,<=,Menor,<,Mayor o igual,>=",
    "Y lógica, AND ,O lógica, OR ,O exclusiva, XOR ,Negación lógica, NOT ,Desplazamiento a la izquierda\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHL ,Desplazamiento a la derecha\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHR ",
    "Máximo común divisor\n\nSe pueden usar uno o más argumentos,GCD(,Mínimo común múltiplo\n\nSe pueden usar uno o más argumentos,LCM(,¿El valor es primo?,IsPrime(,Cantidad de factores primos,NumFact(,menor divisor primo,MinFact(,mayor divisor primo,MaxFact(,Cantidad de divisores,NumDivs(,Suma de divisores,SumDivs(",
    "Primo siguiente,N(,Primo anterior,B(,Cantidad de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,NumDigits(,Suma de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,SumDigits(,Invertir dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,RevDigits(,Concatenar factores primos\n\nPrimer argumento: modo\n0: Primos no repetidos en forma ascendente\n1: Primos no repetidos en forma descendente\n2: Primos repetidos en forma ascendente\n3: Primos repetidos en forma descendente\nSegundo argumento: valor a factorizar,ConcatFact(",
    "Inverso modular\n\nPrimer argumento: valor\nSegundo argumento: módulo,ModInv(,División modular\n\nPrimer argumento: dividendo\nSegundo argumento: divisor\nTercer argumento: módulo,ModDiv(,Exponenciación modular\n\nPrimer argumento: base\nSegundo argumento: exponente\nTercer argumento: módulo,ModPow(,Indicador de Euler,Totient(,Símbolo de Jacobi\n\nPrimer argumento: valor superior\nSegundo argumento: valor inferior,Jacobi(",
    "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partición,P("
  ];
  parens = "Paréntesis izquierdo,(,Paréntesis derecho,),Nueva línea,\u23CE,";
}
else
{
  funcnames =
  [
    "Sum,+,Subtraction,-,Multiplication,*,Division,/,Remainder,%,Power,^,Last answer,ans,Integer square root,sqrt(,Integer root\n\nFirst argument: radicand\nSecond argument: root order,iroot(,Random number\n\nFirst argument: minimum value for random number\nSecond argument: maximum value for random number,Random(,Absolute value,Abs(,Sign,Sign(,Variable,x,Counter,c",
    "Equal,=,Not equal,!=,Greater,>,Not greater,<=,Less,<,Not less,>=",
    "Logic AND, AND ,Logic OR, OR ,Exclusive OR, XOR ,Logic NOT, NOT ,Shift left\n\nLeft operand: value to shift\nRight operand: number of bits, SHL ,Shift right\n\nLeft operand: value to shift\nRight operand: number of bits, SHR ",
    "Greatest Common Divisor\n\nOne or more arguments can be used,GCD(,Least Common Multiple\n\nOne or more arguments can be used,LCM(,The value is prime?,IsPrime(,Number of prime factors,NumFact(,smallest prime divisor,MinFact(,greatest prime divisor,MaxFact(,Number of divisors,NumDivs(,Sum of divisors,SumDivs(",
    "Next prime after,N(,Last prime before,B(,Number of digits\n\nFirst argument: value\nSecond argument: base,NumDigits(,Sum of digits\n\nFirst argument: value\nSecond argument: base,SumDigits(,Reverse digits\n\nFirst argument: value\nSecond argument: base,RevDigits(,Concatenate prime factors\n\nFirst argument: Mode\n0: No repeated primes in ascending order\n1: No repeated primes in descending order\n2: Repeated primes in ascending order\n3: Repeated primes in descending order\nSecond argument: Value to factor,ConcatFact(",
    "Modular inverse\n\nFirst argument: value\nSecond argument: modulus,ModInv(,Modular division\n\nFirst argument: dividend\nSecond argument: divisor\nThird argument: modulus,ModDiv(,Modular power\n\nFirst argument: base\nSecond argument: exponent\nThird argument: modulus,ModPow(,Totient,Totient(,Jacobi symbol\n\nFirst argument: upper value\nSecond argument: lower value,Jacobi(",
    "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partition,P("
  ];
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

function getCalcURLs()
{
  return ["ecmW0000.js",
          "ecm.webmanifest", "ecmc.webmanifest", "ecm-icon-1x.png", "ecm-icon-2x.png", "ecm-icon-4x.png", "ecm-icon-180px.png", "ecm-icon-512px.png", "favicon.ico"];
}

function oneexpr()
{
  btnNext.value = (lang? "Hecho": "Done");
  wzdDescText.innerHTML = (lang? "Paso 1 de 1: Expresión a factorizar": "Step 1 of 1: Expression to factor");
  wzdExamText.innerHTML = "&nbsp;";
  clearWizardTextInput();
  setWizardStep(9);
}

function styleButtons(style1, style2)
{
  btnEval.style.display = style1;
  btnPrime.style.display = style1;
  btnFactor.style.display = style1;
  btnConfig.style.display = style1;
  btnFromFile.style.display = style1;
  btnBlocklyMode.style.display = style1;
  btnOpenWizard.style.display = style1;
  get("functions").style.display = style1;
  get("funcbtns").style.display = style1;
  btnStop.style.display = style2;
  btnMore.style.display = style2;
}

function saveConfig(fromWizard)
{
  if (fromWizard)
  {
    chkHex.checked = chkHexW.checked;
  }
  config = "1" +   // Batch mode
           (chkVerbose.checked? "1" : "0") +
           (chkPretty.checked? "1" : "0") +
           (chkCunningham.checked? "1" : "0") +
           (chkHex.checked? "1" : "0");
  digits = get("digits").value;
  setStorage("ecmConfig", digits+","+config);
}

function showSumSquares()
{
  callWorker("S");  // Indicate worker that user pressed Sum of squares button.
}

function fromWorker(e)
{
  // First character of e is:
  // "1" for intermediate output
  // "2" for ending calculation
  // "4" for sending intermediate data
  // "6" for pausing calculation and showing the Continue button
  // "7" for saving curve number into local storage
  // "8" for saving input expression into local storage
  // "9" for sending data to console.
  // "A" for pausing calculation and showing the Continue button (save file)
  // "B" for sending data to be saved to file and ending calculation.
  // "D" for sending data to div named divisors.
  // "E" for sending data to div named divisors. It includes button More divisors.
  // "K" for showing Blockly errors.
  // "L" for exiting Blockly mode.
  // "M" for loading polynomial factorization application for factorization.
  // "N" for loading polynomial factorization application for evaluation.
  // "S" for sending data to div named sumSquares.
  let firstChar = e.substring(0, 1);
  if (firstChar === "9")
  {
    console.log(e.substring(1));
  }
  else if (firstChar === "8")
  {
    setStorage("ecmFactors", e.substring(1));
    setStorage("ecmCurve", "");
  }
  else if (firstChar === "7")
  {
    setStorage("ecmCurve", e.substring(1));
  }
  else if (firstChar === "D")
  {    // Show divisors.
    get("divisors").innerHTML = e.substring(1);
  }
  else if (firstChar === "E")
  {
    get("divisors").innerHTML = e.substring(1);
    get("showdiv").onclick = function()
    {
      callWorker("D");  // Indicate worker that user pressed Divisors button.
    };
  }
  else if (firstChar === "K")
  {
    get("berror").innerHTML = e.substring(1);
    show("BlocklyErrors");
    hide("BlocklyButtons");
  }
  else if (firstChar === "L")
  {   // Exit Blockly mode.
    history.back();
  }
  else if ((firstChar === "M") || (firstChar === "N"))
  {    // User entered a polynomial. Load calculator to process it.
    loadPolyCalc(firstChar, value.value);
  }
  else if (firstChar === "S")
  {    // Show sum of squares.
    get("sumSquares").innerHTML = e.substring(1);
    get("showSumSq").onclick = showSumSquares;
  }
  else if (firstChar === "T")
  {    // Show sum of squares without button.
    get("sumSquares").innerHTML = e.substring(1);
  }
   else if (firstChar === "4")
  {
    statusDirty = true;
    statusText = e.substring(1);
  }
  else if (firstChar === "5")
  {
    if (e.substring(1, 2) === "1")
    {
      show("sktest");
    }
    else
    {
      hide("sktest");
    }
  }
  else
  {
    resultDirty = true;
    if (firstChar === "2" || firstChar === "B" ||
        firstChar === "6" || firstChar === "A")
    {   // First character passed from web worker is "2".
      statusDirty = true;
      statusText = "";
      styleButtons("inline", "none");  // Enable eval and factor
      if (firstChar === "A" || firstChar === "B")
      {
        tofile = e.substring(1);
        const regexVerbose = /\<span class=\"verbose\">(.*?)\<\/span>/g;
        const regexSpan = /\<span(.*?)>(.*?)\<\/span>/g;
        const regexAbbr = /\<abbr(.*?)>(.*?)\<\/abbr>/g;
        if (config.charAt(1) == "1")
        {    // Discard HTML tags for verbose mode.
          tofile = tofile.replace(regexVerbose, "$1").
                          replace(regexAbbr, "$2").
                          replace(regexSpan, "$2");
        }
        else
        {    // Discard HTML tags for non-verbose mode.
          tofile = tofile.replace(regexVerbose, "").
                          replace(regexSpan, "$2");
        }        
        show("savefile");
        resultText = "";
      }
      else
      {
        resultText = e.substring(1);
      }
      if (firstChar === "6" || firstChar === "A")
      {
        show("cont");
      }
      if (firstChar === "2")
      {
        divisorsDirty = true;
        if (navigator.share)
        {
          show("sharediv");
        }
      }
    }
    else
    {
      resultText = e.substring(1);
    }
  }
}

function comingFromWorker(e)
{
  fromWorker(e.data);
}

function performWork(n, valueText)
{
  let param;
  app = lang + n;
  let charNull = String.fromCharCode(0);
  let helphelp = get("helphelp");
  hide("sharediv");
  if (valueText === "")
  {    // Nothing in input box.
    resultDirty = true;
    resultText = (lang ? "<p>Por favor ingrese una expresión.</p>" :
                         "<p>Please type an expression.</p>");
    return;
  }
  hide("cont");
  hide("help");
  show("helphelp");
  let strHelp = (lang ? "<p class=\"pad\">Aprieta el botón <strong>Ayuda</strong> para obtener ayuda para esta aplicación. Apriétalo de nuevo para retornar a la factorización. También puedes ver <a href=\"/videos/videosEcmc.htm\" target=\"_blank\">videos</a>. Los usuarios con teclado pueden presionar CTRL+ENTER para comenzar la factorización.":
                        "<p class=\"pad\">Press the <strong>Help</strong> button to get help about this application. Press it again to return to the factorization. You can also watch <a href=\"/videos/videosEcm.htm\" target=\"_blank\">videos</a>. Keyboard users can press CTRL+ENTER to start factorization.");
  let strHelp2 = showVersion(lang);
  if (strHelp2 == null)
  {
    return;
  }
  strHelp += strHelp2;
  helphelp.innerHTML = strHelp;
  show("result");
  styleButtons("none", "inline");  // Enable "more" and "stop" buttons
  if (n === 0)
  {
    resultDirty = true;
    resultText = (lang ? "<p>Evaluando la expresión...</p>" :
                         "<p>Evaluating expression...</p>");
  }
  else
  {
    resultDirty = true;
    resultText = (lang ? "<p>Factorizando la expresión...</p>" :
                         "<p>Factoring expression...</p>");
  }
  if (n < -2)
  {
    app += 6;   // Convert to factorization.
  }
  param = digits + "," + app + "," + fromFile + config.substring(1) +
          valueText + charNull + getStorage("ecmFactors");
  if (n === -1 || n === -2)
  {
    param += "," + newCurveOrFactor.value;        // Append new curve number typed by user.
  }
  if (n === -3 || n === -4)
  {
    param += ";" + newCurveOrFactor.value;        // Append new factor typed by user.
  }
  if (!fileContents)
  {
    workerParam = param + charNull;
  }
  else
  {
    callWorker(param + charNull);
  }
}

function dowork(n)
{
  let valueText;
  fromFile = "0";
  if (getFile.value !== "")
  {
    let fileReader = new FileReader();
    fileReader.onload = function(fileLoadedEvent) 
    {
      fromFile = "1";
      valueText = fileLoadedEvent.target.result;
      performWork(n, valueText);
      getFile.value = "";
    };
    fileReader.readAsText(getFile.files[0], "UTF-8");
    value.value = "";
  }
  else
  {
    valueText = value.value.replace(/\u2011/g, "-");
    performWork(n, valueText);
  }
}

function restartFactorization(type)
{
  endWorker();
  dowork(type);
}

function updateVerbose(isVerbose)    
{
  let cssRules = document.styleSheets[0]["cssRules"];
  let index;
  let len = cssRules.length;
  for (index=0; index<len; index++)
  {
    if (cssRules[index >> 0].selectorText === ".verbose")
    {
      cssRules[index >> 0].style["display"] = (isVerbose? "inline": "none");
    }
    if (cssRules[index >> 0].selectorText === ".terse")
    {
      cssRules[index >> 0].style["display"] = (isVerbose? "none": "inline");
    }
  }
}

function fromBlocklyRun(xml)
{
  fromFile = "0";
  performWork(8, xml);
}

function loadScript(scriptUrl)
{
  let myScript = document.createElement("script");
  let xmlhttp = new XMLHttpRequest();
  xmlhttp.open("GET", scriptUrl);
  xmlhttp.onreadystatechange = function()
  {
    if ((xmlhttp.status === 200) && (xmlhttp.readyState === 4))
    {
      scriptsLoaded++;
      myScript.innerHTML = xmlhttp.responseText;
      if (scriptsLoaded === 2)
      {
        document.body.appendChild(script1);
        document.body.appendChild(script2);
        useBlockly(fromBlocklyRun, lang);  // Init Blockly workspace.
      }
    }
  };
  xmlhttp.send();
  return myScript;
}

function initBlockly()
{
  if (blocklyLoaded !== 0)
  {
    useBlockly(null, lang);  // Resize workspace.
    return;
  }
  blocklyLoaded = 1;
  script1 = loadScript("blockly0007.js");
  script2 = loadScript(lang? "es0007.js": "en0007.js");
}

function getFormSendValue()
{
  get("userdata").value = "\n" + value.value + "\n" + divResult.innerHTML + "\n" + get("status").innerHTML;
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
    value.focus();   
  }
  else if (get("wizard").style.display == "block")
  {     // End wizard.
    show("main");
    hide("wizard");
    get("value").focus();
  }
  else if (get("blockmode").style.display == "flex")
  {     // End blockly mode.
    document.activeElement.blur();
    show("main");
    hide("blockmode");
  }
  else if (get("modal-more").style.display == "block")
  {     // End "more" mode.
    hide("modal-more");
  }
  else if (get("modal-config").style.display == "block")
  {     // End configuration mode.
    hide("modal-config");
  }
}

function startUp()
{
  let index;
  value = get("value");
  btnNext = get("next");
  btnEval = get("eval");
  btnPrime = get("prime");
  btnFactor = get("factor");
  btnConfig = get("config");
  btnFromFile = get("fromfile");
  btnBlocklyMode = get("bmode");
  btnMore = get("more");
  btnOpenWizard = get("openwizard");
  btnToFile = get("tofile");
  btnStop = get("stop");
  chkCunningham = get("cunnin");
  chkDecW = get("decW");
  chkHex = get("hex");
  chkHexW = get("hexW");
  chkPretty = get("pretty");
  chkVerbose = get("verbose");
  divResult = get("result");
  getFile = get("getFile");
  newCurveOrFactor = get("curve");
  wzdDescText = get("wzddesc");
  wzdExamText = get("wzdexam");
  wzdInput = get("wzdinput");

  app = lang;
  value.wrap="off";
  btnEval.onclick = function()
  {
    setStorage("ecmFactors","");
    dowork(0);
  };
  btnPrime.onclick = function()
  {
    setStorage("ecmFactors","");
    dowork(6);
  };
  btnFactor.onclick = function()
  {
    setStorage("ecmFactors","");
    dowork(2);
  };
  btnMore.onclick = function()
  {
    get("ncurve").disabled = false;
    history.pushState({id: 4}, "", location.href);
    show("modal-more");
  };
  btnConfig.onclick = function()
  {
    get("digits").value = digits;
    chkVerbose.checked = (config.substring(1,2) === "1");
    chkPretty.checked = (config.substring(2,3) === "1");
    chkCunningham.checked = (config.substring(3,4) === "1");  
    chkHex.checked = (config.substring(4,5) === "1");
    history.pushState({id: 5}, "", location.href);
    show("modal-config");
  };
  btnFromFile.onclick = function()
  {
    getFile.click();
  };
  btnToFile.onclick = downloadResult;
  getFile.onchange = function()
  {
    fileName = getFile.value.replace(/^.*[\\/]/, "");
    if (lang)
    {          // Spanish
      value.value = "Archivo a procesar: " + fileName +
          "\nApriete uno de los botones \"Solo evaluar\", \"Primo\" o \"Factorizar\" para continuar.";
    }
    else
    {          // English
      value.value = "File to process: " + fileName +
          "\nPress \"Only evaluate\", \"Prime\" or \"Factor\" button to continue.";
    }
  };
  btnBlocklyMode.onclick = function()
  {
    hide("main");
    get("blockmode").style.display = "flex";
    history.pushState({id: 3}, "", location.href);
    if (bmodeLoaded === 0)
    {
      initBlockly();
      bmodeLoaded = 1;
    }
  };
  get("exitBlockly").onclick = function()
  {
    history.back();
  };
  btnOpenWizard.onclick = function()
  {
    hide("main");
    show("wizard");
    show("mode");
    get("oneexpr").checked = true;
    history.pushState({id: 2}, "", location.href);
    btnNext.disabled = true;
    wzdInput.value = "";
    wzdInput.focus();
    chkHexW.checked = (config.substring(4,5) === "1");
    chkDecW.checked = (config.substring(4,5) !== "1");
    oneexpr();
  };
  wzdInput.onkeydown = keyDownOnWizard;
  get("btnSentOK").onclick = function()
  {
    history.back();    // End feedback.
  }
  get("btnNotSent").onclick = function()
  {
    history.back();    // End feedback.
  }
  get("oneexpr").onclick = function()
  {
    oneexpr();
  };
  get("loop").onclick = function()
  {
    selectLoop();
  };
  btnNext.onclick = wizardNext;
  wzdInput.oninput = typedOnWizard;
  get("cancel").onclick = function()
  {
    history.back();   // Close configuration mode.
  }
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
    updateVerbose(chkVerbose.checked);
    history.back();   // Close configuration mode.
  };
  get("close-more").onclick = function()
  {
    history.back();   // Close "more" mode.
  };
  get("ncurve").onclick = function()
  {
    history.back();   // Close "more" mode.
    restartFactorization(-2);
  };
  get("nfactor").onclick = function()
  {
    history.back();   // Close "more" mode.
    restartFactorization(-4);
  };
  newCurveOrFactor.oninput = function(event)
  {     // Curve number cannot be greater than 100 000 000.
    get("ncurve").disabled = (newCurveOrFactor.value.length >= 9);
  };
  newCurveOrFactor.onkeydown = function(event)
  {
    let key = event.key;
    let acceptedKeys = ",Backspace,Tab,Right,ArrowRight,Left,ArrowLeft,Cut," +
                       "Control,Meta,Shift,Insert,Delete,Copy,Paste,Home,End," +
                       "0,1,2,3,4,5,6,7,8,9,";
    if (event.ctrlKey || event.metaKey)
    {
      if (key === "c")
      {    // User pressed CTRL-C. Map it to "Copy".
        key = "Copy";
      }
      if (key === "v")
      {    // User pressed CTRL-V. Map it to "Paste".
        key = "Paste";
      }
      if (key === "x")
      {    // User pressed CTRL-X. Map it to "Cut".
        key = "Cut";
      }
    }
    else if (key === "Esc" || key === "Escape")
    {
      history.back();                  // Close "More" mode.
      return;
    }
    if (acceptedKeys.indexOf(","+key+",") < 0)
    {                                  // Key not accepted.
      event.preventDefault();          // Do not propagate this key.  
    }
  };
  btnStop.onclick = function()
  {
    endWorker();
    styleButtons("inline", "none");      // Enable eval and factor
    hide("sktest");    // Hide button if it is present during factorization.
    resultDirty = true;
    resultText += (lang ? "<p>Cálculo detenido por el usuario.</p>" :
                          "<p>Calculation stopped by user</p>");
    statusDirty = true;
    statusText = "";
  };
  value.onkeydown = function (event)
  {
    let keyCode = event.key;
    if (keyCode === "Enter" && event.ctrlKey)
    {
      event.preventDefault();          // Do not propagate Enter key.
      setStorage("ecmFactors","");     // Perform factorization.
      dowork(2);
    }
    return true;
  };
  get("share").onclick = function()
  {
    let shareData =
    {
      title: "Integer Factorization Calculator",
      text: "",
      url: ""
    };
    let tmpHTML = divResult.innerHTML;
    // Convert <sup> and </sup> to exponentiation character.
    tmpHTML = tmpHTML.replace(/<sup>/g, "^");
    tmpHTML = tmpHTML.replace(/<\/sup>/g, "");
    tmpHTML = tmpHTML.replace(/<p>/g, "");
    tmpHTML = tmpHTML.replace(/<\/p>/g, "\n");
    tmpHTML = tmpHTML.replace(/<li>/g, "");
    tmpHTML = tmpHTML.replace(/<\/li>/g, "\n");
    tmpHTML = tmpHTML.replace(/Show divisors/g, "");
    tmpHTML = tmpHTML.replace(/Mostrar divisores/g, "");
    
    // Create a new div element
    let tempDivElement = document.createElement("div");

    // Set the HTML content with the given value
    tempDivElement.innerHTML = tmpHTML;

    // Retrieve the text property of the element 
    shareData.text = tempDivElement.textContent || tempDivElement.innerText || "";
    shareData.url = (window.location.href.split("?")[0]) + "?q=" +
                    encodeURIComponent(value.value);
    navigator.share(shareData).then(function()
    {
    },
    function()
    {
    });
  };
  get("helpbtn").onclick = function()
  {
    let help = get("help");
    let helpStyle = help.style;
    let helphelpStyle = get("helphelp").style;
    let resultStyle = divResult.style;
    if (helpStyle.display === "block" && divResult.innerHTML !== "")     
    {
      helpStyle.display = "none";
      helphelpStyle.display = resultStyle.display = "block";
    }
    else
    {
      hide("sharediv");
      helpStyle.display = "block";
      helphelpStyle.display = resultStyle.display = "none";
    }
  };
  get("skptest").onclick = function()
  {
    hide("sktest");
    restartFactorization(4);
  };
  get("continue").onclick = function()
  {
    hide("cont");
    callWorker("C");  // Indicate worker that user pressed Continue button.
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
  get("wzdfunccat").onchange = function()
  {
    generateFuncButtons("wzdfunccat", "wzdfuncbtns");
  };
  get("value").onfocus = function()
  {
    currentInputBox = get("value");
  };
  get("wzdinput").onfocus = function()
  {
    currentInputBox = get("wzdinput");
  };
  get("comments").oninput = function(_event)
  {
    get("formsend").disabled = (get("comments").value === "");
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
  window.onresize = function(_event)
  {
    let options = {
            "behavior": "auto",
            "block": "center",
            "inline": "center"
        };
    if (document.activeElement === value)
    {  // Center input.
      value.scrollIntoView(options);
    }
    if (document.activeElement === wzdInput)
    {  // Center input.
      wzdInput.scrollIntoView(options);
    }
  };
  setInterval(function()
  {
    if (resultDirty)
    {
      divResult.innerHTML = resultText;
      resultDirty = false;
    }
    if (statusDirty)
    {
      get("status").innerHTML = statusText;
      statusDirty = false;
    }
    if (divisorsDirty)
    {
      let sumSqButton = get("showSumSq");
      if (sumSqButton != null)
      {
        sumSqButton.onclick = showSumSquares;
      }
      if (get("showdiv") != null)
      {
        get("showdiv").onclick = function()
        {
          callWorker("D");  // Indicate worker that user pressed Divisors button.
        };
      }
      divisorsDirty = false;
    }
  }, 100);
  // Generate accordion.
  let acc = document.querySelectorAll("h2");
  let idx, x, y;

  for (idx = 0; idx < acc.length; idx++)
  {
    acc[idx >> 0].addEventListener("click", function()
    {
    // "active" means that panel is being displayed.
      this.children[0].classList.toggle("active");
      let panel = this.nextElementSibling;
      if (panel.style.display === "block")
      {
        panel.style.display = "none";
      }
      else
      {
        panel.style.display = "block";
      }
    });
  }
  get("exprcopy").innerHTML = get("exprorig").innerHTML;
  let c = get("ellCurve");
  let ctx = c.getContext("2d");
  ctx.fillStyle="#FFFFFF";      // White.
  ctx.fillRect(0,0,313,313);    // Clear canvas.
  ctx.fillStyle="#000000";      // Black.
  ctx.strokeStyle="#808000";
  ctx.fillRect(20,0,290,290);   // Clear canvas.
  for (x=0; x<=290; x+=10)      // Draw grid.
  {
    ctx.moveTo(20+x, 0);
    ctx.lineTo(20+x, 290);
    ctx.stroke();
    ctx.moveTo(20, x);
    ctx.lineTo(310, x);
    ctx.stroke();     
  }
  ctx.fillStyle="#00C000";      // Green.
  let ctr;
  for (ctr=0; ctr<points.length; ctr+=2)
  {
    x = points[ctr >> 0];
    y = points[(ctr+1) >> 0];
    ctx.fillRect(20+x*10+1,(28-y)*10+1,9,9);
    if (y !== 0)
    {
      ctx.fillRect(20+x*10+1,(y-1)*10+1,9,9); 
    }
  }
  ctx.fillStyle="#000000";      // Black.
  ctx.font = "15px 'Times New Roman'";
  ctx.fillText("0",5,291);
  ctx.fillText("0",20,308);
  ctx.fillText("28",0,11);
  ctx.fillText("28",297,308);
  ctx.font = "italic "+ctx.font;
  ctx.fillText("y",5,150);
  ctx.fillText("x",160,308);  
  digits = getStorage("ecmConfig");
  if (digits === null || digits === "")
  {
    digits = 6;
    config = "00100";
    setStorage("ecmConfig", digits+","+config);
  }
  else
  {
    index = digits.indexOf(",");
    if (index<0)
    {
      digits = 6;
      config = "00100";
      setStorage("ecmConfig", digits+","+config);
    }
    else
    {
      config = digits.substring(index+1);
      while (config.length < 5)
      {  // Convert legacy configuration.
        config += "0";
      }
      digits = digits.substring(0,index);
      updateVerbose(config.substring(1,2) === "1");
    }
  }
  comingFromPolfact(value);
  registerServiceWorker();
  currentInputBox = get("value");
  generateFuncButtons("funccat", "funcbtns");
  generateFuncButtons("wzdfunccat", "wzdfuncbtns");
  completeFuncButtons("funcbtns");
}
getCalculatorCode("ecmW0000.js", workerParam);
window.addEventListener("load", startUp);
window.addEventListener("popstate", popstate);
window["fromWorker"] = fromWorker;
