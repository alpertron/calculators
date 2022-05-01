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
/* global useBlockly */
/** @define {number} */ var lang = 1;   // Use with Closure compiler.
(function()
{   // This method separates the name space from the Google Analytics code.
var points=[0,6, 2,9, 4,0, 5,6, 7,1, 8,0, 13,9, 14,9, 15,7, 16,7, 17,0, 18,13, 20,5, 22,10, 23,12, 24,6, 27,7];
var wizardStep = 0;
var wizardTextInput;
var worker = 0;
var fileContents = null;
var app;
var blob;
var digits;
var config;
var fromFile;
var tofile;
var fileName;
var workerParam;
var asmjs = typeof(WebAssembly) === "undefined";
var bmodeLoaded = 0;
var statusText = "";
var resultText = "";
var divisorsDirty = false;
var statusDirty = false;
var resultDirty = false;
var calcURLs = ["ecmW0000.js",
               "ecm.webmanifest", "ecmc.webmanifest", "ecm-icon-1x.png", "ecm-icon-2x.png", "ecm-icon-4x.png", "ecm-icon-180px.png", "ecm-icon-512px.png", "favicon.ico"];
var blocklyLoaded = 0;
var scriptsLoaded = 0;
var script1;
var script2;
var funcnames;
var parens;

// DOM resources
var value;
var btnNext;
var btnEval;
var btnPrime;
var btnFactor;
var btnConfig;
var btnFromFile;
var btnBlocklyMode;
var btnOpenWizard;
var btnMore;
var btnToFile;
var btnStop;
var chkCunningham;
var chkDecW;
var chkHex;
var chkHexW;
var chkPretty;
var chkVerbose;
var divResult;
var getFile;
var newCurveOrFactor;
var wzdDescText;
var wzdExamText;
var wzdInput;

if (lang)
{
  funcnames =
  [
    "Suma,+,Resta,-,Multiplicación,*,División,/,Resto,%,Potencia,^,Resultado anterior,ans,Raíz cuadrada entera,sqrt(,Número aleatorio\n\nPrimer argumento: mínimo valor del número aleatorio\nSegundo argumento: máximo valor del número aleatorio,Random(,Valor absoluto,Abs(,Signo,Sign(",
    "Igual,=,Distinto,!=,Mayor,>,Menor o igual,<=,Menor,<,Mayor o igual,>=",
    "Y lógica, AND ,O lógica, OR ,O exclusiva, XOR ,Negación lógica, NOT ,Desplazamiento a la izquierda\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHL ,Desplazamiento a la derecha\n\nOperando izquierdo: valor a desplazar\nOperando derecho: cantidad de bits, SHR ",
    "Máximo común divisor\n\nSe pueden usar uno o más argumentos,GCD(,Mínimo común múltiplo\n\nSe pueden usar uno o más argumentos,LCM(,¿El valor es primo?,IsPrime(,Cantidad de factores primos,NumFact(,menor divisor primo,MinFact(,mayor divisor primo,MaxFact(,Cantidad de divisores,NumDivs(,Suma de divisores,SumDivs(",
    "Primo siguiente,N(,Primo anterior,B(,Cantidad de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,NumDigits(,Suma de dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,SumDigits(,Invertir dígitos\n\nPrimer argumento: valor\nSegundo argumento: base,RevDigits(,Concatenar factores primos\n\nPrimer argumento: modo\n0: Primos no repetidos en forma ascendente\n1: Primos no repetidos en forma descendente\n2: Primos repetidos en forma ascendente\n3: Primos repetidos en forma descendente\nSegundo argumento: valor a factorizar,ConcatFact(",
    "Inverso modular\n\nPrimer argumento: valor\nSegundo argumento: módulo,ModInv(,Exponenciación modular\n\nPrimer argumento: base\nSegundo argumento: exponente\nTercer argumento: módulo,ModPow(,Indicador de Euler,Totient(,Símbolo de Jacobi\n\nPrimer argumento: valor superior\nSegundo argumento: valor inferior,Jacobi(",
    "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partición,P("
  ];
  parens = "Paréntesis izquierdo,(,Paréntesis derecho,),";
}
else
{
  funcnames =
  [
    "Sum,+,Subtraction,-,Multiplication,*,Division,/,Remainder,%,Power,^,Last answer,ans,Integer square root,sqrt(,Random number\n\nFirst argument: minimum value for random number\nSecond argument: maximum value for random number,Random(,Absolute value,Abs(,Sign,Sign(",
    "Equal,=,Not equal,!=,Greater,>,Not greater,<=,Less,<,Not less,>=",
    "Logic AND, AND ,Logic OR, OR ,Exclusive OR, XOR ,Logic NOT, NOT ,Shift left\n\nLeft operand: value to shift\nRight operand: number of bits, SHL ,Shift right\n\nLeft operand: value to shift\nRight operand: number of bits, SHR ",
    "Greatest Common Divisor\n\nOne or more arguments can be used,GCD(,Least Common Multiple\n\nOne or more arguments can be used,LCM(,The value is prime?,IsPrime(,Number of prime factors,NumFact(,smallest prime divisor,MinFact(,greatest prime divisor,MaxFact(,Number of divisors,NumDivs(,Sum of divisors,SumDivs(",
    "Next prime after,N(,Last prime before,B(,Number of digits\n\nFirst argument: value\nSecond argument: base,NumDigits(,Sum of digits\n\nFirst argument: value\nSecond argument: base,SumDigits(,Reverse digits\n\nFirst argument: value\nSecond argument: base,RevDigits(,Concatenate prime factors\n\nFirst argument: Mode\n0: No repeated primes in ascending order\n1: No repeated primes in descending order\n2: Repeated primes in ascending order\n3: Repeated primes in descending order\nSecond argument: Value to factor,ConcatFact(",
    "Modular inverse\n\nFirst argument: value\nSecond argument: modulus,ModInv(,Modular power\n\nFirst argument: base\nSecond argument: exponent\nThird argument: modulus,ModPow(,Totient,Totient(,Jacobi symbol\n\nFirst argument: upper value\nSecond argument: lower value,Jacobi(",
    "Factorial,!,Primorial,#,Fibonacci,F(,Lucas,L(,Partition,P("
  ];
  parens = "Left parenthesis,(,Right parenthesis,),";
}
function get(id)
{
  return document.getElementById(id);
}

function hide(id)
{
  get(id).style.display = "none";
}

function show(id)
{
  get(id).style.display = "block";
}

function oneexpr()
{
  btnNext.value = (lang? "Hecho": "Done");
  wzdDescText.innerHTML = (lang? "Paso 1 de 1: Expresión a factorizar": "Step 1 of 1: Expression to factor");
  wzdExamText.innerHTML = "&nbsp;";
  wizardTextInput = "";
  wizardStep = 9;
}

function setStorage(name, data)
{
  window.localStorage.setItem(name, data);
}

function getStorage(name)
{
  return window.localStorage.getItem(name);
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

function saveConfig()
{    
  config = "1" +   // Batch mode
           (chkVerbose.checked? "1" : "0") +
           (chkPretty.checked? "1" : "0") +
           (chkCunningham.checked? "1" : "0") +
           (chkHex.checked? "1" : "0");
  digits = get("digits").value;
  setStorage("ecmConfig", digits+","+config);
}

function b64decode(str,out)
{
  var ch;
  var idxDest,idxSrc;
  var blocks, leftOver;
  var byte0, byte1, byte2, byte3;
  var conv = new Int8Array(128);
  var len = str.length;
  if (str.charAt(len-1) === "=")
  {
    len--;
  }
  if (str.charAt(len-1) === "=")
  {
    len--;
  }
  blocks=len & (-4);
  for (ch = 65; ch <= 90; ch++)   // A - Z
  {
    conv[ch >> 0] = ch - 65;
  }
  for (ch = 97; ch <= 122; ch++)  // a - z
  {
    conv[ch >> 0] = ch - 71;
  }
  for (ch = 48; ch <= 57; ch++)   // 0 - 9
  {
    conv[ch >> 0] = ch + 4;
  }
  conv[43] = 62;                  // +
  conv[33] = 63;                  // !
  for (idxDest=0,idxSrc=0; idxSrc<blocks; idxDest+=3,idxSrc+=4)
  {
    byte0 = conv[str.charCodeAt(idxSrc)];
    byte1 = conv[str.charCodeAt(idxSrc+1)];
    byte2 = conv[str.charCodeAt(idxSrc+2)];
    byte3 = conv[str.charCodeAt(idxSrc+3)];
    
    out[idxDest >>0 ] = (byte0<<2) + (byte1>>4);
    out[(idxDest+1) >> 0] = (byte1<<4) + (byte2>>2);
    out[(idxDest+2) >> 0] = (byte2<<6) + byte3;
  }
  leftOver = len & 3;
  if (leftOver === 2)
  {
    byte0 = conv[str.charCodeAt(idxSrc)];
    byte1 = conv[str.charCodeAt(idxSrc+1)];
    
    out[idxDest >> 0] = (byte0<<2) + (byte1>>4);
    out[(idxDest+1) >> 0] = byte1<<4;
  }
  else if (leftOver === 3)
  {
    byte0 = conv[str.charCodeAt(idxSrc)];
    byte1 = conv[str.charCodeAt(idxSrc+1)];
    byte2 = conv[str.charCodeAt(idxSrc+2)];
    
    out[idxDest >> 0] = (byte0<<2) + (byte1>>4);
    out[(idxDest+1) >> 0] = (byte1<<4) + (byte2>>2);
    out[(idxDest+2) >> 0] = byte2<<6;
  }
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
    { // First character of e.data is:
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
      var firstChar = e.data.substring(0, 1);
      if (firstChar === "9")
      {
        console.log(e.data.substring(1));
      }
      else if (firstChar === "8")
      {
        setStorage("ecmFactors", e.data.substring(1));
        setStorage("ecmCurve", "");
      }
      else if (firstChar === "7")
      {
        setStorage("ecmCurve", e.data.substring(1));
      }
      else if (firstChar === "D")
      {
        get("divisors").innerHTML = e.data.substring(1);
      }
      else if (firstChar === "E")
      {
        get("divisors").innerHTML = e.data.substring(1);
        get("showdiv").onclick = function()
        {
          callWorker("D");  // Indicate worker that user pressed Divisors button.
        };
      }
      else if (firstChar === "K")
      {
        get("berror").innerHTML = e.data.substring(1);
        show("BlocklyErrors");
        hide("BlocklyButtons");
      }
      else if (firstChar === "L")
      {
        show("main");
        hide("blockmode");
      }
      else if ((firstChar === "M") || (firstChar === "N"))
      {    // User entered a polynomial. Load calculator to process it.
        window.sessionStorage.setItem((firstChar === "M"? "F": "E"),
          value.value);
        window.location.replace(lang? "FACTPOL.HTM": "POLFACT.HTM");
      }
      else if (firstChar === "4")
      {
        statusDirty = true;
        statusText = e.data.substring(1);
      }
      else if (firstChar === "5")
      {
        if (e.data.substring(1, 2) === "1")
        {
          show("skip");
        }
        else
        {
          hide("skip");
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
          hide("modal-more");
          if (firstChar === "A" || firstChar === "B")
          {
            tofile = e.data.substring(1);
            show("savefile");
            resultText = "";
          }
          else
          {
            resultText = e.data.substring(1);
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
          resultText = e.data.substring(1);
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

function performWork(n, valueText)
{
  var param;
  app = lang + n;
  var res = divResult;
  var charNull = String.fromCharCode(0);
  var helphelp = get("helphelp");
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
  helphelp.style.display = "block";
  helphelp.innerHTML = (lang ? "<p class=\"pad\">Aprieta el botón <strong>Ayuda</strong> para obtener ayuda para esta aplicación. Apriétalo de nuevo para retornar a la factorización. También puedes ver <a href=\"/videos/videosEcmc.htm\" target=\"_blank\">videos</a>. Los usuarios con teclado pueden presionar CTRL+ENTER para comenzar la factorización. Esta es la versión "+(asmjs? "asm.js": "WebAssembly")+".</p>":
                               "<p class=\"pad\">Press the <strong>Help</strong> button to get help about this application. Press it again to return to the factorization. You can also watch <a href=\"/videos/videosEcm.htm\" target=\"_blank\">videos</a>. Keyboard users can press CTRL+ENTER to start factorization. This is the "+(asmjs? "asm.js": "WebAssembly")+" version.</p>");
  res.style.display = "block";
  if (typeof(Worker) === "undefined")
  {    // Web workers not supported on this browser.
    resultDirty = true;
    resultText = (lang ? "<p>Esta calculadora necesita Web Workers. Por favor use otro navegador Web.</p>" :
                         "<p>This calculator requires Web Workers. Please use another Web browser.</p>");
    return;
  }
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
  var valueText;
  fromFile = "0";
  if (getFile.value !== "")
  {
    var fileReader = new FileReader();
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
  hide("modal-more");
  if (worker)
  {
    worker.terminate();
  }
  worker = 0;
  dowork(type);
}

function selectLoop()
{   
  btnNext.value = (lang ? "Siguiente": "Next");
  wzdDescText.innerHTML = (lang ? "Paso 1 de 5: Valor inicial de x": "Step 1 of 5: Initial value of x");
  wzdExamText.innerHTML = (lang? "No usar variables <var>x</var> o <var>c</var>. Ejemplo para números de Smith menores que 10000: <code>1</code>": 
                                    "Do not use variables <var>x</var> or <var>c</var>. Example for Smith numbers less than 10000: <code>1</code>");
  wizardStep = 1;
}
  
function wizardNext()
{
  var valueInput = value;
  btnNext.disabled = true;
  switch (++wizardStep)
  {
    case 2:
      wizardTextInput += "x="+wzdInput.value;
      hide("mode");
      wzdDescText.innerHTML = (lang? "Paso 2 de 5: Valor de x para la nueva iteración": "Step 2 of 5: Value of x for new iteration");
      wzdExamText.innerHTML = (lang? "Variables <var>x</var> y/o <var>c</var> requeridas. Ejemplo para números de Smith menores que 10000: <code>x+1</code>":
                                     "Variables <var>x</var> and/or <var>c</var> required. Example for Smith numbers less than 10000: <code>x+1</code>");
      break;
    case 3:
      wizardTextInput += ";x="+wzdInput.value;
      wzdDescText.innerHTML = (lang? "Paso 3 de 5: Condición para finalizar el ciclo": "Step 3 of 5: End loop condition");
      wzdExamText.innerHTML = (lang? "Variables <var>x</var> y/o <var>c</var> requeridas. Ejemplo para números de Smith menores que 10000: <code>x&lt;10000</code>":
                                     "Variables <var>x</var> and/or <var>c</var> required. Example for Smith numbers less than 10000: <code>x&lt;10000</code>");
      break;
    case 4:
      wizardTextInput += ";"+wzdInput.value;
      wzdDescText.innerHTML = (lang? "Paso 4 de 5: Expresión a factorizar": "Step 4 of 5: Expression to factor");
      wzdExamText.innerHTML = (lang? "Variables <var>x</var> y/o <var>c</var> requeridas. Ejemplo para números de Smith menores que 10000: <code>x</code>":
                                     "Variables <var>x</var> and/or <var>c</var> required. Example for Smith numbers less than 10000: <code>x</code>");
      break;
    case 5:
      wizardTextInput += ";"+wzdInput.value;
      btnNext.value = (lang? "Hecho": "Done");
      btnNext.disabled = false;
      wzdDescText.innerHTML = (lang? "Paso 5 de 5: Condición para procesar la expresión": "Step 5 of 5: Process expression condition");
      wzdExamText.innerHTML = (lang? "Variables <var>x</var> y/o <var>c</var> requeridas. Ejemplo para números de Smith menores que 10000: <code>sumdigits(x,10) == sumdigits(concatfact(2,x),10) and not isprime(x)</code>":
                                     "Variables <var>x</var> and/or <var>c</var> required. Example for Smith numbers less than 10000: <code>sumdigits(x,10) == sumdigits(concatfact(2,x),10) and not isprime(x)</code>");
      break;
    case 6:
      if (wzdInput.value !== "")
      {
        wizardTextInput += ";"+wzdInput.value;
      }
      valueInput.value = wizardTextInput;
      wizardStep = 0;
      chkHex.checked = chkHexW.checked;
      saveConfig();
      show("main");
      hide("wizard");
      valueInput.focus();
      break;
    default:
      wizardStep = 0;
      valueInput.value = wzdInput.value;
      chkHex.checked = chkHexW.checked;
      saveConfig();
      show("main");
      hide("wizard");
      valueInput.focus();
      break;
  } 
  if (wizardStep)
  {
    wzdInput.value = "";
    wzdInput.focus();
  }
}

function updateVerbose(isVerbose)    
{
  var cssRules = (document.all) ? document.styleSheets[0]["rules"]: document.styleSheets[0]["cssRules"];
  var index;
  var len = cssRules.length;
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

function endFeedback()
{
  show("main");
  hide("feedback");
  value.focus();   
}

var url = window.location.pathname;
function updateCache(cache)
{
  caches.open("cacheECM").then(function(tempCache)
  {
    tempCache.addAll([url].concat(calcURLs)).then(function()
    {     // Copy cached resources to main cache and delete this one.
      tempCache.matchAll().then(function(responseArr)
      {   // All responses in array responseArr.
        responseArr.forEach(function(responseTempCache, _index, _array)
        {
          cache.put(responseTempCache.url, responseTempCache);
        });
      })
      .finally(function()
      {
        caches.delete("cacheECM");
      });
    });  
  });
}

function fillCache()
{
  // Test whether the HTML is already on the cache.
  caches.open("newCache").then(function(cache)
  {
    cache.match(url).then(function (response)
    {
      if (typeof response === "undefined")
      {     // HTML is not in cache.
        updateCache(cache);
      }
      else
      {     // Response is the HTML contents.
        var date = response.headers.get("last-modified");
            // Request the HTML from the Web server.
            // Use non-standard header to tell Service Worker not to retrieve HTML from cache.
        fetch(url,{headers:{"If-Modified-Since": date, "x-calc": "1"}, cache: "no-store"}).then(function(responseHTML)
        {
          if (responseHTML.status !== 200)
          {
            return;        // HTML could not be retrieved, so go out.
          }
          if (date === responseHTML.headers.get("last-modified"))
          {
            return;        // HTML has not changed, so other files have not been changed. Go out.
          }
          // Read files to new cache.
          // Use temporary cache so if there is any network error, original cache is not changed.
        
          caches.open("cacheECM").then(function(tempCache)
          {                // Do not fetch HTML because it is already fetched.
            tempCache.addAll(calcURLs).then(function()
            {              // Copy cached resources to main cache and delete this one.
              tempCache.matchAll().then(function(responseArr)
              {            // All responses in array responseArr.
                responseArr.forEach(function(responseTempCache, _index, _array)
                {
                  var urlTemp = responseTempCache.url;
                  var indexZero = url.indexOf("00");
                  if (indexZero > 0)
                  {        // There is an old version of this resource on cache to be erased.
                    cache.keys().then(function(keys)
                    {
                      keys.forEach(function(requestCache, _idx, _arr)
                      {    // Traverse cache.
                        if (requestCache.url.substring(0, indexZero+2) === urlTemp.substring(0, indexZero+2) &&
                            requestCache.url.substring(indexZero+2, indexZero+4) !== urlTemp.substring(indexZero+2, indexZero+4) &&
                            requestCache.url.substring(indexZero+4) === urlTemp.substring(indexZero+4))
                        {  // Old version of asset found (different number and same prefix and suffix). Delete it from cache.
                          cache.delete(requestCache);
                        }  
                      });
                      // Put resource into cache after old resource has been erased.
                      cache.put(urlTemp, responseTempCache);
                    });
                  }
                  else
                  {   // Put resource into cache (no old resource into cache). 
                    cache.put(urlTemp, responseTempCache);
                  }
                });
                cache.put(url, responseHTML);
              });
            })
            .finally(function()
            {
              caches.delete("cacheECM");
            });
          })
          .catch (function()     // Cannot fetch HTML.
          {
            updateCache(cache);
          });
        });
      }
    });
  });
}

function fromBlocklyRun(xml)
{
  fromFile = "0";
  performWork(8, xml);
}

function loadScript(scriptUrl)
{
  var myScript = document.createElement("script");
  var xmlhttp = new XMLHttpRequest();
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
  script1 = loadScript("blockly0005.js");
  script2 = loadScript(lang? "es0005.js": "en0005.js");
}

function generateFuncButtons(optionCategory, funcButtons, inputId)
{
  var button;
  var catIndex;
  var funcbtns = get(funcButtons);
  var catnbr = get(optionCategory).selectedIndex;
  var funcname = (parens + funcnames[catnbr]).split(",");
  // Append all buttons to document fragment instead of funcbtns
  // and finally append the fragment to funcbtns to minimize redraws.
  var fragment = document.createDocumentFragment();
  for (catIndex = 0; catIndex < funcname.length/2; catIndex++)
  {
    button = document.createElement("button");
    button.setAttribute("type", "button");        // Indicate this is a button, not submit.
    button.setAttribute("title", funcname[catIndex*2]);  // Text of tooltip.
    button.innerHTML = funcname[catIndex*2 + 1];         // Text of button.
    button.onclick = function()
    {
      var input = get(inputId);
      input.focus();
      var start = input.selectionStart;
      input.value = input.value.substring(0, start) +
                    this.innerText +
                    input.value.substring(input.selectionEnd, input.value.length);
        // Place the caret at the end of the appended text.
      input.selectionStart = start + this.innerText.length;
      input.selectionEnd = input.selectionStart;
    };
    fragment.appendChild(button);
  }
  funcbtns.innerHTML = "";
  funcbtns.appendChild(fragment);
}

function completeFuncButtons(funcButtons, inputId)
{
  var button;
  var catIndex;
  var funcname = (parens + funcnames[0]).split(",");
  var funcbtns = get(funcButtons);
  for (catIndex = 0; catIndex < funcname.length/2; catIndex++)
  {
    button = funcbtns.children[catIndex];
    button.setAttribute("title", funcname[catIndex*2]);  // Text of tooltip.
    button.onclick = function()
    {
      var input = get(inputId);
      input.focus();
      var start = input.selectionStart;
      input.value = input.value.substring(0, start) +
                    this.innerText +
                    input.value.substring(input.selectionEnd, input.value.length);
        // Place the caret at the end of the appended text.
      input.selectionStart = start + this.innerText.length;
      input.selectionEnd = input.selectionStart;
    };
  } 
}

function startUp()
{
  var index, ecmFactor;
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
    show("modal-more");
  };
  btnConfig.onclick = function()
  {
    get("digits").value = digits;
    chkVerbose.checked = (config.substr(1,1) === "1");
    chkPretty.checked = (config.substr(2,1) === "1");
    chkCunningham.checked = (config.substr(3,1) === "1");  
    chkHex.checked = (config.substr(4,1) === "1");
    show("modal-config");
  };
  btnFromFile.onclick = function()
  {
    getFile.click();
  };
  btnToFile.onclick = function()
  {
    hide("savefile");
    var fileBlob = new Blob([tofile], { type: "text/plain" });
    var fileUrl = URL.createObjectURL(fileBlob);
    var a = document.createElement("a");
    a.href = fileUrl;
    a.download = fileName;
    var clickHandler = function()
    {
      setTimeout(function()
      {
        URL.revokeObjectURL(fileUrl);
        this.removeEventListener("click", clickHandler);
      },
      150);
    };
    a.addEventListener("click", clickHandler, false);
    a.click();
  };
  getFile.onchange = function()
  {
    fileName = getFile.value.replace(/^.*[\\\/]/, "");
    if (lang)
    {          // Spanish
      value.value = "Archivo a procesar: " + fileName +
          "\nApriete el botón \"Solo evaluar\", \"Primo\" or \"Factorizar\" para continuar.";
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
    if (bmodeLoaded === 0)
    {
      initBlockly();
      bmodeLoaded = 1;
    }
  };
  get("exitBlockly").onclick = function()
  {
    show("main");
    hide("blockmode");
  };
  btnOpenWizard.onclick = function()
  {
    hide("main");
    show("wizard");
    show("mode");
    get("oneexpr").checked = true;
    btnNext.disabled = true;
    wzdInput.value = "";
    wzdInput.focus();
    chkHexW.checked = (config.substr(4,1) === "1");
    chkDecW.checked = (config.substr(4,1) !== "1");
    oneexpr();
  };
  wzdInput.onkeydown = function (event)
  {
    if (event.keyCode === 10 || event.keyCode === 13)
    {
      event.preventDefault();          // Do not propagate Enter key.
      if (!btnNext.disabled)
      {                                // Next button is not disabled.
        wizardNext();                  // Perform same operation as if the user had pressed Next button.
      }
    }
    if (event.altKey)
    {                                  // User pressed ALT key.
      if (event.keyCode === 80)
      {                                // User pressed ALT-P.
        event.preventDefault();        // Do not propagate key.
        if (get("oneexpr").checked)
        {
          get("oneexpr").checked = false;
          get("loop").checked = true;
          selectLoop();
        }
        else
        {
          get("oneexpr").checked = true;
          get("loop").checked = false;
          oneexpr();
        }
      }
      else if (event.keyCode === 68)
      {                                // User pressed ALT-D.
        event.preventDefault();        // Do not propagate key.
        chkDecW.checked = true;
        chkHexW.checked = false;
      }
      else if (event.keyCode === 72)
      {                                // User pressed ALT-H.
        event.preventDefault();        // Do not propagate key.
        chkDecW.checked = false;
        chkHexW.checked = true;
      }
    }
    return true;
  };
  get("oneexpr").onclick = function()
  {
    oneexpr();
  };
  get("loop").onclick = function()
  {
    selectLoop();
  };
  btnNext.onclick = function()
  {
    wizardNext();
  };
  wzdInput.oninput = function()
  {
    var inputValue = wzdInput.value;
    if (inputValue !== "")
    {         // User typed something on input box.
      if (wizardStep === 1 || wizardStep === 9 ||
          (inputValue.lastIndexOf("x") >= 0 || inputValue.lastIndexOf("c") >= 0 ||
          inputValue.lastIndexOf("X") >= 0 || inputValue.lastIndexOf("C") >= 0))
      {       // At least one x or c. Indicate valid.
        btnNext.disabled = false;
      }
      else
      {
        btnNext.disabled = true;
      }
    }
    else if (wizardStep === 5)
    {         // Last step is optional, so empty input is valid.
      btnNext.disabled = false;
    }
    else
    {         // For required input, empty input is invalid.
      btnNext.disabled = true;
    }
  };
  get("cancel").onclick = function()
  {
    show("main");
    hide("wizard");
  };
  get("close-config").onclick = function()
  {
    hide("modal-config");
  };
  get("cancel-config").onclick = function()
  {
    hide("modal-config");
  };
  get("save-config").onclick = function()
  {
    saveConfig();
    updateVerbose(chkVerbose.checked);
    hide("modal-config");
  };
  get("close-more").onclick = function()
  {
    hide("modal-more");
  };
  get("ncurve").onclick = function()
  {
    restartFactorization(-2);
  };
  get("nfactor").onclick = function()
  {
    restartFactorization(-4);
  };
  newCurveOrFactor.onkeypress = function(event)
  {
    return (event.charCode === 8 || event.charCode === 0) ? null : event.charCode >= 48 && event.charCode <= 57;
  };
  btnStop.onclick = function()
  {
    worker.terminate();
    worker = 0;
    styleButtons("inline", "none");      // Enable eval and factor
    hide("skip");    // Hide button if it is present during factorization.
    resultDirty = true;
    resultText += (lang ? "<p>Cálculo detenido por el usuario.</p>" :
                          "<p>Calculation stopped by user</p>");
    statusDirty = true;
    statusText = "";
  };
  value.onkeydown = function (event)
  {
    if ((event.keyCode === 10 || event.keyCode === 13) && event.ctrlKey)
    {
      event.preventDefault();          // Do not propagate Enter key.
      setStorage("ecmFactors","");     // Perform factorization.
      dowork(2);
    }
    return true;
  };
  get("share").onclick = function()
  {
    var shareData =
    {
      title: "Integer Factorization Calculator",
      text: "",
      url: ""
    };
    var tmpHTML = divResult.innerHTML;
    // Convert <sup> and </sup> to exponentiation character.
    tmpHTML = tmpHTML.replace(/\<sup\>/g, "\^");
    tmpHTML = tmpHTML.replace(/\<\/sup\>/g, "");
    tmpHTML = tmpHTML.replace(/\<p\>/g, "");
    tmpHTML = tmpHTML.replace(/\<\/p\>/g, "\n");
    tmpHTML = tmpHTML.replace(/\<li\>/g, "");
    tmpHTML = tmpHTML.replace(/\<\/li\>/g, "\n");
    tmpHTML = tmpHTML.replace(/Show divisors/g, "");
    tmpHTML = tmpHTML.replace(/New!/g, "");
    tmpHTML = tmpHTML.replace(/Mostrar divisores/g, "");
    tmpHTML = tmpHTML.replace(/¡Nuevo!/g, "");
    
    // Create a new div element
    var tempDivElement = document.createElement("div");

    // Set the HTML content with the given value
    tempDivElement.innerHTML = tmpHTML;

    // Retrieve the text property of the element 
    shareData.text = tempDivElement.textContent || tempDivElement.innerText || "";
    shareData.url = (window.location.href.split("?")[0]) + "?q=" +
                    encodeURIComponent(value.value);
    navigator.share(shareData);
  };
  get("helpbtn").onclick = function()
  {
    var help = get("help");
    var helpStyle = help.style;
    var helphelpStyle = get("helphelp").style;
    var resultStyle = divResult.style;
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
  get("skiptest").onclick = function()
  {
    hide("skip");
    restartFactorization(4);
  };
  get("continue").onclick = function()
  {
    hide("cont");
    callWorker("C");  // Indicate worker that user pressed Continue button.
  };
  get("formlink").onclick = function()
  {
    hide("main");
    show("feedback");
    get("formfeedback").reset();
    get("name").focus();
    return false;   // Do not follow the link.
  };
  get("formcancel").onclick = function()
  {
    endFeedback();
  };
  get("funccat").onchange = function()
  {
    generateFuncButtons("funccat", "funcbtns", "value");
  };
  get("wzdfunccat").onchange = function()
  {
    generateFuncButtons("wzdfunccat", "wzdfuncbtns", "wzdinput");
  };
  get("formsend").onclick = function()
  {
    var userdata = get("userdata");
    if (get("adduserdata").checked)
    {
      userdata.value = "\n" + value.value + "\n" + divResult.innerHTML + "\n" + get("status").innerHTML;
    }
    else
    {
      userdata.value = "";      
    }
    var xhr = new XMLHttpRequest();
    xhr.onreadystatechange = function(_event)
    {
      if (xhr.readyState === 4) 
      {             // XHR finished.
        if (xhr.status === 200)
        {           // PHP page loaded.
          alert(lang?"Comentarios enviados satisfactoriamente.": "Feedback sent successfully.");
        }
        else
        {           // PHP page not loaded.
          alert(lang?"No se pudieron enviar los comentarios.": "Feedback could not be sent.");
        }
        endFeedback();
      }
    };
    xhr.open("POST", (lang? "/enviomail.php": "/sendmail.php"), true);
    xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    var elements = get("formfeedback").elements;
    var contents = "";
    var useAmp = 0;
    for (var i = 0; i < elements.length; i++)
    {
      var element = elements[i >> 0];
      if (element.type === "radio" && !element.checked)
      {
        continue;
      }
      if (element.name)
      {
        if (useAmp)
        {
          contents += "&";
        }
        contents += element.name + "=" + encodeURIComponent(element.value);
        useAmp++;
      }
    }
    xhr.send(contents);
    return false;   // Send form only through JavaScript.
  };
  window.onclick = function(event)
  {
    var modal = get("modal");
    if (event.target === modal)
    {
      modal.style.display = "none";
    }
  };
  window.onresize = function(_event)
  {
    var options = {
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
  var acc = document.querySelectorAll("h2");
  var idx, x, y;

  for (idx = 0; idx < acc.length; idx++)
  {
    acc[idx >> 0].addEventListener("click", function()
    {
    // "active" means that panel is being displayed.
      this.children[0].classList.toggle("active");
      var panel = this.nextElementSibling;
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
  var c = get("ellCurve");
  var ctx = c.getContext("2d");
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
  var ctr;
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
      config = digits.substr(index+1);
      while (config.length < 5)
      {  // Convert legacy configuration.
        config += "0";
      }
      digits = digits.substr(0,index);
      updateVerbose(config.substr(1,1) === "1");
    }
  }
  var fromPolfact = window.sessionStorage.getItem("F");
  if (fromPolfact != null)
  {    // Number to factor coming from polynomial factorization calculator.
    window.sessionStorage.removeItem("F");
    value.value = fromPolfact;
    dowork(-2);    // Perform factorization.
  }
  fromPolfact = window.sessionStorage.getItem("E");
  if (fromPolfact != null)
  {    // Number to factor coming from polynomial factorization calculator.
    window.sessionStorage.removeItem("E");
    value.value = fromPolfact;
    dowork(0);     // Perform evaluation.
  }
  else
  {
    var search = window.location.search;
    if (search.substring(0,3) === "?q=")
    {
      value.value = decodeURIComponent(search.substring(3));
      dowork(-2);
    }
    else
    {
      ecmFactor = getStorage("ecmFactors");
      if (ecmFactor)
      {          // Continue factoring.
        value.value = ecmFactor.slice(0,ecmFactor.indexOf("="));
        newCurveOrFactor.value = getStorage("ecmCurve");
        dowork(-2);
        newCurveOrFactor.value = "";
      }
    }
  }
  if ("serviceWorker" in navigator)
  { // Attempt to register service worker.
    // There is no need to do anything on registration success or failure in this JavaScript module.
    navigator["serviceWorker"]["register"]("calcSW.js").then(
              function() {/* Nothing to do */}, function() {/* Nothing to do */});
    fillCache();
  }
}
completeFuncButtons("funcbtns", "value");
generateFuncButtons("wzdfunccat", "wzdfuncbtns", "wzdinput");

if (asmjs)
{
  var req = new XMLHttpRequest();
  req.open("GET", "ecmW0000.js", true);
  req.responseType = "arraybuffer";
  req.onreadystatechange = function (_aEvt)
  {
    if (req.readyState === 4 && req.status === 200)
    {
      fileContents = /** @type {ArrayBuffer} */ (req.response);
      if (workerParam)
      {
        callWorker(workerParam);
      }
    }
  };
  req.send(null);
}
else
{
  var wasm = get("wasmb64").text;
  while (wasm.charCodeAt(0) < 32)
  {
    wasm = wasm.substring(1);
  }    
  while (wasm.charCodeAt(wasm.length-1) < 32)
  {
    wasm = wasm.substring(0, wasm.length-1);
  }    
  var length = wasm.length*3/4;
  if (wasm.charCodeAt(wasm.length-1) === 61)
  {
    length--;
  }
  if (wasm.charCodeAt(wasm.length-2) === 61)
  {
    length--;
  }
  fileContents=new Int8Array(length);
  b64decode(wasm, fileContents);
  calcURLs.shift();  // Do not fetch Javascript file that will not be used.
}

window.addEventListener("load", startUp);
})();
