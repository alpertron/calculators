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
/** @define {number} */ var lang = 1;   // Use with Closure compiler.
(function(global)
{   // This method separates the name space from the Google Analytics code.
var worker = 0;
var app;
var blob;
var workerParam;
var fileContents = 0;
var asmjs = typeof(WebAssembly) === "undefined";
function get(x)
{
  return document.getElementById(x);
}
function callWorker(param)
{
  if (!worker)
  {
    if (!blob)
    {
      if (asmjs)
      {    // Asm.js
        blob = new Blob([fileContents],{type: 'text/javascript'});
      }
      else
      {    // WebAssembly
        blob = new Blob([get("worker").textContent],{type: 'text/javascript'});
      }
    }   
    worker = new Worker(window.URL.createObjectURL(blob));
    worker.onmessage = function(e)
    { // First character of e.data is '1' for intermediate text
      // and it is '2' for end of calculation.
      get('result').innerHTML = e.data.substring(1);
      if (e.data.substring(0, 1) == '2')
      {   // First character passed from web worker is '2'.
        get('solve').disabled = false;
        get('steps').disabled = false;
        get('stop').disabled = true;
      }
    }
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
  var param;
  var app = lang + n;
  var res = get('result');
  var coefAText = get('coefA').value.trim();
  var coefBText = get('coefB').value.trim();
  var coefCText = get('coefC').value.trim();
  var coefDText = get('coefD').value.trim();
  var coefEText = get('coefE').value.trim();
  var coefFText = get('coefF').value.trim();
  var digitGroup = get('digits').value.trim();
  get('help').style.display = "none";
  res.style.display = "block";
  var missing = "";
  var zero = String.fromCharCode(0);
  if (coefAText == "")
  {
    missing = (lang? "coeficiente <var>a</var>." : "coefficient <var>a</var>.");
  }
  if (coefBText == "")
  {
    missing = (lang? "coeficiente <var>b</var>." : "coefficient <var>b</var>.");
  }
  if (coefCText == "")
  {
    missing = (lang? "coeficiente <var>c</var>." : "coefficient <var>c</var>.");
  }
  if (coefDText == "")
  {
    missing = (lang? "coeficiente <var>d</var>." : "coefficient <var>d</var>.");
  }
  if (coefEText == "")
  {
    missing = (lang? "coeficiente <var>e</var>." : "coefficient <var>e</var>.");
  }
  if (coefFText == "")
  {
    missing = (lang? "coeficiente <var>f</var>." : "coefficient <var>f</var>.");
  }
  if (missing != "")
  {
    res.innerHTML = (lang? "Por favor ingrese un número o expresión para el "+missing :
                               "Please type a number or expression for the "+missing);
    return;
  }
  get('solve').disabled = true;
  get('steps').disabled = true;
  get('stop').disabled = false;
  res.innerHTML = (lang? "Resolviendo la ecuación cuadrática..." :
                             "Solving the quadratic equation...");
  param = digitGroup + ',' + app + ',' + coefAText + zero + coefBText + zero + coefCText + zero +
                                         coefDText + zero + coefEText + zero + coefFText + zero;
  callWorker(param);
}

function moveNext(e, curr, next)
{    
  var nextInput = get(next);
  if ((e.which == 10 || e.which == 13) && curr.value.trim().length > 0)
  {
    e.preventDefault();
    nextInput.focus();
    nextInput.setSelectionRange(0, nextInput.value.length);
  }
}

function endFeedback()
{
  get("main").style.display = "block";
  get("feedback").style.display = "none";
  get("coefA").focus();   
}

function b64decode(str,out)
{
  var ch, idx;
  var idxDest,idxSrc;
  var blocks, left_over;
  var byte0, byte1, byte2, byte3;
  var conv=new Int8Array(128);
  var len=str.length;
  if(str.charAt(len-1)=="=")
  {
    len--;
  }
  if(str.charAt(len-1)=="=")
  {
    len--;
  }
  blocks=len & (-4);
  for (ch = 65; ch <= 90; ch++)   // A - Z
  {
    conv[ch] = ch - 65;
  }
  for (ch = 97; ch <= 122; ch++)  // a - z
  {
    conv[ch] = ch - 71;
  }
  for (ch = 48; ch <= 57; ch++)   // 0 - 9
  {
    conv[ch] = ch + 4;
  }
  conv[43] = 62;                  // +
  conv[33] = 63;                  // !
  for (idxDest=0,idxSrc=0; idxSrc<blocks; idxDest+=3,idxSrc+=4)
  {
    byte0 = conv[str.charCodeAt(idxSrc)];
    byte1 = conv[str.charCodeAt(idxSrc+1)];
    byte2 = conv[str.charCodeAt(idxSrc+2)];
    byte3 = conv[str.charCodeAt(idxSrc+3)];
    
    out[idxDest] = (byte0<<2) + (byte1>>4);
    out[idxDest+1] = (byte1<<4) + (byte2>>2);
    out[idxDest+2] = (byte2<<6) + byte3;
  }
  left_over = len & 3;
  if (left_over == 2)
  {
    byte0 = conv[str.charCodeAt(idxSrc)];
    byte1 = conv[str.charCodeAt(idxSrc+1)];
    
    out[idxDest] = (byte0<<2) + (byte1>>4);
    out[idxDest+1] = byte1<<4;
  }
  else if (left_over == 3)
  {
    byte0 = conv[str.charCodeAt(idxSrc)];
    byte1 = conv[str.charCodeAt(idxSrc+1)];
    byte2 = conv[str.charCodeAt(idxSrc+2)];
    
    out[idxDest] = (byte0<<2) + (byte1>>4);
    out[idxDest+1] = (byte1<<4) + (byte2>>2);
    out[idxDest+2] = byte2<<6;
  }
}

var calcURLs = ["quadW0000.js",
               "quad.webmanifest", "cuad.webmanifest", "quad-icon-1x.png", "quad-icon-2x.png", "quad-icon-4x.png", "quad-icon-180px.png", "quad-icon-512px.png", "favicon.ico"];

if (!asmjs)
{
  calcURLs.shift();  // Do not fetch Javascript file that will not be used.
}
var url = window.location.pathname;
function fillCache()
{
  // Test whether the HTML is already on the cache.
  caches.open("newCache").then(function(cache)
  {
    cache.match(url).then(function (response)
    {
      if (response === undefined)
      {     // HTML is not in cache.
        UpdateCache(cache);
      }
      else
      {     // Response is the HTML contents.
        var date = response.headers.get('last-modified');
            // Request the HTML from the Web server.
            // Use non-standard header to tell Service Worker not to retrieve HTML from cache.
        fetch(url,{headers:{'If-Modified-Since': date, 'x-calc': '1'}, cache: "no-store"}).then(function(responseHTML)
        {
          if (responseHTML.status != 200)
          {
            return;        // HTML could not be retrieved, so go out.
          }
          if (date == responseHTML.headers.get('last-modified'))
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
                responseArr.forEach(function(responseTempCache, index, array)
                {
                  var urlTemp = responseTempCache.url;
                  var indexZero = url.indexOf("00");
                  if (indexZero > 0)
                  {        // There is an old version of this resource on cache to be erased.
                    cache.keys().then(function(keys)
                    {
                      keys.forEach(function(requestCache, index, array)
                      {    // Traverse cache.
                        if (requestCache.url.substring(0, indexZero+2) == urlTemp.substring(0, indexZero+2) &&
                            requestCache.url.substring(indexZero+2, indexZero+4) != urlTemp.substring(indexZero+2, indexZero+4) &&
                            requestCache.url.substring(indexZero+4) == urlTemp.substring(indexZero+4))
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
            UpdateCache(cache);
          });
        })
      }
    });
  });
}

function UpdateCache(cache)
{
  caches.open("cacheECM").then(function(tempCache)
  {
    tempCache.addAll([url].concat(calcURLs)).then(function()
    {     // Copy cached resources to main cache and delete this one.
      tempCache.matchAll().then(function(responseArr)
      {   // All responses in array responseArr.
        responseArr.forEach(function(responseTempCache, index, array)
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

window.onload = function ()
{
  var param;
  get('stop').disabled = true;
  get('solve').onclick = function ()
  {
    dowork(0);
  }
  get('steps').onclick = function ()
  {
    dowork(2);
  }
  get('stop').onclick = function ()
  {
    worker.terminate();
    worker = 0;
    get('solve').disabled = false;
    get('steps').disabled = false;
    get('stop').disabled = true;
    get('result').innerHTML = 
      (lang? "<p>Cálculo detenido por el usuario.</p>" :
                 "<p>Calculation stopped by user</p>");
  }
  get('helpbtn').onclick = function ()
  {
    get('help').style.display = "block";
    get('result').style.display = "none";
  }
  get('coefA').onkeypress = function(e)
  {
    moveNext(e, this, 'coefB');
  } 
  get('coefB').onkeypress = function(e)
  {
    moveNext(e, this, 'coefC');
  }
  get('coefC').onkeypress = function(e)
  {
    moveNext(e, this, 'coefD');
  }
  get('coefD').onkeypress = function(e)
  {
    moveNext(e, this, 'coefE');
  }
  get('coefE').onkeypress = function(e)
  {
    moveNext(e, this, 'coefF');
  }
  get('coefF').onkeypress = function(e)
  {
    if ((e.which == 10 || e.which == 13) && this.value.trim().length > 0)
    {
      e.preventDefault();
      get('coefA').focus();
      dowork(0);
    }
  }
  get("formlink").onclick = function ()
  {
    get("main").style.display = "none";
    get("feedback").style.display = "block";
    get("formfeedback").reset();
    get("name").focus();
    return false;   // Do not follow the link.
  }
  get("formcancel").onclick = function ()
  {
    endFeedback();
  }
  get("formsend").onclick = function()
  {
    var userdata = get("userdata");
    if (get("adduserdata").checked)
    {
      userdata.value = "ax^2 + bxy + cy^2 + dx + ey + f = 0" +
                       "\na = " + get('coefA').value.trim() + "\nb = " + get('coefB').value.trim() +
                       "\nc = " + get('coefC').value.trim() + "\nd = " + get('coefD').value.trim() +
                       "\ne = " + get('coefE').value.trim() + "\nf = " + get('coefF').value.trim();  
    }
    else
    {
      userdata.value = "";      
    }
    var xhr = new XMLHttpRequest();
    xhr.onreadystatechange = function (event)
    {
      if (xhr.readyState == 4) 
      {             // XHR finished.
        if (xhr.status == 200)
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
    xhr.setRequestHeader('Content-type', 'application/x-www-form-urlencoded');
    var elements = get("formfeedback").elements;
    var contents = "";
    var useAmp = 0;
    for (var i = 0; i < elements.length; i++)
    {
      var element = elements[i];
      if (element.type == "radio" && element.checked == false)
      {
        continue;
      }
      if (element.name)
      {
        if (useAmp)
        {
          contents += '&';
        }
        contents += element.name + "=" + encodeURIComponent(element.value);
        useAmp++;
      }
    }
    xhr.send(contents);
    return false;   // Send form only through JavaScript.
  }
  if ('serviceWorker' in navigator)
  { // Attempt to register service worker.
    // There is no need to do anything on registration success or failure in this JavaScript module.
    navigator["serviceWorker"].register('calcSW.js').then(function() {}, function() {});
    fillCache();
  }
}
if (asmjs)
{
  var req = new XMLHttpRequest();
  req.open('GET', "quadW0000.js", true);
  req.responseType = "arraybuffer";
  req.onreadystatechange = function (aEvt)
  {
    if (req.readyState == 4 && req.status == 200)
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
  var wasm = document.getElementById("wasmb64").text;
  while (wasm.charCodeAt(0) < 32)
  {
    wasm = wasm.substring(1);
  }    
  while (wasm.charCodeAt(wasm.length-1) < 32)
  {
    wasm = wasm.substring(0, wasm.length-1);
  }    
  var length = wasm.length*3/4;
  if (wasm.charCodeAt(wasm.length-1)==61)
  {
    length--;
  }
  if (wasm.charCodeAt(wasm.length-2)==61)
  {
    length--;
  }
  fileContents=new Int8Array(length);
  b64decode(wasm, fileContents); 
}
})(this);
