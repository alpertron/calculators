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
var digits;
var config;
var workerParam;
var fileContents = 0;
var asmjs = typeof(WebAssembly) === "undefined";
function get(x)
{
  return document.getElementById(x);
}

function setStorage(name, data)
{
  localStorage.setItem(name, data);
}

function getStorage(name)
{
  return localStorage.getItem(name);
}

function styleButtons(style1, style2)
{
  get("eval").style.display = style1;
  get("factor").style.display = style1;
  get("config").style.display = style1;
  get("stop").style.display = style2;
  get("more").style.display = style2;
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
    { // First character of e.data is "1" for intermediate text
      // and it is "2" for end of calculation.
      // It is "9" for saving expression to factor into Web Storage.
      var firstChar = e.data.substring(0, 1);
      if (firstChar == "8")
      {
//        setStorage("ecmFactors", e.data.substring(1));
//        setStorage("ecmCurve", "");
      }
      else if (firstChar == "7")
      {
//        setStorage("ecmCurve", e.data.substring(1));
      }
      else if (firstChar == "4")
      {
        get("status").innerHTML = e.data.substring(1);
      }
      else
      {
        get("result").innerHTML = e.data.substring(1);
        if (e.data.substring(0, 1) == "2")
        {   // First character passed from web worker is "2".
          get("status").innerHTML = "";
          styleButtons("inline", "none");  // Enable eval and factor
          get("modal-more").style.display = "none";
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

function dowork(n)
{
  app = lang + n;
  var res = get("result");
  var valueText = get("value").value;
  var helphelp = get("helphelp");
  get("help").style.display = "none";
  helphelp.style.display = "block";
  helphelp.innerHTML = (lang? "<p>Aprieta el botón <strong>Ayuda</strong> para obtener ayuda para esta aplicación. Apriétalo de nuevo para retornar a la factorización.</p>":
                              "<p>Press the <strong>Help</strong> button to get help about this application. Press it again to return to the factorization.</p>");
  res.style.display = "block";
  if (valueText == "")
  {
    res.innerHTML = (lang? "Por favor ingrese una expresión." :
                           "Please type an expression.");
    return;
  }
  styleButtons("none", "inline");  // Enable "more" and "stop" buttons
  res.innerHTML = (lang? "Factorizando la expresión..." :
                         "Factoring expression...");
  var param = digits + "," + app + "," + valueText + String.fromCharCode(0);
  callWorker(param);
}

function endFeedback()
{
  get("main").style.display = "block";
  get("feedback").style.display = "none";
  get("value").focus();   
}

var calcURLs = ["gaussianW0000.js",
               "gaussian.webmanifest", "gausiano.webmanifest", "gaussian-icon-1x.png", "gaussian-icon-2x.png", "gaussian-icon-4x.png", "gaussian-icon-180px.png", "gaussian-icon-512px.png", "favicon.ico"];

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
        var date = response.headers.get("last-modified");
            // Request the HTML from the Web server.
            // Use non-standard header to tell Service Worker not to retrieve HTML from cache.
        fetch(url,{headers:{"If-Modified-Since": date, "x-calc": "1"}, cache: "no-store"}).then(function(responseHTML)
        {
          if (responseHTML.status !== 200)
          {
            return;        // HTML could not be retrieved, so go out.
          }
          if (date == responseHTML.headers.get("last-modified"))
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

window.onload = function ()
{
  var param;
  get("eval").onclick = function ()
  {
    dowork(0);
  };
  get("factor").onclick = function ()
  {
    dowork(2);
  };
  get("stop").onclick = function ()
  {
    worker.terminate();
    worker = 0;
    styleButtons("inline", "none");  // Enable eval and factor
    get("result").innerHTML =
      (lang? "<p>Cálculo detenido por el usuario.</p>" :
             "<p>Calculation stopped by user</p>");
    get("status").innerHTML = "";
  };
  get("more").onclick = function ()
  {
    get("modal-more").style.display = "block";
  };
  get("config").onclick = function ()
  {
    get("digits").value = digits;
    get("batch").checked = (config.substr(0,1)=="1");
    get("verbose").checked = (config.substr(1,1)=="1");
    get("pretty").checked = (config.substr(2,1)=="1");
    get("cunnin").checked = (config.substr(3,1)=="1");  
    get("modal-config").style.display = "block";
  };
  get("close-config").onclick = function ()
  {
    get("modal-config").style.display = "none";
  };
  get("cancel-config").onclick = function ()
  {
    get("modal-config").style.display = "none";
  };
  get("save-config").onclick = function ()
  {
    var oldconfig = config;
    config = (get("batch").checked? "1" :"0") +
             (get("verbose").checked? "1" :"0") +
             (get("pretty").checked? "1" :"0") +
             (get("cunnin").checked? "1" :"0");
    digits = get("digits").value;
    setStorage("ecmConfig", digits+","+config);
    get("modal-config").style.display = "none";
  };
  get("close-more").onclick = function ()
  {
    get("modal-more").style.display = "none";
  };
  get("helpbtn").onclick = function ()
  {
    var help = get("help");
    var helpStyle = help.style;
    var helphelpStyle = get("helphelp").style;
    var result = get("result");
    var resultStyle = result.style;
    if (helpStyle.display == "block" && result.innerHTML != "")     
    {
      helpStyle.display = "none";
      helphelpStyle.display = resultStyle.display = "block";
    }
    else
    {
      helpStyle.display = "block";
      helphelpStyle.display = resultStyle.display = "none";
    }
  };
  get("value").onkeydown = function(evt)
  {
    if (evt.ctrlKey && evt.keyCode === 13)
    {
      dowork(2);
      evt.stopPropagation();
    }
  };
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
      userdata.value = "\n" + get("value").value + "\n" + get("result").innerHTML + "\n" + get("status").innerHTML;
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
    xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
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
          contents += "&";
        }
        contents += element.name + "=" + encodeURIComponent(element.value);
        useAmp++;
      }
    }
    xhr.send(contents);
    return false;   // Send form only through JavaScript.
  }
  window.onclick = function(event)
  {
    var modal = get("modal");
    if (event.target == modal)
    {
      modal.style.display = "none";
    }
  };
  digits = getStorage("ecmConfig");
  if (digits == null || digits == "")
  {
    digits = 6;
    config = "0010";
    setStorage("ecmConfig", digits+","+config);
  }
  else
  {
    var index = digits.indexOf(",");
    if (index<0)
    {
      digits = 6;
      config = "0010";
      setStorage("ecmConfig", digits+","+config);
    }
    else
    {
      config = digits.substr(index+1);
      digits = digits.substr(0,index);
    }
  }
  if ("serviceWorker" in navigator)
  { // Attempt to register service worker.
    // There is no need to do anything on registration success or failure in this JavaScript module.
    navigator["serviceWorker"].register("calcSW.js").then(function() {}, function() {});
    fillCache();
  }
}
if (asmjs)
{
  var req = new XMLHttpRequest();
  req.open("GET", "gaussianW0000.js", true);
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
