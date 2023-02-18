/*
    This file is part of Alpertron Calculators.

    Copyright 2023 Dario Alejandro Alpern

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
/* global asmjs */
/* global calcURLs */
/* global callWorker */
/* global endFeedback */
/* global fileContents */
/* global get */
/* global getFormSendValue */
/* global lang */

function initMenubarEvents()
{  
  let menuItems = document.querySelectorAll("[role=\"menubar\"] > li");
  Array.prototype.forEach.call(menuItems, function(el, i)
  {
    el.addEventListener("click", function(event)
    {
      if (this.getAttribute("aria-expanded") === "false")
      {
        this.setAttribute("aria-expanded", "true");
      }
      else
      {
        this.setAttribute("aria-expanded", "false");
      }
      this.firstElementChild.firstElementChild.firstElementChild.focus();
      event.preventDefault();
      return false;
    });
    
    el.addEventListener("keydown", function(event)
    {
      let nextNode;
      if (event.key === "Enter")
      {
        this.click(event);
        return;
      }
      if (event.key === "ArrowRight")
      {
        nextNode = this.nextElementSibling;
        if (nextNode === null)
        {
          nextNode = this.parentNode.firstElementChild;
        }
        nextNode.focus();
        event.preventDefault();
        return;
      }
      if (event.key === "ArrowLeft")
      {
        nextNode = this.previousElementSibling;
        if (nextNode === null)
        {
          nextNode = this.parentNode.lastElementChild;
        }
        nextNode.focus();
        event.preventDefault();
        return;
      }
      if (event.key === "ArrowUp")
      {
        this.setAttribute("aria-expanded", "true");
        this.firstElementChild.lastElementChild.firstElementChild.focus();
        event.preventDefault();
        return;
      }
      if (event.key === "ArrowDown")
      { 
        this.setAttribute("aria-expanded", "true");
        this.firstElementChild.firstElementChild.firstElementChild.focus();
        event.preventDefault();
      }
    });
    
    el.addEventListener("mouseover", function(event)
    {
      el.addEventListener("mouseover", function(event)
      {
        this.setAttribute("aria-expanded", "true");
      });
      el.addEventListener("mouseout", function(event)
      {
        this.setAttribute("aria-expanded", "false");
      });
    });
    
    let submenuItems = el.querySelectorAll("a");
    Array.prototype.forEach.call(submenuItems, function(el, i)
    {
      el.tabIndex = -1;
      el.addEventListener("click", function(event)
      {
        let parent = this.parentNode.parentNode.parentNode;
        parent.setAttribute("aria-expanded", "false");
        window.location = this.getAttribute("href");
        event.stopImmediatePropagation();
        event.preventDefault();
      });
      el.addEventListener("keydown", function(event)
      {
        let next;
        let parent = this.parentNode.parentNode.parentNode;
        if (event.key === "Tab")
        {
          parent.setAttribute("aria-expanded", "false");
          return;
        }
        if (event.key === "Escape")
        {
          parent.setAttribute("aria-expanded", "false");
          parent.focus();
          event.preventDefault();          
          return;
        }
        if (event.key === "Enter")
        {
          parent.setAttribute("aria-expanded", "false");
          window.location = this.getAttribute("href");
          event.stopImmediatePropagation();
          event.preventDefault();
          return;
        }
        if (event.key === "ArrowRight")
        {
          parent.setAttribute("aria-expanded", "false");
          next = parent.nextElementSibling;
          if (next === null)
          {
            next = parent.parentNode.firstElementChild;
          }
          next.setAttribute("aria-expanded", "true");
          next.firstElementChild.firstElementChild.firstElementChild.focus();
          event.stopImmediatePropagation();
          event.preventDefault();
          return;
        }
        if (event.key === "ArrowLeft")
        {
          parent.setAttribute("aria-expanded", "false");
          next = parent.previousElementSibling;
          if (next === null)
          {
            next = parent.parentNode.lastElementChild;
          }
          next.setAttribute("aria-expanded", "true");
          next.firstElementChild.firstElementChild.firstElementChild.focus();
          event.stopImmediatePropagation();
          event.preventDefault();
          return;
        }
        if (event.key === "ArrowUp" || event.key === "ArrowDown")
        {
          if (event.key === "ArrowUp")
          {
            next = this.parentNode.previousElementSibling;
          }
          else
          {
            next = this.parentNode.nextElementSibling;
          }
          if (next === null)
          {
            parent.setAttribute("aria-expanded", "false");
            parent.focus();
          }
          else
          {
            next.firstElementChild.focus();
          }
          event.stopImmediatePropagation();
          event.preventDefault();
        }
      });
    });
  });
}

let url = window.location.pathname;
let cacheName;

function updateCache(cache)
{
  caches.open(cacheName).then(function(tempCache)
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
        caches.delete(cacheName);
      });
    });  
  });
}

function fillCache(cacheAppName)
{
  cacheName = cacheAppName;
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
        let date = response.headers.get("last-modified");
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
        
          caches.open(cacheName).then(function(tempCache)
          {                // Do not fetch HTML because it is already fetched.
            tempCache.addAll(calcURLs).then(function()
            {              // Copy cached resources to main cache and delete this one.
              tempCache.matchAll().then(function(responseArr)
              {            // All responses in array responseArr.
                responseArr.forEach(function(responseTempCache, _index, _array)
                {
                  let urlTemp = responseTempCache.url;
                  let indexZero = url.indexOf("00");
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
              caches.delete(cacheName);
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

function b64decode(str,out)
{
  let ch;
  let idxDest,idxSrc;
  let blocks, leftOver;
  let byte0, byte1, byte2, byte3;
  let conv = new Int8Array(128);
  let len = str.length;
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

let workPar;
function getCalculatorCode(jsFileName, workerParameter)
{
  workPar = workerParameter;
  if (asmjs)
  {
    let req = new XMLHttpRequest();
    req.open("GET", jsFileName, true);
    req.responseType = "arraybuffer";
    req.onreadystatechange = function (_aEvt)
    {
      if (req.readyState === 4 && req.status === 200)
      {
        fileContents = /** @type {ArrayBuffer} */ (req.response);
        if (workPar)
        {
          callWorker(workPar);
        }
      }
    };
    req.send(null);
  }
  else
  {
    let wasm = get("wasmb64").text;
    while (wasm.charCodeAt(0) < 32)
    {
      wasm = wasm.substring(1);
    }    
    while (wasm.charCodeAt(wasm.length-1) < 32)
    {
      wasm = wasm.substring(0, wasm.length-1);
    }    
    let length = wasm.length*3/4;
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
}

function formSend()
{
  if (get("adduserdata").checked)
  {
    getFormSendValue();
  }
  else
  {
    get("userdata").value = "";
  }
  let xhr = new XMLHttpRequest();
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
  let elements = get("formfeedback").elements;
  let contents = "";
  let useAmp = 0;
  for (let i = 0; i < elements.length; i++)
  {
    let element = elements[i >> 0];
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
}