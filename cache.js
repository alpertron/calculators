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
/* global getCalcURLs */
/* global initMenubarEvents */
let url = window.location.pathname;

function updateCache(cache)
{
  caches.open("cacheTEMP").then(function(tempCache)
  {
    tempCache.addAll([url].concat((typeof(WebAssembly) === "undefined")?
               getCalcURLs():getCalcURLs().slice(1))).then(function()
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
        caches.delete("cacheTEMP");
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
        
          caches.open("cacheTEMP").then(function(tempCache)
          {                // Do not fetch HTML because it is already fetched.
            tempCache.addAll((typeof(WebAssembly) === "undefined")?
               getCalcURLs():getCalcURLs().shift).then(function()
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
              caches.delete("cacheTEMP");
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

function registerServiceWorker()
{
  if ("serviceWorker" in navigator)
  { // Attempt to register service worker.
    // There is no need to do anything on registration success or failure in this JavaScript module.
    navigator["serviceWorker"]["register"]("calcSW.js").then(
              function()
              {
                /* Nothing to do */
              },
              function()
              {
                /* Nothing to do */
              });
    fillCache();
  }
  initMenubarEvents();
}
