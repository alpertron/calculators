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

async function updateCache(cache)
{
  try
  {
    const tempCache = await caches.open("cacheTEMP");
    await tempCache.addAll([url].concat((typeof(WebAssembly) === "undefined")?
               getCalcURLs():getCalcURLs().slice(1)));
    const responseArr = await tempCache.matchAll();
    responseArr.forEach(function(responseTempCache, _index, _array)
    {
      cache.put(responseTempCache.url, responseTempCache).
            then(function(){}, function(){});
    });
  } catch (e)
  {
  } finally
  {
    await caches.delete("cacheTEMP");
  }
}

async function fillCache()
{
  try
  {
    // Test whether the HTML is already on the cache.
    const cache = await caches.open("newCache");
    const response = await cache.match(url); 
    if (typeof response === "undefined")
    {     // HTML is not in cache.
      await updateCache(cache);
    }
    else
    {     // Response is the HTML contents.
      let date = response.headers.get("last-modified");
          // Request the HTML from the Web server.
          // Use non-standard header to tell Service Worker not to retrieve HTML from cache.
      let responseHTML = await fetch(url,
                               {headers:{"If-Modified-Since": date, "x-calc": "1"},
                               cache: "no-store"});
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
      
      try
      {
        let tempCache = await caches.open("cacheTEMP");
        // Do not fetch HTML because it is already fetched.
        await tempCache.addAll((typeof(WebAssembly) === "undefined")?
                  getCalcURLs():getCalcURLs().shift);
        // Copy cached resources to main cache and delete this one.
        let responseArr = tempCache.matchAll();
        // All responses in array responseArr.
        responseArr.forEach(async function(responseTempCache, _index, _array)
        {
          let urlTemp = responseTempCache.url;
          let indexZero = url.indexOf("00");
          if (indexZero > 0)
          {        // There is an old version of this resource on cache to be erased.
            let keys = await cache.keys();
            keys.forEach(function(requestCache, _idx, _arr)
            {    // Traverse cache.
              if (requestCache.url.substring(0, indexZero+2) === urlTemp.substring(0, indexZero+2) &&
                  requestCache.url.substring(indexZero+2, indexZero+4) !== urlTemp.substring(indexZero+2, indexZero+4) &&
                  requestCache.url.substring(indexZero+4) === urlTemp.substring(indexZero+4))
              {  // Old version of asset found (different number and same prefix and suffix). Delete it from cache.
                cache.delete(requestCache).then(function(){}, function(){});
              }  
              // Put resource into cache after old resource has been erased.
              cache.put(urlTemp, responseTempCache).then(function(){}, function(){});
            });
          }
          else
          {   // Put resource into cache (no old resource into cache). 
            await cache.put(urlTemp, responseTempCache);
          }
        });
        await cache.put(url, responseHTML);
      } catch (e)    // Cannot fetch HTML.
      {
        await updateCache(cache);
      }
    }
  } finally
  {
    await caches.delete("cacheTEMP");
  }
}

async function registerServiceWorker()
{
  if ("serviceWorker" in navigator)
  { // Attempt to register service worker.
    // There is no need to do anything on registration success or failure in this JavaScript module.
    await navigator["serviceWorker"]["register"]("calcSW.js");
    await fillCache();
  }
  initMenubarEvents();
}
