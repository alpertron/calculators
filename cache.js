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
  { // Retrieve from server all files. Use temporary cache.
    const tempCache = await caches.open("cacheTEMP");
    // Do not retrieve alternate JavaScript code if
    // WebAssembly is enabled in browser. Always include current HTML.
    await tempCache.addAll([url].concat((typeof(WebAssembly) === "undefined")?
               getCalcURLs():getCalcURLs().slice(1)));
    // For each file in temporary cache, copy it to cache passed as parameter.
    const matchesArr = await tempCache.matchAll();
    for (const match of matchesArr)
    {
      await cache.put(match.url, match);
    }
  }
  catch (e)
  {  // Nothing to do on error (missing file, error retrieving file, etc.)
  }
  finally
  {  // In any case, delete the temporary cache.
    await caches.delete("cacheTEMP");
  }
}

async function fillCache()
{
  try
  {
    // Test whether the current HTML is already on the cache.
    const cache = await caches.open("newCache");
    // url has the format "/xxxx.HTML".
    const response = await cache.match(url);
    if (typeof response === "undefined")
    {     // HTML is not in cache.
      await updateCache(cache);
    }
    else
    {     // Response is the HTML contents.
      let date = response.headers.get("last-modified");
          // Request the HTML from the Web server.
          // Use non-standard header to tell Service Worker
          // not to retrieve HTML from cache.
      let responseHTML = await fetch(url,
                               {headers:{"If-Modified-Since": date, "x-calc": "1"},
                               cache: "no-store"});
      if (responseHTML.status !== 200)
      {
        return;       // HTML could not be retrieved, so go out.
      }
      if (date === responseHTML.headers.get("last-modified"))
      {
        return;       // HTML has not changed,
                      // so other files have not been changed. Go out.
      }
      // Read files to new cache.
      // Use temporary cache so if there is any network error,
      // original cache is not changed.
      
      try
      {
        let tempCache = await caches.open("cacheTEMP");
        // Do not fetch HTML because it is already fetched.
        await tempCache.addAll((typeof(WebAssembly) === "undefined")?
                  getCalcURLs():getCalcURLs().slice(1));
        // Copy cached resources to main cache and delete this one.
        let matchesArr = await tempCache.matchAll();
        // For each match...
        for (const match of matchesArr)
        {
          let urlTemp = match.url;
          let indexZero = url.indexOf("00");
          if (indexZero > 0)
          {      // There is an old version of this resource on cache to be erased.
            let keysArr = await cache.keys();
            // For each key...
            for (const key of keysArr)
            {    // Traverse cache.
              if (key.url.startsWith(urlTemp.substring(0, indexZero+2)) &&
                  key.url.substring(indexZero+2, indexZero+4) !== urlTemp.substring(indexZero+2, indexZero+4) &&
                  key.url.endsWith(urlTemp.substring(indexZero+4)))
              {  // Old version of asset found (different number and same prefix
                 // and suffix). Delete it from cache.
                await cache.delete(key);
              }  
            }
          }
             // Put resource into cache. 
          await cache.put(urlTemp, match);
        }
        await cache.put(url, responseHTML);
      }
      catch (e)    // Cannot fetch HTML.
      {
        await updateCache(cache);
      }
    }
  }
  finally
  {
    await caches.delete("cacheTEMP");
  }
}

function registerServiceWorker()
{
  initMenubarEvents();
  if ("serviceWorker" in navigator)
  { // Attempt to register service worker.
    // There is no need to do anything on registration success or failure in this JavaScript module.
    navigator["serviceWorker"]["register"]("calcSW.js").
      then(fillCache, function()
    {        
    });
  }
}
