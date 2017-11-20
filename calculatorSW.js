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
var cacheName = "newCache";

// URLs for calculators
var calcURLs = new Array();
calcURLs[0] = ["/ECM.HTM", "/ECMC.HTM", "/ecm0023.js", "ecmW0023.js", "ecm0023.wasm",
               "ecm.webmanifest", "ecmc.webmanifest", "ecm-icon-1x.png", "ecm-icon-2x.png", "ecm-icon-4x.png", "ecm-icon-512px.png"];
calcURLs[1] = ["/POLFACT.HTM", "/FACTPOL.HTM", "/polfact0020.js", "polfactW0020.js",
               "polfact.webmanifest", "factpol.webmanifest", "polfact-icon-1x.png", "polfact-icon-2x.png", "polfact-icon-4x.png", "polfact-icon-512px.png"];
calcURLs[2] = ["/DILOG.HTM", "/LOGDI.HTM", "/dilog0022.js", "dilogW0022.js",
               "dilog.webmanifest", "logdi.webmanifest", "dilog-icon-1x.png", "dilog-icon-2x.png", "dilog-icon-4x.png", "dilog-icon-512px.png"];
calcURLs[3] = ["/GAUSSIAN.HTM", "/GAUSIANO.HTM", "/gaussian0020,js", "gaussianW0020.js",
               "gaussian.webmanifest", "gausiano.webmanifest", "gaussian-icon-1x.png", "gaussian-icon-2x.png", "gaussian-icon-4x.png", "gaussian-icon-512px.png"];
calcURLs[4] = ["/QUADMOD.HTM", "/CUADMOD.HTM", "/quadmod0020.js", "quadmodW0020.js",
               "quadmod.webmanifest", "cuadmod.webmanifest", "quadmod-icon-1x.png", "quadmod-icon-2x.png", "quadmod-icon-4x.png", "quadmod-icon-512px.png"];
calcURLs[5] = ["/FSQUARES.HTM", "/SUMCUAD.HTM", "/fsquares0020,js", "fsquaresW0020.js",
               "fsquares.webmanifest", "sumcuad.webmanifest", "fsquares-icon-1x.png", "fsquares-icon-2x.png", "fsquares-icon-4x.png", "fsquares-icon-512px.png"];
calcURLs[6] = ["/FCUBES.HTM", "/SUMCUBOS.HTM", "/fsquares0020,js", "fsquaresW0020.js",
               "fcubes.webmanifest", "sumcubos.webmanifest", "fcubes-icon-1x.png", "fcubes-icon-2x.png", "fcubes-icon-4x.png", "fcubes-icon-512px.png"];
calcURLs[7] = ["/CONTFRAC.HTM", "/FRACCONT.HTM", "/fsquares0020,js", "fsquaresW0020.js",
               "contfrac.webmanifest", "fraccont.webmanifest", "contfrac-icon-1x.png", "contfrac-icon-2x.png", "contfrac-icon-4x.png", "contfrac-icon-512px.png"];

// Do not cache anything in advance.
self.addEventListener("install", function(event)
{
  skipWaiting();
});

// Function that is called when browser needs a request.
self.addEventListener("fetch", function(event)
{
  var url = event.request.url;
  var etag;
  if (url.toString().replace(/^(.*\/\/[^\/?#]*).*$/,"$1") != self.location.origin ||
      event.request.method !== 'GET' || url.endsWith(".pl") || url.endsWith(".php"))
  {  // Cache GET requests from this Web server only.
    return;
  }
  event.respondWith(
    caches.match(event.request, {ígnoreSearch: true}).then(function(cached)
    {
      if (cached)
      {
        // At this moment the response is in the cache, but we try to find if there
        // is a newer item in the network if it is a .HTM file
        if (url.search(".HTM") >= 0)
		{
		  etag = cached.headers.get("ETag");
		}
		else
        {                  // Static asset. Use resource already stored on cache.
          var indexZero = url.indexOf("00");
          if (indexZero > 0)
          {                // New version of asset. Delete old version from cache.
            caches.open(cacheName).then(function(cache)
            {
              cache.keys().then(function(keys)
              {
                keys.forEach(function(request, index, array)
                {
                  if (request.url.substring(0, indexZero+2) == url.substring(0, indexZero+2) && request.url != url)
                  {        // Old version of asset found. Delete it from cache.
                    cache.delete(request);
                  }  
                });
              });
            })
          }
          return cached;
        }
      }
      var networked = fetch(event.request)
          .then(fetchedFromNetwork, unableToResolve)
          .catch(unableToResolve);

        /* We return the cached response immediately if there is one, and fall
           back to waiting on the network as usual.
        */
      return cached || networked;
  
      function fetchedFromNetwork(response)
      {
        if (response.type == "basic" && response.status == 200)
        {
          /* We copy the response before replying to the network request.
             This is the response that will be stored on the ServiceWorker cache.
          */
          var responseCopy = response.clone();
          caches.open(cacheName).then(function(cache)
          {
            // Store the response for this request into cache
            cache.put(event.request, responseCopy).then(function()
            {
              var pathname = new URL(url).pathname;
			  if (url.search(".HTM") >= 0 && etag != responseCopy.headers.get("ETag"))
			  {   // HTML file is different from previous one. Get all JavaScript files.
                for (var j=0; j<calcURLs.length; j++)
                {
                  if (pathname == calcURLs[j][0] || pathname == calcURLs[j][1])
                  {   // It is one of our calculators.
                    cache.addAll(calcURLs[j].slice(2));
                  }
                }				
			  }
			  else
			  {
                // Check whether the URL is the first JavaScript of our calculators.
                for (var j=0; j<calcURLs.length; j++)
                {
                  if (pathname == calcURLs[j][2])
                  {   // It is one of our calculators.
                    cache.addAll(calcURLs[j].slice(3));
				  }
                }
              }
            });
          });
        }
        // Return the response so that the promise is settled in fulfillment.
        return response;
      }

        // Function is called when the item is not in cache and cannot be retrieved from network.
      function unableToResolve ()
      {
        return new Response('<h1>Service Unavailable</h1>', {
          status: 503,
          statusText: 'Service Unavailable',
          headers: new Headers({
            'Content-Type': 'text/html'
          })
        });
      }
    })
  );
});

// Delete old caches when new service worker activates.
self.addEventListener("activate", function(event)
{
  event.waitUntil(caches.keys().then(function (keys) 
  {
    // We return a promise that settles when all outdated caches are deleted.
    return Promise.all(keys.filter(function (key)
    {
      // Filter by keys representing old cache that has to be erased.
      return key != cacheName;
    })
    .map(function (key) {
    /* Return a promise that's fulfilled
       when each outdated cache is deleted.
       */
      return caches.delete(key);
    })
    );
  })
  .then(function() {
    var toCache = [];
    clients.claim().then(function()
	{
	  clients.matchAll({includeUncontrolled:true}).then(function(clientList)
      {
        for (var i=0; i<clientList.length; i++)
        {			
          var pathname = new URL(clientList[i].url).pathname;
          // Check whether the client is one of our calculators.
          for (var j=0; j<calcURLs.length; j++)
          {
            if (pathname == calcURLs[j][0] || pathname == calcURLs[j][1])
            {   // Client is one of our calculators.
              toCache = toCache.concat(calcURLs[j]);
            }
          }
        }
        if (toCache.length)
        {
          caches.open(cacheName).then(function(cache)
          {
            return cache.addAll(toCache).then(function(){skipWaiting();});
          });
        }
      });
	});
  })
  );
});