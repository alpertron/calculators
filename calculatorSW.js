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
var version = '0017::';




self.addEventListener("install", function(event)
{  // Add requests to cache during service worker installation.
  event.waitUntil(caches
      .open(version + 'ecm')
      .then(function(cache) {
        return cache.addAll([
          '/ECM.HTM',
          '/ECMC.HTM',
          '/ecm0006.js',
          '/ecmW0006.js',
          '/ecm0006.wasm',
		  '/ecm.json',
		  '/ecm-icon-1x.png',
		  '/ecm-icon-2x.png',
		  '/ecm-icon-4x.png',
        ]);
      })
      .then(function() {
        // Install complete.
      })
  );
});

// Function that is called when browser needs a request.
self.addEventListener("fetch", function(event) {

  if (event.request.method !== 'GET') {
    // Cache GET requests only.
    return;
  }
  event.respondWith(
    caches.match(event.request).then(function(cached) {
		// At this moment the response is in the cache, but we try to find if there
		// is a newer item in the network.
        var networked = fetch(event.request)
          .then(fetchedFromNetwork, unableToResolve)
          .catch(unableToResolve);

        /* We return the cached response immediately if there is one, and fall
           back to waiting on the network as usual.
        */
        return cached || networked;

        function fetchedFromNetwork(response) {
		  if (response.url.startsWith("https://www.alpertron.com.ar")) {
            /* We copy the response before replying to the network request.
               This is the response that will be stored on the ServiceWorker cache.
            */
            var cacheCopy = response.clone();
            caches.open(version + 'pages').then(function add(cache) {
                // Store the response for this request into cache
                cache.put(event.request, cacheCopy);
              })
              .then(function() {
              });
          }
          // Return the response so that the promise is settled in fulfillment.
          return response;
        }

		// Function is called when the item is not in cache and cannot be retrieved from network.
        function unableToResolve () {
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
      // Filter by keys that don't start with the latest version prefix.
      return !key.startsWith(version);
    })
    .map(function (key) {
    /* Return a promise that's fulfilled
       when each outdated cache is deleted.
       */
      return caches.delete(key);
    })
    );
  })
  .then(function() {})
  );
});