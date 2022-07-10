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
// Do not cache anything in advance.
self.addEventListener("install", function(_event)
{
  self.skipWaiting();
});

function fetchedFromNetwork(response)
{
         // Return the response so that the promise is settled in fulfillment.
  return response;
}

         // Function is called when the item is not in cache and cannot be retrieved from network.
function unableToResolve()
{
  return new Response("<h1>Cannot connect to Web server</h1>", {
          status: 503,
          statusText: "Service Unavailable",
          headers: new Headers({"Content-Type": "text/html"})
  });
}

// Function that is called when browser needs a request.
self.addEventListener("fetch", function(event)
{
  const url = event.request.url;
      // Check if special header indicating not to read from cache has arrived.
  const nocache = event.request.headers.get("x-calc");
  const noQueryString;
  if (url.toString().replace(/^(.*\/\/[^\/?#]*).*$/,"$1") !== self.location.origin ||
      event.request.method !== "GET" || url.endsWith(".pl") || url.endsWith(".php"))
  {  // Cache GET requests from this Web server only.
    return;
  }
  // Erase any query information.
  const QmarkOffset = url.indexOf("?");
  if (QmarkOffset < 0)
  {
    noQueryString = url;
  }
  else
  {
    noQueryString = url.substring(0, QmarkOffset);
  }
  event.respondWith(  
    caches.match(noQueryString).then(function(cached)
    {
      if (cached && (nocache !== "1" || !navigator.onLine))
      {        // At this moment the response is in the cache.
        return cached;
      }
               // URL not in cache. Get resource from network.
      return fetch(event.request)
          .then(fetchedFromNetwork, unableToResolve)
          .catch(unableToResolve);
  
    })
  );
});

