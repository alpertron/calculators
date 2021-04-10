/*global mergeInto, HEAPU8, Module*/
mergeInto(LibraryManager.library, 
{
  databack: function(data)
  {
    function pointerStringify(offset)
    {
      var charCache = new Array(128);  // Preallocate the cache for the common single byte chars
      var charFromCodePt = String.fromCodePoint || String.fromCharCode;
      var result = [];

      var codePt, byte1;

      result.length = 0;

      while (HEAPU8[offset >> 0] !== 0)
      {
        byte1 = HEAPU8[(offset++) >> 0];
        if (byte1 <= 0x7F)
        {
          codePt = byte1;
        }
        else if (byte1 <= 0xDF)
        {
          codePt = ((byte1 & 0x1F) << 6) | (HEAPU8[(offset++) >> 0] & 0x3F);
        }
        result.push(charCache[codePt >> 0] || (charCache[codePt >> 0] = charFromCodePt(codePt)));
      }
      return result.join("");
    }
    self.postMessage(pointerStringify(data));
  },
  stamp: function()
  {
    return Math.floor(new Date().getTime() / 1000);
  },
  tenths: function()
  {
    return Math.floor(new Date().getTime() / 100);
  },
  startSkipTest: function()
  {
    self.postMessage("51");   // Show Skip Test button on screen
  },
  endSkipTest: function()
  {
    self.postMessage("52");   // Hide Skip Test button from screen
  },
  getCunn: function(data)
  {
    function pointerStringify(offset)
    {
      var charCache = new Array(128);  // Preallocate the cache for the common single byte chars
      var charFromCodePt = String.fromCodePoint || String.fromCharCode;
      var result = [];

      var codePt, byte1;
 
      result.length = 0;

      while (HEAPU8[offset >> 0] !== 0)
      {
        byte1 = HEAPU8[(offset++) >> 0];
        if (byte1 <= 0x7F)
        {
          codePt = byte1;
        }
        else if (byte1 <= 0xDF)
        {
          codePt = ((byte1 & 0x1F) << 6) | (HEAPU8[(offset++) >> 0] & 0x3F);
        }
        result.push(charCache[codePt >> 0] || (charCache[codePt >> 0] = charFromCodePt(codePt)));
      }
      return result.join("");
    }
    var copyString = Module.cwrap("copyString", "number", ["string"]);
    var req = new XMLHttpRequest();
    req.open("GET", pointerStringify(data), false);
    req.send(null);
    if (req.status === 200)
    {
      copyString(req.responseText);
    }
  }
});
