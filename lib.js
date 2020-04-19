mergeInto(LibraryManager.library, 
{
  databack: function(data)
  {
    function PointerStringify(offset)
    {
      var charCache = new Array(128);  // Preallocate the cache for the common single byte chars
      var charFromCodePt = String.fromCodePoint || String.fromCharCode;
      var result = [];

      var codePt, byte1;

      result.length = 0;

      while (HEAPU8[offset] != 0)
      {
        byte1 = HEAPU8[offset++];
        if (byte1 <= 0x7F)
        {
          codePt = byte1;
        }
        else if (byte1 <= 0xDF)
        {
          codePt = ((byte1 & 0x1F) << 6) | (HEAPU8[offset++] & 0x3F);
        }
        result.push(charCache[codePt] || (charCache[codePt] = charFromCodePt(codePt)));
      }
      return result.join('');
    }
    self.postMessage(PointerStringify(data));
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
    function PointerStringify(offset)
    {
      var charCache = new Array(128);  // Preallocate the cache for the common single byte chars
      var charFromCodePt = String.fromCodePoint || String.fromCharCode;
      var result = [];

      var codePt, byte1;
 
      result.length = 0;

      while (HEAPU8[offset] != 0)
      {
        byte1 = HEAPU8[offset++];
        if (byte1 <= 0x7F)
        {
          codePt = byte1;
        }
        else if (byte1 <= 0xDF)
        {
          codePt = ((byte1 & 0x1F) << 6) | (HEAPU8[offset++] & 0x3F);
        }
        result.push(charCache[codePt] || (charCache[codePt] = charFromCodePt(codePt)));
      }
      return result.join('');
    }
    var copyString = Module.cwrap('copyString', 'number', ['string']);
    var req = new XMLHttpRequest();
    req.open('GET', PointerStringify(data), false);
    req.send(null);
    if (req.status == 200)
    {
      copyString(req.responseText);
    }
  }
});
