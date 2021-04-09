function convertToString(ptr, str)
{
  var dest = ptr;
  var length = str.length;
  var i, t;
  for (i=0; i<length; i++)
  {
    t = str.charCodeAt(i);
    if (t<128)
    {
      HEAPU8[dest++] = t;
    }
    else
    {
      HEAPU8[dest++] = (t >> 6) + 192;
      HEAPU8[dest++] = (t & 63) + 128;
    }
  }
  HEAPU8[dest] = 0;
}

Module =
{
  "preRun": function()
  {
    self.onmessage = function(e)
    {
      convertToString(_getInputStringPtr(), e.data);
      _doWork();
    };
  },
  "noInitialRun": true,
};

