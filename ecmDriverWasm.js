var exports,HEAPU8;

function PtrToString(ptr)
{
  var t=-1;
  var i = 0;
  var str="", outString="";
  do
  {
    for (i=0; i<1024; i++)
    {
      t = HEAPU8[((ptr++)>>0)];
      if (t==0)
      {
        break;
      }
      if (t>=128)
      {
        t = ((t-192)<<6) + HEAPU8[((ptr++)>>0)] - 128;
      }
      str += String.fromCharCode(t);
    }
    outString += str;
    str = "";
  } while (t!=0);
  outString += str + String.fromCharCode(0);
  return outString;
}

function ConvertToString(ptr, str)
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

var env =
{
  "databack": function(data)
  {
    self.postMessage(PtrToString(data));
  },
  "tenths": function()
  {
    return Math.floor(new Date().getTime() / 100);
  },
  "getCunn": function(data)
  {
    var req = new XMLHttpRequest();
    req.open('GET', PtrToString(data), false);
    req.send(null);
    if (req.status == 200)
    {
      ConvertToString(exports.getFactorsAsciiPtr(), req.responseText);
    }
  }
};

var info =
{
  "env": env
};  

self.onmessage = function(e)
{
  var request = new XMLHttpRequest();
  request.open('GET', 'ecm0025.wasm');
  request.responseType = 'arraybuffer';
  request.send();

  request.onload = function()
  {
	if (request.status != 200)
	{
      return;
	}
    var bytes = request.response;
    WebAssembly["instantiate"](bytes, info).then(function(results)
    {
      exports = results["instance"]["exports"];
      HEAPU8 = new Uint8Array(exports["memory"]["buffer"]);
      ConvertToString(exports["getInputStringPtr"](), e.data);
      exports["doWork"]();
	  return;
    });
  };
}

