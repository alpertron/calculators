/*
    ULAM SPIRAL
    Copyright(C) 2018 Dario Alejandro Alpern

    This program is free software: you can redistribute it and/or modify
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

/*global instantiate*/
// In order to reduce the number of files to read from Web server, this 
// Javascript file includes both the Javascript in the main thread and the 
// Javascript that drives WebAssembly on its own Web Worker.
(function(global)
{   // This method separates the name space from the Google Analytics code.
var buffer, globals, env, asm;
var zoom, zoomDone, imgData;
var canvas, zoomin, zoomout, center, start;
var isMouseDown;
var currentX, currentY, pixels;
var wasm;
var HEAP8;
var asmjs = typeof(WebAssembly) === "undefined";
var prevX1stTouch, prevY1stTouch;
var prevX2ndTouch, prevY2ndTouch;
var asmGetInformation;
var asmDrawPartialGraphic;
var asmMoveGraphic;
var asmNbrChanged;
var bitsCanvas;
var information;

//##  asmJS goes here (do not change symbols at the left).

function get(x)
{
  return document.getElementById(x);
}

function ptrToString(ptr)
{
  var t=-1;
  var i = 0;
  var str="", outString="";
  do
  {
    for (i=0; i<1024; i++)
    {
      t = HEAP8[((ptr++)>>0)];
      if (t === 0)
      {
        break;
      }
      if (t>=128)
      {
        t = ((t-192)<<6) + HEAP8[((ptr++)>>0)] - 128;
      }
      str += String.fromCharCode(t);
    }
    outString += str;
    str = "";
  } while (t !== 0);
  outString += str;
  return outString;
}

function getHeight()
{
  return Math.min(canvas.scrollHeight, canvas.height);
}

function getWidth()
{
  return Math.min(canvas.scrollWidth, canvas.width);
}

function drawGraphic(ctx, left, top, width, height)
{
  if (width !== 0 && height !== 0)
  {
    var rowNbr, leftPixelSrc, leftPixelDest;
    var startX = left - (getWidth() >> 1);
    var startY = top - (getHeight() >> 1);
    asmDrawPartialGraphic(startX, startX+width, -(startY+height), -startY);
             // Copy from WebAssembly/asm.js buffer to Canvas double buffer.
    leftPixelSrc = top*8192+left*4;
    leftPixelDest = (top+32)*8192+left*4;
    for (rowNbr=0; rowNbr<height; rowNbr++)
    {
      bitsCanvas.set(pixels.subarray(leftPixelSrc, leftPixelSrc + 4*width), leftPixelDest);
      leftPixelSrc += 8192;
      leftPixelDest += 8192;
    }
    ctx.putImageData(imgData, 0, -32, left, 32+top, width, height);
  }
}

function moveGraphic(deltaX, deltaY)
{      
  var ctx = canvas.getContext("2d");
  var width = getWidth();
  var height = getHeight();
  asmMoveGraphic(deltaX, deltaY);
  if (deltaX >= 0)
  {
    if (deltaY >= 0)
    {    // Move to right and bottom.
      if (width-deltaX > 0 && height-deltaY > 0)
      {  // Width and height must be greater than zero to avoid exception.
        ctx.drawImage(ctx.canvas, 0, 0,                          // Topmost pixel of source.
                      width-deltaX, height-deltaY,               // Width and height of source.
                      deltaX, deltaY,                            // Topmost pixel of destination.
                      width-deltaX, height-deltaY);              // Width and height of destination.
        drawGraphic(ctx, 0, 0, width, deltaY);                   // Draw new points at top of graphic.
        drawGraphic(ctx, 0, deltaY, deltaX, height-deltaY);      // Draw new points at left of graphic.
      }
    }
    else
    {    // Move to right and up.
      if (width-deltaX > 0 && height+deltaY > 0)
      {  // Width and height must be greater than zero to avoid exception.
        ctx.drawImage(ctx.canvas, 0, -deltaY,                    // Topmost pixel of source.
                      width-deltaX, height+deltaY,               // Width and height of source.
                      deltaX, 0,                                 // Topmost pixel of destination.
                      width-deltaX, height+deltaY);              // Width and height of destination.
        drawGraphic(ctx, 0, height+deltaY, width, -deltaY);      // Draw new points at bottom of graphic.
        drawGraphic(ctx, 0, 0, deltaX, height+deltaY);           // Draw new points at left of graphic.
      }
    }
  }
  else if (deltaY >= 0)
  {      // Move to left and bottom.
    if (width+deltaX > 0 && height-deltaY > 0)
    {    // Width and height must be greater than zero to avoid exception.
      ctx.drawImage(ctx.canvas, -deltaX, 0,                      // Topmost pixel of source.
                    width+deltaX, height-deltaY,                 // Width and height of source.
                    0, deltaY,                                   // Topmost pixel of destination.
                    width+deltaX, height-deltaY);                // Width and height of destination.
      drawGraphic(ctx, 0, 0, width, deltaY);                     // Draw new points at top of graphic.
      drawGraphic(ctx, width+deltaX, deltaY, -deltaX, height-deltaY); // Draw new points at right of graphic.
    }
  }
  else
  {      // Move to left and up.
    if (width+deltaX > 0 && height+deltaY > 0)
    {    // Width and height must be greater than zero to avoid exception.
      ctx.drawImage(ctx.canvas, -deltaX, -deltaY,                // Topmost pixel of source.
                    width+deltaX, height+deltaY,                 // Width and height of source.
                    0, 0,                                        // Topmost pixel of destination.
                    width+deltaX, height+deltaY);                // Width and height of destination.
      drawGraphic(ctx, 0, height+deltaY, width, -deltaY);        // Draw new points at bottom of graphic.
      drawGraphic(ctx, width+deltaX, 0, -deltaX, height+deltaY); // Draw new points at right of graphic.
    }
  }
}

function updateGraphic(input, nbr)
{
  var startOffset = 11000000|0;
  var idx;
  var value = input.value;
  if (input !== 0)
  {
    for (idx=0; idx<value.length; idx++)  
    {
      HEAP8[(startOffset+idx) >> 0] = value.charCodeAt(idx);
    }
    HEAP8[(startOffset+idx) >> 0] = 0;
  }
  var width = getWidth();
  var height = getHeight();
  asmNbrChanged(startOffset, nbr, width, height);
  drawGraphic(canvas.getContext("2d"), 0, 0, width, height);
}

function b64decode(str,out)
{
  var ch, idx;
  var idxDest,idxSrc;
  var blocks, leftOver;
  var byte0, byte1, byte2, byte3;
  var conv=new Int8Array(128);
  var len=str.length;
  if (str.charAt(len-1) === "=")
  {
    len--;
  }
  if (str.charAt(len-1) === "=")
  {
    len--;
  }
  blocks=len & (-4);
  for (ch = 65; ch <= 90; ch++)   // A - Z
  {
    conv[ch >> 0] = ch - 65;
  }
  for (ch = 97; ch <= 122; ch++)  // a - z
  {
    conv[ch >> 0] = ch - 71;
  }
  for (ch = 48; ch <= 57; ch++)   // 0 - 9
  {
    conv[ch >> 0] = ch + 4;
  }
  conv[43] = 62;                  // +
  conv[33] = 63;                  // !
  for (idxDest=0,idxSrc=0; idxSrc<blocks; idxDest+=3,idxSrc+=4)
  {
    byte0 = conv[str.charCodeAt(idxSrc)];
    byte1 = conv[str.charCodeAt(idxSrc+1)];
    byte2 = conv[str.charCodeAt(idxSrc+2)];
    byte3 = conv[str.charCodeAt(idxSrc+3)];
    
    out[idxDest >>0] = (byte0<<2) + (byte1>>4);
    out[(idxDest+1) >> 0] = (byte1<<4) + (byte2>>2);
    out[(idxDest+2) >> 0] = (byte2<<6) + byte3;
  }
  leftOver = len & 3;
  if (leftOver === 2)
  {
    byte0 = conv[str.charCodeAt(idxSrc)];
    byte1 = conv[str.charCodeAt(idxSrc+1)];
    
    out[idxDest >> 0] = (byte0<<2) + (byte1>>4);
    out[(idxDest+1) >> 0] = byte1<<4;
  }
  else if (leftOver === 3)
  {
    byte0 = conv[str.charCodeAt(idxSrc)];
    byte1 = conv[str.charCodeAt(idxSrc+1)];
    byte2 = conv[str.charCodeAt(idxSrc+2)];
    
    out[idxDest >> 0] = (byte0<<2) + (byte1>>4);
    out[(idxDest+1) >> 0] = (byte1<<4) + (byte2>>2);
    out[(idxDest+2) >> 0] = byte2<<6;
  }
}

function showInfo(ptrInfo)
{
  var infoText = ptrToString(ptrInfo).split("^");
  if (infoText[1] !== center.value)
  {
    center.value = infoText[1];
  } 
  information.innerHTML = infoText[0];
}

function startLowLevelCode()
{
  var length, bytes;
  var info;
  var getPixels;
  var asmJSbuffer;
  var bufSize = 33554432;
  var minusOne = 0xffffffff;
  imgData = canvas.getContext("2d").createImageData(2048, 4096);  // 32 MB;
  buffer = imgData.data.buffer;
  bitsCanvas = new Uint8Array(buffer);
  if (asmjs)
  {                                      // Asm.js initialization.
    asmJSbuffer = new ArrayBuffer(bufSize);
    HEAP8 = new Uint8Array(asmJSbuffer);    // Reserve 32 MB for asm.js variables and buffers.
    env = {"a": {"buffer": asmJSbuffer},
      "abort": function(q) {},
    };
    // check for imul support, and also for correctness ( https://bugs.webkit.org/show_bug.cgi?id=126345 )
    if (!Math["imul"] || Math["imul"](minusOne, 5) !== -5)
    {
      Math["imul"] = function imul(a, b)
      {
        var ah  = a >>> 16;
        var al = a & 0xffff;
        var bh  = b >>> 16;
        var bl = b & 0xffff;
        return (al*bl + ((ah*bl + al*bh) << 16))|0;
      };
    }
    Math.imul = Math["imul"];
    asm = instantiate(env);  // Link asm.js module.
    asmGetInformation = asm.e;
    asmDrawPartialGraphic = asm.c;
    asmMoveGraphic = asm.f;
    asmNbrChanged = asm.g;
    getPixels = asm.d;
    pixels = HEAP8.subarray(getPixels());
    updateGraphic(center, 1);
  }
  else
  {                                      // WebAssembly initialization.
    wasm = get("wasmb64").text;
    while (wasm.charCodeAt(0) < 32)
    {
      wasm = wasm.substring(1);
    }    
    while (wasm.charCodeAt(wasm.length-1) < 32)
    {
      wasm = wasm.substring(0, wasm.length-1);
    }
    length = wasm.length * 3 / 4;
    if (wasm.charCodeAt(wasm.length - 1) === 61)
    {                                    // Base64 ending equal sign found.
      length--;
    }
    if (wasm.charCodeAt(wasm.length - 2) === 61)
    {                                    // Another base64 ending equal sign found.
      length--;
    }
    bytes = new Int8Array(length);
    // Decode Base64.
    b64decode(wasm, bytes);
    info = {};
    WebAssembly["instantiate"](bytes, info).then(function(results)
    {
      asm = results["instance"]["exports"];
      asmGetInformation = asm["getInformation"];
      asmDrawPartialGraphic = asm["drawPartialGraphic"];
      asmMoveGraphic = asm["moveGraphic"];
      asmNbrChanged = asm["nbrChanged"];
      
      HEAP8 = new Uint8Array(asm["memory"]["buffer"]);
      pixels = HEAP8.subarray(asm["getPixels"]());
      updateGraphic(center, 1);
      return;
    });
  }
}

function checkStart()
{
  return (center.value.length < start.value.length || 
        (center.value.length === start.value.length && center.value < start.value));
}

function zoomIn()
{
  if (zoom < 32)
  {
    zoom <<= 1;
    if (zoom === 32)
    {
      zoomin.disabled = true;
    }
    zoomout.disabled = false;
    updateGraphic(0, 3);
  }
}

function zoomOut()
{
  if (zoom > 1)
  {
    zoom >>= 1;
    if (zoom === 1)
    {
      zoomout.disabled = true;
    }
    zoomin.disabled = false;
    updateGraphic(0, 4);
  }
}

function startUp()
{  
  canvas = get("canvas");
  zoomin = get("zoomin");
  zoomout = get("zoomout");
  center = get("center");
  start = get("start");
  information = get("info");
  zoom = 8;
  zoomDone = 0;
  isMouseDown = false;
  var domRect = canvas.getBoundingClientRect();
  canvas.width = domRect.width;
  canvas.height = domRect.height;
  startLowLevelCode();
  zoomin.onclick = function ()
  {
    zoomIn();
  };
  zoomout.onclick = function ()
  {
    zoomOut();
  };
  canvas.onkeydown = function (e)
  {
    if (checkStart())
    {    // Cannot show spiral when center is less than start value.
      return;
    }
    var evt = e || window.event;
    var key = evt.keyCode;
    switch (key)
    {
      case 37: // Left arrow:
        moveGraphic(4, 0);
        showInfo(asmGetInformation(-1, -1));
        evt.preventDefault();          // Do not propagate this key.
        break; 
      case 38: // Up arrow:
        moveGraphic(0, 4);
        showInfo(asmGetInformation(-1, -1));
        evt.preventDefault();          // Do not propagate this key.
        break; 
      case 39: // Right arrow:
        moveGraphic(-4, 0);
        showInfo(asmGetInformation(-1, -1));
        evt.preventDefault();          // Do not propagate this key.
        break; 
      case 40: // Down arrow:
        moveGraphic(0, -4);
        showInfo(asmGetInformation(-1, -1));
        evt.preventDefault();          // Do not propagate this key.
        break; 
    }      
  };
  canvas.onmousemove = function(e)
  {
    if (checkStart())
    {    // Cannot show spiral when center is less than start value.
      return;
    }
    var evt = e || window.event;
    var newX = evt.clientX;
    var newY = evt.clientY;      // Coordinates relative to window.
    var rect = canvas.getBoundingClientRect();
    newX -= rect.left;
    newY -= rect.top;            // Now coordinates are relative to canvas.
    if (isMouseDown)
    {
      if (newX !== currentX || newY !== currentY)
      {
        moveGraphic(newX - currentX, newY - currentY);
        currentX = newX;
        currentY = newY;
      }
    } 
    else
    {
      if (newX < canvas.width && newY < canvas.height)
      {
        showInfo(asmGetInformation(newX, newY));
      }
    }
  };  
  canvas.onmousedown = function(e)
  {
    var evt = e || window.event;
    currentX = evt.clientX;
    currentY = evt.clientY;      // Coordinates relative to window. 
    var rect = canvas.getBoundingClientRect();
    currentX -= rect.left;
    currentY -= rect.top;        // Now coordinates are relative to canvas.
    isMouseDown = true;
  };
  document.onmouseup = function()
  {
    isMouseDown = false;
  };
  canvas.addEventListener("touchstart", function(e)
  {
    var evt = e || window.event;
    var touches = evt.targetTouches;
    var touch = touches[0];
    prevX1stTouch = Math.round(touch.pageX);
    prevY1stTouch = Math.round(touch.pageY);
    if (touches.length === 2)
    {
      touch = touches[1];
      prevX2ndTouch = Math.round(touch.pageX);
      prevY2ndTouch = Math.round(touch.pageY);
      zoomDone = 0;
    }    
  }, false);  
  canvas.addEventListener("touchmove", function(e)
  {
    var touch2, oldDist, newDist, diffX, diffY;
    var evt = e || window.event;
    var touches = evt.targetTouches;
    var touch1 = touches[0];
    var newX = Math.round(touch1.pageX);
    var newY = Math.round(touch1.pageY);
    if (touches.length === 1)
    {      // Drag gesture.
      if (newX !== prevX1stTouch || newY !== prevY1stTouch)
      {
        moveGraphic(newX - prevX1stTouch, newY - prevY1stTouch);
        prevX1stTouch = newX;
        prevY1stTouch = newY;
        showInfo(asmGetInformation(-1, -1));
      }      
    }
    else if (touches.length === 2 && !zoomDone)
    {      // Pinch (zoom in and zoom out) gesture.
      touch2 = evt.touches[1];
           // Get square of distance between fingers.
      diffX = prevX2ndTouch - prevX1stTouch;
      diffY = prevY2ndTouch - prevY1stTouch;
      oldDist = diffX * diffX + diffY * diffY;
      diffX = Math.round(touch2.pageX) - newX;
      diffY = Math.round(touch2.pageY) - newY;
      newDist = diffX * diffX + diffY * diffY;
      if (newDist > oldDist)
      {    // Zoom in.
        zoomIn();
        zoomDone = 1;
      }
      if (newDist < oldDist)
      {    // Zoom out.
        zoomOut();
        zoomDone = 1;
      }
    }
  }, false);
  center.onkeydown = function(e)
  {
    var evt = e || window.event;
    var key = evt.keyCode;
    if (!evt.ctrlKey && !evt.altKey && !evt.metaKey)
    {                                    // No modifier key pressed.
      if (key !== 8 && key !== 9 && key !== 37 && key !== 39 && key !== 45 && key !== 46)
      {                                  // Not backspace, tab, right or left arrow, insert or delete key.
        if (key >= 0x60 && key <= 0x69)
        {
          key -= 0x30;                   // Convert numpad key to standard digit key.
        }
        if (key < 0x30 || key > 0x39 || center.value.length >= 18)
        {                                // Key is not a digit or number is too large.
          evt.preventDefault();          // Do not propagate this key.
        }
      }
    }
  };
  center.oninput = function()
  {
    var ctx;
    if (checkStart())
    {
      information.innerHTML = get("cannotShow").innerHTML;
      ctx = canvas.getContext("2d");
      ctx.beginPath();
      ctx.rect(0, 0, canvas.width, canvas.height);
      ctx.fillStyle = "black";
      ctx.fill();
    }
    else
    {
      information.innerHTML = "";
      updateGraphic(center, 1);
    }
  };
  start.onkeydown = function(e)
  {
    var evt = e || window.event;
    var key = evt.keyCode;
    if (!evt.ctrlKey && !evt.altKey && !evt.metaKey)
    {                                    // No modifier key pressed.
      if (key !== 8 && key !== 9 && key !== 37 && key !== 39 && key !== 45 && key !== 46)
      {                                  // Not backspace, tab, right or left arrow, insert or delete key.
        if (key >= 0x60 && key <= 0x69)
        {
          key -= 0x30;                   // Convert numpad key to standard digit key.
        }
        if (key < 0x30 || key > 0x39 || start.value.length >= 18)
        {                                // Key is not a digit or number is too large.
          evt.preventDefault();          // Do not propagate this key.
        }
      }
    }
  };
  start.oninput = function()
  {
    showInfo(asmGetInformation(-1, -1));
    updateGraphic(start, 2);
  };
  window.onresize = function()
  {
    var domRect = canvas.getBoundingClientRect();
    canvas.width = domRect.width;
    canvas.height = domRect.height;
    updateGraphic(center, 1);
  };
}

window.addEventListener("load", startUp);
})(this);
