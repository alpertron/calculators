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

/* global b64decode */
/* global get */
/* global initMenubarEvents */
/* global instantiate */
/* global isNotSpecialKey */
// In order to reduce the number of files to read from Web server, this 
// Javascript file includes both the Javascript in the main thread and the 
// Javascript that drives WebAssembly on its own Web Worker.
const none = "none";
const block = "block";
const inline = "inline";
const asmjs = typeof(WebAssembly) === "undefined";
let buffer, env, asm;
let zoom, zoomDone, imgData;
let canvas, zoomin, zoomout, center, start;
let isMouseDown;
let currentX, currentY, pixels;
let wasm;
let HEAP8;
let prevX1stTouch, prevY1stTouch;
let prevX2ndTouch, prevY2ndTouch;
let asmGetInformation;
let asmDrawPartialGraphic;
let asmMoveGraphic;
let asmNbrChanged;
let bitsCanvas;
let information;
let applet;
let animform;
let animate;
let stop;
let xincr;
let yincr;
let sincr;
let incrX;
let incrY;
let incrStart;
let oldX, oldY, oldStart;
let delay;
let doanimate;
let cancelanim;
let interval;

//##  asmJS goes here (do not change symbols at the left).

function ptrToString(ptr)
{
  let t=-1;
  let i = 0;
  let str="", outString="";
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
    let rowNbr, leftPixelSrc, leftPixelDest;
    let startX = left - (getWidth() >> 1);
    let startY = top - (getHeight() >> 1);
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
  let ctx = canvas.getContext("2d");
  let width = getWidth();
  let height = getHeight();
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
  let startOffset = 11000000|0;
  let idx;
  let value = input.value;
  if (input !== 0)
  {
    for (idx=0; idx<value.length; idx++)  
    {
      HEAP8[(startOffset+idx) >> 0] = value.charCodeAt(idx);
    }
    HEAP8[(startOffset+idx) >> 0] = 0;
  }
  let width = getWidth();
  let height = getHeight();
  asmNbrChanged(startOffset, nbr, width, height);
  drawGraphic(canvas.getContext("2d"), 0, 0, width, height);
}

function showInfo(ptrInfo)
{
  let infoText = ptrToString(ptrInfo).split("^");
  if (infoText[1] !== center.value)
  {
    center.value = infoText[1];
  } 
  information.innerHTML = infoText[0];
}

function startLowLevelCode()
{
  let length, bytes;
  let info;
  let getPixels;
  let asmJSbuffer;
  let bufSize = 33554432;
  imgData = canvas.getContext("2d").createImageData(2048, 4096);  // 32 MB;
  buffer = imgData.data.buffer;
  bitsCanvas = new Uint8Array(buffer);
  if (asmjs)
  {                                      // Asm.js initialization.
    asmJSbuffer = new ArrayBuffer(bufSize);
    HEAP8 = new Uint8Array(asmJSbuffer);    // Reserve 32 MB for asm.js variables and buffers.
    env = {"a": {"buffer": asmJSbuffer},
      "abort": function(_q)
               {
                 /* Not used*/
               },
    };
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

function animation()
{
  let currX = oldX - incrX * zoom;
  let currY = oldY + incrY * zoom;
  let currStart = oldStart + incrStart;
  let diffX = Math.floor(currX) - Math.floor(oldX);
  let diffY = Math.floor(currY) - Math.floor(oldY);
  let diffStart = Math.floor(currStart) - Math.floor(oldStart);
  if (diffX !== 0 || diffY !== 0)
  {
    moveGraphic(diffX, diffY);
  }
  if (diffStart !== 0)
  {
    start.value = (parseInt(start.value, 10) + diffStart).toString();
    updateGraphic(start, 2);
  }
  if (diffX !== 0 || diffY !== 0 || diffStart !== 0)
  {
    showInfo(asmGetInformation(-1, -1));
  }
  oldX = currX;
  oldY = currY;
  oldStart = currStart;
}

function keydown(evt)
{
  let target = evt.target || evt.srcElement;
  let key = evt.key;
  if (isNotSpecialKey(evt))
  {
    if (!evt.ctrlKey && !evt.altKey && !evt.metaKey)
    {                                  // No modifier key pressed.
      if ((key >= "0" && key <= "9") || key === "-")
      {                                // Digit key has been pressed.
        if (target.value.length >= 18 ||
            (target.value.charAt(0) !== "-" && target.value.length >= 19))
        {                              // Number is too large.
          evt.preventDefault();        // Do not propagate this key.
        }
      }
      else 
      {                                // Not backspace, tab, right or left arrow, insert or delete key.
        evt.preventDefault();          // Do not propagate this key.  
      }
    }
  }
}

function updateInput(input)
{
  let currString = input.value;
  let countSrc;
  let newString = "";
  let c;
  if (currString.indexOf("-") >= 0)
  {
    newString = "-";
  }
  // Delete any character that is not 0 to 9.
  for (countSrc = 0; countSrc < currString.length; countSrc++)
  {
    c = currString.substring(countSrc, countSrc + 1);
    if (c >= "0" && c <= "9")
    {
      newString += c;
    }
  }
  if (input.value !== newString)
  {
    input.value = newString;
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
  animate = get("animate");
  animform = get("animform");
  stop = get("stop");
  xincr = get("xincr");
  yincr = get("yincr");
  sincr = get("sincr");
  delay = get("delay");
  doanimate = get("doanimate");
  cancelanim = get("cancelanim");
  applet = get("applet");
  start.value = "1";
  center.value = "1";
  zoom = 8;
  zoomDone = 0;
  isMouseDown = false;
  let domRect = canvas.getBoundingClientRect();
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
  canvas.onkeydown = function (evt)
  {
    if (checkStart())
    {    // Cannot show spiral when center is less than start value.
      return;
    }
    switch (evt.key)
    {
      case "Left":        // Left arrow:
      case "ArrowLeft":   // Left arrow:
        moveGraphic(4, 0);
        showInfo(asmGetInformation(-1, -1));
        evt.preventDefault();          // Do not propagate this key.
        break; 
      case "Up":         // Up arrow:
      case "ArrowUp":    // Up arrow:
        moveGraphic(0, 4);
        showInfo(asmGetInformation(-1, -1));
        evt.preventDefault();          // Do not propagate this key.
        break; 
      case "Right":      // Right arrow:
      case "ArrowRight": // Right arrow:
        moveGraphic(-4, 0);
        showInfo(asmGetInformation(-1, -1));
        evt.preventDefault();          // Do not propagate this key.
        break; 
      case "Down":       // Down arrow:
      case "ArrowDown":  // Down arrow:
        moveGraphic(0, -4);
        showInfo(asmGetInformation(-1, -1));
        evt.preventDefault();          // Do not propagate this key.
        break; 
    }      
  };
  canvas.onmousemove = function(evt)
  {
    if (checkStart())
    {    // Cannot show spiral when center is less than start value.
      return;
    }
    let newX = evt.clientX;
    let newY = evt.clientY;      // Coordinates relative to window.
    let rect = canvas.getBoundingClientRect();
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
  canvas.onmousedown = function(evt)
  {
    currentX = evt.clientX;
    currentY = evt.clientY;      // Coordinates relative to window. 
    let rect = canvas.getBoundingClientRect();
    currentX -= rect.left;
    currentY -= rect.top;        // Now coordinates are relative to canvas.
    isMouseDown = true;
  };
  document.onmouseup = function()
  {
    isMouseDown = false;
  };
  canvas.addEventListener("touchstart", function(evt)
  {
    let touches = evt.targetTouches;
    let touch = touches[0];
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
  canvas.addEventListener("touchmove", function(evt)
  {
    let touch2, oldDist, newDist, diffX, diffY;
    let touches = evt.targetTouches;
    let touch1 = touches[0];
    let newX = Math.round(touch1.pageX);
    let newY = Math.round(touch1.pageY);
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
  center.onkeydown = keydown;
  center.oninput = function()
  {
    let ctx;
    updateInput(center);
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
  start.onkeydown = keydown;
  start.oninput = function()
  {
    updateInput(start);
    showInfo(asmGetInformation(-1, -1));
    updateGraphic(start, 2);
  };
  window.onresize = function()
  {
    let newDomRect = canvas.getBoundingClientRect();
    canvas.width = newDomRect.width;
    canvas.height = newDomRect.height;
    updateGraphic(center, 1);
  };
  animate.onclick = function()
  {
    sincr.value = "0";
    xincr.value = "0";
    yincr.value = "0";
    delay.value = "1";
    applet.style.display = none;
    animform.style.display = block;
  };
  doanimate.onclick = function()
  {
    animform.style.display = none;
    applet.style.display = inline;
    animate.style.display = none;
    stop.style.display = inline;
    incrX = parseFloat(xincr.value);
    incrY = parseFloat(yincr.value);
    incrStart = parseFloat(sincr.value);
    oldX = 0;
    oldY = 0;
    oldStart = 0;
    interval = setInterval(animation, parseFloat(delay.value) * 1000);    
  };
  cancelanim.onclick = function()
  {
    animform.style.display = none;
    applet.style.display = block;
  };
  stop.onclick = function()
  {
    stop.style.display = none;
    animate.style.display = inline;
    clearInterval(interval);
  };
  initMenubarEvents();
}

window.addEventListener("load", startUp);
