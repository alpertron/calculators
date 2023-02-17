/*
    GAUSSIAN PRIMES
    Copyright(C) 2018 Dario Alejandro Alpern

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Gaussian Primes is distributed in the hope that it will be useful,
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
(function()
{   // This method separates the name space from the Google Analytics code.
const none = "none";
const block = "block";
const inline = "inline";
const asmjs = typeof(WebAssembly) === "undefined";
let buffer, env, asm;
let zoom, zoomDone, imgData;
let canvas, zoomin, zoomout, centerX, centerY;
let isMouseDown;
let currentX, currentY, pixels;
let wasm;
let HEAP8;
let prevX1stTouch, prevY1stTouch;
let prevX2ndTouch, prevY2ndTouch;
let beforeMinus = "";
let asmGetInformation;
let asmDrawPartialGraphic;
let asmMoveGraphic;
let asmNbrChanged;
let bitsCanvas;
let information;
let animform;
let animate;
let stop;
let xincr;
let yincr;
let incrX;
let incrY;
let oldX, oldY;
let delay;
let doanimate;
let cancelanim;
let applet;
let interval;

//##  asmJS goes here (do not change symbols at the left).

function get(x)
{
  return document.getElementById(x);
}

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

function updateGraphic(nbr)
{
  let idx, value, ctr;
  let startOffset = 11000000|0;
  if (nbr === 1)
  {
    ctr = 0;
    value = centerX.value;
    for (idx=0; idx<value.length; idx++, ctr++)
    {
      HEAP8[startOffset+ctr] = value.charCodeAt(idx);
    }
    HEAP8[startOffset+ctr++] = 0;
    value = centerY.value;
    for (idx=0; idx<value.length; idx++, ctr++)
    {
      HEAP8[startOffset+ctr] = value.charCodeAt(idx);
    }
    HEAP8[startOffset+ctr] = 0;
  }
  let width = getWidth();
  let height = getHeight();
  asmNbrChanged(startOffset, nbr, width, height);
  drawGraphic(canvas.getContext("2d"), 0, 0, width, height);
}

function b64decode(str,out)
{
  let ch;
  let idxDest,idxSrc;
  let blocks, leftOver;
  let byte0, byte1, byte2, byte3;
  let conv=new Int8Array(128);
  let len=str.length;
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
  let infoText = ptrToString(ptrInfo).split("^");
  if (infoText[1] !== centerX.value)
  {
    centerX.value = infoText[1];
  } 
  if (infoText[2] !== centerY.value)
  {
    centerY.value = infoText[2];
  } 
  information.innerHTML = infoText[0];
}

function startLowLevelCode()
{
  let length, bytes;
  let info;
  let getPixels;
  let asmJSbuffer;
  imgData = canvas.getContext("2d").createImageData(2048, 4096);  // 32 MB;
  buffer = imgData.data.buffer;          // Reserve 16 MB.
  bitsCanvas = new Uint8Array(buffer);
  if (asmjs)
  {                                      // Asm.js initialization.
    let bufSize = 33554432;
    asmJSbuffer = new ArrayBuffer(bufSize);
    HEAP8 = new Uint8Array(asmJSbuffer);    // Reserve 32 MB for asm.js variables and buffers.
    env = {"a": {"buffer": asmJSbuffer},
      "abort": function(_q)
               {
                 /* Not used*/
               },
    };
    asm = instantiate(env);  // Link asm.js module.
    asmGetInformation = asm.a;
    asmDrawPartialGraphic = asm.a;
    asmMoveGraphic = asm.a;
    asmNbrChanged = asm.a;
    getPixels = asm.a;
    pixels = HEAP8.subarray(getPixels());
    updateGraphic(1);
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
      updateGraphic(1);
    });
  }
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
    updateGraphic(3);
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
    updateGraphic(4);
  }
}

function isNotSpecialKey(event)
{
  let key = event.key;
  let acceptedKeys = ",Backspace,Tab,Right,ArrowRight,Left,ArrowLeft,Cut," +
                     "Control,Meta,Shift,Insert,Delete,Copy,Paste,Home,End,";
  if (event.ctrlKey || event.metaKey)
  {
    if (key === "c")
    {    // User pressed CTRL-C. Map it to "Copy".
      key = "Copy";
    }
    if (key === "v")
    {    // User pressed CTRL-V. Map it to "Paste".
      key = "Paste";
    }
    if (key === "x")
    {    // User pressed CTRL-X. Map it to "Cut".
      key = "Cut";
    }
  }
  return acceptedKeys.indexOf(","+key+",") < 0;
}

function keydown(evt)
{
  let target = evt.target || evt.srcElement;
  let key = evt.key;
  if (isNotSpecialKey(evt))
  {
    if (!evt.ctrlKey && !evt.altKey && !evt.metaKey)
    {                                  // No modifier key pressed.
      if (key >= "0" && key <= "9")
      {                                // Digit key has been pressed.
        if (target.value.length >= 10 ||
            (target.value.charAt(0) !== "-" && target.value.length >= 9))
        {                              // Number is too large.
          evt.preventDefault();        // Do not propagate this key.
        }
      }
      else if (key === "-")
      {                                // Key minus has been pressed.
        if (target.value.indexOf("-") >= 0)
        {                              // There is already a minus sign.
          evt.preventDefault();        // Do not propagate this key.
        }
        else
        {
          beforeMinus = target.value;
        }
      }
      else 
      {                                // Not backspace, tab, right or left arrow, insert or delete key.
        evt.preventDefault();          // Do not propagate this key.  
      }
    }
  }
}
    
function animation()
{
  let currX = oldX - incrX * zoom;
  let currY = oldY + incrY * zoom;
  let diffX = Math.floor(currX) - Math.floor(oldX);
  let diffY = Math.floor(currY) - Math.floor(oldY);
  if (diffX !== 0 || diffY !== 0)
  {
    moveGraphic(diffX, diffY);
    showInfo(asmGetInformation(-1, -1));
  }
  oldX = currX;
  oldY = currY;
}

function startUp()
{  
  canvas = get("canvas");
  zoomin = get("zoomin");
  zoomout = get("zoomout");
  centerX = get("centerX");
  centerY = get("centerY");
  information = get("info");
  animate = get("animate");
  animform = get("animform");
  stop = get("stop");
  xincr = get("xincr");
  yincr = get("yincr");
  delay = get("delay");
  doanimate = get("doanimate");
  cancelanim = get("cancelanim");
  applet = get("applet");
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
  centerX.onkeydown = function(e)
  {
    keydown(e);
  };
  centerX.oninput = function()
  {
    if (beforeMinus !== "" && centerX.value !== "-" + beforeMinus)
    {     // Minus sign is not in the first character. Move it to first character.
      centerX.value = "-" + beforeMinus;
    }
    beforeMinus = "";
    get("info").innerHTML = "";	
    updateGraphic(1);
  };
  centerY.onkeydown = function(e)
  {
    keydown(e);
  };
  centerY.oninput = function()
  {
    if (beforeMinus !== "" && centerY.value !== "-" + beforeMinus)
    {     // Minus sign is not in the first character. Move it to first character.
      centerY.value = "-" + beforeMinus;
    }
    beforeMinus = "";
    get("info").innerHTML = "";
    updateGraphic(1);
  };
  window.onresize = function()
  {
    let newDomRect = canvas.getBoundingClientRect();
    canvas.width = newDomRect.width;
    canvas.height = newDomRect.height;
    updateGraphic(1);
  };
  animate.onclick = function()
  {
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
    oldX = 0;
    oldY = 0;
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
  let menuItems = document.querySelectorAll("[role=\"menubar\"] > li");
  Array.prototype.forEach.call(menuItems, function(el, i)
  {
    el.addEventListener("click", function(event)
    {
      if (this.getAttribute("aria-expanded") === "false")
      {
        this.setAttribute("aria-expanded", "true");
      }
      else
      {
        this.setAttribute("aria-expanded", "false");
      }
      this.firstElementChild.firstElementChild.firstElementChild.focus();
      event.preventDefault();
      return false;
    });
    
    el.addEventListener("keydown", function(event)
    {
      let nextNode;
      if (event.key === "Enter")
      {
        this.click(event);
        return;
      }
      if (event.key === "ArrowRight")
      {
        nextNode = this.nextElementSibling;
        if (nextNode === null)
        {
          nextNode = this.parentNode.firstElementChild;
        }
        nextNode.focus();
        event.preventDefault();
        return;
      }
      if (event.key === "ArrowLeft")
      {
        nextNode = this.previousElementSibling;
        if (nextNode === null)
        {
          nextNode = this.parentNode.lastElementChild;
        }
        nextNode.focus();
        event.preventDefault();
        return;
      }
      if (event.key === "ArrowUp")
      {
        this.setAttribute("aria-expanded", "true");
        this.firstElementChild.lastElementChild.firstElementChild.focus();
        event.preventDefault();
        return;
      }
      if (event.key === "ArrowDown")
      { 
        this.setAttribute("aria-expanded", "true");
        this.firstElementChild.firstElementChild.firstElementChild.focus();
        event.preventDefault();
        return;
      }
    });
    
    el.addEventListener("mouseover", function(event)
    {
      el.addEventListener("mouseover", function(event)
      {
        this.setAttribute("aria-expanded", "true");
      });
      el.addEventListener("mouseout", function(event)
      {
        this.setAttribute("aria-expanded", "false");
      });
    });
    
    let submenuItems = el.querySelectorAll("a");
    Array.prototype.forEach.call(submenuItems, function(el, i)
    {
      el.tabIndex = -1;
      el.addEventListener("click", function(event)
      {
        let parent = this.parentNode.parentNode.parentNode;
        parent.setAttribute("aria-expanded", "false");
        window.location = this.getAttribute("href");
        event.stopImmediatePropagation();
        event.preventDefault();
      });
      el.addEventListener("keydown", function(event)
      {
        let next;
        let parent = this.parentNode.parentNode.parentNode;
        if (event.key === "Tab")
        {
          parent.setAttribute("aria-expanded", "false");
          return;
        }
        if (event.key === "Escape")
        {
          parent.setAttribute("aria-expanded", "false");
          parent.focus();
          event.preventDefault();          
          return;
        }
        if (event.key === "Enter")
        {
          parent.setAttribute("aria-expanded", "false");
          window.location = this.getAttribute("href");
          event.stopImmediatePropagation();
          event.preventDefault();
          return;
        }
        if (event.key === "ArrowRight")
        {
          parent.setAttribute("aria-expanded", "false");
          next = parent.nextElementSibling;
          if (next === null)
          {
            next = parent.parentNode.firstElementChild;
          }
          next.setAttribute("aria-expanded", "true");
          next.firstElementChild.firstElementChild.firstElementChild.focus();
          event.stopImmediatePropagation();
          event.preventDefault();
          return;
        }
        if (event.key === "ArrowLeft")
        {
          parent.setAttribute("aria-expanded", "false");
          next = parent.previousElementSibling;
          if (next === null)
          {
            next = parent.parentNode.lastElementChild;
          }
          next.setAttribute("aria-expanded", "true");
          next.firstElementChild.firstElementChild.firstElementChild.focus();
          event.stopImmediatePropagation();
          event.preventDefault();
          return;
        }
        if (event.key === "ArrowUp" || event.key === "ArrowDown")
        {
          if (event.key === "ArrowUp")
          {
            next = this.parentNode.previousElementSibling;
          }
          else
          {
            next = this.parentNode.nextElementSibling;
          }
          if (next === null)
          {
            parent.setAttribute("aria-expanded", "false");
            parent.focus();
          }
          else
          {
            next.firstElementChild.focus();
          }
          event.stopImmediatePropagation();
          event.preventDefault();
          return;
        }
      });
    });
  });
}

window.addEventListener("load", startUp);
})();
