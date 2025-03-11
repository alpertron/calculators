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
/* global cantShowCanvas */
/* global commonGraphicEvents */
/* global copyStr */
/* global drawGraphic */
/* global get */
/* global getHeight */
/* global getWidth */
/* global initMenubarEvents */
/* global instantiate */
/* global isNotSpecialKey */
/* global moveGraphic */
/* global ptrToString */
/* global setNewDimensionsForCanvas */
/* global startLowLevelCode */
/* global startOffset */
// In order to reduce the number of files to read from Web server, this 
// Javascript file includes both the Javascript in the main thread and the 
// Javascript that drives WebAssembly on its own Web Worker.
const none = "none";
const block = "block";
const inline = "inline";
let fileContents = 0;
let center = null;
let buffer, env, asm;
let zoom, zoomDone, imgData;
let /** HTMLCanvasElement */ canvas;
let zoomin, zoomout, centerX, centerY;
let isMouseDown;
let currentX, currentY, pixels;
let HEAPU8;
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

function updateGraphic(dummy, nbr)
{
  if (nbr === 1)
  {
    copyStr(startOffset, centerX.value + "\x01" + centerY.value);
  }
  let width = getWidth();
  let height = getHeight();
  asmNbrChanged(startOffset, nbr, width, height);
  drawGraphic(canvas.getContext("2d"), 0, 0, width, height);
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

// Perl script will change the fields to match asm.js.
function startLowLevelCodeCallback(asm)
{
  asmGetInformation = asm.a;
  asmDrawPartialGraphic = asm.a;
  asmMoveGraphic = asm.a;
  asmNbrChanged = asm.a;
  let getPixels = asm.a;
  return getPixels;
}

function cantShowCanvas()
{
  return false;
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
  commonGraphicEvents();
  centerX.onkeydown = function(evt)
  {
    if (evt.key === "Up" || evt.key === "ArrowUp")
    {
      canvas.focus();
      evt.preventDefault();          // Do not propagate this key.
    }
    keydown(evt);
  };
  centerX.oninput = function()
  {
    if (beforeMinus !== "" && centerX.value !== "-" + beforeMinus)
    {     // Minus sign is not in the first character. Move it to first character.
      centerX.value = "-" + beforeMinus;
    }
    beforeMinus = "";
    get("info").innerHTML = "";	
    updateGraphic(0, 1);
  };
  centerY.onkeydown = function(evt)
  {
    if (evt.key === "Up" || evt.key === "ArrowUp")
    {
      canvas.focus();
      evt.preventDefault();          // Do not propagate this key.
    }
    keydown(evt);
  };
  centerY.oninput = function()
  {
    if (beforeMinus !== "" && centerY.value !== "-" + beforeMinus)
    {     // Minus sign is not in the first character. Move it to first character.
      centerY.value = "-" + beforeMinus;
    }
    beforeMinus = "";
    get("info").innerHTML = "";
    updateGraphic(0, 1);
  };
  window.onresize = function()
  {
    setNewDimensionsForCanvas();
    updateGraphic(0, 1);
  };
  animate.onclick = function()
  {
    xincr.value = "0";
    yincr.value = "0";
    delay.value = "1";
    applet.style.display = none;
    animform.style.display = block;
    xincr.focus();
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
    setNewDimensionsForCanvas();
    updateGraphic(0, 1);
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
  if (get("left") == null)
  {
    setNewDimensionsForCanvas();
    updateGraphic(0, 1);
  }
  xincr.onkeydown = function(evt)
  {
    if (evt.key === "Down" || evt.key === "ArrowDown")
    {
      yincr.focus();
      evt.preventDefault();          // Do not propagate this key.
    }
  };
  yincr.onkeydown = function(evt)
  {
    if (evt.key === "Up" || evt.key === "ArrowUp")
    {
      xincr.focus();
      evt.preventDefault();          // Do not propagate this key.
    }
    if (evt.key === "Down" || evt.key === "ArrowDown")
    {
      delay.focus();
      evt.preventDefault();          // Do not propagate this key.
    }
  };
  delay.onkeydown = function(evt)
  {
    if (evt.key === "Up" || evt.key === "ArrowUp")
    {
      yincr.focus();
      evt.preventDefault();          // Do not propagate this key.
    }
    if (evt.key === "Down" || evt.key === "ArrowDown")
    {
      doanimate.focus();
      evt.preventDefault();          // Do not propagate this key.
    }
  };
}

window.addEventListener("load", startUp);
