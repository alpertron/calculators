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

/* global cantShowCanvas */
/* global commonGraphicEvents */
/* global drawGraphic */
/* global get */
/* global getHeight */
/* global getWidth */
/* global initMenubarEvents */
/* global instantiate */
/* global isNotSpecialKey */
/* global moveGraphic */
/* global ptrToString */
/* global startLowLevelCode */
// In order to reduce the number of files to read from Web server, this 
// Javascript file includes both the Javascript in the main thread and the 
// Javascript that drives WebAssembly on its own Web Worker.
const none = "none";
const block = "block";
const inline = "inline";
let buffer, env, asm;
let zoom, zoomDone, imgData;
let canvas, zoomin, zoomout, center, start;
let isMouseDown;
let currentX, currentY, pixels;
let HEAPU8;
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

function updateGraphic(input, nbr)
{
  let startOffset = 11000000|0;
  let idx;
  let value = input.value;
  if (input !== 0)
  {
    for (idx=0; idx<value.length; idx++)  
    {
      HEAPU8[(startOffset+idx) >> 0] = value.charCodeAt(idx);
    }
    HEAPU8[(startOffset+idx) >> 0] = 0;
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

// Perl script will change the fields to match asm.js.
function startLowLevelCodeCallback(asm)
{
  asmGetInformation = asm.e;
  asmDrawPartialGraphic = asm.c;
  asmMoveGraphic = asm.f;
  asmNbrChanged = asm.g;
  let getPixels = asm.a;
  return getPixels;
}

function cantShowCanvas()
{
  return (center.value.length < start.value.length || 
        (center.value.length === start.value.length && center.value < start.value));
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
  commonGraphicEvents();
  center.onkeydown = keydown;
  center.oninput = function()
  {
    let ctx;
    updateInput(center);
    if (cantShowCanvas())
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
