/*
    This file is part of Alpertron Calculators.

    Copyright 2023 Dario Alejandro Alpern

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
/* global Android */
/* global asmGetInformation */
/* global asmDrawPartialGraphic */
/* global asmMoveGraphic */
/* global asmNbrChanged */
/* global bitsCanvas */
/* global buffer */
/* global get */
/* global getPixels */
/* global imgData */
/* global pixels */
/* global updateGraphic */
let startOffset;
function startLowLevelCode(type)
{
  let bufSize = 33554432;
  startOffset = 0;
  imgData = canvas.getContext("2d").createImageData(2048, 4096);  // 32 MB;
  buffer = imgData.data.buffer;
  bitsCanvas = new Uint8Array(buffer);
  asmGetInformation = function(x, y)
  {
    return Android.getInformation(x, y);
  }
  asmDrawPartialGraphic = function(xminDisp, xmaxDisp, yminDisp, ymaxDisp)
  {
    Android.drawPartialGraphic(xminDisp, xmaxDisp, yminDisp, ymaxDisp);
  }
  asmMoveGraphic = function(deltaX, deltaY)
  {
    Android.moveGraphic(deltaX, deltaY);
  }
  asmNbrChanged = function(value, inputBoxNbr, newWidth, newHeight)
  {
    return Android.nbrChanged(value, inputBoxNbr, newWidth, newHeight);
  }
  updateGraphic(center, 1);
}

function updateImageData(top, left, width, height)
{
  let rowNbr;
  let leftPixelSrc = top*8192+left*4;
  let leftPixelDest = (top+32)*8192+left*4;
  let rowSize = width * 4;
  let base64 = Android.getPixelArrayData(leftPixelSrc, rowSize, height);
  let base64RowSize = 4*Math.ceil(rowSize/3);
  let uint8Array = new Uint8Array(rowSize);
  for (rowNbr=0; rowNbr<height; rowNbr++)
  {
    let binaryString = atob(base64.substring(base64RowSize*rowNbr, base64RowSize*(rowNbr+1)));
    for (let i = 0; i < binaryString.length; i++)
    {
      uint8Array[i] = binaryString.charCodeAt(i);
    }
    bitsCanvas.set(uint8Array.subarray(0, 4*width), leftPixelDest);
    leftPixelSrc += 8192;
    leftPixelDest += 8192;
  }
}

function copyStr(dummy, str)
{
  startOffset = str;
}

function drawPartialGraphic(xminDisp, xmaxDisp, yminDisp, ymaxDisp)
{
  asmDrawPartialGraphic(xminDisp, xmaxDisp, yminDisp, ymaxDisp);
}

function setNewDimensionsForCanvas()
{
  let width = window.innerWidth;
  let height = window.innerHeight - get("bottom").clientHeight - 50;
  canvas.width = width;
  canvas.height = height;
}
