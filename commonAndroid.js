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
/* global fileContents */
/* global get */
/* global processDKey */
function setStorage(name, data)
{
  Android.setStorage(name, data);
}

function getStorage(name)
{
  return Android.getStorage(name);
}

function getCalculatorCode(fileName, workerParam)
{
  fileContents = 1;
}

function registerServiceWorker()
{
  get("appver").value = Android.getAppVer();
}

function endCalculation()
{
  Android.endCalculation();
}

function getFocusableItems()
{
  let focusableItems = [];
  let buttons = document.getElementsByTagName('button');
  for (let i=0; i<buttons.length; i++)
  {
    if (buttons[i].style.visibility && !buttons[i].disabled)
    {
      focusableItems.push(buttons[i]);
    }
  }
  let inputs = document.getElementsByTagName('input');
  for (let i=0; i<inputs.length; i++)
  {
    if (inputs[i].style.visibility && !buttons[i].disabled)
    {
      focusableItems.push(inputs[i]);
    }
  }
  return focusableItems;
}

function processKey(keyCode)
{
  let focusableItems = getFocusableItems();
  let focusController = this.getFocusController();
  if (keyCode == 19)
  {           // DPAD up
    focusController.moveFocus({x: 0, y: 1});
  }
  else if (keyCode == 20)
  {           // DPAD down
    focusController.moveFocus({x: 0, y: -1});
  }
  else if (keyCode == 21)
  {           // DPAD left
    focusController.moveFocus({x: -1, y: 0});
  }
  else if (keyCode == 22)
  {           // DPAD right
    focusController.moveFocus({x: 1, y: 0});
  }
  else if (keyCode == 23)
  {          // DPAD center
    if (focusController.getCurrentlyFocusedItem())
    {
      focusController.getCurrentlyFocusedItem().onItemClick();
    }
  }
}

window.processDKey = processKey;