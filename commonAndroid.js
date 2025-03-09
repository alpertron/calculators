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
let elementInFocus = null;
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

function isVisible(element)
{
  return element !== null && element.offsetParent !== null &&
         !element.disabled && !element.hidden;
}

function getNextVisibleElement(discardElem)
{
  const elements = document.querySelectorAll("textarea, input, button, p");
  let prevElement = null;
  let firstElement = null;
  let activeElementFound = false;
  let elementsLength = elements.length;
  for (let index=0; index<elementsLength; index++)
  {
    let element = elements[index];
    if (element === elementInFocus)
    {
      if (discardElem)
      {
        return prevElement;
      }
      activeElementFound = true;
    }
    if (isVisible(element))
    {
      if (firstElement == null)
      {
        firstElement = element;
      }
      if (activeElementFound)
      {
        return element;
      }
      prevElement = element;
    }
  }
  if (activeElementFound)
  {
    return prevElement;
  }
  return firstElement;
}

function showFocus(elem)
{
  elementInFocus = elem;
  if (elem == null)
  {   // No element can be focused.
    document.body.focus();
  }
  else
  {
    elem.setAttribute("tabindex", "0");
    elem.setAttribute("readonly", true);
    elem.focus();
    setTimeout(() =>
    {
      elem.removeAttribute("readonly");  // Remove readonly after focusing
    }, 100);
  }
}

function newFocus(event)
{
  if (event.target !== document.body)
  {
    elementInFocus = event.target;
  }
}

document.addEventListener("keydown", function(event)
{
  if (event.keyCode >= 37 && event.keyCode <= 40)
  {    // Directional keys.
    if (!document.hasFocus() || !isVisible(elementInFocus))
    {
      showFocus(getNextVisibleElement(true));
    }
    else if (document.activeElement === document.body)
    {
      showFocus(elementInFocus);
    }
    else
    {
      elementInFocus = document.activeElement;
    }
  }
  if (event.keyCode === 13)
  {    // Enter
    newFocus(event);
  }
});

window.addEventListener("load", function(event)
{
  const elements = document.querySelectorAll("textarea, input, button, p");
  let elementsLength = elements.length;
  for (let index=0; index<elementsLength; index++)
  {
    let element = elements[index];
    element.addEventListener("click", newFocus);
  }
});

function onShowDivisors()
{
  elementInFocus = get("showDiv");
  showFocus(getNextVisibleElement(true));
}

function onShowSumSquares()
{
  elementInFocus = get("showSumSq");
  showFocus(getNextVisibleElement(true));
}

function setFocusTo(newFocusObj)
{
  elementInFocus = newFocusObj;
}
