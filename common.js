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
function get(id)
{
  return document.getElementById(id);
}

function hide(id)
{
  get(id).style.display = "none";
}

function show(id)
{
  get(id).style.display = "block";
}

function topMenuClick(event)
{
  if (event.target.getAttribute("aria-expanded") === "false")
  {
    event.target.setAttribute("aria-expanded", "true");
  }
  else
  {
    event.target.setAttribute("aria-expanded", "false");
  }
  // The target is <span>. The next sibling is <ul>.
  event.target.nextElementSibling.firstElementChild.firstElementChild.focus();
  event.preventDefault();
  return false;
}

function topMenuKeyDown(event)
{
  let currentMenu = event.target;
  let currentLi = currentMenu.parentNode;
  let subMenu = currentMenu.nextElementSibling;
  if (event.key === "Enter")
  {
    currentMenu.click(event);
    return;
  }
  if (event.key === "ArrowRight")
  {
    let nextLi = currentLi.nextElementSibling;
    if (nextLi === null)
    {   // Current element is the last one.
      nextLi = currentLi.parentNode.firstElementChild;
    }
    nextLi.firstElementChild.focus();
    event.preventDefault();
    return;
  }
  if (event.key === "ArrowLeft")
  {
    let prevLi = currentLi.previousElementSibling;
    if (prevLi === null)
    {   // Current element is the first one.
      prevLi = currentLi.parentNode.lastElementChild;
    }
    prevLi.firstElementChild.focus();
    event.preventDefault();
    return;
  }
  if (event.key === "ArrowUp")
  {
    currentMenu.setAttribute("aria-expanded", "true");
    subMenu.lastElementChild.firstElementChild.focus();
    event.preventDefault();
    return;
  }
  if (event.key === "ArrowDown")
  { 
    currentMenu.setAttribute("aria-expanded", "true");
    subMenu.firstElementChild.firstElementChild.focus();
    event.preventDefault();
  }
}

function topMenuMouseEnter(event)
{   
  // The target is <li>, but the attribute is in the <span> tag below it.
  event.target.firstElementChild.setAttribute("aria-expanded", "true");
}

function topMenuMouseLeave(event)
{
  // The target is <li>, but the attribute is in the <span> tag below it.
  event.target.firstElementChild.setAttribute("aria-expanded", "false");
}

function subMenuClick(event)
{
  let parent = event.target.parentNode.parentNode.parentNode;
  parent.setAttribute("aria-expanded", "false");
  window.location = event.target.getAttribute("href");
  event.stopImmediatePropagation();
  event.preventDefault();
}

function subMenuKeyDown(event)
{
  let currentLi = event.target.parentNode.parentNode.parentNode;
  let currentMenu = currentLi.firstElementChild;
  if (event.key === "Tab")
  {
    currentMenu.setAttribute("aria-expanded", "false");
    return;
  }
  if (event.key === "Escape")
  {
    currentMenu.setAttribute("aria-expanded", "false");
    currentMenu.focus();
    event.preventDefault();          
    return;
  }
  if (event.key === "Enter")
  {
    currentMenu.setAttribute("aria-expanded", "false");
    window.location = event.target.getAttribute("href");
    event.stopImmediatePropagation();
    event.preventDefault();
    return;
  }
  if (event.key === "ArrowRight")
  {
    currentMenu.setAttribute("aria-expanded", "false");
    let nextLi = currentLi.nextElementSibling;
    if (nextLi === null)
    {
      nextLi = currentLi.parentNode.firstElementChild;
    }
    currentMenu = nextLi.firstElementChild;
    currentMenu.setAttribute("aria-expanded", "true");
    currentMenu.nextElementSibling.firstElementChild.firstElementChild.focus();
    event.stopImmediatePropagation();
    event.preventDefault();
    return;
  }
  if (event.key === "ArrowLeft")
  {
    currentMenu.setAttribute("aria-expanded", "false");
    let prevLi = currentLi.previousElementSibling;
    if (prevLi === null)
    {
      prevLi = currentLi.parentNode.lastElementChild;
    }
    currentMenu = prevLi.firstElementChild;
    currentMenu.setAttribute("aria-expanded", "true");
    currentMenu.nextElementSibling.firstElementChild.firstElementChild.focus();
    event.stopImmediatePropagation();
    event.preventDefault();
    return;
  }
  if (event.key === "ArrowUp" || event.key === "ArrowDown")
  {
    let nextMenuItem;
    if (event.key === "ArrowUp")
    {
      nextMenuItem = event.target.parentNode.previousElementSibling;
    }
    else
    {
      nextMenuItem = event.target.parentNode.nextElementSibling;
    }
    if (nextMenuItem === null)
    {
      currentMenu.setAttribute("aria-expanded", "false");
      currentMenu.focus();
    }
    else
    {
      nextMenuItem.firstElementChild.focus();
    }
    event.stopImmediatePropagation();
    event.preventDefault();
  }
}

function initMenubarEvents()
{  
  let menuItems = document.querySelectorAll("[role=\"menubar\"] > li");
  Array.prototype.forEach.call(menuItems, function(el, i)
  {
    // Capture click and keydown events in the <span> element.
    el.firstElementChild.addEventListener("click", topMenuClick);
    el.firstElementChild.addEventListener("keydown", topMenuKeyDown);
    el.addEventListener("mouseenter", topMenuMouseEnter);
    el.addEventListener("mouseleave", topMenuMouseLeave);
    // Get all <a> elements in the list <li>.
    let submenuItems = el.querySelectorAll("a");
    Array.prototype.forEach.call(submenuItems, function(el, i)
    {
      el.tabIndex = -1;
      el.addEventListener("click", subMenuClick);
      el.addEventListener("keydown", subMenuKeyDown);
    });
  });
}

function b64decode(str,out)
{
  let ch;
  let idxDest,idxSrc;
  let blocks, leftOver;
  let byte0, byte1, byte2, byte3;
  let conv = new Int8Array(128);
  let len = str.length;
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
    byte0 = conv[str.codePointAt(idxSrc)];
    byte1 = conv[str.codePointAt(idxSrc+1)];
    byte2 = conv[str.codePointAt(idxSrc+2)];
    byte3 = conv[str.codePointAt(idxSrc+3)];
    
    out[idxDest >>0 ] = (byte0<<2) + (byte1>>4);
    out[(idxDest+1) >> 0] = (byte1<<4) + (byte2>>2);
    out[(idxDest+2) >> 0] = (byte2<<6) + byte3;
  }
  leftOver = len & 3;
  if (leftOver === 2)
  {
    byte0 = conv[str.codePointAt(idxSrc)];
    byte1 = conv[str.codePointAt(idxSrc+1)];
    
    out[idxDest >> 0] = (byte0<<2) + (byte1>>4);
    out[(idxDest+1) >> 0] = byte1<<4;
  }
  else if (leftOver === 3)
  {
    byte0 = conv[str.codePointAt(idxSrc)];
    byte1 = conv[str.codePointAt(idxSrc+1)];
    byte2 = conv[str.codePointAt(idxSrc+2)];
    
    out[idxDest >> 0] = (byte0<<2) + (byte1>>4);
    out[(idxDest+1) >> 0] = (byte1<<4) + (byte2>>2);
    out[(idxDest+2) >> 0] = byte2<<6;
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

function initGraphicFunctionPointers()
{     // Not used if not Android.
}
