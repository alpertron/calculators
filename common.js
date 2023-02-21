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
/* global asmjs */
/* global endFeedback */
/* global fileContents */
/* global get */
/* global lang */

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

function initMenubarEvents()
{  
  let menuItems = document.querySelectorAll("[role=\"menubar\"] > li");
  Array.prototype.forEach.call(menuItems, function(el, i)
  {
    el.addEventListener("click", function(event)
    {
      if (event.target.getAttribute("aria-expanded") === "false")
      {
        event.target.setAttribute("aria-expanded", "true");
      }
      else
      {
        event.target.setAttribute("aria-expanded", "false");
      }
      event.target.firstElementChild.firstElementChild.firstElementChild.focus();
      event.preventDefault();
      return false;
    });
    
    el.addEventListener("keydown", function(event)
    {
      let nextNode;
      if (event.key === "Enter")
      {
        event.target.click(event);
        return;
      }
      if (event.key === "ArrowRight")
      {
        nextNode = event.target.nextElementSibling;
        if (nextNode === null)
        {
          nextNode = event.target.parentNode.firstElementChild;
        }
        nextNode.focus();
        event.preventDefault();
        return;
      }
      if (event.key === "ArrowLeft")
      {
        nextNode = event.target.previousElementSibling;
        if (nextNode === null)
        {
          nextNode = event.target.parentNode.lastElementChild;
        }
        nextNode.focus();
        event.preventDefault();
        return;
      }
      if (event.key === "ArrowUp")
      {
        event.target.setAttribute("aria-expanded", "true");
        event.target.firstElementChild.lastElementChild.firstElementChild.focus();
        event.preventDefault();
        return;
      }
      if (event.key === "ArrowDown")
      { 
        event.target.setAttribute("aria-expanded", "true");
        event.target.firstElementChild.firstElementChild.firstElementChild.focus();
        event.preventDefault();
      }
    });
    
    el.addEventListener("mouseenter", function(event)
    {
      event.target.setAttribute("aria-expanded", "true");
    });
   
    el.addEventListener("mouseleave", function(event)
    {
      event.target.setAttribute("aria-expanded", "false");
    });
    
    let submenuItems = el.querySelectorAll("a");
    Array.prototype.forEach.call(submenuItems, function(el, i)
    {
      el.tabIndex = -1;
      el.addEventListener("click", function(event)
      {
        let parent = event.target.parentNode.parentNode.parentNode;
        parent.setAttribute("aria-expanded", "false");
        window.location = event.target.getAttribute("href");
        event.stopImmediatePropagation();
        event.preventDefault();
      });
      el.addEventListener("keydown", function(event)
      {
        let next;
        let parent = event.target.parentNode.parentNode.parentNode;
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
          window.location = event.target.getAttribute("href");
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
            next = event.target.parentNode.previousElementSibling;
          }
          else
          {
            next = event.target.parentNode.nextElementSibling;
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
        }
      });
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
    byte0 = conv[str.charCodeAt(idxSrc)];
    byte1 = conv[str.charCodeAt(idxSrc+1)];
    byte2 = conv[str.charCodeAt(idxSrc+2)];
    byte3 = conv[str.charCodeAt(idxSrc+3)];
    
    out[idxDest >>0 ] = (byte0<<2) + (byte1>>4);
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
