/*
    This file is part of Alpertron Calculators.

    Copyright 2015 Dario Alejandro Alpern

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
/* global get */
/* global getFormSendValue */
/* global hide */
/* global show */
function formSend()
{
  if (get("adduserdata").checked)
  {
    getFormSendValue();
  }
  else
  {
    get("userdata").value = "";
  }
  let xhr = new XMLHttpRequest();
  xhr.onreadystatechange = function(_event)
  {
    if (xhr.readyState === 4) 
    {             // XHR finished.
      if (xhr.status === 200)
      {           // PHP page loaded.
        hide("sending");
        show("sentOK");
      }
      else
      {           // PHP page not loaded.
        hide("sending");
        show("notSent");
      }
    }
  };
  xhr.open("POST", "/" + get("sendmail").textContent + ".php", true);
  xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
  let elements = get("formfeedback").elements;
  let contents = "";
  let useAmp = 0;
  for (let i = 0; i < elements.length; i++)
  {
    let element = elements[i >> 0];
    if (element.type === "radio" && !element.checked)
    {
      continue;
    }
    if (element.name)
    {
      if (useAmp)
      {
        contents += "&";
      }
      contents += element.name + "=" + encodeURIComponent(element.value);
      useAmp++;
    }
  }
  xhr.send(contents);
  hide("feedback");
  show("sending");
  return false;   // Send form only through JavaScript.
}

function clickFormLink(event)
{    
  hide("sentOK");
  hide("notSent");
  hide("main");
  show("feedback");
  get("formfeedback").reset();
  get("formsend").disabled = true;
  get("name").focus();
  history.pushState({id: 1}, "", location.href);
  event.preventDefault();   // Do not follow the link.
}
