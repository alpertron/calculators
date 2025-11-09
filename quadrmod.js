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
/* global callWorker */
/* global changeInputmode */
/* global clickFormLink */
/* global completeFuncButtons */
/* global endCalculation */
/* global endWorker */
/* global formSend */
/* global generateFuncButtons */
/* global get */
/* global getCalculatorCode */
/* global getConfig */
/* global hide */
/* global registerServiceWorker */
/* global setStorage */
/* global show */
let fileContents = 0;
let currentInputBox;
let funcnames;
let parens;
let verboseValue;
let prettyValue;
let CunninghamValue;
let hexValue;
let digits;
let config;

function getFuncNames()
{
  return funcnames;
}

function getParens()
{
  return parens;
}

function fromWorker(e)
{
  // First character of e is "1" for intermediate text
  // and it is "2" for end of calculation.
  get("result").innerHTML = e.substring(1);
  if (e.substring(0, 1) === "2")
  {   // First character passed from web worker is "2".
    get("solve").disabled = false;
    get("stop").disabled = true;
    endCalculation();
  }
}

function saveConfig()
{
  config = "1" +   // Batch mode
           verboseValue +
           prettyValue +
           CunninghamValue +
           hexValue +
           (get("kbd")[1].selected? "1" : "0");
  digits = get("digits").value.trim();
  setStorage("ecmConfig", digits+","+config);
  changeInputmode(get("kbd")[1].selected);
}

function comingFromWorker(e)
{
  fromWorker(e.data);
}

function dowork()
{
  let param;
  let res = get("result");
  let quadrText = get("quad").value.trim();
  let linText = get("lin").value.trim();
  let constText = get("const").value.trim();
  let modText = get("mod").value.trim();
  let digitGroup = digits;
  hide("help");
  show("result");
  let missing = "";
  if (quadrText === "")
  {
    missing = get("missingQuad").textContent;
  }
  if (linText === "")
  {
    missing = get("missingLin").textContent;
  }
  if (constText === "")
  {
    missing = get("missingCons").textContent;
  }
  if (modText === "")
  {
    missing = get("missingMod").textContent;
  }
  if (missing !== "")
  {
    res.innerHTML = get("missing").textContent.replace("{missing}", missing);
    return;
  }
  get("solve").disabled = true;
  get("stop").disabled = false;
  res.innerHTML = get("solving").textContent;
  param = digitGroup + "," + "," + quadrText + String.fromCharCode(0) + linText +
    String.fromCharCode(0) + constText +String.fromCharCode(0) + modText + String.fromCharCode(0);
  callWorker(param);
}

function getCalcURLs()
{
  return [addLangToFilename("quadmodW0000.js"),
          "quadmod.webmanifest", "cuadmod.webmanifest", "quadmod-icon-1x.png",
          "quadmod-icon-2x.png", "quadmod-icon-4x.png", "quadmod-icon-180px.png",
          "quadmod-icon-512px.png", "favicon.ico"];
}

function getFormSendValue()
{
  get("userdata").value = "ax^2 + bx + c = 0 (mod n)" + 
                     "\na = " + get("quad").value + "\nb = " + get("lin").value +
                     "\nc = " + get("const").value + "\nn = " + get("mod").value;
}

function popstate(event)
{
  if (get("feedback").style.display === "block" ||
      get("sentOK").style.display === "block" ||
      get("notSent").style.display === "block")
  {
    show("main");
    hide("feedback");
    hide("sentOK");
    hide("notSent");
    get("quad").focus();   
  }
  else if (get("modal-config").style.display === "block")
  {     // End configuration mode.
    hide("modal-config");
  }
}

function startUp()
{
  funcnames =
  [
    get("btn1").textContent,
    get("btn2").textContent,
    get("btn3").textContent,
    get("btn4").textContent,
    get("btn5").textContent,
    get("btn6").textContent,
    get("btn7").textContent,
    get("btn8").textContent
  ];
  parens = get("parens").textContent;
  get("btnSentOK").onclick = function()
  {
    history.back();
  };
  get("btnNotSent").onclick = function()
  {
    history.back();
  };
  get("stop").disabled = true;
  get("solve").onclick = dowork;
  get("stop").onclick = function()
  {
    endWorker();
    get("solve").disabled = false;
    get("stop").disabled = true;
    get("result").innerHTML = get("stopped").innerHTML;
  };
  get("helpbtn").onclick = function()
  {
    show("help");
    hide("result");
  };
  get("quad").onfocus = function()
  {
    currentInputBox = get("quad");
  };
  get("lin").onfocus = function()
  {
    currentInputBox = get("lin");
  };
  get("const").onfocus = function()
  {
    currentInputBox = get("const");
  };
  get("mod").onfocus = function()
  {
    currentInputBox = get("mod");
  };
  get("funccat").onchange = function()
  {
    generateFuncButtons("funccat", "funcbtns");
  };
  get("clrinput").onclick = function()
  {
    get("quad").value = "";
    get("lin").value = "";
    get("const").value = "";
    get("mod").value = "";
    get("quad").focus();
  };
  get("config").onclick = function()
  {
    verboseValue = config.charAt(1);
    prettyValue = config.charAt(2);
    CunninghamValue = config.charAt(3);
    hexValue = config.charAt(4);
    get("digits").value = digits;
    history.pushState({id: 5}, "", location.href);
    get("kbd")[+config.charAt(5)].selected = "selected";
    show("modal-config");
  };
  get("close-config").onclick = function()
  {
    history.back();   // Close configuration mode.
  };
  get("cancel-config").onclick = function()
  {
    history.back();   // Close configuration mode.
  };
  get("save-config").onclick = function()
  {
    saveConfig();
    history.back();   // Close configuration mode.
  };
  getConfig();
  get("formlink").onclick = clickFormLink;
  get("formcancel").onclick = function()
  {
    history.back();
  };
  get("comments").oninput = function(_event)
  {
    get("formsend").disabled = (get("comments").value === "");
  };
  get("formsend").onclick = formSend;
  currentInputBox = get("quad");
  registerServiceWorker();
  completeFuncButtons("funcbtns");
}

window.addEventListener("load", startUp);
window.addEventListener("popstate", popstate);
getCalculatorCode("quadmodW0000.js", false);
window["fromWorker"] = fromWorker;
