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
const exprTextEs = "Por favor ingrese un número o expresión para ";
const exprTextEn = "Please type a number or expression for the ";
let fileContents = 0;
let result, dlog, stop, base, pow, mod, digits, helpbtn, formlink;
let formfeedback, formcancel, formsend, userdata;
let currentInputBox;
let funcnames;
let parens;
let verboseValue;
let prettyValue;
let CunninghamValue;
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
  result.innerHTML = e.substring(1);
  if (e.substring(0, 1) === "2")
  {   // First character passed from web worker is "2".
    dlog.disabled = false;
    stop.disabled = true;
    endCalculation();
  }
}

function saveConfig()
{
  config = "1" +   // Batch mode
           verboseValue +
           prettyValue +
           CunninghamValue +
           (get("hex").checked? "1" : "0") +
           (get("kbd")[1].selected? "1" : "0");
  digits = get("digits").value.trim();
  setStorage("ecmConfig", digits+","+config);
  changeInputmode(get("kbd")[1].selected);
}

function comingFromWorker(e)
{
  fromWorker(e.data);
}

function dowork(n)
{
  const app = n + (config.charAt(4) === "1"? 16: 0);
  const baseText = base.value;
  const powText = pow.value;
  const modText = mod.value;
  hide("help");
  show("result");
  if (baseText === "")
  {
    result.innerHTML = get("missingBase").textContent;
    return;
  }
  if (powText === "")
  {
    result.innerHTML = get("missingPower").textContent;
    return;
  }
  if (modText === "")
  {
    result.innerHTML = get("missingMod").textContent;
    return;
  }
  dlog.disabled = true;
  stop.disabled = false;
  result.innerHTML = get("solving").innerHTML;
  const param = digits + "," + app + "," + baseText + String.fromCharCode(0) + powText +
  String.fromCharCode(0) + modText + String.fromCharCode(0);
  callWorker(param);
}

function getCalcURLs()
{
  return [addLangToFilename("dilogW0000.js"),
          "dilog.webmanifest", "logdi.webmanifest", "dilog-icon-1x.png",
          "dilog-icon-2x.png", "dilog-icon-4x.png", "dilog-icon-180px.png",
          "dilog-icon-512px.png", "favicon.ico"];
}

function getFormSendValue()
{
  userdata.value = get("formatUserData").textContent.replace("{base}", base.value).
          replace("{power}", pow.value).replace("{mod}", mod.value);
}

function popstate(event)
{
  if (get("feedback").style.display === "block" ||
      get("sentOK").style.display === "block" ||
      get("notSent").style.display === "block")
  {         // End feedback.
    show("main");
    hide("feedback");
    hide("sentOK");
    hide("notSent");
    base.focus();   
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
  result = get("result");
  dlog = get("dlog");
  stop = get("stop");
  base = get("base");
  pow = get("pow");
  mod = get("mod");
  helpbtn = get("helpbtn");
  formlink = get("formlink");
  formfeedback = get("formfeedback");
  formcancel = get("formcancel");
  formsend = get("formsend");
  userdata = get("userdata");
  stop.disabled = true;
  dlog.onclick = function()
  {
    dowork(0);
  };
  stop.onclick = function()
  {
    endWorker();
    dlog.disabled = false;
    stop.disabled = true;
    result.innerHTML = get("stopped").innerHTML;
  };
  get("btnSentOK").onclick = function()
  {
    history.back();
  };
  get("btnNotSent").onclick = function()
  {
    history.back();
  };
  helpbtn.onclick = function()
  {
    show("help");
    hide("result");
  };
  base.onfocus = function()
  {
    currentInputBox = base;
  };
  pow.onfocus = function()
  {
    currentInputBox = pow;
  };
  mod.onfocus = function()
  {
    currentInputBox = mod;
  };
  get("funccat").onchange = function()
  {
    generateFuncButtons("funccat", "funcbtns");
  };
  formlink.onclick = clickFormLink;
  get("comments").oninput = function(_event)
  {
    get("formsend").disabled = (get("comments").value === "");
  };
  get("clrinput").onclick = function()
  {
    base.value = "";
    pow.value = "";
    mod.value = "";
    base.focus();
  };
  get("config").onclick = function()
  {
    verboseValue = config.charAt(1);
    prettyValue = config.charAt(2);
    CunninghamValue = config.charAt(3);
    get("digits").value = digits;
    get("hex").checked = (config.charAt(4) === "1");
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
  formcancel.onclick = function()
  {
    history.back();
  };
  formsend.onclick = formSend;
  currentInputBox = base;
  generateFuncButtons("funccat", "funcbtns");
  registerServiceWorker();
}
getCalculatorCode("dilogW0000.js", false);
window.addEventListener("load", startUp);
window.addEventListener("popstate", popstate);
window["fromWorker"] = fromWorker;
