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
let coefA;
let coefB;
let coefC;
let coefD;
let coefE;
let coefF;

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
    get("steps").disabled = false;
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

function dowork(n)
{
  let param;
  let app = n;
  let res = get("result");
  let coefAText = coefA.value.trim();
  let coefBText = coefB.value.trim();
  let coefCText = coefC.value.trim();
  let coefDText = coefD.value.trim();
  let coefEText = coefE.value.trim();
  let coefFText = coefF.value.trim();
  let digitGroup = get("digits").value.trim();
  hide("help");
  show("result");
  let missing = "";
  let zero = String.fromCharCode(0);
  if (coefAText === "")
  {
    missing = "<var>a</var>";
  }
  if (coefBText === "")
  {
    missing = "<var>b</var>";
  }
  if (coefCText === "")
  {
    missing = "<var>c</var>";
  }
  if (coefDText === "")
  {
    missing = "<var>d</var>";
  }
  if (coefEText === "")
  {
    missing = "<var>e</var>";
  }
  if (coefFText === "")
  {
    missing = "<var>f</var>";
  }
  if (missing !== "")
  {
    res.innerHTML = get("missing").textContent.replace("{missing}", missing);
    return;
  }
  get("solve").disabled = true;
  get("steps").disabled = true;
  get("stop").disabled = false;
  res.innerHTML = get("solving").textContent;
  param = digitGroup + "," + app + "," + coefAText + zero + coefBText + zero + coefCText + zero +
                                         coefDText + zero + coefEText + zero + coefFText + zero;
  callWorker(param);
}

function moveNext(e, curr, next)
{    
  let nextInput = get(next);
  if (e.key === "Enter" && curr.value.trim().length > 0)
  {
    e.preventDefault();
    nextInput.focus();
    nextInput.setSelectionRange(0, nextInput.value.length);
  }
}

function getCalcURLs()
{
  return [addLangToFilename("quadW0000.js"),
          "quad.webmanifest", "cuad.webmanifest", "quad-icon-1x.png",
          "quad-icon-2x.png", "quad-icon-4x.png", "quad-icon-180px.png",
          "quad-icon-512px.png", "favicon.ico"];
}

function getFormSendValue()
{
  get("userdata").value = "ax^2 + bxy + cy^2 + dx + ey + f = 0" +
                     "\na = " + coefA.value.trim() + "\nb = " + coefB.value.trim() +
                     "\nc = " + coefC.value.trim() + "\nd = " + coefD.value.trim() +
                     "\ne = " + coefE.value.trim() + "\nf = " + coefF.value.trim();  
}

function popstate(event)
{
  if (get("feedback").style.display === "block" ||
      get("sentOK").style.display === "block" ||
      get("notSent").style.display === "block")
  {           // End feedback.
    show("main");
    hide("feedback");
    hide("sentOK");
    hide("notSent");
    coefA.focus();   
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
  coefA = get("coefA");
  coefB = get("coefB");
  coefC = get("coefC");
  coefD = get("coefD");
  coefE = get("coefE");
  coefF = get("coefF");
  get("btnSentOK").onclick = function()
  {
    history.back();
  };
  get("btnNotSent").onclick = function()
  {
    history.back();
  };
  get("stop").disabled = true;
  get("solve").onclick = function()
  {
    dowork(0);
  };
  get("steps").onclick = function()
  {
    dowork(2);
  };
  get("stop").onclick = function()
  {
    endWorker();
    get("solve").disabled = false;
    get("steps").disabled = false;
    get("stop").disabled = true;
    get("result").innerHTML = get("stopped").textContent;
  };
  get("helpbtn").onclick = function()
  {
    show("help");
    hide("result");
  };
  get("clrinput").onclick = function()
  {
    coefA.value = "";
    coefB.value = "";
    coefC.value = "";
    coefD.value = "";
    coefE.value = "";
    coefF.value = "";
    coefA.focus();
  };
  coefA.onkeydown = function(e)
  {
    moveNext(e, this, "coefB");
  };
  coefB.onkeydown = function(e)
  {
    moveNext(e, this, "coefC");
  };
  coefC.onkeydown = function(e)
  {
    moveNext(e, this, "coefD");
  };
  coefD.onkeydown = function(e)
  {
    moveNext(e, this, "coefE");
  };
  coefE.onkeydown = function(e)
  {
    moveNext(e, this, "coefF");
  };
  coefF.onkeydown = function(e)
  {
    if (e.key === "Enter" && e.target.value.trim().length > 0)
    {
      e.preventDefault();
      coefA.focus();
      dowork(0);
    }
  };
  coefA.onfocus = function()
  {
    currentInputBox = coefA;
  };
  coefB.onfocus = function()
  {
    currentInputBox = coefB;
  };
  coefC.onfocus = function()
  {
    currentInputBox = coefC;
  };
  coefD.onfocus = function()
  {
    currentInputBox = coefD;
  };
  coefE.onfocus = function()
  {
    currentInputBox = coefE;
  };
  coefF.onfocus = function()
  {
    currentInputBox = coefF;
  };
  get("funccat").onchange = function()
  {
    generateFuncButtons("funccat", "funcbtns");
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
  currentInputBox = coefA;
  generateFuncButtons("funccat", "funcbtns");
  registerServiceWorker();
}
window.addEventListener("load", startUp);
window.addEventListener("popstate", popstate);
getCalculatorCode("quadW0000.js", false);
window["fromWorker"] = fromWorker;
