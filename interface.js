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
/* global clearWizardTextInput */
/* global clickFormLink */
/* global endCalculation */
/* global endWorker */
/* global get */
/* global getConfig */
/* global formSend */
/* global generateFuncButtons */
/* global getCalculatorCode */
/* global hide */
/* global keyDownOnWizard */
/* global registerServiceWorker */
/* global selectLoop */
/* global setStorage */
/* global setWizardStep */
/* global show */
/* global typedOnWizard */
/* global wizardNext */
/** @define {number} */ const app = 0;   // Use with Closure compiler.
let fileContents = 0;
let currentInputBox;
let funcnames;
let parens;
let value;
let verboseValue;
let prettyValue;
let CunninghamValue;
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

function styleButtons(style1, style2)
{
  get("calc").style.display = style1;
  if (get("calc5") !== null)
  {
    get("calc5").style.display = style1;
  }
  if (get("calcN") !== null)
  {
    get("calcN").style.display = style1;
  }
  if ((app !== 4) && (app !== 5))
  {    // Continued fraction applet does not use wizard.
    get("openwizard").style.display = style1;
  }
  if (get("stop") !== null)
  {
    get("stop").style.display = style2;
  }
}

let calcURLs;

if (app < 2)
{
  calcURLs = [addLangToFilename("fsquaresW0000.js"),
              "fsquares.webmanifest", "sumcuad.webmanifest", "fsquares-icon-1x.png",
              "fsquares-icon-2x.png", "fsquares-icon-4x.png", "fsquares-icon-180px.png",
              "fsquares-icon-512px.png", "favicon.ico"];
}
else if (app < 4)
{
  calcURLs = [addLangToFilename("fcubesW0000.js"),
              "fcubes.webmanifest", "sumcubos.webmanifest", "fcubes-icon-1x.png",
              "fcubes-icon-2x.png", "fcubes-icon-4x.png", "fcubes-icon-180px.png",
              "fcubes-icon-512px.png", "favicon.ico"];
}
else if (app < 6)
{
  calcURLs = [addLangToFilename("contfracW0000.js"),
              "contfrac.webmanifest", "fraccont.webmanifest", "contfrac-icon-1x.png",
              "contfrac-icon-2x.png", "contfrac-icon-4x.png", "contfrac-icon-180px.png",
              "contfrac-icon-512px.png", "favicon.ico"];
}
else
{
  calcURLs = [addLangToFilename("tsqcubesW0000.js"),
              "tsqcubes.webmanifest", "tcuadcub.webmanifest", "tsqcubes-icon-1x.png",
              "tsqcubes-icon-2x.png", "tsqcubes-icon-4x.png", "tsqcubes-icon-180px.png",
              "tsqcubes-icon-512px.png", "favicon.ico"];
}

function getCalcURLs()
{
  return calcURLs;
}

function fromWorker(e)
{
  // First character of e is:
  // "1" for intermediate output
  // "2" for end calculation
  // "4" for sending data to status line
  // "6" for pausing calculation and showing the Continue button
  let firstChar = e.substring(0, 1);
  if (firstChar === "4")
  {
    get("status").innerHTML = e.substring(1);
  }
  else if (firstChar === "9")
  {
    console.log(e.substring(1)+ "\n");
  }
  else
  {
    get("result").innerHTML = e.substring(1);
    if (firstChar === "2" || firstChar === "6")
    {   // First character passed from web worker is "2" or "6".
      get("status").innerHTML = "";
      styleButtons("inline", "none");  // Enable buttons that must be enabled when applet is not running
      if (firstChar === "6")
      {
        show("cont");
      }
      else
      {
        endCalculation();
      }
    }
  }
}

function comingFromWorker(e)
{
  fromWorker(e.data);
}

function saveConfig(fromWizard)
{
  if (fromWizard)
  {
    get("hex").checked = get("hexW").checked;
  }
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

function performCalc(from)
{
  let valueA, valueB, valueC, digitGroup;
  let res = get("result");
  let id;
  show("result");
  valueA = get("num").value;
  if (valueA === "")
  {
    if ((app === 4) || (app === 5))
    {
      res.innerHTML = get("missingNum").textContent;
    }
    else
    {
      res.innerHTML = get("missing").textContent;
    }
    return;
  }
  if ((app === 4) || (app === 5))
  {
    valueB = get("delta").value;
    if (valueB === "")
    {
      res.innerHTML = get("missingSqr").textContent;
      return;
    }
    valueC = get("den").value;
    if (valueC === "")
    {
      res.innerHTML = get("missingDen").textContent;
      return;
    }
  }
  digitGroup = get("digits").value;
  hide("help");
  if (app === 0)    // Closure compiler cannot optimize switch, so a series of "if" instructions is used.
  {
    res.innerHTML = "Computing sum of squares...";
  }
  else if (app === 1)
  {
    res.innerHTML = "Calculando suma de cuadrados...";
  }
  else if (app === 4)
  {
    res.innerHTML = "Computing continued fraction expansion...";
  }
  else if (app === 5)
  {
    res.innerHTML = "Calculando desarrollo en fracciones continuas...";
  }
  res.scrollIntoView({behavior: "smooth"});
  let param = "";
  if ((app === 6) || (app === 7))
  {         // Sum of two squares and a power.
    param = from + ",";
  }
  let options = app + ((config.charAt(4) === "1")? 64: 0);
  if ((app === 4) || (app === 5))
  {         // Continued fractions.
    options += get("converg").checked? 32: 0;
  }
  param += digitGroup + "," + options + "," + valueA + String.fromCharCode(0);
  if ((app === 4) || (app === 5))
  {         // Continued fractions.
    param += valueB + String.fromCharCode(0) + valueC + String.fromCharCode(0);
  }
  else
  {
    styleButtons("none", "inline");  // Enable "stop" button
  }
  hide("cont");
  let helphelp = get("helphelp");
  const version = (typeof(WebAssembly) === "undefined"? "nowebassy": "webassy");
  helphelp.innerHTML = "<p>" + get("firstLine").innerHTML + " " +
                       get(version).innerHTML + "</p>";
  callWorker(param);
  show("helphelp");
  helphelp.innerHTML = get("firstLine").innerHTML;
}

function oneexpr()
{
  get("next").value = get("done").textContent;
  get("wzddesc").innerHTML = get("wzddesc1").textContent;
  get("wzdexam").innerHTML = "&nbsp;";
  clearWizardTextInput();
  setWizardStep(9);
}

function getFormSendValue()
{
  let userdata = get("userdata");
  if ((app !== 4) && (app !== 5))
  {     // Not continued fraction application.
    userdata.value = "\n" + get("num").value + "\n" + get("result").innerHTML + "\n" + get("status").innerHTML;
  }
  else
  {     // Continued fraction application.
    userdata.value = "\nnum = " + get("num").value + "\ndelta = " + get("delta").value + "\nden = " + get("den").value;         
  }
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
    get("num").focus();
  }
  else if (get("wizard") != null &&
           get("wizard").style.display === "block")
  {         // End wizard.
    show("main");
    hide("wizard");
    get("num").focus();
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
  let param;
  value = get("num");
  if ((app !== 4) && (app !== 5))
  {    // Not continued fraction.
    funcnames.push(get("btn9").textContent); 
    get("num").onkeydown = function(e)
    {
      let digitGroup = get("digits").value;
      let res = get("result");
      show("result");
      let input = get("num").value;
      let keyCode = e.key;
      if (keyCode === "Enter")
      {  // Used pressed Enter key
        if (input === "")
        {
          res.innerHTML = get("missing").innerHTML;
          return;
        }
        if (app === 0)
        {
          res.innerHTML = get("computing").textContent;
        }
        param = digitGroup + "," + app + "," + input + String.fromCharCode(0);
        styleButtons("none", "inline");  // Enable "stop" button
        callWorker(param);
      }
    };
  }
  get("btnSentOK").onclick = function()
  {
    history.back();
  };
  get("btnNotSent").onclick = function()
  {
    history.back();
  };
  get("calc").onclick = function()
  {
    performCalc(3);
  };
  if (get("calc5") !== null)
  {
    get("calc5").onclick = function()
    {
      performCalc(5);
    };
  }
  if (get("calcN") !== null)
  {
    get("calcN").onclick = function()
    {
      let expon = get("expon").value;
      if (!/^\d+$/.test(expon) || parseInt(expon, 10)%2 === 0)
      {
        get("result").innerHTML = get("expodd").innerHTML;
      }
      else
      {
        performCalc(expon);
      }
    };
  }
  if ((app !== 4) && (app !== 5))
  {    // Continued fraction applet does not use wizard.
    get("openwizard").onclick = function()
    {
      get("exprwiz").innerHTML = get("expr").innerHTML;
      hide("main");
      show("wizard");
      get("wzdupper").style.display = "grid";
      get("oneexpr").checked = true;
      get("next").disabled = true;
      get("decW").checked = !get("hexW").checked;
      get("wzdinput").value = "";
      get("wzdinput").focus();
      history.pushState({id: 2}, "", location.href);
      oneexpr();
    };
    get("wzdinput").onkeydown = keyDownOnWizard;
    get("oneexpr").onclick = function()
    {
      oneexpr();
    };
    get("loop").onclick = function()
    {
      selectLoop();
    };
    get("next").onclick = wizardNext;
    get("wzdinput").oninput = typedOnWizard;
    get("cancel").onclick = function()
    {
      history.back();
    };
  }
  if (get("stop") !== null)
  {
    get("stop").onclick = function()
    {
      endWorker();
      styleButtons("inline", "none");  // Enable buttons that have to be enabled when applet is not running.
      get("result").innerHTML = get("stopped").innerHTML;
      get("status").innerHTML = "";
    };
  }
  get("continue").onclick = function()
  {
    hide("cont");
    callWorker("C");  // Indicate worker that user pressed Continue button.
  };
  get("num").onkeydown = function(event)
  {
    let keyCode = event.key;
    if (keyCode === "Enter" && event.ctrlKey)
    {
      event.preventDefault();          // Do not propagate Enter key.
      performCalc(0);                  // Perform calculation.
    }
    return true;
  };
  get("helpbtn").onclick = function()
  {
    if (get("help").style.display === "block" && get("result").innerHTML !== "")
    {
      hide("help");
      show("helphelp");
      show("result");
    }
    else
    {
      show("help");
      hide("helphelp");
      hide("result");
    }
  };
  get("formlink").onclick = clickFormLink;
  get("formcancel").onclick = function()
  {
    history.back();
  };
  get("funccat").onchange = function()
  {
    generateFuncButtons("funccat", "funcbtns");
  };
  if (get("wzdfunccat"))
  {
    get("wzdfunccat").onchange = function()
    {
      generateFuncButtons("wzdfunccat", "wzdfuncbtns");
    };
  }
  get("clrinput").onclick = function()
  {
    if ((app === 4) || (app === 5))
    {
      get("delta").value = "";
      get("den").value = "";
    }
    get("num").value = "";
    get("num").focus();
  };
  get("num").onfocus = function()
  {
    currentInputBox = get("num");
  };
  if (get("wzdinput"))
  {
    get("wzdinput").onfocus = function()
    {
      currentInputBox = get("wzdinput");
    };
  }
  else
  {
    get("delta").onfocus = function()
    {
      currentInputBox = get("delta");
    };
    get("den").onfocus = function()
    {
      currentInputBox = get("den");
    };
  }
  get("comments").oninput = function(_event)
  {
    get("formsend").disabled = (get("comments").value === "");
  };
  get("formsend").onclick = formSend;
  currentInputBox = get("num");
  generateFuncButtons("funccat", "funcbtns");
  if (get("wzdfunccat"))
  {
    generateFuncButtons("wzdfunccat", "wzdfuncbtns");
  }
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
    saveConfig(false);
    history.back();   // Close configuration mode.
  };
  getConfig();
  registerServiceWorker();
}
if (app < 2)
{
  getCalculatorCode("fsquaresW0000.js", false);
}
else if (app < 4)
{
  getCalculatorCode("fcubesW0000.js", false);
}
else if (app < 6)
{
  getCalculatorCode("contfracW0000.js", false);
}
else
{
  getCalculatorCode("tsqcubesW0000.js", false);
}
window.addEventListener("load", startUp);
window.addEventListener("popstate", popstate);
window["fromWorker"] = fromWorker;
