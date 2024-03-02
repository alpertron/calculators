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
/* global clickFormLink */
/* global endCalculation */
/* global endWorker */
/* global formSend */
/* global get */
/* global getCalculatorCode */
/* global hide */
/* global registerServiceWorker */
/* global show */
/** @define {number} */ const lang = 1;   // Use with Closure compiler.
let busy = false;
let workerParam;
let fileContents = 0;

function fromWorker(e)
{
  // First character of e is "1" for intermediate text
  // and it is "2" for end of calculation.
  let firstChar = e.substring(0, 1);
  if ((firstChar === "M") || (firstChar === "N"))
  {    // User entered a number. Load calculator to process it.
    window.sessionStorage.setItem((firstChar === "M"? "F": "E"),
      get("poly").value);
    window.location.replace(lang? "ECMC.HTM": "ECM.HTM");
    return;
  }
  let result = get("result");
  result.innerHTML = e.substring(1);
  if (e.substring(0, 1) === "2")
  {   // First character passed from web worker is "2".
    get("eval").disabled = false;
    get("factor").disabled = false;
    get("stop").disabled = true;
    busy = false;
    result.setAttribute("aria-live", "polite");
    endCalculation();
  }
  else if (!busy)
  {
    busy = true;
    result.setAttribute("aria-live", "off");
  }
}

function comingFromWorker(e)
{
  fromWorker(e.data);
}

function dowork(n)
{
  let app = lang + n + (get("out").value.charCodeAt(0)-48)*4;
  let res = get("result");
  let polyText = get("poly").value;
  let modText = get("mod").value;
  let digitGroup = get("digits").value;
  hide("help");
  show("result");
  if (polyText === "")
  {
    res.innerHTML = (lang? "Por favor ingrese una expresión para el polinomio a evaluar." :
                           "Please type an expression for the polynomial to evaluate.");
    return;
  }
  if (modText === "")
  {
    res.innerHTML = (lang? "Por favor ingrese un número o expresión para el módulo." :
                           "Please type a number or expression for the modulus.");
    return;
  }
  get("eval").disabled = true;
  get("factor").disabled = true;
  get("stop").disabled = false;
  res.innerHTML = (lang? "Factorizando el polinomio..." :
                         "Factoring polynomial...");
  let param = digitGroup + "," + app + "," + modText + String.fromCharCode(0) + polyText +
  String.fromCharCode(0);
  if (!fileContents)
  {
    workerParam = param;
  }
  else
  {
    callWorker(param);
  }
}

function getCalcURLs()
{
  return ["polfactW0000.js",
          "polfact.webmanifest", "factpol.webmanifest", "polfact-icon-1x.png", "polfact-icon-2x.png", "polfact-icon-4x.png", "polfact-icon-180px.png", "polfact-icon-512px.png", "favicon.ico"];
}

function getFormSendValue()
{
  get("userdata").value = get("poly").value + " (mod " + get("mod").value + ")";
}

function popstate(event)
{
  if (get("feedback").style.display == "block" ||
      get("sentOK").style.display == "block" ||
      get("notSent").style.display == "block")
  {        // End feedback.
    show("main");
    hide("feedback");
    hide("sentOK");
    hide("notSent");
    get("poly").focus();   
  }
}

function startUp()
{
  get("btnSentOK").onclick = function()
  {
    history.back();
  }
  get("btnNotSent").onclick = function()
  {
    history.back();
  }
  get("stop").disabled = true;
  get("eval").onclick = function ()
  {
    dowork(2);
  };
  get("factor").onclick = function ()
  {
    dowork(0);
  };
  get("stop").onclick = function ()
  {
    endWorker();
    get("eval").disabled = false;
    get("factor").disabled = false;
    get("stop").disabled = true;
    get("result").innerHTML = 
      (lang? "<p>Factorización detenida por el usuario.</p>" :
             "<p>Factorization stopped by user</p>");
  };
  get("helpbtn").onclick = function ()
  {
    show("help");
    hide("result");
  };
  get("poly").oninput = function ()
  {
    let input = get("poly");
    let loc = input.selectionStart;
    if (input.value.substring(loc-1, loc) === ".")
    {
      input.value = input.value.replace(".", "x^");
      setTimeout(function()
      {    // Workaround for Firefox on Android. Delete extra space.
        input.value = input.value.replace("x ^", "x^");
        input.selectionStart = loc + 2;
        input.selectionEnd = loc + 2;
      }, 30);   
    }
  };
  get("formlink").onclick = clickFormLink;
  get("formcancel").onclick = function ()
  {
    history.back();
  };
  get("comments").oninput = function(_event)
  {
    get("formsend").disabled = (get("comments").value === "");
  };
  get("formsend").onclick = formSend;

  // Generate accordion.
  let acc = document.querySelectorAll("h2");
  let idx;

  for (idx = 0; idx < acc.length; idx++)
  {
    acc[idx >> 0].addEventListener("click", function()
    {
    // "active" means that panel is being displayed.
      this.children[0].classList.toggle("active");
      let panel = this.nextElementSibling;
      if (panel.style.display === "block")
      {
        panel.style.display = "none";
      }
      else
      {
        panel.style.display = "block";
      }
    });
  }
  let fromEcm = window.sessionStorage.getItem("F");
  let polyTextArea = get("poly");
  if (fromEcm != null)
  {    // Polynomial to factor coming from integer factorization calculator.
    window.sessionStorage.removeItem("F");
    polyTextArea.value = fromEcm;
    dowork(0);    // Factor polynomial.
  }
  fromEcm = window.sessionStorage.getItem("E");
  if (fromEcm != null)
  {    // Polynomial to factor coming from integer factorization calculator.
    window.sessionStorage.removeItem("E");
    polyTextArea.value = fromEcm;
    dowork(2);    // Evaluate polynomial.
  }
  let search = window.location.search;
  if (search.startsWith("?q="))
  {    // Process query string converting Tex to ASCII math.
    search = decodeURIComponent(search.substring(3));
    let fracNesting = new Array(300);
    let outBuffer = "";
    let nestingIndex = 0;
    let queryLength = search.length;
    let index = 0;
    while (index < queryLength)
    {
      let character = search.substring(index, index+1);
      if (character === "{")
      {        // Open paren in Tex.
        outBuffer += "(";
        nestingIndex++;
        index++;
      }
      else if (character === "}")
      {        // Close paren in Tex.
        outBuffer += ")";
        nestingIndex--;
        if (fracNesting[+nestingIndex] != null)
        {      // End of numerator. Insert slash.
          outBuffer += "/";
          fracNesting[+nestingIndex] = null;
        }
        index++;
      }
      else if (search.substring(index, index+5) === "\\frac")
      {        // Start of fraction. Indicate nesting index.
        fracNesting[+nestingIndex] = 1;
        index += 5;
      }
      else
      {        // Other character. Copy it as is.
        outBuffer += character;
        index++;
      }
    }
    polyTextArea.value = outBuffer;
    dowork(0);    // Factor polynomial.
  }
  registerServiceWorker();
};
getCalculatorCode("polfactW0000.js", workerParam);
window.addEventListener("load", startUp);
window.addEventListener("popstate", popstate);
window["fromWorker"] = fromWorker;
