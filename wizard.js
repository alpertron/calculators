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
/* global get */
/* global hide */
/* global oneexpr */
/* global saveConfig */
/* global show */
/* global value */
let wizardStep = 0;
let wizardTextInput;

function setWizardStep(stepNbr)
{
  wizardStep = stepNbr;
}

function clearWizardTextInput()
{
  wizardTextInput = "";
}

function selectLoop()
{   
  get("next").value = get("next").innerHTML;
  get("wzddesc").innerHTML = get("wzddesc1").innerHTML;
  get("wzdexam").innerHTML = get("novars").innerHTML;
  wizardStep = 1;
}
  
function wizardNext()
{
  let nextBtn = get("next");
  let wzdDescText = get("wzddesc");
  let wzdExamText = get("wzdexam");
  let wzdInput = get("wzdinput");
  let valueInput = value;
  let textExample = get("varsreq").innerHTML;
  nextBtn.disabled = true;
  switch (++wizardStep)
  {
    case 2:
      wizardTextInput += "x="+wzdInput.value;
      hide("wzdupper");
      wzdDescText.innerHTML = get("wzddesc2").innerHTML;
      wzdExamText.innerHTML = textExample + "<code>x+1</code>";
      break;
    case 3:
      wizardTextInput += ";x="+wzdInput.value;
      wzdDescText.innerHTML = get("wzddesc3").innerHTML;
      wzdExamText.innerHTML = textExample + "<code>x&lt;10000</code>";
      break;
    case 4:
      wizardTextInput += ";"+wzdInput.value;
      wzdDescText.innerHTML = get("wzddesc4").innerHTML;
      wzdExamText.innerHTML = textExample + "<code>x</code>";
      break;
    case 5:
      wizardTextInput += ";"+wzdInput.value;
      nextBtn.value = get("done").innerHTML;
      nextBtn.disabled = false;
      wzdDescText.innerHTML = get("wzddesc5").innerHTML;
      wzdExamText.innerHTML = textExample + "<code>sumdigits(x,10) == sumdigits(concatfact(2,x),10) and not isprime(x)</code>";
      break;
    case 6:
      if (wzdInput.value !== "")
      {
        wizardTextInput += ";"+wzdInput.value;
      }
      valueInput.value = wizardTextInput;
      wizardStep = 0;
      saveConfig(true);
      show("main");
      hide("wizard");
      history.back();
      valueInput.focus();
      break;
    default:
      wizardStep = 0;
      saveConfig(true);
      valueInput.value = wzdInput.value;
      show("main");
      hide("wizard");
      history.back();
      valueInput.focus();
      break;
  } 
  if (wizardStep)
  {
    wzdInput.value = "";
    wzdInput.focus();
  }
}

function typedOnWizard()
{
  let inputValue = get("wzdinput").value;
  let nextBtn = get("next");
  if (inputValue !== "")
  {         // User typed something on input box.
    if (wizardStep === 1 || wizardStep === 9 || (inputValue.lastIndexOf("x") >= 0 || inputValue.lastIndexOf("c") >= 0 ||
        inputValue.lastIndexOf("X") >= 0 || inputValue.lastIndexOf("C") >= 0))
    {       // At least one x or c. Indicate valid.
      nextBtn.disabled = false;
    }
    else
    {
      nextBtn.disabled = true;
    }
  }
  else if (wizardStep === 5)
  {         // Last step is optional, so empty input is valid.
    nextBtn.disabled = false;
  }
  else
  {         // For required input, empty input is invalid.
    nextBtn.disabled = true;
  }
}

function keyDownOnWizard(event)
{
  let keyCode = event.key;
  if (keyCode === "Enter")
  {
    if (!get("next").disabled)
    {                                // Next button is not disabled.
      wizardNext();                  // Perform same operation as if the user had pressed Next button.
    }
    event.stopPropagation();         // Do not propagate key.
    event.preventDefault();
  }
  if (keyCode === "Escape" || keyCode === "Esc")
  {
    show("main");
    hide("wizard");
  }
  if (event.altKey)
  {                                  // User pressed ALT key.
    if (keyCode === "P")
    {                                // User pressed ALT-P.
      event.stopPropagation();       // Do not propagate key.
      event.preventDefault();
      if (get("oneexpr").checked)
      {
        get("oneexpr").checked = false;
        get("loop").checked = true;
        selectLoop();
      }
      else
      {
        get("oneexpr").checked = true;
        get("loop").checked = false;
        oneexpr();
      }
    }
    else if (keyCode === "D")
    {                                // User pressed ALT-D.
      event.stopPropagation();       // Do not propagate key.
      event.preventDefault();
      get("decW").checked = true;
      get("hexW").checked = false;
    }
    else if (keyCode === "H")
    {                                // User pressed ALT-H.
      event.stopPropagation();       // Do not propagate key.
      event.preventDefault();
      get("decW").checked = false;
      get("hexW").checked = true;
    }
  }
  return true;
}
