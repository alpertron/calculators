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
/* global lang */
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
  get("next").value = (lang ? "Siguiente": "Next");
  get("wzddesc").innerHTML = (lang ? "Paso 1 de 5: Valor inicial de x": "Step 1 of 5: Initial value of x");
  get("wzdexam").innerHTML = (lang? "No usar variables <var>x</var> o <var>c</var>. Ejemplo para números de Smith menores que 10000: <code>1</code>": 
                                       "Do not use variables <var>x</var> or <var>c</var>. Example for Smith numbers less than 10000: <code>1</code>");
  wizardStep = 1;
}
  
function wizardNext()
{
  let nextBtn = get("next");
  let wzdDescText = get("wzddesc");
  let wzdExamText = get("wzdexam");
  let wzdInput = get("wzdinput");
  let valueInput = value;
  let textExample = (lang? "Variables <var>x</var> y/o <var>c</var> requeridas. Ejemplo para números de Smith menores que 10000: <code>":
                           "Variables <var>x</var> and/or <var>c</var> required. Example for Smith numbers less than 10000: <code>");
  nextBtn.disabled = true;
  switch (++wizardStep)
  {
    case 2:
      wizardTextInput += "x="+wzdInput.value;
      hide("mode");
      wzdDescText.innerHTML = (lang? "Paso 2 de 5: Valor de x para la nueva iteración": "Step 2 of 5: Value of x for new iteration");
      wzdExamText.innerHTML = textExample + "x+1</code>";
      break;
    case 3:
      wizardTextInput += ";x="+wzdInput.value;
      wzdDescText.innerHTML = (lang? "Paso 3 de 5: Condición para finalizar el ciclo": "Step 3 of 5: End loop condition");
      wzdExamText.innerHTML = textExample + "x&lt;10000</code>";
      break;
    case 4:
      wizardTextInput += ";"+wzdInput.value;
      wzdDescText.innerHTML = (lang? "Paso 4 de 5: Expresión a factorizar": "Step 4 of 5: Expression to factor");
      wzdExamText.innerHTML = textExample + "x</code>";
      break;
    case 5:
      wizardTextInput += ";"+wzdInput.value;
      nextBtn.value = (lang? "Hecho": "Done");
      nextBtn.disabled = false;
      wzdDescText.innerHTML = (lang? "Paso 5 de 5: Condición para procesar la expresión": "Step 5 of 5: Process expression condition");
      wzdExamText.innerHTML = textExample + "sumdigits(x,10) == sumdigits(concatfact(2,x),10) and not isprime(x)</code>";
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
