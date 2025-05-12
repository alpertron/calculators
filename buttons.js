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
/* global currentInputBox */
/* global getFuncNames */
/* global get */
/* global getParens */
function buttonClick(event)
{
  event.preventDefault();  // Prevent button from getting focus.
  let input = currentInputBox;
  let start = input.selectionStart;
  let chars = event.target.getAttribute("totalchars");
  if (chars === "\u23CE")
  {
    chars = "\n";
  }
  let newValue = input.value.substring(0, start) +
                 chars +
                 input.value.substring(input.selectionEnd);
    // When appending a function, place the cursor just after
    // the opening paren. Otherwise, place the cursor at the
    // end of the inserted characters.
  let offset = chars.indexOf("(") + 1;
  if (offset === 0)
  { // Not a function. Place cursor at the end of these characters.
    offset = chars.length;
  }
  input.focus();
  input.value = newValue;
  setTimeout(() => {  // Required for Android TV.
    input.selectionEnd = start + offset;
    input.selectionStart = input.selectionEnd;
  }, window["isTV"]? 100: 0);
}

function setButtonInfo(button, funcname, catIndex)
{
  let btnName = funcname[catIndex*2 + 1];
  button.setAttribute("title", funcname[catIndex*2]);  // Text of tooltip.
  let nbrArguments = parseInt(btnName.slice(-1), 10);
  if (nbrArguments >= 1 && nbrArguments <= 9)
  {
    let text = btnName.slice(0, -1) + "(";
    // Set text of button.
    button.innerHTML = text;
    // Text to be displayed when the button is pressed.
    button.setAttribute("totalchars", text + ",".repeat(nbrArguments-1) + ")");
  }
  else
  {
    // Set text of button.
    button.innerHTML = btnName;
    // Text to be displayed when the button is pressed.
    button.setAttribute("totalchars", btnName);
  }
  button.onclick = buttonClick;
}
 
function generateFuncButtons(optionCategory, funcButtons)
{
  let button;
  let catIndex;
  let funcbtns = get(funcButtons);
  let catnbr = get(optionCategory).selectedIndex;
  let funcname = (getParens() + getFuncNames()[+catnbr]).split(",");
  // Append all buttons to document fragment instead of funcbtns
  // and finally append the fragment to funcbtns to minimize redraws.
  let fragment = document.createDocumentFragment();
  let buttonSpan = document.createElement("span");
  fragment.appendChild(buttonSpan);
  for (catIndex = 0; catIndex < funcname.length/2; catIndex++)
  {
    if (funcButtons === "wzdfuncbtns" && catIndex === 2)
    {   // Do not show Enter button on wizard.
      continue;
    }
    let btnName = funcname[catIndex*2 + 1];
    if (btnName == "")
    {
      buttonSpan = document.createElement("span");
      fragment.appendChild(buttonSpan);      
    }
    else
    {
      button = document.createElement("button");
      button.setAttribute("type", "button");        // Indicate this is a button, not submit.
      setButtonInfo(button, funcname, catIndex);
      button.classList.add("funcbtn");
      buttonSpan.appendChild(button);
    }
  }
  funcbtns.innerHTML = "";
  funcbtns.appendChild(fragment);
}

// Add tooltips and click handlers to buttons already defined in HTML.
function completeFuncButtons(funcButtons)
{
  let button;
  let catIndex;
  let funcname = (getParens() + getFuncNames()[0]).split(",");
  let funcBtns = get(funcButtons);
  let spanIndex = 0;
  let spanBtns = funcBtns.children[+spanIndex];
  if (spanBtns.tagName !== "SPAN")
  {
    spanBtns = funcBtns;
  }
  let btnIndex = 0;
  for (catIndex = 0; catIndex < funcname.length/2; catIndex++)
  {
    if (funcname[catIndex*2] === "")
    {
      if (spanBtns.tagName === "SPAN")
      {
        spanIndex++;
        btnIndex = 0;
        spanBtns = funcBtns.children[+spanIndex];
      }
    }
    else
    {
      button = spanBtns.children[+btnIndex];
      setButtonInfo(button, funcname, catIndex);
      btnIndex++;
    }
  } 
}
