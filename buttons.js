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
  let input = currentInputBox;
  input.focus();
  let start = input.selectionStart;
  let chars = event.target.getAttribute("totalchars");
  if (chars === "\u23CE")
  {
    chars = "\n";
  }
  input.value = input.value.substring(0, start) +
                chars +
                input.value.substring(input.selectionEnd);
    // Place the caret at the end of the appended text.
  let offset = chars.indexOf("(") + 1;
  if (offset == 0)
  {
    offset = chars.length;
  }
  input.selectionStart = start + offset;
  input.selectionEnd = input.selectionStart;
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
  for (catIndex = 0; catIndex < funcname.length/2; catIndex++)
  {
    button = document.createElement("button");
    button.setAttribute("type", "button");        // Indicate this is a button, not submit.
    button.setAttribute("title", funcname[catIndex*2]);  // Text of tooltip.
    let btnName = funcname[catIndex*2 + 1];
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
    button.classList.add("funcbtn");
    button.onclick = buttonClick;
    fragment.appendChild(button);
  }
  funcbtns.innerHTML = "";
  funcbtns.appendChild(fragment);
}

function completeFuncButtons(funcButtons)
{
  let button;
  let catIndex;
  let funcname = (getParens() + getFuncNames()[0]).split(",");
  let funcbtns = get(funcButtons);
  for (catIndex = 0; catIndex < funcname.length/2; catIndex++)
  {
    button = funcbtns.children[+catIndex];
    button.setAttribute("title", funcname[catIndex*2]);  // Text of tooltip.
    button.onclick = buttonClick;
  } 
}
