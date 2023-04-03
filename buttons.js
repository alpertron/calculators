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
  input.value = input.value.substring(0, start) +
                event.target.innerText +
                input.value.substring(input.selectionEnd);
    // Place the caret at the end of the appended text.
  input.selectionStart = start + event.target.innerText.length;
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
    button.innerHTML = funcname[catIndex*2 + 1];         // Text of button.
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
