"use strict";
/* global goog */
/* global Blockly */
/** @type {function(string)} */
let fromBlocklyRun;
/** @type {number} */
goog.provide("BigIntField");
goog.require("Blockly");
goog.require("Blockly.Comment");
goog.require("Blockly.VerticalFlyout");
goog.require("Blockly.FlyoutButton");
goog.require("Blockly.Toolbox");
goog.require("Blockly.Trashcan");
goog.require("Blockly.VariablesDynamic");
goog.require("Blockly.ZoomControls");
goog.require("Blockly.ShortcutItems");
goog.require("Blockly.ContextMenuItems");
goog.require("Blockly.Mutator");
goog.require("Blockly.Warning");
goog.require("Blockly.FieldDropdown");
goog.require("Blockly.FieldLabelSerializable");
goog.require("Blockly.FieldTextInput");
goog.require("Blockly.FieldVariable");
goog.require("Blockly.geras.Renderer");
goog.require("Blockly.Themes.Classic");
goog.require("Blockly.Constants.Logic");
goog.require("Blockly.Constants.Loops");
goog.require("Blockly.Constants.Variables");
goog.require("Blockly.Constants.VariablesDynamic");
goog.require("Blockly.Xml");

/** @type {Function} */
let blocklyResize;
let workspace;
let BigIntField;

function get(id)
{
  return document.getElementById(id);
}

function hide(id)
{
  get(id).style.display = "none";
}

function show(id)
{
  get(id).style.display = "block";
}

function BigIntValidator(newValue)
{
  /** @type {number} */
  let count;
  /** @type {boolean} */
  let insideDigits = false;
  /** @type {boolean} */
  let insideHex = false;
  /** @type {boolean} */
  let isMinus = false;
  /** @type {string} */
  let output = "";
  for (count=0; count<newValue.length; count++)
  {
    /** @type {string} */
    let c = newValue.charAt(count);
    if (c === " ")
    {
      continue;
    }
    if (insideHex)
    {
      if ((c < "0" || c > "9") && (c < "A" || c > "F") && (c < "a" || c > "f"))
      {              // Not a hex character.
        return null;
      }
      output += c;
    }
    else if (insideDigits)
    {
      if (c < "0" || c > "9")
      {              // Not a digit.
        return null;
      }
      output += c;
    }
    else
    {
      if (c >= "0" && c <= "9") 
      {
        if (isMinus)
        {
          output = "-";
        }
        if (c === "0" && (newValue.charAt(count+1) === "X" || newValue.charAt(count+1) === "x"))
        {
          insideHex = true;
          count++;      // Discard hex prefix.
        }
        else
        {
          insideDigits = true;
          output += c;  // Append first digit.
        }
      }
      else if (c === "-")
      {
        isMinus = !isMinus;
      }
      else if (c !== "+")
      {
        return null;
      }
    }
  }
  if (!insideDigits)
  {
    return null;
  }
  return output;
}

/** @param {function(string)} callback */
function useBlockly(callback)
{
  if (callback == null)
  {
    blocklyResize(null);
    Blockly.svgResize(workspace);
    return;
  }
  fromBlocklyRun = callback;
  /** @type {!Array<!Object>} */
  let blocksUncompressed = new Array();
  /** @type {number} */
  let index;
  /** @type {BlockInfo} */
  let destArray;
  /** @type {Array<string>} */
  let defineBlocks;
  /** @type {number} */
  let uncompressedIndex = 0;
  /** @type {number} */
  let groupNbr = 65;
  /** @type {number} */
  let itemNbr = 65;
  /** @type {string} */
  let ecmToolbar = "<xml>" +
    "<category name=\"" +
    get("flowctrl").textContent +
    "\" colour=\"230\">" + 
      "{controls_repeat_ext[TIMES]}" +
      "{controls_if}" +
      "{controls_for[FROM][TO][BY]}" +
      "{controls_whileUntil}" +
    "</category>" +
    "<category name=\"" +
    get("basicmath").textContent +
    "\" colour=\"355\">" +
      "<block type=\"M\"><field name=\"1\">5</field></block>";

  const div = get("defineblock");
  const defineBlocks = Array.from(div.children).map(el => el.textContent.trim());
  for (index=0; index<defineBlocks.length; index++)
  {
    destArray = {};
    /** @type {Array<string>} */
    let oneBlock = defineBlocks[+index].split(";");
    /** @type {number} */
    let nbr = +oneBlock[0];
    /** @type {string} */
    let message = oneBlock[1];
    if (nbr >= 4000)
    {
      ecmToolbar += "</category><category name=\"" + message + "\" colour=\"" + (nbr - 4000) + "\">";
      groupNbr++;
      itemNbr = 65;
      continue;
    }
    ecmToolbar += "{" + String.fromCharCode(groupNbr, itemNbr);
    destArray.type = String.fromCharCode(groupNbr, itemNbr);
    if (message.startsWith("concat"))
    {
      const div = get("options");
      const options = Array.from(div.querySelectorAll("span"))
        .map((span, index) => [span.textContent.trim(), index.toString()]);
      destArray.args0 = [{"type": "input_value", "name": "1"},
                         {"type": "field_dropdown", "name": "2",
                          "options": options}];
      ecmToolbar += "[1][2]}";
    }
    else if (message.indexOf("%3") >= 0)
    {
      destArray.args0 = [{"type": "input_value", "name": "1"},
                         {"type": "input_value", "name": "2"},
                         {"type": "input_value", "name": "3"}];
      ecmToolbar += "[1][2][3]}";
    }
    else if (message.indexOf("%2") >= 0)
    {
      destArray.args0 = [{"type": "input_value", "name": "1"},
                         {"type": "input_value", "name": "2"}];
      ecmToolbar += "[1][2]}";
    }
    else
    {
      destArray.args0 = [{"type": "input_value", "name": "1"}];
      ecmToolbar += "[1]}";
    }
    destArray.message0 = message;
    if (oneBlock[0] >= 1000)
    {       // Statement.
      destArray.previousStatement = null;
      destArray.nextStatement = null;
      destArray.output = undefined;
    }
    else
    {
      destArray.previousStatement = undefined;
      destArray.nextStatement = undefined;
      destArray.output = null;
    }
    destArray.inputsInline = true;
    destArray.colour = oneBlock[0] % 1000;
    blocksUncompressed[+uncompressedIndex] = destArray;
    uncompressedIndex++;
    itemNbr++;
  }
  Blockly.defineBlocksWithJsonArray(blocksUncompressed);
  ecmToolbar += "</category>" +
    "<category name=\"Variables\" custom=\"VARIABLE\" colour=\"330\"></category>" +
    "</xml>";
  let blocklyArea = get("blocklyArea");
  let blocklyDiv = get("blocklyDiv");
  /** @type {string} */
  let myToolbar = ecmToolbar.replace(/\{(\w+)([[\]\w]*)}/g, "<block type=\"$1\">$2</block>");
  myToolbar = myToolbar.replace(/\[(\w+)]/g, "<value name=\"$1\"><shadow type=\"M\"><field name=\"1\">5</field></shadow></value>");
  Blockly.Blocks["M"] =
  {
    init:
    /** @this {Blockly.Block} */
    function()
    {
      /** @type {Blockly.FieldTextInput} */
      let field = new Blockly.FieldTextInput("5", BigIntValidator);
      field.setTooltip(get("inputint").textContent);
      this.appendDummyInput().appendField(field, "1");
      this.setOutput(true, null);
      this.setColour(355);
    },
  };
  workspace = Blockly.inject(blocklyDiv,
      {"toolbox": myToolbar,
      "zoom": {"controls": true, "wheel": true},
      "media": "/",
      "css": false});
      
  blocklyResize = function(_e)
  {
    // Compute the absolute coordinates and dimensions of blocklyArea.
    let element = blocklyArea;
    /** @type {number} */
    let x = 0;
    /** @type {number} */
    let y = 0;
    do
    {
      x += element.offsetLeft;
      y += element.offsetTop;
      element = element.offsetParent;
    } while (element);
    // Position blocklyDiv over blocklyArea.
    blocklyDiv.style.left = x + "px";
    blocklyDiv.style.top = y + "px";
    blocklyDiv.style.width = blocklyArea.offsetWidth + "px";
    blocklyDiv.style.height = (window.innerHeight - y) + "px";
    Blockly.svgResize(workspace);
  };
  window.addEventListener("resize", blocklyResize, false);
  get("deleteBlocks").onclick = function()
  {
    /** @type {number} */
    let count = workspace.getAllBlocks().length;
    if (count < 2 || window.confirm("Delete all " + count + " blocks?"))
    {
      workspace.clear();
    }
  };
  get("bload").onclick = function()
  {
    /** @type {string} */
    let filename = get("bfilename").value.trim();
    /** @type {string|null} */
    let contents = window.localStorage.getItem("blockly"+filename);
    if (contents == null)
    {
      alert("File "+ filename +" not found");
      return;
    }
    workspace.clear();
    let xml = Blockly.Xml.textToDom(contents);
    Blockly.Xml.domToWorkspace(xml, workspace);
  };
  get("bsave").onclick = function()
  {
    /** @type {string} */
    let filename = get("bfilename").value.trim();
    let xml = Blockly.Xml.workspaceToDom(workspace);
    let contents = Blockly.Xml.domToText(xml);
    window.localStorage.setItem("blockly"+filename, contents);
  };
  get("runBlockly").onclick = function()
  {
    let xml = Blockly.Xml.workspaceToDom(workspace);
    /** @type {string} */
    let contents = Blockly.Xml.domToText(xml);
    console.log(contents);
    fromBlocklyRun(contents);
  };
  get("berrcont").onclick = function()
  {
    hide("BlocklyErrors");
    show("BlocklyButtons");
  };
  blocklyResize(null);
  Blockly.svgResize(workspace);
}

window["useBlockly"] = useBlockly;
