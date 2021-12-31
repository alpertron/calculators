"use strict";
var goog;
var fromBlocklyRun;
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
goog.require("Blockly.Blocks.procedures");

var blocklyResize;
var workspace;
var BigIntField;
function get(id)
{
  return document.getElementById(id);
}

function setStorage(name, data)
{
  window.localStorage.setItem(name, data);
}

function getStorage(name)
{
  return window.localStorage.getItem(name);
}

function BigIntValidator(newValue)
{
  var count;
  var insideDigits = false;
  var insideHex = false;
  var isMinus = false;
  var output = "";
  for (count=0; count<newValue.length; count++)
  {
    var c = newValue.charAt(count);
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
    if (insideDigits)
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
          count++;    // Discard hex prefix.
          continue;
        }
        insideDigits = true;
        output += c;
      }
      else if (c === "+")
      {
        continue;
      }
      else if (c === "-")
      {
        isMinus = !isMinus;
      }
      else
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

function useBlockly(callback)
{
  if (callback == null)
  {
    blocklyResize(null);
    Blockly.svgResize(workspace);
    return;
  }
  fromBlocklyRun = callback;
  var blocksUncompressed = new Array();
  var index;
  var destArray;
  var uncompressedIndex = 0;
  var groupNbr = 65;
  var itemNbr = 65;
  var ecmToolbar = "<xml>" +
    "<category name=\"Flow Control\" colour=\"230\">" + 
      "{controls_repeat_ext[TIMES]}" +
      "{controls_if}" +
      "{controls_for[FROM][TO][BY]}" +
      "{controls_whileUntil}" +
    "</category>" +
    "<category name=\"Basic Math\" colour=\"355\">" +
      "<block type=\"M\"><field name=\"1\">5</field></block>";
 
  var defineBlocks =
  [
    "355;%1 + %2",
    "355;%1 - %2",
    "355;%1 × %2",
    "355;%1 ÷ %2",
    "355;remainder of %1 ÷ %2",
    "355;%1 to the %2",
    "355;square root of %1",
    "355;absolute value of %1",
    "355;random integer from %1 to %2",
    "4263;Comparisons",
    "263;%1 = %2",
    "263;%1 ≠ %2",
    "263;%1 \x3E %2",
    "263;%1 ≤ %2",
    "263;%1 \x3C %2",
    "263;%1 ≥ %2",
    "4203;Logic",
    "203;%1 AND %2",
    "203;%1 OR %2",
    "203;%1 XOR %2",
    "203;NOT %1",
    "203;shift left %1 by %2 bits",
    "203;shift right %1 by %2 bits",
    "4143;Divisibility",
    "143;gcd of %1 and %2",
    "143;lcm of %1 and %2",
    "143;is %1 prime",
    "143;number of prime factors of %1",
    "143;smallest prime factor of %1",
    "143;greatest prime factor of %1",
    "143;number of divisors of %1",
    "143;sum of divisors of %1",
    "4026;Recreational Math",
    "26;next prime after %1",
    "26;last prime before %1",
    "26;number of digits of %1 in base %2",
    "26;sum of digits of %1 in base %2",
    "26;reverse digits of %1 in base %2",
    "26;concatenate prime factors of %1 using mode %2",
    "4077;Number Theory",
    "77;inverse of %1 modulo %2",
    "77;%1 to the %2 modulo %3",
    "77;totient of %1",
    "77;Jacobi symbol of %1 over %2",
    "4055;Other",
    "55;factorial of %1",
    "55;factorial of order %1 of %2",
    "55;primorial of %1",
    "55;Element %1 of Fibonacci sequence",
    "55;Element %1 of Lucas sequence",
    "55;partition of %1",
    "4170;Output",
    "1170;print %1",
    "1170;print %1 in hex",
    "1170;print prime factors of %1",
    "1170;print prime factors of %1 in hex",
  ];
  for (index=0; index<defineBlocks.length; index++)
  {
    destArray = new Array();
    var oneBlock = defineBlocks[index].split(";");
    var nbr = +oneBlock[0];
    var message = oneBlock[1];
    if (nbr >= 4000)
    {
      ecmToolbar += "</category><category name=\"" + message + "\" colour=\"" + (nbr - 4000) + '">';
      groupNbr++;
      itemNbr = 65;
      continue;
    }
    ecmToolbar += "{" + String.fromCharCode(groupNbr, itemNbr);
    /** @suppress {checkTypes} */
    destArray["type"] = String.fromCharCode(groupNbr, itemNbr);
    if (message.indexOf("%concat") >= 0)
    {
      /** @suppress {checkTypes} */
      destArray["args0"] = [{"type": "input_value", "name": "1"},
                            {"type": "dropdown", "name": "2",
                             "options": [["0", "0"], ["1", "1"], ["2", "2"], ["3", "3"]]}];
      ecmToolbar += "[1]}";
    }
    else if (message.indexOf("%3") >= 0)
    {
      /** @suppress {checkTypes} */
      destArray["args0"] = [{"type": "input_value", "name": "1"},
                            {"type": "input_value", "name": "2"},
                            {"type": "input_value", "name": "3"}];
      ecmToolbar += "[1][2][3]}";
    }
    else if (message.indexOf("%2") >= 0)
    {
      /** @suppress {checkTypes} */
      destArray["args0"] = [{"type": "input_value", "name": "1"},
                            {"type": "input_value", "name": "2"}];
      ecmToolbar += "[1][2]}";
    }
    else
    {
      /** @suppress {checkTypes} */
      destArray["args0"] = [{"type": "input_value", "name": "1"}];
      ecmToolbar += "[1]}";
    }
    /** @suppress {checkTypes} */
    destArray["message0"] = message;
    if (oneBlock[0] >= 1000)
    {       // Statement.
      /** @suppress {checkTypes} */
      destArray["previousStatement"] = null;
      /** @suppress {checkTypes} */
      destArray["nextStatement"] = null;
    }
    else
    {
      /** @suppress {checkTypes} */
      destArray["output"] = null;
    }
    /** @suppress {checkTypes} */
    destArray["inputsInline"] = true;
    /** @suppress {checkTypes} */
    destArray["colour"] = oneBlock[0] % 1000;
    blocksUncompressed[uncompressedIndex] = destArray;
    uncompressedIndex++;
    itemNbr++;
  }
  Blockly.defineBlocksWithJsonArray(blocksUncompressed);
  ecmToolbar += "</category>" +
    "<category name=\"Variables\" custom=\"VARIABLE\"></category>" +
    "<category name=\"Functions\" custom=\"PROCEDURE\"></category>" +
    "</xml>";
  var blocklyArea = get("blocklyArea");
  var blocklyDiv = get("blocklyDiv");
  var myToolbar = ecmToolbar.replace(/\{(\w+)([\[\]\w]*)}/g, "<block type=\"$1\">$2</block>");
  myToolbar = myToolbar.replace(/\[(\w+)]/g, '<value name="$1"><shadow type="M"><field name="1">5</field></shadow></value>');
  Blockly.Blocks["M"] =
  {
    "init": function()
    {
      var field = new Blockly.FieldTextInput("5", BigIntValidator);
      field.setTooltip("Enter an integer. Precede it by 0x to enter a hex number.");
      this.appendDummyInput().appendField(field, "1");
      this.setOutput(true, null);
      this.setColour(355);
    },
  };
  workspace = Blockly.inject(blocklyDiv,
      {"toolbox": myToolbar,
      "zoom": {"controls": true, "wheel": true,}});
      
  blocklyResize = function(e) {
    // Compute the absolute coordinates and dimensions of blocklyArea.
    var element = blocklyArea;
    var x = 0;
    var y = 0;
    do {
      x += element.offsetLeft;
      y += element.offsetTop;
      element = element.offsetParent;
    } while (element);
    // Position blocklyDiv over blocklyArea.
    blocklyDiv.style.left = x + 'px';
    blocklyDiv.style.top = y + 'px';
    blocklyDiv.style.width = blocklyArea.offsetWidth + 'px';
    blocklyDiv.style.height = (window.innerHeight - y) + 'px';
    Blockly.svgResize(workspace);
  };
  window.addEventListener('resize', blocklyResize, false);
  get("deleteBlocks").onclick = function()
  {
    var count = workspace.getAllBlocks().length;
    if (count < 2 || window.confirm('Delete all ' + count + ' blocks?'))
    {
      workspace.clear();
    };
  };
  get("bload").onclick = function()
  {
    var filename = get("bfilename").value.trim();
    var contents = getStorage("blockly"+filename);
    if (contents == null)
    {
      alert("File "+ filename +" not found");
      return;
    }
    workspace.clear();
    var xml = Blockly.Xml.textToDom(contents);
    Blockly.Xml.domToWorkspace(xml, workspace);
  };
  get("bsave").onclick = function()
  {
    var filename = get("bfilename").value.trim();
    var xml = Blockly.Xml.workspaceToDom(workspace);
    var contents = Blockly.Xml.domToText(xml);
    setStorage("blockly"+filename, contents);
  };
  get("runBlockly").onclick = function()
  {
    var xml = Blockly.Xml.workspaceToDom(workspace);
    var contents = Blockly.Xml.domToText(xml);
    alert("Shortly this application will support running Blockly!!!");
//    console.log(contents);
    fromBlocklyRun(contents);
  };
  blocklyResize(null);
  Blockly.svgResize(workspace);
}

window["useBlockly"] = useBlockly;
