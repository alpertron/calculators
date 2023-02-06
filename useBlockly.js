"use strict";
/* global goog */
/* global Blockly */
/** @type {function(string)} */
let fromBlocklyRun;
/** @type {number} */
let language;
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
/** @param {number} lang */
function useBlockly(callback, lang)
{
  if (callback == null)
  {
    blocklyResize(null);
    Blockly.svgResize(workspace);
    return;
  }
  language = lang;
  fromBlocklyRun = callback;
  /** @type {!Array<!Object>} */
  let blocksUncompressed = new Array();
  /** @type {number} */
  let index;
  /** @type {Array} */
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
    (lang? "Control de flujo": "Flow Control") +
    "\" colour=\"230\">" + 
      "{controls_repeat_ext[TIMES]}" +
      "{controls_if}" +
      "{controls_for[FROM][TO][BY]}" +
      "{controls_whileUntil}" +
    "</category>" +
    "<category name=\"" +
    (lang? "Matemática básica": "Basic Math") +
    "\" colour=\"355\">" +
      "<block type=\"M\"><field name=\"1\">5</field></block>";
 
  /** @type {Array<string>} */
  let defineBlocksEn =
  [
    "355;%1 + %2",
    "355;%1 - %2",
    "355;%1 × %2",
    "355;%1 ÷ %2",
    "355;remainder of %1 ÷ %2",
    "355;%1 to the %2",
    "355;square root of %1",
    "355;%2-th root of %1",
    "355;random integer from %1 to %2",
    "355;absolute value of %1",
    "355;sign of %1",
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
    "203;%1 shifted left by %2 bits",
    "203;%1 shifted right by %2 bits",
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
    "26;concatenate prime factors of %1 %2",
    "4077;Number Theory",
    "77;inverse of %1 modulo %2",
    "77;%1 divided by %2 modulo %3",
    "77;%1 to the %2 modulo %3",
    "77;totient of %1",
    "77;Jacobi symbol of %1 over %2",
    "4055;Other",
    "55;factorial of %1",
    "55;factorial of order %1 of %2",
    "55;primorial of %1",
    "55;element %1 of Fibonacci sequence",
    "55;element %1 of Lucas sequence",
    "55;partition of %1",
    "4170;Output",
    "1170;print %1",
    "1170;print %1 in hex",
    "1170;print prime factors of %1",
    "1170;print prime factors of %1 in hex",
    "1170;print is %1 prime",
    "1170;print is %1 prime in hex",
  ];
  
  /** @type {Array<string>} */
  let defineBlocksEs =
  [
    "355;%1 + %2",
    "355;%1 - %2",
    "355;%1 × %2",
    "355;%1 ÷ %2",
    "355;resto de %1 ÷ %2",
    "355;%1 a la %2",
    "355;raíz cuadrada de %1",
    "355;raíz %1-ésima de %2",
    "355;entero aleatorio entre %1 y %2",
    "355;valor absoluto de %1",
    "355;signo de %1",
    "4263;Comparaciones",
    "263;%1 = %2",
    "263;%1 ≠ %2",
    "263;%1 \x3E %2",
    "263;%1 ≤ %2",
    "263;%1 \x3C %2",
    "263;%1 ≥ %2",
    "4203;Lógica",
    "203;%1 AND %2",
    "203;%1 OR %2",
    "203;%1 XOR %2",
    "203;NOT %1",
    "203;%1 desplazado a la izquierda %2 bits",
    "203;%1 desplazado a la derecha %2 bits",
    "4143;Divisibilidad",
    "143;mcd de %1 y %2",
    "143;mcm de %1 y %2",
    "143;es %1 primo",
    "143;cantidad de factores primos de %1",
    "143;menor factor primo de %1",
    "143;mayor factor primo de %1",
    "143;cantidad de divisores de %1",
    "143;suma de los divisores de %1",
    "4026;Matemática recreativa",
    "26;siguiente primo después de %1",
    "26;último primo antes de %1",
    "26;cantidad de dígitos de %1 en base %2",
    "26;suma de dígitos de %1 en base %2",
    "26;inversión de dígitos de %1 en base %2",
    "26;concatenación de factores primos de %1 %2",
    "4077;Teoría de números",
    "77;inverso de %1 módulo %2",
    "77;%1 dividido %2 módulo %3",
    "77;%1 a la %2 módulo %3",
    "77;indicador de Euler de %1",
    "77;símbolo de Jacobi de %1 sobre %2",
    "4055;Otros",
    "55;factorial de %1",
    "55;factorial de orden %1 de %2",
    "55;primorial de %1",
    "55;elemento %1 de la secuencia de Fibonacci",
    "55;elemento %1 de la secuencia de Lucas",
    "55;particiones of %1",
    "4170;Salida",
    "1170;mostrar %1",
    "1170;mostrar %1 en hexa",
    "1170;mostrar factores primos de %1",
    "1170;mostrar factores primos de %1 en hexa",
    "1170;mostrar si %1 es primo",
    "1170;mostrar si %1 es primo en hexa",
  ];

  if (lang)
  {
    defineBlocks = defineBlocksEs;
  }
  else
  {
    defineBlocks = defineBlocksEn;
  }
  for (index=0; index<defineBlocks.length; index++)
  {
    destArray = new Array();
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
    /** @suppress {checkTypes} */
    destArray["type"] = String.fromCharCode(groupNbr, itemNbr);
    if (message.indexOf("concat") === 0)
    {
      let options;
      if (lang)
      {
        options = [["no repetidos en orden ascendente", "0"],
                   ["no repetidos en orden descendente", "1"],
                   ["repetidos en orden ascendente", "2"],
                   ["repetidos en orden descendente", "3"]];
      }
      else
      {
        options = [["not repeated in ascending order", "0"],
                   ["not repeated in descending order", "1"],
                   ["repeated in ascending order", "2"],
                   ["repeated in descending order", "3"]];
      }  
      /** @suppress {checkTypes} */
      destArray["args0"] = [{"type": "input_value", "name": "1"},
                            {"type": "field_dropdown", "name": "2",
                             "options": options}];
      ecmToolbar += "[1][2]}";
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
  let myToolbar = ecmToolbar.replace(/\{(\w+)([\[\]\w]*)}/g, "<block type=\"$1\">$2</block>");
  myToolbar = myToolbar.replace(/\[(\w+)]/g, "<value name=\"$1\"><shadow type=\"M\"><field name=\"1\">5</field></shadow></value>");
  Blockly.Blocks["M"] =
  {
    init:
    /** @this {Blockly.Block} */
    function()
    {
      /** @type {Blockly.FieldTextInput} */
      let field = new Blockly.FieldTextInput("5", BigIntValidator);
      field.setTooltip(lang?"Ingrese un número entero. Los números hexadecimales deben tener el prefijo 0x":
        "Enter an integer. Precede it by 0x to enter a hex number.");
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
      
  blocklyResize = function(_e) {
    // Compute the absolute coordinates and dimensions of blocklyArea.
    let element = blocklyArea;
    /** @type {number} */
    let x = 0;
    /** @type {number} */
    let y = 0;
    do {
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
    let contents = getStorage("blockly"+filename);
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
    setStorage("blockly"+filename, contents);
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
