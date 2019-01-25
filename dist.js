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

// In order to reduce the number of files to read from Web server, this 
// Javascript file includes both the Javascript in the main thread and the 
// Javascript that drives WebAssembly on its own Web Worker.
(function(global)
{   // This method separates the name space from the Google Analytics code.
var testNbr, questionNbr, score, withCountries;
var parcInterp, finalInterp;
var kmText, degreeText, currentText, finalText, pointText;
var direction = [];
var ordinal = [];
var cityData = [];
var cityName = [];
var countryName = [];
var Latitude = [];
var Longitude = [];
var cityIdx = [];
var Rad = Math.PI/(180 * 60); // Minutes to radians
var cityFrom, cityTo;
function get(x)
{
  return document.getElementById(x);
}

function ClearScreen()
{
  var donotshow = document.getElementsByClassName("donotshow");
  var index;
  for (index=0; index<donotshow.length; index++)
  {
    donotshow[index].style.display = "none";
  }
  get("mainmenu").style.display = "none";
  get("help").style.display = "none";
}

function ShowMainMenu()
{
  ClearScreen();
  get("mainmenu").style.display = "block";
  get("help").style.display = "block";
}

function SetTextToClass(text, className)
{
  var elems = document.getElementsByClassName(className);
  var index;
  for (index=0; index<elems.length; index++)
  {
    elems[index].innerHTML = text;
  }
}

function getCityName(cityIndex, withCountries)
{
  if (withCountries)
  {
    return cityName[cityIndex]+" ("+countryName[cityIndex]+")";
  }
  else
  {
    return cityName[cityIndex];
  }
}

function StartTestType1()
{
  ClearScreen();
  testNbr = 1;
  questionNbr = 1;
  score = 0;
  Test1(false);
}

function StartTestType2()
{
  ClearScreen();
  testNbr = 2;
  questionNbr = 1;
  score = 0;
  Test2();
}

function Test1(complete)
{
  var cityNameFrom, cityNameTo;
  var LatFrom, LatTo;
  var LongFrom, LongTo;
  var playerDist, trueDist;
  var partialScore;
  var test11, test12;
  var dist1;
  var show, index;
  ClearScreen();
  if (questionNbr != 0)
  {      // Inside test 1.
    if (complete == false)
    {
      cityFrom = Math.floor(Math.random()*cityData.length);
      do
      {
        cityTo = Math.floor(Math.random()*cityData.length);
      } while (cityFrom == cityTo);
    }
    get("ord1").innerHTML = ordinal[questionNbr - 1];
    get("score1").innerHTML = score + " " + pointText + (score == 1?"":"s");
  }
  else
  {      // Query distance between two cities.
    cityFrom = parseInt(get("cityFrom").value, 10) - 1;
    cityTo = parseInt(get("cityTo").value, 10) - 1;
  }
  cityNameFrom = getCityName(cityFrom, withCountries || questionNbr == 0);
  cityNameTo = getCityName(cityTo, withCountries || questionNbr == 0);
  LatFrom = Latitude[cityFrom];
  LongFrom = Longitude[cityFrom];
  LatTo = Latitude[cityTo];
  LongTo = Longitude[cityTo];
  SetTextToClass(cityNameFrom, "cityctry_from");
  SetTextToClass(cityNameTo, "cityctry_to");
  show = document.getElementsByClassName(questionNbr == 0?"findDist": "notFindDist");
  for (index=0; index<show.length; index++)
  {
    show[index].style.display = "block";
  }
  show = document.getElementsByClassName(questionNbr == 0?"notFindDist": "findDist");
  for (index=0; index<show.length; index++)
  {
    show[index].style.display = "none";
  }
  if (complete)
  {
    playerDist = parseInt(get("dist1").value);
    trueDist = getDistance(LatFrom, LongFrom, LatTo, LongTo);
    get("dist12_1").innerHTML = playerDist + " " + kmText + (playerDist==1?"":"s");
    get("dist12_2").innerHTML = get("dist12_3").innerHTML = trueDist + " " + kmText + (trueDist==1?"":"s");
    cityNameFrom = getCityName(cityFrom, false);  // Do not append country name.
    cityNameTo = getCityName(cityTo, false);
    SetTextToClass(cityNameFrom, "city_from");
    SetTextToClass(cityNameTo, "city_to");
    get("dirdeg1").innerHTML = getDirection(LatFrom, LongFrom, LatTo, LongTo);
    get("dirdeg2").innerHTML = getDirection(LatTo, LongTo, LatFrom, LongFrom);
    get("sun_from_0621").innerHTML = getDaytime(-LatFrom, LongFrom);
    get("sun_from_1221").innerHTML = getDaytime(LatFrom, LongFrom);
    get("sun_to_0621").innerHTML = getDaytime(-LatTo, LongFrom);
    get("sun_to_1221").innerHTML = getDaytime(LatTo, LongFrom);
    if (questionNbr != 0)
    {
      partialScore = Math.round(100 - 100 * Math.abs(Math.log(trueDist/playerDist)));
      if (partialScore < 0) partialScore = 0;
      if (partialScore == 100) Interpretation = parcInterp[0];
      else if (partialScore > 94) Interpretation = parcInterp[1];
      else if (partialScore > 89) Interpretation = parcInterp[2];
      else if (partialScore > 79) Interpretation = parcInterp[3];
      else if (partialScore > 69) Interpretation = parcInterp[4];
      else if (partialScore > 59) Interpretation = parcInterp[5];
      else if (partialScore > 39) Interpretation = parcInterp[6];
      else if (partialScore > 19) Interpretation = parcInterp[7];
      else Interpretation = parcInterp[8];
      get("score2").innerHTML = partialScore;
      score += partialScore;
      get("parcInter").innerHTML = Interpretation;
      get("score3").innerHTML = score;
      if (questionNbr == 10)
      {
        if (score > 949) Interpretation = finalInterp[0];
        else if (score > 899) Interpretation = finalInterp[1];
        else if (score > 799) Interpretation = finalInterp[2];
        else if (score > 699) Interpretation = finalInterp[3];
        else if (score > 599) Interpretation = finalInterp[4];
        else if (score > 399) Interpretation = finalInterp[5];
        else Interpretation = finalInterp[6];
        get("finalInter").innerHTML = Interpretation;
        get("finalInter").style.display = "block";
        get("scoreType").innerHTML = finalText;
      }
      else
      {
        get("finalInter").style.display = "none";
        get("scoreType").innerHTML = currentText;
      }
      questionNbr++;
      test11 = document.getElementsByClassName("test1_1");
      for (index=0; index<test11.length; index++)
      {
        test11[index].style.display="none";
      }
      test12 = document.getElementsByClassName("test1_2");
      for (index=0; index<test12.length; index++)
      {
        test12[index].style.display="block";
      }
      get("nextq").style.display = (questionNbr>10? "none": "inline");
    }
  }
  else
  {
    test11 = document.getElementsByClassName("test1_1");
    for (index=0; index<test11.length; index++)
    {
      test11[index].style.display="block";
    }
    test12 = document.getElementsByClassName("test1_2");
    for (index=0; index<test12.length; index++)
    {
      test12[index].style.display="none";
    }
  }
  get("test1").style.display="block";
  if (complete == false && questionNbr != 0)
  {
    dist1 = get("dist1");
    dist1.focus();
    dist1.scrollIntoView();
    dist1.value = "";
  }
}

function getDistance(Lat1, Lon1, Lat2, Lon2)
{
  var q = Math.cos(Lat1*Rad)*Math.cos(Lat2*Rad)*Math.cos((Lon1-Lon2)*Rad)+Math.sin(Lat1*Rad)*Math.sin(Lat2*Rad);
  if (q >= 1)
  {
    return 0;
  }
  if (q < -1)
  {
    return 20012;
  }
  return Math.round(6370*Math.acos(q));
}

function getDirection(Lat1, Lon1, Lat2, Lon2)
{
  var Angle;
  var AngAux;
  var Difference;

  if (Lon1 == Lon2)
  {
    Angle = (Lat1>Lat2)?180:0;
  }
  else
  {
    if (Lat1 == 90*60 || Lat2 == -90*60)
    {
      Angle = 180;
    }
    else
    {
      if (Lat1 == -90*60 || Lat2 == 90*60)
      {
        Angle = 0;
      }
      else
      {
        Difference = Lon2 - Lon1;
        Angle = 90 - (Math.round(180/Math.PI*Math.atan((Math.tan(Lat2*Rad)*Math.cos(Lat1*Rad)-Math.sin(Lat1*Rad)*Math.cos(Difference*Rad))/Math.sin(Difference*Rad))));
        Difference = (Difference + 360*60) % (360*60);
        if (Difference > 180*60)
        {
          Angle = (Angle + 180) % 360;
        }
      }
    }
  }
  AngAux = (Angle + 348) % 360;
  return Angle + " " + degreeText + (Angle==1?"":"s") + " (" + direction[Math.floor(AngAux/22.5)] + ")";
}

function getDaytime(Lat, Lon)
{
  var q;
  var minutesDaytime;

  if (Math.abs(Lat) == 90*60)
  {
    q = Lat;
  }
  else
  {    // Use generalized equation (all values in degrees): (sin 0.83 + sin latitude * sin d) / (cos latitude * cos d)
       // where d = Sun declination = 23.437 degrees for solstices.
    q = (-0.01449 + Math.sin(Lat*Rad)*0.39774) / (Math.cos(Lat*Rad) * 0.9175);
  }
  if (q >= 1)
  {
    minutesDaytime = 0;
  }
  else if (q <= -1)
  {
    minutesDaytime = 1440;
  }
  else
  {
    minutesDaytime = Math.round(458.3662*Math.acos(q));
  }
  return Math.floor(minutesDaytime/60) + "h " + (minutesDaytime%60) + "m";
}

function Test2()
{
  var idx, innerIdx, randomValue;
  
  ClearScreen();
  get("ord2").innerHTML = ordinal[questionNbr - 1];
  for (idx=0; idx<6; idx++)
  {
    do
    {
      randomValue = Math.floor(Math.random()*cityData.length);
      for (innerIdx=0; innerIdx<idx; innerIdx++)
      {
        if (cityIdx[innerIdx] == randomValue)
        {
          break;    
        }
      }
    } while (idx != innerIdx);
    cityIdx[idx] = randomValue;
    get("city"+(idx+1)).innerHTML = getCityName(randomValue, withCountries);
  }
  get("score4").innerHTML = score + " " + pointText + (score == 1?"":"s");
  get("order").value = "";
  get("test2").style.display="block";
}

function ShowResultsTest2(playerInput)
{
  var combination, position, parcScore;
  var arrayDist = new Array(25);
  var arrayDistComb = new Array(24);
  var arrayComb = new Array(24);
  var i = 0;
  var j, k, u, v, temp;
  var text = "<ol>";
  for (j=0; j<5; j++)
  {
    for (k=1; k<=5; k++)
    {
      arrayDist[i++] = getDistance(Latitude[cityIdx[j]], Longitude[cityIdx[j]], Latitude[cityIdx[k]], Longitude[cityIdx[k]]);
    }
  }
  for (j=0; j<24; j++)
  {
    arrayComb[j] = "123412431324134214231432213421432314234124132431312431423214324134123421412341324213423143124321".substring(j*4,j*4+4);
    combination = "0" + arrayComb[j] + "5";
    arrayDistComb[j] = 0;
    for (k=0; k<5; k++)
    {
      u = parseInt(combination.substring(k,k+1), 10);
      v = parseInt(combination.substring(k+1,k+2), 10);
      arrayDistComb[j] += arrayDist[u*5+v-1];
    }
    for (k=j; k>0; k--)
    {
      if (arrayDistComb[k] > arrayDistComb[k-1])
      {
        break;
      }
      i = arrayDistComb[k];
      arrayDistComb[k] = arrayDistComb[k-1];
      arrayDistComb[k-1] = i;
      temp = arrayComb[k];
      arrayComb[k] = arrayComb[k-1];
      arrayComb[k-1] = temp;
    }
  }
  position = -1;
  for (j=0; j<24; j++)
  {
    combination = arrayComb[j];
    if (combination == playerInput)
    {
      position = j;
      text += '<li class="red">';
    }
    else
    {
      text += '<li>';
    }
    text += getCityName(cityIdx[0], withCountries) + " &rarr; ";
    for (k=0; k<4; k++)
    {
      text += getCityName(cityIdx[parseInt(combination.substring(k, k+1), 10)], withCountries) + " &rarr; ";
    }
	parcScore = Math.round((1-j/23)*(1-j/23)*100);
    text += getCityName(cityIdx[5], withCountries) + ": "+ arrayDistComb[j] + " km (" + parcScore + " " + pointText +(parcScore == 1?"":"s") + ")";
  }
  text += "</ol>";
  get("list").innerHTML = text;
  parcScore = Math.round((1-position/23)*(1-position/23)*100);
  score += parcScore;
  questionNbr++;
  get("option").innerHTML = position + 1;
  get("score5").innerHTML = parcScore + " " + pointText + (parcScore == 1?"":"s");
  get("score6").innerHTML = score;
  get("nextq2").style.display = (questionNbr>10? "none": "inline");
  get("test2").style.display = "none";
  get("test2_2").style.display = "block";
}

function ListCities(onlyList)
{
  var listCitiesHTML = "<ol>";
  var idx;
  var degminLat, signLat;
  var degminLon, signLon;
  ClearScreen();
  for (idx=0; idx<cityData.length; idx++)
  {
    signLat = 0;
    degminLat = Latitude[idx];
    if (degminLat < 0)
    {
      degminLat = -degminLat;
      signLat = -1;
    }
    signLon = 0;
    degminLon = Longitude[idx];
    if (degminLon < 0)
    {
      degminLon = -degminLon;
      signLon = -1;
    }
    listCitiesHTML += "<li>" + getCityName(idx, true) + ": " +
                      Math.floor(degminLat/60)+"&deg; " + (degminLat%60)+"' " + (signLat? southText: northText) + ", " +
                      Math.floor(degminLon/60)+"&deg; " + (degminLon%60)+"' " + (signLon? westText: eastText) + "</li>";
  }
  listCitiesHTML += "</ol>";
  get("allCities").innerHTML = listCitiesHTML;
  get("lc").style.display = (onlyList?"block": "none");
  get("fc").style.display = (onlyList?"none": "block");
  get("cityFrom").value = "";
  get("cityTo").value = "";
  get("listCities").style.display = "block";
  if (onlyList == false)
  {
    get("cityFrom").focus();            // Enter first city.
  }
}

function grayFindDistButton()
{
  var cityFrom = get("cityFrom").value;
  var cityTo = get("cityTo").value;
  var findDist = get("findDist");
  if (cityFrom.length > 0 && parseInt(cityFrom,10) >= 1 && parseInt(cityFrom,10) <= cityData.length &&
      cityTo.length > 0 && parseInt(cityTo,10) >= 1 && parseInt(cityTo,10) <= cityData.length)
  {
    findDist.disabled = false;
  }
  else
  {
    findDist.disabled = true;
  }
}

function startUp()
{
  var index, degmin, s;
  var degLat, minLat;
  var degLon, minLon;
  var data = get("cities").innerHTML;
  cityData = data.split("\n");
  data = get("direction").innerHTML;
  direction = data.split(",");
  data = get("ordinal").innerHTML;
  ordinal = data.split(",");
  data = get("parcInterp").innerHTML;
  parcInterp = data.split("\n");
  data = get("finalInterp").innerHTML;
  finalInterp = data.split("\n");
  kmText = get("km").innerHTML;
  degreeText = get("degree").innerHTML;
  currentText = get("current").innerHTML;
  finalText = get("final").innerHTML;
  pointText = get("point").innerHTML;
  northText = get("North").innerHTML;
  southText = get("South").innerHTML;
  eastText = get("East").innerHTML;
  westText = get("West").innerHTML;
  for (index=0; index<cityData.length; index++)
  {
    s = cityData[index];
    cityName[index] = s.substring(0, s.indexOf(","));
    countryName[index] = s.substring(s.indexOf(",")+1, s.indexOf(":"));
    degmin = parseInt(s.substring(s.indexOf(":")+1, s.indexOf(";")), 10);
    degLat = Math.floor(Math.abs(degmin) / 100);
    minLat = Math.abs(degmin) % 100;
    Latitude[index] = (60 * degLat + minLat) * (degmin>0? 1:-1);       // Convert degrees and minutes to minutes.
    degmin = parseInt(s.substring(s.indexOf(";")+1, s.length), 10);
    degLon = Math.floor(Math.abs(degmin) / 100);
    minLon = Math.abs(degmin) % 100;
    Longitude[index] = (60 * degLon + minLon) * (degmin>0? 1:-1);      // Convert degrees and minutes to minutes.
  }
  get("1with").onclick = function ()
  {
    withCountries = true;
    StartTestType1();
  };
  get("1without").onclick = function ()
  {
    withCountries = false;
    StartTestType1();
  };
  get("2with").onclick = function ()
  {
    withCountries = true;
    StartTestType2();
  };
  get("2without").onclick = function ()
  {
    withCountries = false;
    StartTestType2();
  };
  get("listcity").onclick = function ()
  {
    ClearScreen();
    ListCities(true);
  };
  get("finddist").onclick = function ()
  {
    ClearScreen();
    get("findDist").disabled = true;
    ListCities(false);
  };
  get("dist1").onkeydown = function (event)
  {
    var key = event.keyCode;
    if (key == 10 || key == 13)
    {
      event.preventDefault();          // Do not propagate Enter key.
      if (get("dist1").value.length > 0)
      { 
        Test1(true);                   // Second part of test 1.
      }
    }
    if (key != 8 && key != 9 && key != 37 && key != 39 && key != 45)
    {                                  // Not backspace, tab, right or left arrow or insert key.
      if (key >= 0x60 && key <= 0x69)
      {
        key -= 0x30;                   // Convert numpad key to standard digit key.
      }
      if (key < 0x30 || key > 0x39 || get("dist1").value.length == 6)
      {                                // Key is not a digit or number is too large.
        event.preventDefault();        // Do not propagate this key.
      }
    }
  };
  get("order").onkeydown = function (event)
  {
    var key = event.keyCode;
    var value = get("order").value;
    if ((key == 10 || key == 13) && value.length == 4)
    {
      event.preventDefault();          // Do not propagate Enter key.
      ShowResultsTest2(value);         // Second part of test 2.
    }
    if (key != 8 && key != 9 && key != 37 && key != 39 && key != 45)
    {                                  // Not backspace, tab, right or left arrow or insert key.
      if (key >= 0x60 && key <= 0x69)
      {
        key -= 0x30;                   // Convert numpad key to standard digit key.
      }
      if (key < 0x31 || key > 0x34 || value.indexOf(String.fromCharCode(key)) >= 0)
      {                                // Key is not a digit 1 to 4 or digit is already used.
        event.preventDefault();        // Do not propagate this key.
      }
    }
  };
  get("nextq").onclick = function ()
  {
    if (questionNbr > 0 && questionNbr < 11)
    {
      Test1(false);
    }
    else
    {
      ShowMainMenu();
    }
  };
  get("exitt").onclick = function ()
  {
    ShowMainMenu();
  };
  get("nextq2").onclick = function ()
  {
    if (questionNbr < 11)
    {
      Test2();
    }
    else
    {
      ShowMainMenu();
    }
  };
  get("exitt2").onclick = function ()
  {
    ShowMainMenu();
  };
  get("exitt3").onclick = function ()
  {
    ShowMainMenu();
  };
  get("cityFrom").onkeydown = function (event)
  {
    var key = event.keyCode;
    var value = get("cityFrom").value;
    if ((key == 10 || key == 13) && parseInt(value,10) >= 1 && parseInt(value,10) <= cityData.length)
    {
      event.preventDefault();          // Do not propagate Enter key.
      get("cityTo").focus();           // Enter second city.
    }
    if (key != 8 && key != 9 && key != 37 && key != 39 && key != 45)
    {                                  // Not backspace, tab, right or left arrow or insert key.
      if (key >= 0x60 && key <= 0x69)
      {
        key -= 0x30;                   // Convert numpad key to standard digit key.
      }
      if (key < 0x30 || key > 0x39 || value.length == 3)
      {                                // Key is not a digit or 3 digits is already used.
        event.preventDefault();        // Do not propagate this key.
      }
    }
  };
  get("cityFrom").oninput = grayFindDistButton;
  get("cityTo").onkeydown = function (event)
  {
    var key = event.keyCode;
    var value = get("cityTo").value;
    if ((key == 10 || key == 13) && parseInt(value,10) >= 1 && parseInt(value,10) <= cityData.length)
    {
      event.preventDefault();          // Do not propagate Enter key.
      questionNbr = 0;
      Test1(true);                     // Show distance between cities.
    }
    if (key != 8 && key != 9 && key != 37 && key != 39 && key != 45)
    {                                  // Not backspace, right or left arrow or insert key.
      if (key >= 0x60 && key <= 0x69)
      {
        key -= 0x30;                   // Convert numpad key to standard digit key.
      }
      if (key < 0x30 || key > 0x39 || value.length == 3)
      {                                // Key is not a digit or 3 digits is already used.
        event.preventDefault();        // Do not propagate this key.
      }
    }
  };
  get("cityTo").oninput = grayFindDistButton;
  get("findDist").onclick = function ()
  {
    questionNbr = 0;
    Test1(true);                     // Show distance between cities.
  };
}

addEventListener("load", startUp);

})(this);

