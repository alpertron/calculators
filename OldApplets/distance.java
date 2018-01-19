/* Distance between cities */
/* Written by Dario Alejandro Alpern (August 1998) */
import java.applet.*;
import java.awt.*;

public final class Distance extends Applet {

  private final float Rad = (float)Math.PI/(float)(180 * 60); /* Minutos a radianes */
  private Checkbox cb1, cb2, cb3, cb4, cb5, cb6;
  private Button btnMenuPrincipal;
  private List Localidad1, Localidad2, Localidad3;
  private Button btnGetDistance;
  private Button btnVerCiudades;
  private Button btnSiguientePregunta;
  private boolean WithCountries;
  private boolean EnMenuPrincipal = false;
  private int Pregunta = 0;
  private int Ciud1, Ciud2;
  private GridBagLayout gridbag;
  private GridBagConstraints c = new GridBagConstraints();
  private TextField textField;
  private int KmJugador, test, Puntos;
  private int[] Ciud = new int[6];
  private final String[] Rumbos = {"NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW","N"};
  private final String[] Ordinal = {"First", "Second", "Third", "Fourth", "Fifth", "Sixth", "Seventh", "Eighth", "Ninth", "Last"};
  private final String[] dataCity = {
    "Acapulco,Mexico:1651;-9954",
    "Alexandria,Egypt:3112;2957",
    "Amman,Jordan:3157;3557",
    "Amsterdam,Netherlands:5223;454",
    "Anchorage,Alaska:6113;-14954",
    "Ankara,Turkey:3956;3252",
    "Asuncion,Paraguay:-2518;-5737",
    "Atenas,Grecia:3759;2344",
    "Auckland,New Zealand:-3652;17444",
    "Baghdad,Iraq:3318;4424",
    "Bangkok,Thailand:1345;10031",
    "Barcelona,Spain:4123;211",
    "Beijing,China:3955;11623",
    "Beirut,Lebanon:3352;3532",
    "Berlin,Germany:5232;1326",
    "Bogota,Colombia:437;-7405",
    "Bombay,India:1856;7249",
    "Brasilia,Brazil:-1547;-4754",
    "Brussels,Belgium:5051;421",
    "Bucharest,Romania:4427;2606",
    "Buenos Aires,Argentina:-3437;-5823",
    "Cairo,Egypt:3004;3114",
    "Calcutta,India:2230;8818",
    "Canberra,Australia:-3519;14904",
    "Cape Town,South Africa:-3354;1825",
    "Caracas,Venezuela:1030;-6654",
    "Chicago,USA:4152;-8737",
    "Copenhaguen,Denmark:5540;1234",
    "Cordoba,Argentina:-3124;-6411",
    "Dakar,Senegal:1442;-1727",
    "Damascus,Syria:3332;3620",
    "Dublin,Ireland:5321;-616",
    "Godthaab,Greenland:6411;-5144",
    "Guadalajara,Mexico:2042;-10319",
    "Guatemala,Guatemala:1437;-9031",
    "Havana,Cuba:2306;-8222",
    "Helsinki,Finland:6010;2456",
    "Hong Kong,China:2216;11412",
    "Honolulu,Hawaii:2119;-15750",
    "Istanbul,Turkey:4100;2858",
    "Jerusalem,Israel:3148;3512",
    "Karachi,Pakistán:2451;6701",
    "Kathmandu,Nepal:2743;8519",
    "Khartoum,Sudan:1538;3232",
    "Kiev,Ukraine:5025;3034",
    "La Paz,Bolivia:-1629;-6808",
    "Lagos,Nigeria:627;326",
    "Lima,Peru:-1205;-7703",
    "Lisbon,Portugal:3844;-908",
    "London,United Kingdom:5129;-8",
    "Los Angeles,USA:3403;-11814",
    "Luxembourg,Luxembourg:4937;608",
    "Madrid,Spain:4025;-341",
    "Managua,Nicaragua:1208;-8616",
    "Mar del Plata,Argentina:-3800;-5733",
    "Mecca,Saudi Arabia:2126;3949",
    "Mendoza,Argentina:-3254;-6850",
    "Mexico,Mexico:1922;-9911",
    "Miami,USA:2546;-8012",
    "Milan,Italy:4529;911",
    "Monrovia,Liberia:618;-1047",
    "Monte Carlo,Monaco:4344;726",
    "Montevideo,Uruguay:-3454;-5609",
    "Montreal,Canada:4530;-7336",
    "Moscow,Russia:5545;3738",
    "New York,USA:4043;-7401",
    "Nicosia,Cyprus:3510;3322",
    "North Pole,No country:9000;0",
    "Oslo,Norway:5955;1046",
    "Ottawa,Canada:4525;-7543",
    "Panama,Panama:901;-7930",
    "Paris,France:4852;221",
    "Prague,Czechoslovakia:5005;1425",
    "Puebla,Mexico:1903;-9812",
    "Quito,Ecuador:-13;-7831",
    "Reykjavik,Iceland:6408;-2154",
    "Rio de Janeiro,Brazil:-2253;-4314",
    "Rome,Italy:4154;1230",
    "Rosario,Argentina:-3258;-6039",
    "Saint Petersburg,Russia:5956;3019",
    "Sao Paulo,Brazil:-2331;-4638",
    "San Francisco,USA:3746;-12225",
    "San Jose,Costa Rica:956;-8405",
    "San Salvador,El Salvador:1342;-8912",
    "Santiago,Chile:-3327;-7039",
    "Seoul,South Korea:3732;12659",
    "Shanghai,China:3114;12129",
    "Sofia,Bulgaria:4241;2320",
    "South Pole,No country:-9000;0",
    "Stockholm,Sweden:5919;1804",
    "Sydney,Australia:-3353;15112",
    "Taipei,Taiwan:2503;12131",
    "Tegucigalpa,Honduras:1406;-8712",
    "Tokyo,Japan:3542;13946",
    "Toronto,Canada:4339;-7922",
    "Tripoli,Libya:3253;1317",
    "Ulaanbaatar,Mongolia:4755;10656",
    "Ushuaia,Argentina:-5448;-6818",
    "Vancouver,Canada:4915;-12307",
    "Vienna,Austria:4812;1623",
    "Vladivostok,Russia:4308;13154",
    "Warsaw,Poland:5214;2100",
    "Washington,USA:3854;-7702",
    "Yaounde,Cameroon:352;1131",
    "Zurich,Switzerland:4723;832"
    };
  private final String[] City = new String[dataCity.length];
  private final String[] Country = new String[dataCity.length];
  private final int[] Latitud = new int[dataCity.length];
  private final int[] Longitud = new int[dataCity.length];

  private void ShowMainMenu() {
    setLayout(new GridLayout(0,1));
    setFont(new Font("Helvetica", Font.PLAIN, 20));
    add(new Label("MAIN MENU",Label.CENTER));
    add(cb1);
    add(cb2);
    add(cb3);
    add(cb4);
    add(cb5);
    add(cb6);
    add(btnMenuPrincipal);
    validate();
    EnMenuPrincipal = true;
    Pregunta = 0;
    }                        /* fin método ShowMainMenu */

  private void EraseApplet() {
    removeAll();
    EnMenuPrincipal = false;
    }                        /* fin método EraseApplet */

  public boolean handleEvent(Event e) {
    if (e.id == Event.ACTION_EVENT && e.target == btnMenuPrincipal ||
        e.key == Event.ENTER && EnMenuPrincipal) {
      if (cb1.getState()) {  /* Comenzar con test tipo 1 (con países) */
        WithCountries = true;
        StartTestTipo1();
        }
      if (cb2.getState()) {  /* Comenzar con test tipo 1 (sin países) */
        WithCountries = false;
        StartTestTipo1();
        }
      if (cb3.getState()) {  /* Comenzar con test tipo 2 (con países) */
        EraseApplet();
        WithCountries = true;
        StartTestTipo2();
        }
      if (cb4.getState()) {  /* Comenzar con test tipo 2 (sin países) */
        EraseApplet();
        WithCountries = false;
        StartTestTipo2();
        }
      if (cb5.getState()) {  /* Ver listado de ciudades */
        EraseApplet();
        VerCiudades();
        }
      if (cb6.getState()) {  /* Averiguar distancia entre ciudades */
        EraseApplet();
        GetDistance();
        }
      }
    if (e.id == Event.ACTION_EVENT && e.target == btnVerCiudades ||
        (e.key == Event.ESCAPE && EnMenuPrincipal == false)) {
      EraseApplet();
      ShowMainMenu();
      }
    if (e.id == Event.ACTION_EVENT && e.target == btnGetDistance &&
        Localidad1.getSelectedIndex() >= 0 &&
        Localidad2.getSelectedIndex() >= 0) {
      Test1(true);
      }
    if (e.id == Event.ACTION_EVENT && e.target == btnSiguientePregunta) {
      if (Pregunta > 0 && Pregunta < 11) {
        if (test == 1) {
          Test1(false);
          }
        else {
          Test2();
          }
        }
      else {
        EraseApplet();
        ShowMainMenu();
        }
      }
    if (e.target == textField && e.id == Event.KEY_PRESS) {
      if (e.key == Event.ENTER) {
        if (test == 1) {
          if (textField.getText().length() != 0) {
            KmJugador = Integer.parseInt(textField.getText());
            Test1(true);
            }
          }
        else {                    /* test 2 */
          if (textField.getText().length() == 4) {
            MostrarTablaTest2(textField.getText());
            }
          }
        }
      if (e.key >= 48 && e.key <= 57) {
        if (test == 1) {
          if (textField.getText().length() == 6) {
            return true;
            }
          }
        else {     /* test 2 */
          if (e.key > 52 || e.key == 48 || 
              textField.getText().indexOf(e.key-48+"") != -1) {
            return true;
            }          
          }
        }  
      else if (e.key >= 32) {
        return true;
        }
      }
    return super.handleEvent(e);
    }                     /* Fin método handleEvent */

  public void init() {
    String s;
    CheckboxGroup cbg;
    int i,j,GraLat,MinLat,GraLon,MinLon;

    Localidad1 = new List();
    Localidad2 = new List();
    Localidad3 = new List();
    for(i=0; i<dataCity.length; i++) {
      s = dataCity[i];
      City[i] = s.substring(0, s.indexOf(","));
      Country[i] = s.substring(s.indexOf(",")+1, s.indexOf(":"));
      j = Integer.parseInt(s.substring(s.indexOf(":")+1, s.indexOf(";")));
      GraLat = Math.abs(j) / 100;
      MinLat = Math.abs(j) % 100;
      Latitud[i] = 60 * (j/100) + (j%100);     /* Convertir a minutos */
      j = Integer.parseInt(s.substring(s.indexOf(";")+1, s.length()));
      GraLon = Math.abs(j) / 100;
      MinLon = Math.abs(j) % 100;
      Longitud[i] = 60 * (j/100) + (j%100);    /* Convertir a minutos */
      s = GetCity(i, true);
      Localidad1.addItem(s);
      Localidad2.addItem(s);
      Localidad3.addItem((i+1) + ") " + s + ": "+
         GraLat + "° " + MinLat + "' " + (Latitud[i] < 0?"S, ":"N, ") +
         GraLon + "° " + MinLon + "' " + (Latitud[i] < 0?"W":"E"));
      }
    String s1 = "Start test number 1";
    String s2 = "Start test number 2";
    String s3 = " (with countries)";
    String s4 = " (without countries)";
    cbg = new CheckboxGroup();
    cb1 = new Checkbox(s1+s3, cbg, true);
    cb2 = new Checkbox(s1+s4, cbg, false);
    cb3 = new Checkbox(s2+s3, cbg, false);
    cb4 = new Checkbox(s2+s4, cbg, false);
    cb5 = new Checkbox("Show list of cities", cbg, false);
    cb6 = new Checkbox("Get distance between cities", cbg, false);
    btnMenuPrincipal = new Button("Click here to continue");
    removeAll();
    ShowMainMenu();
    }           /* end init method*/

  private String GetCity(int NumeroCiudad, boolean WithCountries) {
    if (WithCountries) {
      return City[NumeroCiudad]+" ("+Country[NumeroCiudad]+")";
      }
    else {
      return City[NumeroCiudad];
      }
    }           /* end GetCity method */

  private void GetDistance() {
    gridbag = new GridBagLayout();
    c = new GridBagConstraints();
    btnGetDistance = new Button("Get Distance");
    setLayout(gridbag);
    c.weightx = 1.0;
    c.gridwidth = GridBagConstraints.RELATIVE;  /* Penúltimo en la fila */
    c.fill = GridBagConstraints.BOTH;
    c.weighty = 1.0;
    gridbag.setConstraints(Localidad1, c);
    add(Localidad1);
    c.gridwidth = GridBagConstraints.REMAINDER;  /* Ultimo en la fila */
    gridbag.setConstraints(Localidad2, c);
    add(Localidad2);
    c.fill = GridBagConstraints.CENTER;
    c.weighty = 0.0;
    gridbag.setConstraints(btnGetDistance, c);
    add(btnGetDistance);
    validate();   
    }             /* fin método GetDistance */

  private void StartTestTipo1() {
    test = 1;
    Pregunta = 1;
    Puntos = 0;
    Test1(false);
    }             /* fin método StartTestTipo1 */

  private void StartTestTipo2() {
    test = 2;
    Pregunta = 1;
    Puntos = 0;
    Test2();
    }             /* fin método StartTestTipo2 */

  private void VerCiudades() {
    gridbag = new GridBagLayout();
    c = new GridBagConstraints();
    btnVerCiudades = new Button("Click here to return");
    setLayout(gridbag);
    c.weightx = 1.0;
    c.gridwidth = GridBagConstraints.REMAINDER;  /* Ultimo en la fila */
    c.fill = GridBagConstraints.BOTH;
    c.weighty = 1.0;
    gridbag.setConstraints(Localidad3, c);
    add(Localidad3);
    c.fill = GridBagConstraints.CENTER;
    c.weighty = 0.0;
    gridbag.setConstraints(btnVerCiudades, c);
    add(btnVerCiudades);
    validate();   
    }              /* fin método VerCiudades */

  private void makeLabel(String Caption, Color color) {
    Label label = new Label(Caption);
    gridbag.setConstraints(label, c);
    add(label);
    label.setForeground(color);
    }              /* fin método makeLabel */

  private int Distancia(int Lat1, int Lon1, int Lat2, int Lon2) {
    float q = (float) (Math.cos(Lat1*Rad)*Math.cos(Lat2*Rad)*Math.cos((Lon1-Lon2)*Rad)+Math.sin(Lat1*Rad)*Math.sin(Lat2*Rad));
    if (q >= 1) return 0;
    if (q < -1) return 20012;
    return (int)(Math.round(6370*Math.acos(q)));
    }              /* fin método Distancia */

  private String ObtenerRumbo(int Lat1, int Lon1, int Lat2, int Lon2) {
    int Angulo;
    int AngAux;
    float Diferencia;

    if (Lon1 == Lon2) Angulo = (Lat1>Lat2)?180:0;
    else {
      if (Lat1 == 90*60 || Lat2 == (float)(-90*60)) Angulo = 180;
      else {
        if (Lat1 == -90*60 || Lat2 == 90*60) Angulo = 0;
        else {
          Diferencia = Lon2 - Lon1;
          Angulo = 90 - (int)(Math.round(180/Math.PI*Math.atan((Math.tan(Lat2*Rad)*Math.cos(Lat1*Rad)-Math.sin(Lat1*Rad)*Math.cos(Diferencia*Rad))/Math.sin(Diferencia*Rad))));
          Diferencia = (Diferencia + 360*60) % (360*60);
          if (Diferencia > 180*60) {
            Angulo = (Angulo + 180) % 360;
            }
          }
        }
      }
    AngAux = (Angulo + 348) % 360;
    return Angulo + " degree" + (Angulo==1?"":"s") + " (" + Rumbos[(int)Math.floor((float)AngAux/22.5)] + ")";
    }              /* fin método ObtenerRumbo */

  private String HorasSol(int Lat, int Lon) {
    double q;
    int MinutosSol;

    if (Math.abs(Lat) == 90*60) q = Lat;
    else q = .4338 * Math.tan(Lat*Rad);
    if (q >= 1) MinutosSol = 0;
    else if (q <= -1) MinutosSol = 1440;
    else MinutosSol = (int)Math.round(458.3662*Math.acos(q));
    return (int)(MinutosSol/60) + "h " + (MinutosSol%60) + "m";
    }              /* fin método HorasSol */

  private void Test1(boolean Completo) {
    int Distcia;
    int PuntajeParcial;
    String Interpretacion;
    String Ci1, Ci2;
    int Lat1, Lat2, Lon1, Lon2;

    EraseApplet();
    setFont(new Font("Helvetica", Font.PLAIN, 12));
    gridbag = new GridBagLayout();
    c = new GridBagConstraints();
    c.gridwidth = GridBagConstraints.REMAINDER;  /* Ultimo en la fila */
    c.weightx = 1.0;
    c.fill = GridBagConstraints.BOTH;
    setLayout(gridbag);
    if (Pregunta != 0) {          /* En test tipo 1 */
      if (Completo == false) {    /* Generate cities */
        Ciud1 = (int)Math.floor(Math.random()*dataCity.length);
        do {
          Ciud2 = (int)Math.floor(Math.random()*dataCity.length);
          } while (Ciud1 == Ciud2);
        setFont(new Font("Helvetica", Font.BOLD, 14));
        makeLabel("Your current score is: "+Puntos+" point"+(Puntos==1?"":"s"),Color.red);
        }
      }
    else {                        /* Distance between cities */
      Ciud1 = Localidad1.getSelectedIndex();
      Ciud2 = Localidad2.getSelectedIndex();
      }
    Ci1 = GetCity(Ciud1, WithCountries || Pregunta == 0);
    Ci2 = GetCity(Ciud2, WithCountries || Pregunta == 0);
    Lat1 = Latitud[Ciud1];
    Lon1 = Longitud[Ciud1];
    Lat2 = Latitud[Ciud2];
    Lon2 = Longitud[Ciud2];
    if (Pregunta != 0) {
      makeLabel(Ordinal[Pregunta - 1]+" question: What is the distance between the following cities?", Color.black);
      }
    else {
      makeLabel("The cities are:", Color.black);
      }
    makeLabel(Ci1+" and "+Ci2, Color.black);
    makeLabel(" ", Color.black);
    if (Completo == false) {
      c.weightx = 0.0;
      c.gridwidth = GridBagConstraints.RELATIVE;
      makeLabel("Distance in kilometers:", Color.black);
      c.gridwidth = GridBagConstraints.REMAINDER;
      c.gridheight = GridBagConstraints.REMAINDER;
      textField = new TextField(6);
      c.weightx = 1.0;
      c.fill = GridBagConstraints.NONE;
      gridbag.setConstraints(textField, c);
      add(textField);
      c.fill = GridBagConstraints.BOTH;
      }
    else {
      Distcia = Distancia(Lat1, Lon1, Lat2, Lon2);
      if (Pregunta != 0) {
        makeLabel("You think that the distance is "+KmJugador+" kilometer"+(KmJugador==1?"":"s")+", but the exact value is "+Distcia+" kilometer"+(Distcia==1?"":"s")+".", Color.red);
        }
      else {
        makeLabel("The distance between places is "+Distcia+" kilometer"+(Distcia==1?"":"s")+".", Color.red);
        }
      makeLabel("I can also tell you that:", Color.black);
      makeLabel("To go from "+GetCity(Ciud1,false)+" to "+GetCity(Ciud2,false)+": "+ObtenerRumbo(Lat1, Lon1, Lat2, Lon2), Color.black);
      makeLabel("To go from "+GetCity(Ciud2,false)+" to "+GetCity(Ciud1,false)+": "+ObtenerRumbo(Lat2, Lon2, Lat1, Lon1), Color.black);
      makeLabel("where: 0º = North; 90º = East; 180º = South; 270º = West.", Color.black);
      c.gridwidth = GridBagConstraints.RELATIVE;
      makeLabel("Sun visible in " + Ci1 + ":", Color.black);
      c.gridwidth = GridBagConstraints.REMAINDER;
      makeLabel("Sun visible in " + Ci2 + ":", Color.black);
      c.gridwidth = GridBagConstraints.RELATIVE;
      makeLabel("June 21st: "+HorasSol(-Lat1, Lon1), Color.black);
      c.gridwidth = GridBagConstraints.REMAINDER;
      makeLabel("June 21st: "+HorasSol(-Lat2, Lon2), Color.black);
      c.gridwidth = GridBagConstraints.RELATIVE;
      makeLabel("December 21st: "+HorasSol(Lat1, Lon1), Color.black);
      c.gridwidth = GridBagConstraints.REMAINDER;
      makeLabel("December 21st: "+HorasSol(Lat2, Lon2), Color.black);
      if (Pregunta != 0) {
        PuntajeParcial = (int)Math.round(100 - 100 * Math.abs(Math.log((double)Distcia/(double)KmJugador)));
        if (PuntajeParcial < 0) PuntajeParcial = 0;
        if (PuntajeParcial == 100) Interpretacion = "Unbelievable!";
        else if (PuntajeParcial > 94) Interpretacion = "Excellent!";
        else if (PuntajeParcial > 89) Interpretacion = "very good";
        else if (PuntajeParcial > 79) Interpretacion = "ok.";
        else if (PuntajeParcial > 69) Interpretacion = "not so bad, regular, bah";
        else if (PuntajeParcial > 59) Interpretacion = "hhhmmmm...";
        else if (PuntajeParcial > 39) Interpretacion = "poor";
        else if (PuntajeParcial > 19) Interpretacion = "Yeeech!!";
        else Interpretacion = "Puajjj!!";
        Puntos += PuntajeParcial;
        makeLabel("This time you score "+PuntajeParcial+" ("+Interpretacion+"), so your "+(Pregunta==10?"final":"actual")+" score is "+Puntos+".", Color.red);
        if (Pregunta == 10) {
          if (Puntos > 949) Interpretacion = "I should tell you something against my principles: Congratulations!";
          else if (Puntos > 899) Interpretacion = "Almost perfect. ALMOST. Next time you will perform worse.";
          else if (Puntos > 799) Interpretacion = "Apparently you know geography (or you are lucky).";
          else if (Puntos > 699) Interpretacion = "Continue playing, this game is for you";
          else if (Puntos > 599) Interpretacion = "Obviously, you need to practice (several months, ha, ha).";
          else if (Puntos > 399) Interpretacion = "It seems that you hate geography.";
          else Interpretacion = "You will perform better playing 'Doom'";
          c.gridwidth = GridBagConstraints.RELATIVE;
          makeLabel(Interpretacion, Color.red);
          }
        Pregunta++;
        }
      BotonSiguientePregunta();
      }
    validate();
    if (Completo == false) {
      textField.requestFocus();
      }
    else {
      btnSiguientePregunta.requestFocus();
      }
    }          /* fin método Test1 */

  private void BotonSiguientePregunta() {
    c.gridwidth = GridBagConstraints.REMAINDER;
    c.fill = GridBagConstraints.NONE;
    btnSiguientePregunta = new Button((Pregunta>10 || Pregunta==0)?"Main menu":"Next question");
    gridbag.setConstraints(btnSiguientePregunta, c);
    add(btnSiguientePregunta);
    }          /* fin método BotonSiguientePregunta */

  private void Test2() {
    int j, k, azar;

    EraseApplet();
    setFont(new Font("Helvetica", Font.PLAIN, 12));
    gridbag = new GridBagLayout();
    c = new GridBagConstraints();
    c.gridwidth = GridBagConstraints.REMAINDER;  /* Ultimo en la fila */
    c.weightx = 1.0;
    c.fill = GridBagConstraints.BOTH;
    setLayout(gridbag);
    for (j=0; j<6; j++) {
      do {
        azar = (int)Math.floor(Math.random()*dataCity.length);
        for (k=0; k<j; k++) {
          if (Ciud[k] == azar) break;
          }
        } while (k != j);
      Ciud[j] = azar;
      }
    setFont(new Font("Helvetica", Font.BOLD, 14));
    makeLabel("Your current score is: "+Puntos+" point"+(Puntos==1?"":"s"),Color.red);
    makeLabel(Ordinal[Pregunta - 1]+" question: How do you travel from", Color.black);
    makeLabel(GetCity(Ciud[0], WithCountries) + " to " + GetCity(Ciud[5], WithCountries), Color.black);
    makeLabel("if you must pass through the following four cities so the total distance is minimum?", Color.black);
    makeLabel(" ", Color.black);
    for (j=1; j<5; j++) {
      makeLabel(j + ") " + GetCity(Ciud[j], WithCountries), Color.black);
      }
    makeLabel(" ", Color.black);
    c.weightx = 0.0;
    c.gridwidth = GridBagConstraints.RELATIVE;
    makeLabel("Ej.: If you think the order is 1, 2, 3, 4, please type 1234", Color.black);
    c.gridwidth = GridBagConstraints.REMAINDER;
    c.gridheight = GridBagConstraints.REMAINDER;
    textField = new TextField(6);
    c.weightx = 1.0;
    c.fill = GridBagConstraints.NONE;
    gridbag.setConstraints(textField, c);
    add(textField);
    c.fill = GridBagConstraints.BOTH;
    validate();
    textField.requestFocus();
    }          /* fin método Test2 */

  private void MostrarTablaTest2(String orden) {
    int lugar, i, j, k, u, v;
    String combinacion;
    int[] arrayDist = new int[25];
    int[] arrayDistComb = new int[24];
    String[] arrayComb = new String[24];
    List distancias;
    String temp;

    EraseApplet();
    setFont(new Font("Helvetica", Font.PLAIN, 10));
    gridbag = new GridBagLayout();
    c = new GridBagConstraints();
    c.gridwidth = GridBagConstraints.REMAINDER;  /* Ultimo en la fila */
    c.weightx = 1.0;
    c.weighty = 1.0;
    c.fill = GridBagConstraints.BOTH;
    setLayout(gridbag);
    i = 0;
    for (j=0; j<5; j++) {
      for (k=1; k<=5; k++) {
        arrayDist[i++] = Distancia(Latitud[Ciud[j]], Longitud[Ciud[j]], Latitud[Ciud[k]], Longitud[Ciud[k]]);
        }
      }
    for (j=0; j<24; j++) {
      arrayComb[j] = "123412431324134214231432213421432314234124132431312431423214324134123421412341324213423143124321".substring(j*4,j*4+4);
      combinacion = "0" + arrayComb[j] + "5";
      arrayDistComb[j] = 0;
      for (k=0; k<5; k++) {
        u = Integer.parseInt(combinacion.substring(k,k+1));
        v = Integer.parseInt(combinacion.substring(k+1,k+2));
        arrayDistComb[j] += arrayDist[u*5+v-1];
        }
      for (k=j; k>0; k--) {
        if (arrayDistComb[k] > arrayDistComb[k-1]) break;
        i = arrayDistComb[k];
        arrayDistComb[k] = arrayDistComb[k-1];
        arrayDistComb[k-1] = i;
        temp = arrayComb[k];
        arrayComb[k] = arrayComb[k-1];
        arrayComb[k-1] = temp;
        }
      }
    distancias = new List();
    lugar = -1;
    setFont(new Font("Courier", Font.PLAIN, 12));
    for (j=0; j<24; j++) {
      combinacion = arrayComb[j];
      if (combinacion.compareTo(orden) == 0) {
        lugar = j;
        }
      distancias.addItem((lugar==j?"-> ":"   ")+(j+1)+") "+GetCity(Ciud[0], WithCountries)+" -> "+ GetCity(Ciud[Integer.parseInt(combinacion.substring(0,1))], WithCountries)+" -> " + GetCity(Ciud[Integer.parseInt(combinacion.substring(1,2))], WithCountries) + " -> ");
      distancias.addItem("   " + GetCity(Ciud[Integer.parseInt(combinacion.substring(2,3))], WithCountries)+" -> "+ GetCity(Ciud[Integer.parseInt(combinacion.substring(3,4))], WithCountries)+" -> "+ GetCity(Ciud[5], WithCountries)+": "+ arrayDistComb[j] + " km (" + Math.round(100-100*j/23) + " points)");
      if (j<23)  distancias.addItem(" ");
      }
    gridbag.setConstraints(distancias, c);
    add(distancias);
    c.weighty = 0.0;
    i = Math.round(100-100*lugar/23);
    Puntos += i;
    makeLabel("You have chosen the option " + (lugar+1) + " so you deserve " + i + " points.", Color.red);
    c.gridwidth = GridBagConstraints.RELATIVE;
    makeLabel("Your "+(Pregunta==10?"final":"actual")+" score is "+Puntos+".", Color.red);
    Pregunta++;
    BotonSiguientePregunta();
    validate();
    distancias.makeVisible(lugar*3);
    }          /* fin método MostrarTablaTest2 */

  }            /* fin clase Dist */
