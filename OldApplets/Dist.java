/* Distancia entre ciudades */
/* Hecho por Dario Alejandro Alpern en agosto de 1998 */
import java.applet.*;
import java.awt.*;

public final class Dist extends Applet {

  private final float Rad = (float)Math.PI/(float)(180 * 60); /* Minutos a radianes */
  private Checkbox cb1, cb2, cb3, cb4, cb5, cb6;
  private Button btnMenuPrincipal;
  private List Localidad1, Localidad2, Localidad3;
  private Button btnAveriguarDistancia;
  private Button btnVerCiudades;
  private Button btnSiguientePregunta;
  private boolean ConPaises;
  private boolean EnMenuPrincipal = false;
  private int Pregunta = 0;
  private int Ciud1, Ciud2;
  private GridBagLayout gridbag;
  private GridBagConstraints c = new GridBagConstraints();
  private TextField textField;
  private int KmJugador, test, Puntos;
  private int[] Ciud = new int[6];
  private final String[] Rumbos = {"NNE","NE","ENE","E","ESE","SE","SSE","S","SSO","SO","OSO","O","ONO","NO","NNO","N"};
  private final String[] Ordinal = {"Primera", "Segunda", "Tercera", "Cuarta", "Quinta", "Sexta", "Séptima", "Octava", "Novena", "Última"};
  private final String[] datosCiud = {
    "Acapulco,México:1651;-9954",
    "Alejandría,Egipto:3112;2957",
    "Amman,Jordania:3157;3557",
    "Amsterdam,Holanda:5223;454",
    "Anchorage,Alaska:6113;-14954",
    "Ankara,Turquía:3956;3252",
    "Asunción,Paraguay:-2518;-5737",
    "Atenas,Grecia:3759;2344",
    "Auckland,Nueva Zelanda:-3652;17444",
    "Bagdad,Irak:3318;4424",
    "Bangkok,Thailandia:1345;10031",
    "Barcelona,España:4123;211",
    "Beijing,China:3955;11623",
    "Beirut,Líbano:3352;3532",
    "Berlín,Alemania:5232;1326",
    "Bogotá,Colombia:437;-7405",
    "Bombay,India:1856;7249",
    "Brasilia,Brasil:-1547;-4754",
    "Bruselas,Bélgica:5051;421",
    "Bucarest,Rumania:4427;2606",
    "Buenos Aires,Argentina:-3437;-5823",
    "Calcuta,India:2230;8818",
    "Canberra,Australia:-3519;14904",
    "Caracas,Venezuela:1030;-6654",
    "Chicago,EEUU:4152;-8737",
    "Ciudad del Cabo,Sudáfrica:-3354;1825",
    "Copenhague,Dinamarca:5540;1234",
    "Córdoba,Argentina:-3124;-6411",
    "Dakar,Senegal:1442;-1727",
    "Damasco,Siria:3332;3620",
    "Dublín,Irlanda:5321;-616",
    "El Cairo,Egipto:3004;3114",
    "Estambul,Turquía:4100;2858",
    "Estocolmo,Suecia:5919;1804",
    "Godthaab,Groenlandia:6411;-5144",
    "Guadalajara,México:2042;-10319",
    "Guatemala,Guatemala:1437;-9031",
    "Helsinki,Finlandia:6010;2456",
    "Hong Kong,China:2216;11412",
    "Honolulú,Hawaii:2119;-15750",
    "Jerusalén,Israel:3148;3512",
    "Karachi,Pakistán:2451;6701",
    "Katmandú,Nepal:2743;8519",
    "Khartoum,Sudán:1538;3232",
    "Kiev,Ucrania:5025;3034",
    "La Habana,Cuba:2306;-8222",
    "La Meca,Arabia Saudí:2126;3949",
    "La Paz,Bolivia:-1629;-6808",
    "Lagos,Nigeria:627;326",
    "Lima,Perú:-1205;-7703",
    "Lisboa,Portugal:3844;-908",
    "Londres,Reino Unido:5129;-8",
    "Los Ángeles,EEUU:3403;-11814",
    "Luxemburgo,Luxemburgo:4937;608",
    "Madrid,España:4025;-341",
    "Managua,Nicaragua:1208;-8616",
    "Mar del Plata,Argentina:-3800;-5733",
    "Mendoza,Argentina:-3254;-6850",
    "México,México:1922;-9911",
    "Miami,EEUU:2546;-8012",
    "Milán,Italia:4529;911",
    "Monrovia,Liberia:618;-1047",
    "Montecarlo,Mónaco:4344;726",
    "Montevideo,Uruguay:-3454;-5609",
    "Montreal,Canadá:4530;-7336",
    "Moscú,Rusia:5545;3738",
    "Nicosia,Chipre:3510;3322",
    "Nueva York,EEUU:4043;-7401",
    "Oslo,Noruega:5955;1046",
    "Ottawa,Canadá:4525;-7543",
    "Panamá,Panamá:901;-7930",
    "París,Francia:4852;221",
    "Polo Norte,Sin país:9000;0",
    "Polo Sur,Sin país:-9000;0",
    "Praga,República Checa:5005;1425",
    "Puebla,México:1903;-9812",
    "Quito,Ecuador:-13;-7831",
    "Reykjavik,Islandia:6408;-2154",
    "Río de Janeiro,Brasil:-2253;-4314",
    "Roma,Italia:4154;1230",
    "Rosario,Argentina:-3258;-6039",
    "San Francisco,EEUU:3746;-12225",
    "San José,Costa Rica:956;-8405",
    "San Pablo,Brasil:-2331;-4638",
    "San Petersburgo,Rusia:5956;3019",
    "San Salvador,El Salvador:1342;-8912",
    "Santiago,Chile:-3327;-7039",
    "Seúl,Corea del Sur:3732;12659",
    "Shanghai,China:3114;12129",
    "Sofía,Bulgaria:4241;2320",
    "Sydney,Australia:-3353;15112",
    "Taipei,Taiwán:2503;12131",
    "Tegucigalpa,Honduras:1406;-8712",
    "Tokio,Japón:3542;13946",
    "Toronto,Canadá:4339;-7922",
    "Trípoli,Libia:3253;1317",
    "Ulan Bator,Mongolia:4755;10656",
    "Ushuaia,Argentina:-5448;-6818",
    "Vancouver,Canadá:4915;-12307",
    "Varsovia,Polonia:5214;2100",
    "Viena,Austria:4812;1623",
    "Vladivostok,Rusia:4308;13154",
    "Washington,EEUU:3854;-7702",
    "Yaoundé,Camerún:352;1131",
    "Zurich,Suiza:4723;832"
    };
  private final String[] Ciudad = new String[datosCiud.length];
  private final String[] Pais = new String[datosCiud.length];
  private final int[] Latitud = new int[datosCiud.length];
  private final int[] Longitud = new int[datosCiud.length];

  private void MostrarMenuPrincipal() {
    setLayout(new GridLayout(0,1));
    setFont(new Font("Helvetica", Font.PLAIN, 20));
    add(new Label("MENÚ PRINCIPAL",Label.CENTER));
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
    }                        /* fin método MostrarMenuPrincipal */

  private void LimpiarApplet() {
    removeAll();
    EnMenuPrincipal = false;
    }                        /* fin método LimpiarApplet */

  public boolean handleEvent(Event e) {
    if (e.id == Event.ACTION_EVENT && e.target == btnMenuPrincipal ||
        e.key == Event.ENTER && EnMenuPrincipal) {
      if (cb1.getState()) {  /* Comenzar con test tipo 1 (con países) */
        ConPaises = true;
        IniciarTestTipo1();
        }
      if (cb2.getState()) {  /* Comenzar con test tipo 1 (sin países) */
        ConPaises = false;
        IniciarTestTipo1();
        }
      if (cb3.getState()) {  /* Comenzar con test tipo 2 (con países) */
        LimpiarApplet();
        ConPaises = true;
        IniciarTestTipo2();
        }
      if (cb4.getState()) {  /* Comenzar con test tipo 2 (sin países) */
        LimpiarApplet();
        ConPaises = false;
        IniciarTestTipo2();
        }
      if (cb5.getState()) {  /* Ver listado de ciudades */
        LimpiarApplet();
        VerCiudades();
        }
      if (cb6.getState()) {  /* Averiguar distancia entre ciudades */
        LimpiarApplet();
        AveriguarDistancia();
        }
      }
    if (e.id == Event.ACTION_EVENT && e.target == btnVerCiudades ||
        (e.key == Event.ESCAPE && EnMenuPrincipal == false)) {
      LimpiarApplet();
      MostrarMenuPrincipal();
      }
    if (e.id == Event.ACTION_EVENT && e.target == btnAveriguarDistancia &&
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
        LimpiarApplet();
        MostrarMenuPrincipal();
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
    for(i=0; i<datosCiud.length; i++) {
      s = datosCiud[i];
      Ciudad[i] = s.substring(0, s.indexOf(","));
      Pais[i] = s.substring(s.indexOf(",")+1, s.indexOf(":"));
      j = Integer.parseInt(s.substring(s.indexOf(":")+1, s.indexOf(";")));
      GraLat = Math.abs(j) / 100;
      MinLat = Math.abs(j) % 100;
      Latitud[i] = 60 * (j/100) + (j%100);     /* Convertir a minutos */
      j = Integer.parseInt(s.substring(s.indexOf(";")+1, s.length()));
      GraLon = Math.abs(j) / 100;
      MinLon = Math.abs(j) % 100;
      Longitud[i] = 60 * (j/100) + (j%100);    /* Convertir a minutos */
      s = ObtenerCiudad(i, true);
      Localidad1.addItem(s);
      Localidad2.addItem(s);
      Localidad3.addItem((i+1) + ") " + s + ": "+
         GraLat + "° " + MinLat + "' lat " + (Latitud[i] < 0?"S, ":"N, ") +
         GraLon + "° " + MinLon + "' long " + (Latitud[i] < 0?"O":"E"));
      }
    String s1 = "Comenzar con el test tipo 1";
    String s2 = "Comenzar con el test tipo 2";
    String s3 = " (con países)";
    String s4 = " (sin países)";
    cbg = new CheckboxGroup();
    cb1 = new Checkbox(s1+s3, cbg, true);
    cb2 = new Checkbox(s1+s4, cbg, false);
    cb3 = new Checkbox(s2+s3, cbg, false);
    cb4 = new Checkbox(s2+s4, cbg, false);
    cb5 = new Checkbox("Ver el listado de ciudades", cbg, false);
    cb6 = new Checkbox("Averiguar la distancia entre ciudades", cbg, false);
    btnMenuPrincipal = new Button("Apriete aquí para seguir");
    removeAll();
    MostrarMenuPrincipal();
    }           /* fin método init */

  private String ObtenerCiudad(int NumeroCiudad, boolean ConPaises) {
    if (ConPaises) {
      return Ciudad[NumeroCiudad]+" ("+Pais[NumeroCiudad]+")";
      }
    else {
      return Ciudad[NumeroCiudad];
      }
    }           /* fin método ObtenerCiudad */

  private void AveriguarDistancia() {
    gridbag = new GridBagLayout();
    c = new GridBagConstraints();
    btnAveriguarDistancia = new Button("Averiguar Distancia");
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
    gridbag.setConstraints(btnAveriguarDistancia, c);
    add(btnAveriguarDistancia);
    validate();   
    }             /* fin método AveriguarDistancia */

  private void IniciarTestTipo1() {
    test = 1;
    Pregunta = 1;
    Puntos = 0;
    Test1(false);
    }             /* fin método IniciarTestTipo1 */

  private void IniciarTestTipo2() {
    test = 2;
    Pregunta = 1;
    Puntos = 0;
    Test2();
    }             /* fin método IniciarTestTipo2 */

  private void VerCiudades() {
    gridbag = new GridBagLayout();
    c = new GridBagConstraints();
    btnVerCiudades = new Button("Apriete aquí para volver");
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
    return Angulo + " grado" + (Angulo==1?"":"s") + " (" + Rumbos[(int)Math.floor((float)AngAux/22.5)] + ")";
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

    LimpiarApplet();
    setFont(new Font("Helvetica", Font.PLAIN, 12));
    gridbag = new GridBagLayout();
    c = new GridBagConstraints();
    c.gridwidth = GridBagConstraints.REMAINDER;  /* Ultimo en la fila */
    c.weightx = 1.0;
    c.fill = GridBagConstraints.BOTH;
    setLayout(gridbag);
    if (Pregunta != 0) {          /* En test tipo 1 */
      if (Completo == false) {    /* Generar las ciudades */
        Ciud1 = (int)Math.floor(Math.random()*datosCiud.length);
        do {
          Ciud2 = (int)Math.floor(Math.random()*datosCiud.length);
          } while (Ciud1 == Ciud2);
        setFont(new Font("Helvetica", Font.BOLD, 14));
        makeLabel("Su puntaje actual es: "+Puntos+" punto"+(Puntos==1?"":"s"),Color.red);
        }
      }
    else {                        /* Distancia entre dos ciudades */
      Ciud1 = Localidad1.getSelectedIndex();
      Ciud2 = Localidad2.getSelectedIndex();
      }
    Ci1 = ObtenerCiudad(Ciud1, ConPaises || Pregunta == 0);
    Ci2 = ObtenerCiudad(Ciud2, ConPaises || Pregunta == 0);
    Lat1 = Latitud[Ciud1];
    Lon1 = Longitud[Ciud1];
    Lat2 = Latitud[Ciud2];
    Lon2 = Longitud[Ciud2];
    if (Pregunta != 0) {
      makeLabel(Ordinal[Pregunta - 1]+" pregunta: ¿Cuál es la distancia entre las siguientes ciudades?", Color.black);
      }
    else {
      makeLabel("Las ciudades son:", Color.black);
      }
    makeLabel(Ci1+" y "+Ci2, Color.black);
    makeLabel(" ", Color.black);
    if (Completo == false) {
      c.weightx = 0.0;
      c.gridwidth = GridBagConstraints.RELATIVE;
      makeLabel("Distancia en kilómetros:", Color.black);
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
        makeLabel("Según ud. la distancia es de "+KmJugador+" kilómetro"+(KmJugador==1?"":"s")+", pero en realidad es de "+Distcia+" kilómetro"+(Distcia==1?"":"s")+".", Color.red);
        }
      else {
        makeLabel("La distancia entre ambos lugares es de "+Distcia+" kilómetro"+(Distcia==1?"":"s")+".", Color.red);
        }
      makeLabel("Además le informo que:", Color.black);
      makeLabel("Para ir de "+ObtenerCiudad(Ciud1,false)+" a "+ObtenerCiudad(Ciud2,false)+": "+ObtenerRumbo(Lat1, Lon1, Lat2, Lon2), Color.black);
      makeLabel("Para ir de "+ObtenerCiudad(Ciud2,false)+" a "+ObtenerCiudad(Ciud1,false)+": "+ObtenerRumbo(Lat2, Lon2, Lat1, Lon1), Color.black);
      makeLabel("donde: 0º = norte; 90º = este; 180º = sur; 270º = oeste.", Color.black);
      c.gridwidth = GridBagConstraints.RELATIVE;
      makeLabel("Sol sobre el horizonte de " + Ci1 + ":", Color.black);
      c.gridwidth = GridBagConstraints.REMAINDER;
      makeLabel("Sol sobre el horizonte de " + Ci2 + ":", Color.black);
      c.gridwidth = GridBagConstraints.RELATIVE;
      makeLabel("El 21 de junio: "+HorasSol(-Lat1, Lon1), Color.black);
      c.gridwidth = GridBagConstraints.REMAINDER;
      makeLabel("El 21 de junio: "+HorasSol(-Lat2, Lon2), Color.black);
      c.gridwidth = GridBagConstraints.RELATIVE;
      makeLabel("El 21 de diciembre: "+HorasSol(Lat1, Lon1), Color.black);
      c.gridwidth = GridBagConstraints.REMAINDER;
      makeLabel("El 21 de diciembre: "+HorasSol(Lat2, Lon2), Color.black);
      if (Pregunta != 0) {
        PuntajeParcial = (int)Math.round(100 - 100 * Math.abs(Math.log((double)Distcia/(double)KmJugador)));
        if (PuntajeParcial < 0) PuntajeParcial = 0;
        if (PuntajeParcial == 100) Interpretacion = "¡impresionante!";
        else if (PuntajeParcial > 94) Interpretacion = "¡excelente!";
        else if (PuntajeParcial > 89) Interpretacion = "¡muy bien!";
        else if (PuntajeParcial > 79) Interpretacion = "bien, che.";
        else if (PuntajeParcial > 69) Interpretacion = "No está mal: regular, bah";
        else if (PuntajeParcial > 59) Interpretacion = "mmmm...";
        else if (PuntajeParcial > 39) Interpretacion = "ay, ay, ay, ¡qué mediocre!";
        else if (PuntajeParcial > 19) Interpretacion = "¡¡Yeeech!!";
        else Interpretacion = "¡¡Puajjj!!";
        Puntos += PuntajeParcial;
        makeLabel("Su puntaje en esta ronda es "+PuntajeParcial+" ("+Interpretacion+"), por lo tanto su puntaje "+(Pregunta==10?"final":"actual")+" es "+Puntos+".", Color.red);
        if (Pregunta == 10) {
          if (Puntos > 949) Interpretacion = "Aunque vaya en contra de mis principios, lo felicito.";
          else if (Puntos > 899) Interpretacion = "Casi perfecto. CASI. La próxima vez no le irá tan bien.";
          else if (Puntos > 799) Interpretacion = "Se podría decir que sabe del tema. Pero todavía le falta.";
          else if (Puntos > 699) Interpretacion = "Siga jugando, va por el buen camino.";
          else if (Puntos > 599) Interpretacion = "Evidentemente, le falta práctica (varios meses de práctica).";
          else if (Puntos > 399) Interpretacion = "Parece que su fuerte no es la geografía.";
          else Interpretacion = "Me parece que ud. tendría que jugar al 'Doom'.";
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
    btnSiguientePregunta = new Button((Pregunta>10 || Pregunta==0)?"Menú principal":"Siguiente pregunta");
    gridbag.setConstraints(btnSiguientePregunta, c);
    add(btnSiguientePregunta);
    }          /* fin método BotonSiguientePregunta */

  private void Test2() {
    int j, k, azar;

    LimpiarApplet();
    setFont(new Font("Helvetica", Font.PLAIN, 12));
    gridbag = new GridBagLayout();
    c = new GridBagConstraints();
    c.gridwidth = GridBagConstraints.REMAINDER;  /* Ultimo en la fila */
    c.weightx = 1.0;
    c.fill = GridBagConstraints.BOTH;
    setLayout(gridbag);
    for (j=0; j<6; j++) {
      do {
        azar = (int)Math.floor(Math.random()*datosCiud.length);
        for (k=0; k<j; k++) {
          if (Ciud[k] == azar) break;
          }
        } while (k != j);
      Ciud[j] = azar;
      }
    setFont(new Font("Helvetica", Font.BOLD, 14));
    makeLabel("Su puntaje actual es: "+Puntos+" punto"+(Puntos==1?"":"s"),Color.red);
    makeLabel(Ordinal[Pregunta - 1]+" pregunta: ¿En qué orden se debe pasar por las siguientes", Color.black);
    makeLabel("cuatro ciudades para lograr la mínima distancia partiendo de", Color.black);
    makeLabel(ObtenerCiudad(Ciud[0], ConPaises) + " y llegando a " + ObtenerCiudad(Ciud[5], ConPaises) + "?", Color.black);
    makeLabel(" ", Color.black);
    for (j=1; j<5; j++) {
      makeLabel(j + ") " + ObtenerCiudad(Ciud[j], ConPaises), Color.black);
      }
    makeLabel(" ", Color.black);
    c.weightx = 0.0;
    c.gridwidth = GridBagConstraints.RELATIVE;
    makeLabel("Ej.: Si cree que el orden es 1, 2, 3, 4, digite 1234", Color.black);
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

    LimpiarApplet();
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
      distancias.addItem((lugar==j?"-> ":"   ")+(j+1)+") "+ObtenerCiudad(Ciud[0], ConPaises)+" -> "+ ObtenerCiudad(Ciud[Integer.parseInt(combinacion.substring(0,1))], ConPaises)+" -> " + ObtenerCiudad(Ciud[Integer.parseInt(combinacion.substring(1,2))], ConPaises) + " -> ");
      distancias.addItem("   " + ObtenerCiudad(Ciud[Integer.parseInt(combinacion.substring(2,3))], ConPaises)+" -> "+ ObtenerCiudad(Ciud[Integer.parseInt(combinacion.substring(3,4))], ConPaises)+" -> "+ ObtenerCiudad(Ciud[5], ConPaises)+": "+ arrayDistComb[j] + " km (" + Math.round(100-100*j/23) + " puntos)");
      if (j<23)  distancias.addItem(" ");
      }
    gridbag.setConstraints(distancias, c);
    add(distancias);
    c.weighty = 0.0;
    i = Math.round(100-100*lugar/23);
    Puntos += i;
    makeLabel("Ud. eligió la opción " + (lugar+1) + " por lo que le corresponden " + i + " puntos.", Color.red);
    c.gridwidth = GridBagConstraints.RELATIVE;
    makeLabel("Su puntaje "+(Pregunta==10?"final":"actual")+" es "+Puntos+".", Color.red);
    Pregunta++;
    BotonSiguientePregunta();
    validate();
    distancias.makeVisible(lugar*3);
    }          /* fin método MostrarTablaTest2 */

  }            /* fin clase Dist */
