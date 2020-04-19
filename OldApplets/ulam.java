/* Espiral de Ulam. Hecho por Dario Alpern el 28 de octubre de 2000 */
/* Segunda version: 14 de febrero de 2003 */

import java.applet.Applet;
import java.awt.*;

public final class ulam extends Applet {
  ulamCanvas ulamSpiral;
  Image offScreenImage;
  int controlType = 0;
  Graphics offScreenGraphics;
  int width, height;
  int xCenter = 0, yCenter = 0, thickness = 8;
  boolean repaintUlamSpiral = true;
  protected long startNbr = 1;
  Label labelPosition;
  Label labelNWSE;
  Label labelSWNE;
  TextField textCenter;
  TextField textStartNbr;
  FramedArea framedArea;
  Controls controls;
  Image controlImage;
  Image controlImage2;
  boolean showAlgebraic = false;
  boolean showDelta = false;
  static final int multiple[][] = new int[25][97];
  static final int primes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,
                               43,47,53,59,61,67,71,73,79,83,89,97};

  static {
    int i,j,k;
    for (i=0; i<primes.length; i++) {
      k = primes[i];
      for (j=k/2+1; j>=0; j--) {
        multiple[i][j*j%k] = 1;
        }
      }
    }

  public void init() {
    controlImage = getImage(getCodeBase(), "botones.gif");
    controlImage2 = getImage(getCodeBase(), "botones2.gif");

    setLayout(null);

    setBackground(Color.lightGray);
    framedArea = new FramedArea(this);
    labelPosition = new Label("");
    labelNWSE = new Label("");
    labelSWNE = new Label("");
    textCenter = new TextField(15);
    textStartNbr = new TextField("1",15);
    controls = new Controls(this);
    width = size().width;
    height = size().height;
    framedArea.reshape(0, 0, width, height-96);
    labelPosition.reshape(0, height-64, width, 24);
    labelPosition.setAlignment(Label.CENTER);
    labelNWSE.reshape(0, height-44, width, 24);
    labelNWSE.setAlignment(Label.CENTER);
    labelSWNE.reshape(0, height-24, width, 24);
    labelSWNE.setAlignment(Label.CENTER);
    textCenter.reshape(width-139-297-139-13, height-92, 139, 24);
    controls.reshape(width-297-139-13, height-92, 297, 24);
    textStartNbr.reshape(width-139-13, height-92, 139, 24);
    add(framedArea);
    add(labelPosition);
    add(labelNWSE);
    add(labelSWNE);
    add(textCenter);
    add(controls);
    add(textStartNbr);
    validate();
    height -= 96;
    offScreenImage = createImage(width, height);
    offScreenGraphics = offScreenImage.getGraphics();
    offScreenGraphics.setColor(Color.black);
    drawUlamSpiral();
    framedArea.canvas.requestFocus();
    }  /* end method init */

  private void drawUlamSpiral() {
    drawPartialUlamSpiral(xCenter - width/2/thickness-1,
                          xCenter + width/2/thickness+1,
                          yCenter - height/2/thickness-1,
                          yCenter + height/2/thickness+1);
    }

  private void drawPartialUlamSpiral(int xmin, int xmax, int ymin, int ymax) {
    int x, y;
    long n=0;
    for (x = xmin; x <= xmax; x++) {
      for (y = ymin; y <= ymax; y++) {
        if (x>=0 && x>=y && x>=-y) {
          n = 4*(long)x*(long)x+3*(long)x+(long)y+startNbr;
          }
        else {
          if (y>=0 && y>x && y>-x) {
            n = 4*(long)y*(long)y-3*(long)y-(long)x+startNbr;
            }
          else {
            if (x<=0 && x<=y && x<=-y) {
              n = 4*(long)x*(long)x+(long)x-(long)y+startNbr;
              }
            else {
              n = 4*(long)y*(long)y-(long)y+(long)x+startNbr;
              }
            }
          }
        setPoint(n, x, y);
        }  /* end for y */
      }    /* end for x */
    if (repaintUlamSpiral) {
      ulamSpiral.repaint();
      if (xCenter>=0 && xCenter>=yCenter && xCenter>=-yCenter) {
        n = 4*(long)xCenter*(long)xCenter+3*(long)xCenter+(long)yCenter+startNbr;
        }
      else {
        if (yCenter>=0 && yCenter>xCenter && yCenter>-xCenter) {
          n = 4*(long)yCenter*(long)yCenter-3*(long)yCenter-(long)xCenter+startNbr;
          }
        else {
          if (xCenter<=0 && xCenter<=yCenter && xCenter<=-yCenter) {
            n = 4*(long)xCenter*(long)xCenter+(long)xCenter-(long)yCenter+startNbr;
            }
          else {
            n = 4*(long)yCenter*(long)yCenter-(long)yCenter+(long)xCenter+startNbr;
            }
          }
        }
      textCenter.setText(String.valueOf(n));
      }
    }      /* end method drawUlamSpiral */

  private void setPoint(long n, int x, int y) {
    int xPhysical, yPhysical, absx, absy;
    Color algebraicColor = showAlgebraic? Color.blue:Color.black;
    long t;
    xPhysical = width/2 + (x-xCenter)*thickness;
    yPhysical = height/2 - (y-yCenter)*thickness;
    if (n!=2 && n%2==0) {
      offScreenGraphics.setColor(algebraicColor);
      }
    else if (isPrime(n)) {
      offScreenGraphics.setColor(Color.green);
      }
    else {
      if (x > y && x > -y) {
        t = x;
        }
      else if (x < y && x < -y) {
        t = -x;
        }
      else if (y > 0) {
        t = y;
        }
      else {
        t = -y;
        }
      if (x + y >= 0) {
        if (x - y >= 0) {    /* Right quadrant */
          if (algebraicFactor(4, 4, n-4*t*t-4*t) ||
              algebraicFactor(4, 2, n-4*t*t-2*t)) {
            offScreenGraphics.setColor(algebraicColor);
            }
          else {
            offScreenGraphics.setColor(Color.black);
            }
          }
        else {                             /* Upper quadrant */
          if (algebraicFactor(4, -4, n-4*t*t+4*t) ||
              algebraicFactor(4, -2, n-4*t*t+2*t)) {
            offScreenGraphics.setColor(algebraicColor);
            }
          else {
            offScreenGraphics.setColor(Color.black);
            }
          }
        }
      else {
        if (x - y >= 0) {    /* Lower quadrant */
          if (algebraicFactor(4, 0, n-4*t*t) ||
              algebraicFactor(4, 2, n-4*t*t-2*t)) {
            offScreenGraphics.setColor(algebraicColor);
            }
          else {
            offScreenGraphics.setColor(Color.black);
            }
          }
        else {                             /* Left quadrant */
          if (algebraicFactor(4, 0, n-4*t*t) ||
              algebraicFactor(4, -2, n-4*t*t+2*t)) {
            offScreenGraphics.setColor(algebraicColor);
            }
          else {
            offScreenGraphics.setColor(Color.black);
            }
          }
        }
      }
    if (thickness < 3) {
      offScreenGraphics.fillRect(xPhysical, yPhysical, thickness, thickness);
      }
    else {
      offScreenGraphics.fillRect(xPhysical, yPhysical, thickness, thickness);
      offScreenGraphics.setColor(Color.white);
      absx = (x>0?x:-x);
      absy = (y>0?y:-y);
      if (absx >= absy && !(x==-y && y<0)) {
        offScreenGraphics.drawLine(xPhysical, yPhysical,
                    xPhysical, yPhysical+thickness-1);
        }
      if (absx < absy && !(x==y-1 && y>0) || absx==absy && y<=0) {
        offScreenGraphics.drawLine(xPhysical, yPhysical+thickness-1,
                    xPhysical+thickness-1, yPhysical+thickness-1);
        }
      }
    }     /* end method setPoint */

  private boolean algebraicFactor(long a, long b, long c) {
    long delta = b*b - 4*a*c;
    double sqdelta = Math.sqrt(delta);
    long t1 = (long)(((double)(-b) + sqdelta) / (double)a);
    long t2 = (long)(((double)(-b) - sqdelta) / (double)a);
    return (a*t1*t1+2*b*t1+4*c == 0 && a*t2*t2+2*b*t2+4*c == 0);
                    /* Twice the roots are integer numbers */
    }

// Verificar si el n£mero es primo utilizando el m‚todo de potenciaci¢n que
// figura en http://www.utm.edu/research/primes/prove2.html (tercer par grafo)
// de Chris K. Caldwell:
//
// Si tenemos que verificar el n£mero impar "n", hacemos n-1=(2^s)d donde "d"
// es impar y s>0: decimos que n es a-SPRP si a^d=1 (mod n) o bien
// (a^d)^(2^r)=-1 para alg£n valor r tal que 0<=r<s.
//
// Entonces:
// - Si n < 9.080.191 es 31 y 73-SPRP entonces n es primo.
// - Si n < 4.759.123.141 es 2, 7 y 61-SPRP entonces n es primo.
// - Si n < 1.000.000.000.000 es 2, 13, 23 y 1662803-SPRP entonces n es primo.
// - Si n < 2.152.302.898.747 es 2, 3, 5, 7 y 11-SPRP entonces n es primo.
// - Si n < 3.474.749.660.383 es 2, 3, 5, 7, 11 y 13-SPRP entonces n es primo.
// - Si n < 341.550.071.728.321 es 2, 3, 5, 7, 11, 13 y 17-SPRP entonces n es
//   primo.

  private boolean isPrime(long n) {
    long div;

    if (n==1) {return false;}
    if (n%2==0) {if (n==2) {return true;} else {return false;}}
    if (n%3==0) {if (n==3) {return true;} else {return false;}}
    if (n%5==0) {if (n==5) {return true;} else {return false;}}
    if (n%7==0) {if (n==7) {return true;} else {return false;}}
    if (n%11==0) {if (n==11) {return true;} else {return false;}}
    if (n%13==0) {if (n==13) {return true;} else {return false;}}
    if (n%17==0) {if (n==17) {return true;} else {return false;}}
    if (n%19==0) {if (n==19) {return true;} else {return false;}}
    if (n%23==0) {if (n==23) {return true;} else {return false;}}
    if (n%29==0) {if (n==29) {return true;} else {return false;}}
    if (n%31==0) {if (n==31) {return true;} else {return false;}}
    if (n%37==0) {if (n==37) {return true;} else {return false;}}
    if (n%41==0) {if (n==41) {return true;} else {return false;}}
    if (n%43==0) {if (n==43) {return true;} else {return false;}}
    if (n%47==0) {if (n==47) {return true;} else {return false;}}
    if (n%53==0) {if (n==53) {return true;} else {return false;}}
    if (n%59==0) {if (n==59) {return true;} else {return false;}}
    if (n%61==0) {if (n==61) {return true;} else {return false;}}
    if (n%67==0) {if (n==67) {return true;} else {return false;}}
    if (n%71==0) {if (n==71) {return true;} else {return false;}}
    if (n%73==0) {if (n==73) {return true;} else {return false;}}
    if (n%79==0) {if (n==79) {return true;} else {return false;}}
    if (n%83==0) {if (n==83) {return true;} else {return false;}}
    if (n%89==0) {if (n==89) {return true;} else {return false;}}
    if (n%97==0) {if (n==97) {return true;} else {return false;}}
    if (n < 10000) {return true;}
    if (n < 9080191) {
      if (StrongPseudoprime(n, 31)==false) {return false;}
      if (StrongPseudoprime(n, 73)==false) {return false;}
      return true;
      }
    if (StrongPseudoprime(n, 2)==false) {return false;}
    if (n < (long)475912314*(long)10) {
      if (StrongPseudoprime(n, 7)==false) {return false;}
      if (StrongPseudoprime(n, 61)==false) {return false;}
      return true;
      }
    if (n < (long)100000000*(long)10000) {
      if (StrongPseudoprime(n, 13)==false) {return false;}
      if (StrongPseudoprime(n, 23)==false) {return false;}
      if (StrongPseudoprime(n, 1662803)==false) {return false;}
      return true;
      }
    if (StrongPseudoprime(n, 3)==false) {return false;}
    if (StrongPseudoprime(n, 5)==false) {return false;}
    if (StrongPseudoprime(n, 7)==false) {return false;}
    if (StrongPseudoprime(n, 11)==false) {return false;}
    if (n < (long)215230289*(long)10000) {return true;}
    if (StrongPseudoprime(n, 13)==false) {return false;}
    if (n < (long)347474966*(long)10000) {return true;}
    if (StrongPseudoprime(n, 17)==false) {return false;}
    return true;
    }           /* end method isPrime */

  private boolean StrongPseudoprime(long n, long base) {
    long exp, prod;
    int r, s;
    exp = n-1;
    s = 0;
    while ((exp & 1) == 0) {
      exp >>= 1;
      s++;
      }
    prod = 1;
    if (n < 2000000000) {
      while (true) {
        if ((exp & 1) == 1) {
          prod = (base*prod)%n;
          }
        exp >>= 1;
        if (exp == 0) {break;}
        base = (base*base)%n;
        }
      if (prod == 1) {return true;}
      for (r=0; r<s; r++) {
        if (prod == n-1) {return true;}
        prod = (prod*prod)%n;
        }
      return false;
      }
    while (true) {
      if ((exp & 1) == 1) {
        prod = base*prod-n*(long)(((double)base*(double)prod)/(double)n);
        if (prod>=n) {prod-=n;}
        if (prod<0) {prod+=n;}
        }
      exp >>= 1;
      if (exp == 0) {break;}
      base = base*base-n*(long)(((double)base*(double)base)/(double)n);
      if (base>=n) {base-=n;}
      if (base<0) {base+=n;}
      }
    if (prod == 1) {return true;}
    for (r=0; r<s; r++) {
      if (prod == n-1) {return true;}
      prod = prod*prod-n*(long)(((double)prod*(double)prod)/(double)n);
      if (prod>=n) {prod-=n;}
      if (prod<0) {prod+=n;}
      }
    return false;
    }   /* end method StrongPseudoprime */

    protected void ChangeDir() {
      controlType = 1-controlType;
      controls.repaint();
      }

    protected void Left() {
      xCenter -= 32 / thickness;
             /* Copiar grafico hacia la derecha */
      offScreenGraphics.copyArea(0, 0, width - 32, height, 32, 0);

             /* Rellenar parte izquierda del grafico */
      drawPartialUlamSpiral(xCenter - width/2/thickness-1,
                            xCenter - (width/2-32)/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter + height/2/thickness+1);
      }

    protected void Up() {
      yCenter += 32 / thickness;
            /* Copiar grafico hacia abajo */
      offScreenGraphics.copyArea(0, 0, width, height - 32, 0, 32);

            /* Rellenar parte superior del grafico */
      drawPartialUlamSpiral(xCenter - width/2/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter + (height/2-32)/thickness-1,
                            yCenter + height/2/thickness+1);
      }

    protected void Down() {
      yCenter -= 32 / thickness;
            /* Copiar grafico hacia arriba */
      offScreenGraphics.copyArea(0, 32, width, height - 32, 0, -32);

            /* Rellenar parte inferior del grafico */
      drawPartialUlamSpiral(xCenter - width/2/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter - (height/2-32)/thickness+1);
      }

    protected void Right() {
      xCenter += 32 / thickness;
            /* Copiar grafico hacia la izquierda */
      offScreenGraphics.copyArea(32, 0, width - 32, height, -32, 0);

            /* Rellenar parte derecha del grafico */
      drawPartialUlamSpiral(xCenter + (width/2-32)/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter + height/2/thickness+1);
      }

    protected void UpLeft() {
      xCenter -= 32 / thickness;
      yCenter += 32 / thickness;
            /* Copiar grafico hacia abajo a la derecha */
      offScreenGraphics.copyArea(0, 0, width - 32, height - 32, 32, 32);
      repaintUlamSpiral = false;

            /* Rellenar parte superior del grafico */
      drawPartialUlamSpiral(xCenter - width/2/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter + (height/2-32)/thickness-1,
                            yCenter + height/2/thickness+1);
      repaintUlamSpiral = true;

             /* Rellenar parte izquierda del grafico */
      drawPartialUlamSpiral(xCenter - width/2/thickness-1,
                            xCenter - (width/2-32)/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter + (height/2-32)/thickness+1);
      }

    protected void UpRight() {
      xCenter += 32 / thickness;
      yCenter += 32 / thickness;
            /* Copiar grafico hacia abajo a la izquierda */
      offScreenGraphics.copyArea(32, 0, width - 32, height - 32, -32, 32);
      repaintUlamSpiral = false;

            /* Rellenar parte superior del grafico */
      drawPartialUlamSpiral(xCenter - width/2/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter + (height/2-32)/thickness-1,
                            yCenter + height/2/thickness+1);
      repaintUlamSpiral = true;

            /* Rellenar parte derecha del grafico */
      drawPartialUlamSpiral(xCenter + (width/2-32)/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter + (height/2-32)/thickness+1);
      }

    protected void DownLeft() {
      xCenter -= 32 / thickness;
      yCenter -= 32 / thickness;
            /* Copiar grafico hacia arriba a la derecha */
      offScreenGraphics.copyArea(0, 32, width - 32, height - 32, 32, -32);
      repaintUlamSpiral = false;

            /* Rellenar parte inferior del grafico */
      drawPartialUlamSpiral(xCenter - width/2/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter - (height/2-32)/thickness+1);
      repaintUlamSpiral = true;

             /* Rellenar parte izquierda del grafico */
      drawPartialUlamSpiral(xCenter - width/2/thickness-1,
                            xCenter - (width/2-32)/thickness+1,
                            yCenter - (height/2-32)/thickness-1,
                            yCenter + height/2/thickness+1);
      }

    protected void DownRight() {
      xCenter += 32 / thickness;
      yCenter -= 32 / thickness;
            /* Copiar grafico hacia arriba a la izquierda */
      offScreenGraphics.copyArea(32, 32, width - 32, height - 32, -32, -32);
      repaintUlamSpiral = false;

            /* Rellenar parte inferior del grafico */
      drawPartialUlamSpiral(xCenter - width/2/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter - (height/2-32)/thickness+1);
      repaintUlamSpiral = true;

            /* Rellenar parte derecha del grafico */
      drawPartialUlamSpiral(xCenter + (width/2-32)/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter - (height/2-32)/thickness-1,
                            yCenter + height/2/thickness+1);
      }

    protected void In() {
      if (thickness < 32) {
        thickness = thickness * 2;
        drawUlamSpiral();
        }
      }

    protected void Out() {
      if (thickness > 1) {
        thickness = thickness / 2;
        drawUlamSpiral();
        }
      }

    protected void Delta() {
      if (showDelta) showDelta = false;
      else showDelta = true;
      return;
      }

    protected void Algebraic() {
      if (showAlgebraic) showAlgebraic = false;
      else showAlgebraic = true;
      drawUlamSpiral();
      return;
      }

    protected void Increment() {
      if (startNbr < 9999999999999L) {
        startNbr++;
        textStartNbr.setText(String.valueOf(startNbr));
        long n = Long.parseLong(textCenter.getText());
        if (n < startNbr) {
          textCenter.setText(String.valueOf(startNbr));
          }  
        changeCenter();
        }
      }

    protected void Decrement() {
      if (startNbr > 1) {
        startNbr--;
        textStartNbr.setText(String.valueOf(startNbr));
        changeCenter();
        }
      }

    private void changeCenter() {
      long n = Long.parseLong(textCenter.getText()) - startNbr + 1;
      long a = ((long)Math.sqrt((double)(n-1))+1)/2;
      if (n < 4*a*a-2*a+2) {
        xCenter = (int)(4*a*a - 3*a + 1 - n);
        yCenter = (int)a;
        }
      else if (n <= 4*a*a) {
        xCenter = (int)(-a);
        yCenter = (int)(4*a*a - a + 1 - n);
        }
      else if (n <= 4*a*a+2*a) {
        xCenter = (int)(n - 4*a*a - a - 1);
        yCenter = (int)(-a);
        }
      else {
        xCenter = (int)a;
        yCenter = (int)(n - 4*a*a - 3*a - 1);
        }
      drawUlamSpiral();
      }

    public boolean handleEvent(Event e) {
      long a, n;

      if (e.id == Event.KEY_PRESS &&
           (e.target == textCenter || e.target == textStartNbr)) {
        if (e.key >= '0' && e.key <= '9') {
          if (((TextField)(e.target)).getText().length() == 14) {
            return true;
            }
          }
        else {
          if (e.key == Event.ENTER) {
            if (e.target == textCenter) {
              if (textCenter.getText().length() == 0) {
                n = 1;
                }
              else {
                n = Long.parseLong(textCenter.getText());
                if (n < startNbr) {
                  n = startNbr;
                  }
                }
              }
            else if (e.target == textStartNbr) {
              n = Long.parseLong(textCenter.getText());
              if (textStartNbr.getText().length() == 0) {
                startNbr = 1;
                }
              else {
                startNbr = Long.parseLong(textStartNbr.getText());
                if (n < startNbr) {
                  textCenter.setText(String.valueOf(startNbr));
                  }
                }
              }
            else {
              return true;
              }
            changeCenter();
            return true;
            }
          if (e.key != Event.BACK_SPACE && e.key != Event.DELETE) {
            return true;
            }
          }
        }
      return super.handleEvent(e);
      }
  } /* end applet ulam */

class FramedArea extends Panel {
   protected ulamCanvas canvas;
   public FramedArea(ulam controller) {
      super();

      //Set layout to one that makes its contents as big as possible.
      setLayout(new GridLayout(1,0));

      add(canvas = new ulamCanvas(controller));
      validate();
      }

   public Insets insets() {
      return new Insets(4,4,5,5);
      }

   public void paint(Graphics g) {
      Dimension d = size();
      Color bg = getBackground();
 
      g.setColor(bg);
      g.draw3DRect(0, 0, d.width - 1, d.height - 1, true);
      g.draw3DRect(3, 3, d.width - 7, d.height - 7, false);
      }
   }      /* end class FramedArea */

class ulamCanvas extends Canvas {
    ulam controller;

    public ulamCanvas(ulam controller) {
       super();
       this.controller = controller;
       controller.ulamSpiral = this;
       }

    public boolean keyDown(Event event, int key) {
       if (key == Event.DOWN) {
          controller.Down();
          return true;
          }
       else if (key == Event.UP) {
          controller.Up();
          return true;
          }
       else if (key == Event.LEFT) {
          controller.Left();
          return true;
          }
       else if (key == Event.RIGHT) {
          controller.Right();
          return true;
          }
       return false;
       }

    public boolean mouseExit(Event event, int x, int y) {
       controller.labelPosition.setText("");
       controller.labelSWNE.setText("");
       controller.labelNWSE.setText("");
       return true;
       }

    public boolean mouseMove(Event event, int x, int y) {
       long n, r, t;
       int xLogical, yLogical;
       String info;

       xLogical = x-controller.width/2;
       if (xLogical>=0) {
         xLogical = xLogical / controller.thickness;
         }
       else {
         xLogical = (xLogical-controller.thickness+1) / controller.thickness;
         }
       xLogical += controller.xCenter;
       yLogical = controller.height/2-y;
       if (yLogical>=0) {
         yLogical = yLogical / controller.thickness;
         }
       else {
         yLogical = (yLogical-controller.thickness+1) / controller.thickness;
         }
       yLogical += controller.yCenter+1;
       if (xLogical>=0 && xLogical>=yLogical && xLogical>=-yLogical) {
         n = 4*(long)xLogical*(long)xLogical+3*(long)xLogical+(long)yLogical+
             controller.startNbr;
         }
       else {
         if (yLogical>=0 && yLogical>xLogical && yLogical>-xLogical) {
           n = 4*(long)yLogical*(long)yLogical-3*(long)yLogical-(long)xLogical+
               controller.startNbr;
           }
         else {
           if (xLogical<=0 && xLogical<=yLogical && xLogical<=-yLogical) {
             n = 4*(long)xLogical*(long)xLogical+(long)xLogical-(long)yLogical+
                 controller.startNbr;
             }
           else {
             n = 4*(long)yLogical*(long)yLogical-(long)yLogical+(long)xLogical+
                 controller.startNbr;
             }
           }
         }
       if (xLogical > yLogical && xLogical > -yLogical) {
         t = xLogical;
         }
       else if (xLogical < yLogical && xLogical < -yLogical) {
         t = -xLogical;
         }
       else if (yLogical > 0) {
         t = yLogical;
         }
       else {
         t = -yLogical;
         }
       if (xLogical + yLogical >= 0) {
         if (xLogical - yLogical >= 0) {    /* Right quadrant */
           info = "SW-NE: 4t^2 + 4t";
           ShowLabel(controller.labelSWNE, info, 4, 4, n-4*t*t-4*t);
           info = "NW-SE: 4t^2 + 2t";
           ShowLabel(controller.labelNWSE, info, 4, 2, n-4*t*t-2*t);
           }
         else {                             /* Upper quadrant */
           info = "SW-NE: 4t^2 - 4t";
           ShowLabel(controller.labelSWNE, info, 4, -4, n-4*t*t+4*t);
           info = "NW-SE: 4t^2 - 2t";
           ShowLabel(controller.labelNWSE, info, 4, -2, n-4*t*t+2*t);
           }
         }
       else {
         if (xLogical - yLogical >= 0) {    /* Lower quadrant */
           info = "SW-NE: 4t^2";
           ShowLabel(controller.labelSWNE, info, 4, 0, n-4*t*t);
           info = "NW-SE: 4t^2 + 2t";
           ShowLabel(controller.labelNWSE, info, 4, 2, n-4*t*t-2*t);
           }
         else {                             /* Left quadrant */
           info = "SW-NE: 4t^2";
           ShowLabel(controller.labelSWNE, info, 4, 0, n-4*t*t);
           info = "NW-SE: 4t^2 - 2t";
           ShowLabel(controller.labelNWSE, info, 4, -2, n-4*t*t+2*t);
           }
         }

       info += " (t="+t+")";
       controller.labelPosition.setText("x="+xLogical+", y="+yLogical+", n="+n+", t="+t);
       return false;
       }   

    private void ShowLabel(Label label, String text, long a, long b, long c) {
       int i,p;
       long delta, t1, t2;
       double sqdelta;
       StringBuffer text2 = new StringBuffer(text);
       boolean firstTime = true;
       if (c > 0) {
         text2 = text2.append(" + ").append(c);
         }
       else if (c < 0) {
         text2 = text2.append(" - ").append(-c);
         }
       delta = b*b - 4*a*c;
       if (controller.showDelta) {
         text2 = text2.append(" (delta=").append(delta).append(")");
         }
       if (b != 0 || c != 0) {
         sqdelta = Math.sqrt(delta);
         t1 = (long)(((double)(-b) + sqdelta) / (double)a);
         t2 = (long)(((double)(-b) - sqdelta) / (double)a);
         if (a*t1*t1+2*b*t1+4*c == 0 && a*t2*t2+2*b*t2+4*c == 0) {
                    /* Twice the roots are integer numbers */
           text2 = text2.append(" = ").append(show(t1)).
                         append(" ").append(show(t2));
           }
         else {
           if (c%2 == 0) {
             if (b%4 == 0 && c%4 == 0) {
               text2 = text2.append(" = 4 (2t^2");
               a/=4; b/=4; c/=4;
               }
             else {
               text2 = text2.append(" = 2 (t^2");
               a/=2; b/=2; c/=2;
               }
             if (b != 0) {
               if (b<0) {
                 text2 = text2.append(" - ");
                 if (b != -1) {
                   text2 = text2.append(-b);
                   }
                 }
               else {
                 text2 = text2.append(" + ");
                 if (b != 1) {
                   text2 = text2.append(b);
                   }
                 }
               text2 = text2.append("t");
               }
             if (c<0) {
               text2 = text2.append(" - ").append(-c);
               }
             else {
               text2 = text2.append(" + ").append(c);
               }
             text2 = text2.append(")");
             }
           else {
             for (i=0; i<ulam.primes.length; i++) {
               p = ulam.primes[i];
               if (ulam.multiple[i][(int)(delta%p+p)%p] == 0) {
                 if (firstTime) {
                   firstTime = false;
                   text2 = text2.append(" (").append(p);
                   }
                 else {
                   text2 = text2.append(", ").append(p);
                   }
                 }
               }         /* end for */
             if (firstTime == false) {
               text2.append(")");
               }
             }           /* end if */
           }             /* end if */
         }               /* end if */
       label.setText(text2.toString());
       }

    private String show(long root) {
      if (root > 0) {
        return "(2t - "+root+")";
        }
      if (root < 0) {
        return "(2t + "+(-root)+")";
        }
      return "2t";
      }

    public void update(Graphics g) {
       paint(g);
       }

    public void paint(Graphics g) {
       if (controller.offScreenImage == null) {
         g.setColor(Color.black);
         g.fillRect(0, 0, controller.width, controller.height);
         }
       else {
         g.drawImage(controller.offScreenImage, 0, 0, this);
         }
       }
    }          /* end class ulamCanvas */


class Controls extends Canvas implements Runnable {
    ulam controller;
    int buttonNbr;
    Thread actionThread;

    public Controls(ulam controller) {
       super();
       this.controller = controller;
       }

    public void update(Graphics g) {
       paint(g);
       }

    public void paint (Graphics g) {
       if (controller.controlType == 0) {
         g.drawImage(controller.controlImage, 0, 0, this);
         }
       else {
         g.drawImage(controller.controlImage2, 0, 0, this);
         }
       }

    public boolean mouseExit(Event event, int x, int y) {
       if (buttonNbr != 0) {
         buttonNbr = 0;
         actionThread.interrupt();
         }
       controller.framedArea.canvas.requestFocus();
       return true;
       }

    public boolean mouseUp(Event event, int x, int y) {
       if (buttonNbr != 0) {
         buttonNbr = 0;
         actionThread.interrupt();
         }
       controller.framedArea.canvas.requestFocus();
       return true;
       }

    public boolean lostFocus(Event evt, Object what) {
       if (buttonNbr != 0) {
         buttonNbr = 0;
         actionThread.interrupt();
         }
       controller.framedArea.canvas.requestFocus();
       return true;
       }

    public boolean mouseDown(Event event, int x, int y) {
       if (y>1 && y<22 && x%27>1 && x%27<24) {    /* Inside button */
         buttonPressed(x/27+1);
         return true;
         }
       buttonNbr = 0;
       return false;
       }

    public void buttonPressed(int buttonNbr) {
       if (actionThread != null && actionThread.isAlive()) {
         return;
         }
       if (this.buttonNbr == 0) {
         this.buttonNbr = buttonNbr;
         actionThread = new Thread(this);  /* Start button action thread */
         actionThread.start();
         }
       else {
         this.buttonNbr = buttonNbr;
         }
       }

    public void run() {
      executeAction(buttonNbr);
      try {
        Thread.sleep(1000);
        } catch (InterruptedException ie) {};
      while (buttonNbr != 0) {
        executeAction(buttonNbr);
        try {
          Thread.sleep(100);
          } catch (InterruptedException ie) {};
        }
      }

    private void executeAction(int buttonNbr) {
      switch (buttonNbr) {
        case 1:
           if (controller.controlType == 0) {
             controller.Left();
             }
           else {
             controller.UpLeft();
             }
           break;
        case 2:
           if (controller.controlType == 0) {
             controller.Up();
             }
           else {
             controller.DownRight();
             }
           break;
        case 3:
           if (controller.controlType == 0) {
             controller.Down();
             }
           else {
             controller.UpRight();
             }
           break;
        case 4:
           if (controller.controlType == 0) {
             controller.Right();
             }
           else {
             controller.DownLeft();
             }
           break;
        case 5:
           controller.ChangeDir();
           break;
        case 6:
           controller.In();
           break;
        case 7:
           controller.Out();
           break;
        case 8:
           controller.Delta();
           break;
        case 9:
           controller.Algebraic();
           break;
        case 10:
           controller.Increment();
           break;
        case 11:
           controller.Decrement();
           break;
        }
      }
    }           /* end class Controls */
