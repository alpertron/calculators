/* Gauss primes. Written by Dario Alpern on June 15th, 2003 */

import java.applet.Applet;
import java.awt.*;

public final class gausspr extends Applet {
  gaussprCanvas gaussPrimesPicture;
  Image offScreenImage;
  int controlType = 0;
  Graphics offScreenGraphics;
  int width, height;
  int xCenter = 0, yCenter = 0, thickness = 8;
  boolean repaintGaussPrimesPicture = true;
  protected long startNbr = 1;
  Label labelPosition;
  Label labelPlus;
  Label labelI;
  TextField textCenterReal;
  TextField textCenterImag;
  gaussprFramedArea framedArea;
  gaussprControls controls;
  Image controlImage;
  Image controlImage2;

  public void init() {
    controlImage = getImage(getCodeBase(), "botones.gif");
    controlImage2 = getImage(getCodeBase(), "botones2.gif");

    setLayout(null);

    setBackground(Color.lightGray);
    framedArea = new gaussprFramedArea(this);
    labelPosition = new Label("", Label.CENTER);
    labelPlus = new Label("+", Label.CENTER);
    labelI = new Label("i", Label.CENTER);
    textCenterReal = new TextField("0", 7);
    textCenterImag = new TextField("0", 7);
    controls = new gaussprControls(this);
    width = size().width;
    height = size().height;
    framedArea.reshape(0, 0, width, height-34);
    textCenterReal.reshape(10, height-30, 89, 24);
    labelPlus.reshape(100, height-30, 19, 24);
    textCenterImag.reshape(120, height-30, 89, 24);
    labelI.reshape(210, height-30, 19, 24);
    controls.reshape(width-27*7-139-13, height-30, 27*7, 24);
    labelPosition.reshape(width-139-13, height-30, 139, 24);
    add(framedArea);
    add(labelPosition);
    add(labelPlus);
    add(labelI);
    add(textCenterReal);
    add(textCenterImag);
    add(controls);
    validate();
    height -= 34;
    offScreenImage = createImage(width, height);
    offScreenGraphics = offScreenImage.getGraphics();
    offScreenGraphics.setColor(Color.black);
    drawGaussPrimesPicture();
    framedArea.canvas.requestFocus();
    }  /* end method init */

  private void drawGaussPrimesPicture() {
    drawPartialGaussPrimesPicture(xCenter - width/2/thickness-1,
                          xCenter + width/2/thickness+1,
                          yCenter - height/2/thickness-1,
                          yCenter + height/2/thickness+1);
    }

  private void drawPartialGaussPrimesPicture(int xmin, int xmax, int ymin, int ymax) {
    int x, y;
    for (x = xmin; x <= xmax; x++) {
      for (y = ymin; y <= ymax; y++) {
        setPoint(x, y);
        }  /* end for y */
      }    /* end for x */
    if (repaintGaussPrimesPicture) {
      gaussPrimesPicture.repaint();
      }
    }      /* end method drawGaussPrimesPicture */

  private void setPoint(int x, int y) {
    int xPhysical, yPhysical, absx, absy;
    long t, delta;
    xPhysical = width/2 + (x-xCenter)*thickness;
    yPhysical = height/2 - (y-yCenter)*thickness;
    if (x==0) {
      if (y>=0) {
        offScreenGraphics.setColor((y&3)==3 && isPrime(y)?Color.green:Color.black);
        }
      else {
        offScreenGraphics.setColor((y&3)==1 && isPrime(-y)?Color.green:Color.black);
        }
      }
    else if (y==0) {
      if (x>=0) {
        offScreenGraphics.setColor((x&3)==3 && isPrime(x)?Color.green:Color.black);
        }
      else {
        offScreenGraphics.setColor((x&3)==1 && isPrime(-x)?Color.green:Color.black);
        }
      }
    else {
      delta = (long)x*(long)x+(long)y*(long)y;
      offScreenGraphics.setColor(isPrime(delta)?Color.green:Color.black);
      }
    offScreenGraphics.fillRect(xPhysical, yPhysical, thickness, thickness);
    if (thickness >= 3) {
      if (x==0) {
        offScreenGraphics.setColor(Color.white);
        offScreenGraphics.drawLine(xPhysical+thickness/2, yPhysical,
                    xPhysical+thickness/2, yPhysical+thickness-1);
        }
      else if ((x%10) == 0) {
        offScreenGraphics.setColor(Color.white);
        offScreenGraphics.drawLine(xPhysical+thickness/2, yPhysical+thickness/4,
                    xPhysical+thickness/2, yPhysical+3*thickness/4);
        }
      if (y==0) {
        offScreenGraphics.setColor(Color.white);
        offScreenGraphics.drawLine(xPhysical, yPhysical+thickness/2,
                    xPhysical+thickness-1, yPhysical+thickness/2);
        }
      else if ((y%10) == 0) {
        offScreenGraphics.setColor(Color.white);
        offScreenGraphics.drawLine(xPhysical+thickness/4, yPhysical+thickness/2,
                    xPhysical+3*thickness/4, yPhysical+thickness/2);
        }
      }
    }     /* end method setPoint */

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

    if (n==0 || n==1) {return false;}
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
      textCenterReal.setText(xCenter+"");
             /* Copiar grafico hacia la derecha */
      offScreenGraphics.copyArea(0, 0, width - 32, height, 32, 0);

             /* Rellenar parte izquierda del grafico */
      drawPartialGaussPrimesPicture(xCenter - width/2/thickness-1,
                            xCenter - (width/2-32)/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter + height/2/thickness+1);
      }

    protected void Up() {
      yCenter += 32 / thickness;
      textCenterImag.setText(yCenter+"");
            /* Copiar grafico hacia abajo */
      offScreenGraphics.copyArea(0, 0, width, height - 32, 0, 32);

            /* Rellenar parte superior del grafico */
      drawPartialGaussPrimesPicture(xCenter - width/2/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter + (height/2-32)/thickness-1,
                            yCenter + height/2/thickness+1);
      }

    protected void Down() {
      yCenter -= 32 / thickness;
      textCenterImag.setText(yCenter+"");
            /* Copiar grafico hacia arriba */
      offScreenGraphics.copyArea(0, 32, width, height - 32, 0, -32);

            /* Rellenar parte inferior del grafico */
      drawPartialGaussPrimesPicture(xCenter - width/2/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter - (height/2-32)/thickness+1);
      }

    protected void Right() {
      xCenter += 32 / thickness;
      textCenterReal.setText(xCenter+"");
            /* Copiar grafico hacia la izquierda */
      offScreenGraphics.copyArea(32, 0, width - 32, height, -32, 0);

            /* Rellenar parte derecha del grafico */
      drawPartialGaussPrimesPicture(xCenter + (width/2-32)/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter + height/2/thickness+1);
      }

    protected void UpLeft() {
      xCenter -= 32 / thickness;
      yCenter += 32 / thickness;
      textCenterReal.setText(xCenter+"");
      textCenterImag.setText(yCenter+"");
            /* Copiar grafico hacia abajo a la derecha */
      offScreenGraphics.copyArea(0, 0, width - 32, height - 32, 32, 32);
      repaintGaussPrimesPicture = false;

            /* Rellenar parte superior del grafico */
      drawPartialGaussPrimesPicture(xCenter - width/2/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter + (height/2-32)/thickness-1,
                            yCenter + height/2/thickness+1);
      repaintGaussPrimesPicture = true;

             /* Rellenar parte izquierda del grafico */
      drawPartialGaussPrimesPicture(xCenter - width/2/thickness-1,
                            xCenter - (width/2-32)/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter + (height/2-32)/thickness+1);
      }

    protected void UpRight() {
      xCenter += 32 / thickness;
      yCenter += 32 / thickness;
      textCenterReal.setText(xCenter+"");
      textCenterImag.setText(yCenter+"");
            /* Copiar grafico hacia abajo a la izquierda */
      offScreenGraphics.copyArea(32, 0, width - 32, height - 32, -32, 32);
      repaintGaussPrimesPicture = false;

            /* Rellenar parte superior del grafico */
      drawPartialGaussPrimesPicture(xCenter - width/2/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter + (height/2-32)/thickness-1,
                            yCenter + height/2/thickness+1);
      repaintGaussPrimesPicture = true;

            /* Rellenar parte derecha del grafico */
      drawPartialGaussPrimesPicture(xCenter + (width/2-32)/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter + (height/2-32)/thickness+1);
      }

    protected void DownLeft() {
      xCenter -= 32 / thickness;
      yCenter -= 32 / thickness;
      textCenterReal.setText(xCenter+"");
      textCenterImag.setText(yCenter+"");
            /* Copiar grafico hacia arriba a la derecha */
      offScreenGraphics.copyArea(0, 32, width - 32, height - 32, 32, -32);
      repaintGaussPrimesPicture = false;

            /* Rellenar parte inferior del grafico */
      drawPartialGaussPrimesPicture(xCenter - width/2/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter - (height/2-32)/thickness+1);
      repaintGaussPrimesPicture = true;

             /* Rellenar parte izquierda del grafico */
      drawPartialGaussPrimesPicture(xCenter - width/2/thickness-1,
                            xCenter - (width/2-32)/thickness+1,
                            yCenter - (height/2-32)/thickness-1,
                            yCenter + height/2/thickness+1);
      }

    protected void DownRight() {
      xCenter += 32 / thickness;
      yCenter -= 32 / thickness;
      textCenterReal.setText(xCenter+"");
      textCenterImag.setText(yCenter+"");
            /* Copiar grafico hacia arriba a la izquierda */
      offScreenGraphics.copyArea(32, 32, width - 32, height - 32, -32, -32);
      repaintGaussPrimesPicture = false;

            /* Rellenar parte inferior del grafico */
      drawPartialGaussPrimesPicture(xCenter - width/2/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter - height/2/thickness-1,
                            yCenter - (height/2-32)/thickness+1);
      repaintGaussPrimesPicture = true;

            /* Rellenar parte derecha del grafico */
      drawPartialGaussPrimesPicture(xCenter + (width/2-32)/thickness-1,
                            xCenter + width/2/thickness+1,
                            yCenter - (height/2-32)/thickness-1,
                            yCenter + height/2/thickness+1);
      }

    protected void In() {
      if (thickness < 32) {
        thickness = thickness * 2;
        drawGaussPrimesPicture();
        }
      }

    protected void Out() {
      if (thickness > 1) {
        thickness = thickness / 2;
        drawGaussPrimesPicture();
        }
      }

    private void changeCenter() {
      try {
        xCenter = Integer.parseInt(textCenterReal.getText().trim());
        } catch (Exception e) {xCenter = 0;}
      try {
        yCenter = Integer.parseInt(textCenterImag.getText().trim());
        } catch (Exception e) {yCenter = 0;}
      drawGaussPrimesPicture();
      }

    public boolean handleEvent(Event e) {
      long a, n;
      String content;

      if (e.id == Event.KEY_PRESS &&
           (e.target == textCenterReal || e.target == textCenterImag)) {
        if (e.key >= '0' && e.key <= '9') {
          content = ((TextField)(e.target)).getText();
          if (content.length() == 7) {
            if (content.charAt(0) != '-') {
              return true;
              }
            }
          else if (content.length() == 8) {
            return true;
            }
          }
        else {
          if (e.key == Event.ENTER) {
            changeCenter();
            framedArea.canvas.requestFocus();
            return true;
            }
          else if (e.key == Event.TAB) {
            if (e.target == textCenterReal) {
              textCenterImag.requestFocus();
              }
            else {
              textCenterReal.requestFocus();
              }
            }
          else if (e.key == '-') {
            content = ((TextField)(e.target)).getText();
            if (content.length() != 0) {
              return true;
              }
            }
          else if (e.key != Event.BACK_SPACE && e.key != Event.DELETE) {
            return true;
            }
          }
        }
      return super.handleEvent(e);
      }
  } /* end applet gausspr */

class gaussprFramedArea extends Panel {
   protected gaussprCanvas canvas;
   public gaussprFramedArea(gausspr controller) {
      super();

      //Set layout to one that makes its contents as big as possible.
      setLayout(new GridLayout(1,0));

      add(canvas = new gaussprCanvas(controller));
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
   }      /* end class gaussprFramedArea */

class gaussprCanvas extends Canvas {
    gausspr controller;

    public gaussprCanvas(gausspr controller) {
       super();
       this.controller = controller;
       controller.gaussPrimesPicture = this;
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
       return true;
       }

    public boolean mouseMove(Event event, int x, int y) {
       long n, r, t;
       int xLogical, yLogical;

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
       if (yLogical == 0) {
         controller.labelPosition.setText(xLogical+"");
         }
       else if (yLogical > 0) {
         if (xLogical == 0) {
           controller.labelPosition.setText(yLogical+"i");
           }
         else {
           controller.labelPosition.setText(xLogical+" + "+yLogical+"i");
           }
         }
       else {
         if (xLogical == 0) {
           controller.labelPosition.setText(yLogical+"i");
           }
         else {
           controller.labelPosition.setText(xLogical+" - "+(-yLogical)+"i");
           }
         }
       return false;
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
    }          /* end class gaussprCanvas */


class gaussprControls extends Canvas implements Runnable {
    gausspr controller;
    int buttonNbr;
    Thread actionThread;

    public gaussprControls(gausspr controller) {
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
        }
      }
    }           /* end class gaussprControls */
