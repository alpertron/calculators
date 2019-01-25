/*
    GAUSSIAN PRIMES
    Copyright(C) 2018 Dario Alejandro Alpern

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Gaussian Primes is distributed in the hope that it will be useful,
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
var buffer, globals, env, asm;
var zoom, zoomDone, imgData;
var canvas, zoomin, zoomout, centerX, centerY;
var isMouseDown;
var currentX, currentY, pixels;
var HEAP8, HEAPU8;
var asmjs = typeof(WebAssembly) === "undefined";
var ShowInformation;
var prevX1stTouch, prevY1stTouch;
var prevX2ndTouch, prevY2ndTouch;
var beforeMinus = "";
var bitsCanvas;
// Webassembly module encoded in Base64.
var wasm="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iHwiB0IfiEL/////D4MgADUCBCIJIAV+fCAIIAdC/////weDIAkgBn58IganbEH/////B3GtIgcgCn58IAcgC34gBnxCH4h8IgunQf////8HcSEBIAtCH4inIQACQCALIANBAWqsQh+GWg0AIAMgAEcNAiABIARJDQILIAAgA2sgASAEayIBQR91akH/////B3EhACABQf////8HcSEBDAELIAEoAgAhASAAKAIAIQACQAJAQQAoAhAiA0H//wFKDQAgASAAbCADbyEBDAELIAGsIACsfiADrIGnIQELQQAhAAsgAiAANgIEIAIgATYCAAvrCAkBfwN+DH8BfgF/AX4CfwF+BX9BAEEAKAIEQRBrIh42AgQCQEEALQAkDQBBAEEBOgAkC0EAKAIsIhsgAWpBACgCDCIYdUEAKAIoIgFqrCEHIBsgAGogGHUgAWoiF6whFkEAKAIwIgBBAEEAKAI0ayIBIANrIBh1a6whBiAAIAEgAmsgGHVrIgSsIQUgHkEIakEEaiEVAkADQCAWIAdVDQEgFqciCEEAIBdrIBZCf1UbIglBA3EhCiAIQQpvIQ0gFiAWfiIZQh+IpyELIBmnIQwgBSEZIAQhGgJAA0AgGSAGVQ0BQQAoArgEQQJtIBmnIg9BACgCMGsgGHRrQQAoAjRqIRAgCEEAKAIoayAYdEEAKAK0BEECbWpBACgCLGshDgJAAkACQCAIRQ0AIA9FDQEgFSAZIBl+IhRCH4inIAtqQf////8HcTYCACAeIBSnIAxqQf////8HcTYCCEGAgIN4QYCAgHggHkEIahABGyEDDAILIBVBADYCACAeIA9BACAaayAZQn9VGyIANgIIQYCAgHghAyAAQQNxQQNHDQFBgICDeEGAgIB4IB5BCGoQARshAwwBCyAVQQA2AgAgHiAJNgIIQYCAgHghAyAKQQNHDQBBgICDeEGAgIB4IB5BCGoQARshAwtBACgCuAQiEkEBQQAoAgwiGHQiACAQaiIBIAEgEkobIRNBACgCtAQiESAAIA5qIgAgACARShshAiAOQQAgDkEAShshHSAQQQAgEEEAShsiHCEbAkADQCAbIBNODQEgG0ELdCAdakECdEGABWohASAdIQACQANAIAAgAk4NASABIAM2AgAgAEEBaiEAIAFBBGohAQwACwsgG0EBaiEbDAALCwJAIBhBAkgNAAJAAkAgCEUNACANDQFBASAYQX9qdCIAIA5qIgFBAEgNASABIBFODQEgEkEBIBhBfmp0IBBqIhsgAGoiACAAIBJKGyEDIBtBACAbQQBKGyIAQQt0IAFqQQJ0QYAFaiEBA0AgACADTg0CIAFBwIGDfjYCACAAQQFqIQAgAUGAwABqIQEMAAsLQQEgGEF/anQgDmoiAEEASA0AIAAgEU4NACAAIBxBC3RqQQJ0QYAFaiEAA0AgHCATTg0BIABBwIGDfjYCACAcQQFqIRwgAEGAwABqIQAMAAsLAkAgD0UNACAPQQpvDQFBASAYQX9qdCIAIBBqIgFBAEgNASABIBJODQEgEUEBIBhBfmp0IA5qIgMgAGoiACAAIBFKGyECIANBACADQQBKGyIAIAFBC3RqQQJ0QYAFaiEBA0AgACACTg0CIAFBwIGDfjYCACAAQQFqIQAgAUEEaiEBDAALC0EBIBhBf2p0IBBqIgBBAEgNACAAIBJODQAgAEELdCAdakECdEGABWohAANAIB0gAk4NASAAQcCBg342AgAgHUEBaiEdIABBBGohAAwACwsgGkEBaiEaIBlCAXwhGQwACwsgF0EBaiEXIBZCAXwhFgwACwtBACAeQRBqNgIEC4kBAQN/AkAgAUF/Sg0AIABBLToAACAAQQFqIQBBACABayEBC0EAIQRBgJTr3AMhAwNAAkAgASADbSICIARyRQ0AIAAgAkEwajoAAEEBIQQgAEEBaiEAIAEgAiADbGshAQsgA0ETSiECIANBCm0hAyACDQALIABBADoAASAAIAFBMGo6AAAgAEEBagvUAQEEf0EAQQA6AEACQCAAQQBIDQBBACgCMCEDQQAoAjQhBEEAKAK4BCEFQcAAQQAoAiwgAGpBACgCtARBfm1qQQAoAgwiAnVBACgCKGoQBSIAQSA6AAACQAJAIAMgBCABayAFQQJtaiACdWpBAWoiAUEASA0AIABBAWpBq8AAOwAAIABBA2ogARAFIQAMAQsgAEEBakGtwAA7AAAgAEEDakEAIAFrEAUhAAsgAEHpADsAAAtBwARBACgCKBAFGkHgBEEAKAIwEAUaQcAAQcAEQeAEEAALVQECf0EAQQAoAihBACgCLCAAayICQQAoAgwiAHVqNgIoQQBBASAAdEF/aiIDIAJxNgIsQQBBACgCNCABaiIBIANxNgI0QQBBACgCMCABIAB1ajYCMAvaAQBBACADNgK4BEEAIAI2ArQEAkACQAJAIAFBAUcNAEEAIAAQCTYCKEEAQQA2AiwgAEF/aiEBA0AgAUEBaiIBLQAADQALQQAhAkEAIAFBAWoQCTYCMAwBC0EAKAIMIQICQCABQQNHDQAgAkEFRg0CQQAgAkEBajYCDEEAQQAoAixBAXQ2AixBACgCNEEBdCECDAELIAJFDQFBACACQX9qNgIMQQBBACgCLEEBdTYCLEEAKAI0QQF1IQILQQAgAjYCNEHABEEAKAIoEAUaQeAEQQAoAjAQBRoLQQALTwEDfyAAQQFqIAAgAC0AACIBQS1GGyEAQQAhAwJAA0AgACwAACICRQ0BIABBAWohACADQQpsIAJqQVBqIQMMAAsLQQAgA2sgAyABQS1GGwsfAQF/IABBf2ohAQNAIAFBAWoiAS0AAA0ACyABIABrCwUAQYAFCwueAQUAQQQLBGBFDwEAQQwLBAMAAAAAQaCFgAgLGQIDBQcLDRETFx0fJSkrLzU7PUNHSU9TWWEAQcCFgAgLUAAAAAAAAAAA/wcAAAAAAADV9RQAAAAAALFxggEAAAAAx32hPwEAAAA77j8f6gMAAN8cOAdSBgAAwciyUkZtAgDByLJSRm0CAAAALHY1SERTAEGQhoAICwoCAwUHCw0RExcA";
var initializer = [3,0,0,0,0,0,0,0,0,0,0,0,255,7,0,0,0,0,0,0,213,245,20,0,0,0,0,0,177,113,130,1,0,0,0,0,199,125,161,63,1,0,0,0,59,238,63,31,234,3,0,0,223,28,56,7,82,6,0,0,193,200,178,82,70,109,2,0,193,200,178,82,70,109,2,0,0,0,44,118,53,72,68,83,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,2,3,5,7,11,13,17,19,23,0];
var myAsmJs = (function(global, env, buffer) {
 "use asm";
 var a = new global.Int8Array(buffer);
 var c = new global.Int32Array(buffer);
 var l = env.STACKTOP | 0;
 var D = global.Math.floor;
 var R = global.Math.imul;
 var _ = env._showInfo;
 function ha(b) {
  b = b | 0;
  var d = 0, e = 0, f = 0, g = 0, h = 0, i = 0, j = 0, k = 0, m = 0, n = 0, o = 0, p = 0, q = 0, r = 0, s = 0, t = 0, u = 0, v = 0, w = 0, x = 0, y = 0, z = 0, A = 0;
  y = l;
  l = l + 32 | 0;
  t = y + 16 | 0;
  u = y + 8 | 0;
  v = y;
  w = c[b >> 2] | 0;
  c[32] = w;
  s = c[b + 4 >> 2] | 0;
  c[33] = s;
  e = (s | 0) == 0;
  a : do if (e) {
   switch (w | 0) {
   case 1:
    {
     b = 0;
     break a;
    }
   case 2:
    break;
   default:
    {
     x = 4;
     break a;
    }
   }
   b = 1;
  } else x = 4; while (0);
  b : do if ((x | 0) == 4) if (!(w & 1)) b = 0; else {
   b = 1;
   do {
    d = a[92 + b >> 0] | 0;
    if (e) {
     if ((w | 0) == (d | 0)) {
      b = 1;
      break b;
     }
     if (!((w >>> 0) % (d >>> 0) | 0)) {
      b = 0;
      break b;
     }
    } else if (!((((R(2147483648 % (d >>> 0) | 0, (s >>> 0) % (d >>> 0) | 0) | 0) + ((w >>> 0) % (d >>> 0) | 0) | 0) >>> 0) % (d >>> 0) | 0)) {
     b = 0;
     break b;
    }
    b = b + 1 | 0;
   } while (b >>> 0 < 25);
   ga();
   c[t >> 2] = c[35];
   r = t + 4 | 0;
   c[r >> 2] = c[36];
   q = w & 2147483646;
   m = (q | 0) == 0;
   p = m ? 31 : 0;
   q = m ? s : q;
   m = (q & 65535 | 0) == 0;
   p = m ? p + 16 | 0 : p;
   q = m ? q >>> 16 : q;
   m = (q & 255 | 0) == 0;
   p = m ? p + 8 | 0 : p;
   q = m ? q >>> 8 : q;
   m = (q & 15 | 0) == 0;
   p = m ? p + 4 | 0 : p;
   q = m ? q >>> 4 : q;
   m = (q & 3 | 0) == 0;
   p = ((m ? q >>> 2 : q) & 1 ^ 1) + (m ? p + 2 | 0 : p) | 0;
   m = e ? 0 : 31;
   q = (e ^ 1) & 1;
   k = e ? w : s;
   j = k >>> 0 > 65535;
   m = j ? m + 16 | 0 : m;
   k = j ? k >>> 16 : k;
   j = (k & 65280 | 0) == 0;
   m = j ? m : m + 8 | 0;
   k = j ? k : k >>> 8;
   j = (k & 240 | 0) == 0;
   m = j ? m : m + 4 | 0;
   k = j ? k : k >>> 4;
   j = (k & 12 | 0) == 0;
   m = ((j ? k : k >>> 2) >>> 1 & 1) + (j ? m : m + 2 | 0) | 0;
   j = 1 << ((m | 0) % 31 | 0);
   k = u + 4 | 0;
   h = m + -1 | 0;
   m = (m | 0) > (p | 0);
   n = v + 4 | 0;
   o = p + -1 | 0;
   d = 1;
   g = 0;
   i = 0;
   while (1) {
    b = c[12 + ((g | 1) << 2) >> 2] | 0;
    if ((b | 0) >= (s | 0)) {
     if ((b | 0) != (s | 0)) {
      b = 1;
      break b;
     }
     if ((c[12 + (g << 2) >> 2] | 0) >= (w | 0)) {
      b = 1;
      break b;
     }
    }
    b = a[117 + i >> 0] | 0;
    do {
     fa(t, 140, t);
     d = d + 1 | 0;
    } while (d >>> 0 < b >>> 0);
    e = c[t >> 2] | 0;
    c[u >> 2] = e;
    c[k >> 2] = c[r >> 2];
    if (m) {
     b = q;
     e = h;
     f = j;
     while (1) {
      z = f >>> 1;
      A = (z | 0) == 0;
      b = (A << 31 >> 31) + b | 0;
      f = A ? 1073741824 : z;
      ca(u, u, u);
      if (c[128 + (b << 2) >> 2] & f | 0) ca(u, t, u);
      if ((e | 0) <= (p | 0)) break; else e = e + -1 | 0;
     }
     b = o;
     e = c[u >> 2] | 0;
    } else b = h;
    if (!((e | 0) == (c[35] | 0) ? (c[k >> 2] | 0) == (c[36] | 0) : 0)) x = 23;
    c : do if ((x | 0) == 23 ? (x = 0, fa(u, 140, v), c[v >> 2] | c[n >> 2] | 0) : 0) {
     if ((b | 0) <= -1) {
      b = 0;
      break b;
     }
     while (1) {
      ca(u, u, u);
      if ((c[u >> 2] | 0) == (c[35] | 0) ? (c[k >> 2] | 0) == (c[36] | 0) : 0) {
       b = 0;
       break b;
      }
      fa(u, 140, v);
      if (!(c[v >> 2] | c[n >> 2])) break c;
      if ((b | 0) > 0) b = b + -1 | 0; else {
       b = 0;
       break b;
      }
     }
    } while (0);
    g = g + 2 | 0;
    i = i + 1 | 0;
   }
  } while (0);
  l = y;
  return b | 0;
 }
 function ja(a, b) {
  a = a | 0;
  b = b | 0;
  var d = 0, e = 0, f = 0, g = 0, h = 0, i = 0, j = 0, k = 0, m = 0, n = 0, o = 0, p = 0, q = 0, r = 0, s = 0, t = 0, u = 0, v = 0, w = 0, x = 0, y = 0, z = 0, A = 0;
  A = l;
  l = l + 32 | 0;
  d = A + 16 | 0;
  e = A + 8 | 0;
  f = A;
  w = c[2] | 0;
  v = (a - (c[39] | 0) << w) + ((c[38] | 0) / 2 | 0) - (c[40] | 0) | 0;
  w = ((c[41] | 0) / 2 | 0) - (b - (c[42] | 0) << w) + (c[43] | 0) | 0;
  m = (a | 0) == 0;
  do if (m) {
   s = (b | 0) > -1 ? b : 0 - b | 0;
   c[f >> 2] = s;
   c[f + 4 >> 2] = 0;
   if ((s & 3 | 0) == 3) {
    i = (ha(f) | 0) == 0;
    i = i ? -16777216 : -16728064;
   } else i = -16777216;
  } else {
   if (b | 0) {
    ia(a, a, d);
    ia(b, b, e);
    ea(d, e, f);
    i = (ha(f) | 0) == 0;
    i = i ? -16777216 : -16728064;
    break;
   }
   s = (a | 0) > -1 ? a : 0 - a | 0;
   c[f >> 2] = s;
   c[f + 4 >> 2] = 0;
   if ((s & 3 | 0) == 3) {
    i = (ha(f) | 0) == 0;
    i = i ? -16777216 : -16728064;
   } else i = -16777216;
  } while (0);
  r = (v | 0) > 0 ? v : 0;
  d = (w | 0) > 0 ? w : 0;
  e = c[2] | 0;
  k = 1 << e;
  s = k + v | 0;
  j = c[38] | 0;
  s = (s | 0) > (j | 0) ? j : s;
  k = k + w | 0;
  j = c[41] | 0;
  k = (k | 0) > (j | 0) ? j : k;
  j = (d | 0) < (k | 0);
  if (j) {
   h = (r | 0) < (s | 0);
   g = d;
   do {
    if (h) {
     e = r;
     f = 262144 + (g << 11 << 2) + (r << 2) | 0;
     while (1) {
      c[f >> 2] = i;
      e = e + 1 | 0;
      if ((e | 0) >= (s | 0)) break; else f = f + 4 | 0;
     }
    }
    g = g + 1 | 0;
   } while ((g | 0) < (k | 0));
   e = c[2] | 0;
  }
  a : do if ((e | 0) > 1) {
   if (m) {
    e = (1 << e + -1) + v | 0;
    if (!((e | 0) >= (c[38] | 0) | (e | 0) < 0 | j ^ 1)) {
     e = 262144 + (d << 11 << 2) + (e << 2) | 0;
     while (1) {
      c[e >> 2] = -4144960;
      d = d + 1 | 0;
      if ((d | 0) >= (k | 0)) break; else e = e + 8192 | 0;
     }
    }
   } else if ((((a | 0) % 10 | 0 | 0) == 0 ? (n = 1 << e + -1, o = n + v | 0, (o | 0) > -1 & (o | 0) < (c[38] | 0)) : 0) ? (p = (1 << e + -2) + w | 0, q = p + n | 0, p = (p | 0) > 0 ? p : 0, n = c[41] | 0, q = (q | 0) > (n | 0) ? n : q, (p | 0) < (q | 0)) : 0) {
    e = 262144 + (p << 11 << 2) + (o << 2) | 0;
    d = p;
    while (1) {
     c[e >> 2] = -4144960;
     d = d + 1 | 0;
     if ((d | 0) >= (q | 0)) break; else e = e + 8192 | 0;
    }
   }
   if (!b) {
    d = (1 << (c[2] | 0) + -1) + w | 0;
    if (!((d | 0) > -1 & (d | 0) < (c[41] | 0) & (r | 0) < (s | 0))) break;
    e = r;
    d = 262144 + (d << 11 << 2) + (r << 2) | 0;
    while (1) {
     c[d >> 2] = -4144960;
     e = e + 1 | 0;
     if ((e | 0) >= (s | 0)) break a; else d = d + 4 | 0;
    }
   }
   if ((((b | 0) % 10 | 0 | 0) == 0 ? (t = c[2] | 0, u = 1 << t + -1, x = u + w | 0, (x | 0) > -1 & (x | 0) < (c[41] | 0)) : 0) ? (y = (1 << t + -2) + v | 0, z = y + u | 0, y = (y | 0) > 0 ? y : 0, w = c[38] | 0, z = (z | 0) > (w | 0) ? w : z, (y | 0) < (z | 0)) : 0) {
    e = y;
    d = 262144 + (x << 11 << 2) + (y << 2) | 0;
    while (1) {
     c[d >> 2] = -4144960;
     e = e + 1 | 0;
     if ((e | 0) >= (z | 0)) break; else d = d + 4 | 0;
    }
   }
  } while (0);
  l = A;
  return;
 }
 function ca(a, b, d) {
  a = a | 0;
  b = b | 0;
  d = d | 0;
  var e = 0, f = 0, g = 0, h = 0, i = 0, j = 0.0, k = 0.0, l = 0.0, m = 0.0, n = 0.0, o = 0, p = 0, q = 0.0, r = 0.0;
  g = c[33] | 0;
  e = c[a >> 2] | 0;
  f = c[b >> 2] | 0;
  h = c[32] | 0;
  if (!g) {
   da(e, f, d, h);
   a = 0;
  } else {
   o = c[b + 4 >> 2] | 0;
   j = +(e | 0);
   p = R(f, e) | 0;
   q = +(f | 0);
   b = c[34] | 0;
   p = (R(p, b) | 0) & 2147483647;
   l = +(p | 0);
   n = +(h | 0);
   r = +D(+((j * q + n * l) * 4.656612873077393e-10 + .5));
   p = (~~r >>> 0) + (R(o, e) | 0) + (R(p, g) | 0) & 2147483647;
   k = +(g | 0);
   m = +(o | 0);
   i = ~~((r + (j * m + k * l) + (p >>> 0 < 1073741824 ? 536870912.0 : -536870912.0)) * 4.656612873077393e-10) >>> 0;
   e = c[a + 4 >> 2] | 0;
   l = +(e | 0);
   b = (R(p + (R(e, f) | 0) | 0, b) | 0) & 2147483647;
   j = +(b | 0);
   n = +D(+((q * l + +(p >>> 0) + n * j) * 4.656612873077393e-10 + .5));
   b = (~~n >>> 0) + (R(e, o) | 0) + (R(b, g) | 0) + i | 0;
   e = b & 2147483647;
   a = ~~((n + (m * l + k * j + +(i >>> 0)) + (e >>> 0 < 1073741824 ? 536870912.0 : -536870912.0)) * 4.656612873077393e-10) >>> 0;
   if (a >>> 0 <= g >>> 0 ? e >>> 0 < h >>> 0 | (a | 0) != (g | 0) : 0) b = e; else {
    b = b - h & 2147483647;
    a = (e - h >> 31) - g + a & 2147483647;
   }
   c[d >> 2] = b;
  }
  c[d + 4 >> 2] = a;
  return;
 }
 function ba(a) {
  a = a | 0;
  var b = 0, d = 0, e = 0.0, f = 0, g = 0, h = 0, i = 0, j = 0.0, k = 0;
  f = c[32] | 0;
  k = a + 8 | 0;
  d = c[a >> 2] | 0;
  h = ~~+D(+((+(c[k >> 2] | 0) * 2147483648.0 + +(c[a + 4 >> 2] | 0) + +(d | 0) * 4.656612873077393e-10) / (+(c[33] | 0) + +(f | 0) * 4.656612873077393e-10) + .5)) >>> 0;
  h = (h | 0) < 0 ? 2147483647 : h;
  j = +(h | 0);
  i = 0;
  b = 0;
  e = 0.0;
  while (1) {
   g = d + b - (R(f, h) | 0) & 2147483647;
   e = e + (+(b | 0) + (+(d | 0) - j * +(f | 0)));
   d = e < 0.0;
   b = ~~+D(+(((g >>> 0 < 1073741824 ? 536870912.0 : -536870912.0) + (d ? e + 4611686018427388000.0 : e)) * 4.656612873077393e-10));
   c[a + (i << 2) >> 2] = g;
   g = i + 1 | 0;
   if ((g | 0) == 2) break;
   i = g;
   e = d ? -2147483648.0 : 0.0;
   f = c[128 + (g << 2) >> 2] | 0;
   d = c[a + (g << 2) >> 2] | 0;
  }
  i = (c[k >> 2] | 0) + b & 2147483647;
  c[k >> 2] = i;
  if (i | 0) {
   b = 0;
   d = 0;
   while (1) {
    i = a + (d << 2) | 0;
    b = (c[i >> 2] | 0) + b + (c[128 + (d << 2) >> 2] | 0) | 0;
    c[i >> 2] = b & 2147483647;
    d = d + 1 | 0;
    if ((d | 0) == 2) break; else b = b >>> 31;
   }
   c[k >> 2] = 0;
  }
  return;
 }
 function ma(b, d) {
  b = b | 0;
  d = d | 0;
  var e = 0, f = 0;
  a[177] = 0;
  if ((b | 0) > -1) {
   e = c[2] | 0;
   f = (c[42] | 0) + 1 + ((c[43] | 0) - d + ((c[41] | 0) / 2 | 0) >> e) | 0;
   d = la(177, ((c[40] | 0) + b + ((c[38] | 0) / -2 | 0) >> e) + (c[39] | 0) | 0) | 0;
   b = d + 1 | 0;
   a[d >> 0] = 32;
   e = d + 2 | 0;
   if ((f | 0) > -1) {
    a[b >> 0] = 43;
    a[e >> 0] = 32;
    d = la(d + 3 | 0, f) | 0;
   } else {
    a[b >> 0] = 45;
    a[e >> 0] = 32;
    d = la(d + 3 | 0, 0 - f | 0) | 0;
   }
   a[d >> 0] = 105;
   a[d + 1 >> 0] = 0;
  }
  la(677, c[39] | 0) | 0;
  la(707, c[42] | 0) | 0;
  _(177, 677, 707);
  return;
 }
 function oa(a, b, d, e) {
  a = a | 0;
  b = b | 0;
  d = d | 0;
  e = e | 0;
  var f = 0;
  c[38] = d;
  c[41] = e;
  do if ((b | 0) != 1) {
   a = c[2] | 0;
   if ((b | 0) == 3) {
    if ((a | 0) == 5) break;
    c[2] = a + 1;
    c[40] = c[40] << 1;
    a = c[43] << 1;
    f = 8;
    break;
   } else {
    if (!a) break;
    c[2] = a + -1;
    c[40] = c[40] >> 1;
    a = c[43] >> 1;
    f = 8;
    break;
   }
  } else {
   c[39] = pa(a) | 0;
   c[40] = 0;
   c[42] = pa(a + ((qa(a) | 0) + 1) | 0) | 0;
   a = 0;
   f = 8;
  } while (0);
  if ((f | 0) == 8) {
   c[43] = a;
   la(677, c[39] | 0) | 0;
   la(707, c[42] | 0) | 0;
  }
  return 0;
 }
 function ka(b, d, e, f) {
  b = b | 0;
  d = d | 0;
  e = e | 0;
  f = f | 0;
  var g = 0, h = 0, i = 0, j = 0;
  if (!(a[176] | 0)) aa();
  h = c[39] | 0;
  j = c[40] | 0;
  i = c[2] | 0;
  b = (j + b >> i) + h | 0;
  h = (j + d >> i) + h | 0;
  j = c[42] | 0;
  d = 0 - (c[43] | 0) | 0;
  g = j - (d - e >> i) | 0;
  e = j - (d - f >> i) | 0;
  if ((b | 0) <= (h | 0)) {
   f = (g | 0) > (e | 0);
   d = b;
   while (1) {
    if (!f) {
     b = g;
     while (1) {
      ja(d, b);
      if ((b | 0) < (e | 0)) b = b + 1 | 0; else break;
     }
    }
    if ((d | 0) < (h | 0)) d = d + 1 | 0; else break;
   }
  }
  return;
 }
 function la(b, c) {
  b = b | 0;
  c = c | 0;
  var d = 0, e = 0, f = 0;
  if ((c | 0) < 0) {
   a[b >> 0] = 45;
   e = 0;
   b = b + 1 | 0;
   f = 1e9;
   c = 0 - c | 0;
  } else {
   e = 0;
   f = 1e9;
  }
  while (1) {
   d = (c | 0) / (f | 0) | 0;
   if (!(d | e)) d = e; else {
    c = c - (R(d, f) | 0) | 0;
    a[b >> 0] = d + 48;
    d = 1;
    b = b + 1 | 0;
   }
   if ((f | 0) > 19) {
    e = d;
    f = (f | 0) / 10 | 0;
   } else break;
  }
  f = b + 1 | 0;
  a[b >> 0] = c + 48;
  a[f >> 0] = 0;
  return f | 0;
 }
 function fa(a, b, d) {
  a = a | 0;
  b = b | 0;
  d = d | 0;
  var e = 0, f = 0, g = 0, h = 0;
  f = c[32] | 0;
  g = c[33] | 0;
  e = (c[b >> 2] | 0) + (c[a >> 2] | 0) | 0;
  h = e & 2147483647;
  a = (e >>> 31) + (c[a + 4 >> 2] | 0) + (c[b + 4 >> 2] | 0) | 0;
  if (a >>> 0 <= g >>> 0 ? (h | 0) < (f | 0) | (a | 0) != (g | 0) : 0) b = h; else {
   b = e - f & 2147483647;
   a = a - g + (h - f >> 31) | 0;
  }
  c[d >> 2] = b;
  c[d + 4 >> 2] = a & 2147483647;
  return;
 }
 function pa(b) {
  b = b | 0;
  var c = 0, d = 0, e = 0;
  e = (a[b >> 0] | 0) == 45;
  b = e ? b + 1 | 0 : b;
  c = a[b >> 0] | 0;
  if (!(c << 24 >> 24)) b = 0; else {
   d = b;
   b = 0;
   do {
    d = d + 1 | 0;
    b = (b * 10 | 0) + -48 + (c << 24 >> 24) | 0;
    c = a[d >> 0] | 0;
   } while (c << 24 >> 24 != 0);
  }
  return (e ? 0 - b | 0 : b) | 0;
 }
 function ga() {
  var a = 0, b = 0;
  c[36] = 0;
  if (!(c[33] | 0)) c[35] = 1; else {
   b = c[32] | 0;
   a = R(2 - (R(b, b) | 0) | 0, b) | 0;
   a = R(2 - (R(a, b) | 0) | 0, a) | 0;
   a = R(2 - (R(a, b) | 0) | 0, a) | 0;
   c[34] = (R(2 - (R(a, b) | 0) | 0, 0 - a | 0) | 0) & 2147483647;
   c[37] = 1;
   c[35] = 0;
   ba(140);
  }
  return;
 }
 function da(a, b, d, e) {
  a = a | 0;
  b = b | 0;
  d = d | 0;
  e = e | 0;
  var f = 0;
  f = R(b, a) | 0;
  if ((e | 0) < 32768) a = (f | 0) % (e | 0) | 0; else {
   a = f - (R(~~+D(+(+(a | 0) * +(b | 0) / +(e | 0) + .5)), e) | 0) | 0;
   a = ((a | 0) < 0 ? e : 0) + a | 0;
  }
  c[d >> 2] = a;
  return;
 }
 function na(a, b) {
  a = a | 0;
  b = b | 0;
  var d = 0, e = 0;
  e = (c[40] | 0) - a | 0;
  d = c[2] | 0;
  c[39] = (c[39] | 0) + (e >> d);
  a = (1 << d) + -1 | 0;
  c[40] = e & a;
  b = (c[43] | 0) + b | 0;
  c[42] = (c[42] | 0) + (b >> d);
  c[43] = b & a;
  return;
 }
 function ia(a, b, d) {
  a = a | 0;
  b = b | 0;
  d = d | 0;
  var e = 0;
  e = (R(b, a) | 0) & 2147483647;
  c[d >> 2] = e;
  c[d + 4 >> 2] = ~~((+(a | 0) * +(b | 0) + (e >>> 0 < 1073741824 ? 536870912.0 : -536870912.0)) * 4.656612873077393e-10) >>> 0;
  return;
 }
 function ea(a, b, d) {
  a = a | 0;
  b = b | 0;
  d = d | 0;
  var e = 0;
  e = (c[b >> 2] | 0) + (c[a >> 2] | 0) | 0;
  c[d >> 2] = e & 2147483647;
  c[d + 4 >> 2] = (e >>> 31) + (c[a + 4 >> 2] | 0) + (c[b + 4 >> 2] | 0) & 2147483647;
  return;
 }
 function qa(b) {
  b = b | 0;
  var c = 0;
  c = b;
  while (1) if (!(a[c >> 0] | 0)) break; else c = c + 1 | 0;
  return c - b | 0;
 }
 function aa() {
  a[176] = 1;
  return;
 }
 return {
  _moveGraphic: na,
  _nbrChanged: oa,
  _ShowInformation: ma,
  _drawPartialGraphic: ka
 };
});

function get(x)
{
  return document.getElementById(x);
}

function PtrToString(ptr)
{
  var t=-1;
  var i = 0;
  var str="", outString="";
  do
  {
    for (i=0; i<1024; i++)
    {
      t = HEAP8[((ptr++)>>0)];
      if (t==0)
      {
        break;
      }
      if (t>=128)
      {
        t = ((t-192)<<6) + HEAP8[((ptr++)>>0)] - 128;
      }
      str += String.fromCharCode(t);
    }
    outString += str;
    str = "";
  } while (t!=0);
  outString += str;
  return outString;
}

function getHeight()
{
  return Math.min(canvas.scrollHeight, canvas.height);
}

function getWidth()
{
  return Math.min(canvas.scrollWidth, canvas.width);
}

function moveGraphic(deltaX, deltaY)
{      
  var ctx = canvas.getContext("2d");
  var width = getWidth();
  var height = getHeight();
  asm["_moveGraphic"](deltaX, deltaY);
  if (deltaX >= 0)
  {
    if (deltaY >= 0)
    {    // Move to right and bottom.
      if (width-deltaX > 0 && height-deltaY > 0)
      {  // Width and height must be greater than zero to avoid exception.
        ctx.drawImage(ctx.canvas, 0, 0,                          // Topmost pixel of source.
                      width-deltaX, height-deltaY,               // Width and height of source.
                      deltaX, deltaY,                            // Topmost pixel of destination.
                      width-deltaX, height-deltaY);              // Width and height of destination.
        drawGraphic(ctx, 0, 0, width, deltaY);                   // Draw new points at top of graphic.
        drawGraphic(ctx, 0, deltaY, deltaX, height-deltaY);      // Draw new points at left of graphic.
      }
    }
    else
    {    // Move to right and up.
      if (width-deltaX > 0 && height+deltaY > 0)
      {  // Width and height must be greater than zero to avoid exception.
        ctx.drawImage(ctx.canvas, 0, -deltaY,                    // Topmost pixel of source.
                      width-deltaX, height+deltaY,               // Width and height of source.
                      deltaX, 0,                                 // Topmost pixel of destination.
                      width-deltaX, height+deltaY);              // Width and height of destination.
        drawGraphic(ctx, 0, height+deltaY, width, -deltaY);      // Draw new points at bottom of graphic.
        drawGraphic(ctx, 0, 0, deltaX, height+deltaY);           // Draw new points at left of graphic.
      }
    }
  }
  else if (deltaY >= 0)
  {      // Move to left and bottom.
    if (width+deltaX > 0 && height-deltaY > 0)
    {    // Width and height must be greater than zero to avoid exception.
      ctx.drawImage(ctx.canvas, -deltaX, 0,                      // Topmost pixel of source.
                    width+deltaX, height-deltaY,                 // Width and height of source.
                    0, deltaY,                                   // Topmost pixel of destination.
                    width+deltaX, height-deltaY);                // Width and height of destination.
      drawGraphic(ctx, 0, 0, width, deltaY);                     // Draw new points at top of graphic.
      drawGraphic(ctx, width+deltaX, deltaY, -deltaX, height-deltaY); // Draw new points at right of graphic.
    }
  }
  else
  {      // Move to left and up.
    if (width+deltaX > 0 && height+deltaY > 0)
    {    // Width and height must be greater than zero to avoid exception.
      ctx.drawImage(ctx.canvas, -deltaX, -deltaY,                // Topmost pixel of source.
                    width+deltaX, height+deltaY,                 // Width and height of source.
                    0, 0,                                        // Topmost pixel of destination.
                    width+deltaX, height+deltaY);                // Width and height of destination.
      drawGraphic(ctx, 0, height+deltaY, width, -deltaY);        // Draw new points at bottom of graphic.
      drawGraphic(ctx, width+deltaX, 0, -deltaX, height+deltaY); // Draw new points at right of graphic.
    }
  }
}

function drawGraphic(ctx, left, top, width, height)
{
  if (width != 0 && height != 0)
  {
    var rowNbr, leftPixelSrc, leftPixelDest;
    var startX = left - (getWidth() >> 1);
    var startY = top - (getHeight() >> 1);
    asm["_drawPartialGraphic"](startX, startX+width, -(startY+height), -startY);
    if (!asmjs)
    {         // Using WebAssembly: copy from WebAssembly buffer to Canvas double buffer.
      leftPixelSrc = top*8192+left*4;
      leftPixelDest = (top+32)*8192+left*4;
      for (rowNbr=0; rowNbr<height; rowNbr++)
      {
        bitsCanvas.set(pixels.subarray(leftPixelSrc, leftPixelSrc + 4*width), leftPixelDest);
        leftPixelSrc += 8192;
        leftPixelDest += 8192;
      }
    }
    ctx.putImageData(imgData, 0, -32, left, 32+top, width, height);
  }
}

function updateGraphic(nbr)
{
  var idx, value, ctr;
  if (nbr == 1)
  {
    ctr = 0;
    value = centerX.value;
    for (idx=0; idx<value.length; idx++, ctr++)
    {
      HEAP8[11000000+ctr] = value.charCodeAt(idx);
    }
    HEAP8[11000000+ctr++] = 0;
    value = centerY.value;
    for (idx=0; idx<value.length; idx++, ctr++)
    {
      HEAP8[11000000+ctr] = value.charCodeAt(idx);
    }
    HEAP8[11000000+ctr] = 0;
  }
  var width = getWidth();
  var height = getHeight();
  asm["_nbrChanged"](11000000, nbr, width, height);
  drawGraphic(canvas.getContext("2d"), 0, 0, width, height);
}

function charDecode(ch)
{
  if (ch >= 65 && ch <= 90)
  {                            // Character between A and Z.
    ch -= 65;                  // Convert to range 0 to 25.
  } 
  else if (ch >= 97 && ch <= 122)
  {                            // Character between a and z.
    ch -= 71;                  // Convert to range 26 to 51.
  }
  else if (ch >= 48 && ch <= 57)
  {                            // Character between 0 and 9.
    ch += 4;                   // Convert to range 52 to 61.
  }
  else if (ch == 43)
  {                            // Character is a plus sign.
    ch = 62;                   // Convert to code 62.
  }
  else if (ch == 47)
  {                            // Character is a slash.
    ch = 63;                   // Convert to code 63.
  }
  return ch;
}

function b64decode(str, out)
{
  var ch;
  var idxDest, idxSrc;
  var blocks,left_over;
  var len = str.length;
  // Ignore 
  if (str.charAt(len-1) == '=')
  {
    len--;
  }
  if (str.charAt(len-1) == '=')
  {
    len--;
  }
  blocks = len & -4;
  for (idxDest=0,idxSrc=0; idxSrc < blocks; idxDest += 3,idxSrc += 4)
  {
    out[idxDest] = (charDecode(str.charCodeAt(idxSrc)) << 2) + ((charDecode(str.charCodeAt(idxSrc+1)) & 0x30) >> 4);
    out[idxDest+1] = (charDecode(str.charCodeAt(idxSrc+1)) << 4) + (charDecode(str.charCodeAt(idxSrc+2)) >> 2);
    out[idxDest+2] = (charDecode(str.charCodeAt(idxSrc+2)) << 6) + charDecode(str.charCodeAt(idxSrc+3));
  }
  left_over = len & 3;
  if (left_over == 2)
  {
    out[idxDest] = (charDecode(str.charCodeAt(idxSrc)) << 2) + ((charDecode(str.charCodeAt(idxSrc+1)) & 0x30) >> 4);
    out[idxDest+1] = (charDecode(str.charCodeAt(idxSrc+1)) << 4);
  }
  else if (left_over == 3)
  {
    out[idxDest] = (charDecode(str.charCodeAt(idxSrc)) << 2) + ((charDecode(str.charCodeAt(idxSrc+1)) & 0x30) >> 4);
    out[idxDest+1] = (charDecode(str.charCodeAt(idxSrc+1)) << 4) + (charDecode(str.charCodeAt(idxSrc+2)) >> 2);
    out[idxDest+2] = charDecode(str.charCodeAt(idxSrc+2)) << 6;
  }
}

function showInfo(bottomText, centerXText, centerYText)
{
  var bottom = PtrToString(bottomText);
  var centerXValue = PtrToString(centerXText);
  if (centerXValue != centerX.value)
  {
    centerX.value = centerXValue;
  } 
  var centerYValue = PtrToString(centerYText);
  if (centerYValue != centerY.value)
  {
    centerY.value = centerYValue;
  } 
  info.innerHTML = bottom;
}

function startLowLevelCode()
{
  var length, bytes;
  var info;
  imgData = canvas.getContext("2d").createImageData(2048, 4096);  // 32 MB;
  buffer = imgData.data.buffer;          // Reserve 16 MB.
  if (asmjs)
  {                                      // Asm.js initialization.
    HEAP8 = new Uint8Array(buffer);
    HEAP8.fill(0);                       // Initialize entire array to zero.
    HEAP8.set(initializer, 8);           // Initialized data starts from offset 8.
    global = {"Math": Math, "Int8Array": Int8Array, "Int32Array": Int32Array};
    env = {"STACKTOP": 12000000, "STACKMAX": 13000000,
      abortStackOverflow: function(q) {},
      _showInfo: showInfo};
    // check for imul support, and also for correctness ( https://bugs.webkit.org/show_bug.cgi?id=126345 )
    if (!Math['imul'] || Math['imul'](0xffffffff, 5) !== -5) Math['imul'] = function imul(a, b)
    {
      var ah  = a >>> 16;
      var al = a & 0xffff;
      var bh  = b >>> 16;
      var bl = b & 0xffff;
      return (al*bl + ((ah*bl + al*bh) << 16))|0;
    };
    Math.imul = Math['imul'];
    asm = myAsmJs(global, env, buffer);  // Link asm.js module.
    ShowInformation = asm["_ShowInformation"];
  }
  else
  {                                      // WebAssembly initialization.
    length = wasm.length * 3 / 4;
    if (wasm.charCodeAt(wasm.length - 1) == 61)
    {                                    // Base64 ending equal sign found.
      length--;
    }
    if (wasm.charCodeAt(wasm.length - 2) == 61)
    {                                    // Another base64 ending equal sign found.
      length--;
    }
    bytes = new Int8Array(length);
    // Decode Base64.
    b64decode(wasm, bytes);
    info = {"env": {"_showInfo": showInfo}};
    WebAssembly["instantiate"](bytes, info).then(function(results)
    {
      asm = results["instance"]["exports"];
      ShowInformation = asm["_ShowInformation"];
      HEAP8 = new Uint8Array(asm["memory"]["buffer"]);
      bitsCanvas = new Uint8Array(buffer);
      pixels = HEAP8.subarray(asm["_getPixels"]());
      updateGraphic(1);
      return;
    });
  }
}

function zoomIn()
{
  if (zoom < 32)
  {
    zoom <<= 1;
    if (zoom == 32)
    {
      zoomin.disabled = true;
    }
    zoomout.disabled = false;
    updateGraphic(3);
  }
}

function zoomOut()
{
  if (zoom > 1)
  {
    zoom >>= 1;
    if (zoom == 1)
    {
      zoomout.disabled = true;
    }
    zoomin.disabled = false;
    updateGraphic(4);
  }
}

function keydown(e) 
{
  var evt = e || window.event;
  var target = evt.target || evt.srcElement;
  var key = evt.keyCode;
  if (evt.ctrlKey == false && evt.altKey == false && evt.metaKey == false)
  {                                  // No modifier key pressed.
    if (key >= 0x60 && key <= 0x69)
    {
      key -= 0x30;                   // Convert numpad key to standard digit key.
    }
    if (key >= 0x30 && key <= 0x39)
    {                                // Digit key has been pressed.
      if (target.value.length >= 10 || (target.value.charAt(0) != '-' && target.value.length >= 9))
      {                              // Number is too large.
        evt.preventDefault();        // Do not propagate this key.
      }
    }
    else if (key == 109 || key == 189)
    {                                // Key minus has been pressed.
      if (target.value.indexOf("-") >= 0)
      {                              // There is already a minus sign.
        evt.preventDefault();        // Do not propagate this key.
      }
	  else
	  {
		beforeMinus = target.value;
	  }
    }
    else if (key != 8 && key != 9 && key != 37 && key != 39 && key != 45 && key != 46)
    {                                // Not backspace, tab, right or left arrow, insert or delete key.
      evt.preventDefault();          // Do not propagate this key.  
    }
  }
}
    
function startUp()
{  
  canvas = get("canvas");
  zoomin = get("zoomin");
  zoomout = get("zoomout");
  centerX = get("centerX");
  centerY = get("centerY");
  zoom = 8;
  zoomDone = 0;
  isMouseDown = false;
  var domRect = canvas.getBoundingClientRect();
  canvas.width = domRect.width;
  canvas.height = domRect.height;
  startLowLevelCode();
  zoomin.onclick = function ()
  {
    zoomIn();
  };
  zoomout.onclick = function ()
  {
    zoomOut();
  };
  canvas.onkeydown = function (e)
  {
    var evt = e || window.event;
    var key = evt.keyCode;
    switch (key)
    {
      case 37: // Left arrow:
        moveGraphic(4, 0);
        ShowInformation(-1, -1);
        evt.preventDefault();          // Do not propagate this key.
        break; 
      case 38: // Up arrow:
        moveGraphic(0, 4);
        ShowInformation(-1, -1);
        evt.preventDefault();          // Do not propagate this key.
        break; 
      case 39: // Right arrow:
        moveGraphic(-4, 0);
        ShowInformation(-1, -1);
        evt.preventDefault();          // Do not propagate this key.
        break; 
      case 40: // Down arrow:
        moveGraphic(0, -4);
        ShowInformation(-1, -1);
        evt.preventDefault();          // Do not propagate this key.
        break; 
    }      
  };
  canvas.onmousemove = function(e)
  {
    var evt = e || window.event;
    var newX = evt.clientX;
    var newY = evt.clientY;      // Coordinates relative to window.
    var rect = canvas.getBoundingClientRect();
    newX -= rect.left;
    newY -= rect.top;            // Now coordinates are relative to canvas.
    if (isMouseDown)
    {
      if (newX != currentX || newY != currentY)
      {
        moveGraphic(newX - currentX, newY - currentY);
        currentX = newX;
        currentY = newY;
      }
    } 
    else
    {
      if (newX < canvas.width && newY < canvas.height)
      {
        asm["_ShowInformation"](newX, newY);
      }
    }
  };  
  canvas.onmousedown = function(e)
  {
    var evt = e || window.event;
    currentX = evt.clientX;
    currentY = evt.clientY;      // Coordinates relative to window. 
    var rect = canvas.getBoundingClientRect();
    currentX -= rect.left;
    currentY -= rect.top;        // Now coordinates are relative to canvas.
    isMouseDown = true;
  };
  document.onmouseup = function()
  {
    isMouseDown = false;
  };
  canvas.addEventListener("touchstart", function(e)
  {
    var evt = e || window.event;
    var touches = evt.targetTouches;
    touch = touches[0];
    prevX1stTouch = Math.round(touch.pageX);
    prevY1stTouch = Math.round(touch.pageY);
    if (touches.length == 2)
    {
      touch = touches[1];
      prevX2ndTouch = Math.round(touch.pageX);
      prevY2ndTouch = Math.round(touch.pageY);
      zoomDone = 0;
    }    
  }, false);  
  canvas.addEventListener("touchmove", function(e)
  {
    var touch2, oldDist, newDist, diffX, diffY;
    var evt = e || window.event;
    var touches = evt.targetTouches;
    var touch1 = touches[0];
    var newX = Math.round(touch1.pageX);
    var newY = Math.round(touch1.pageY);
    if (touches.length == 1)
    {      // Drag gesture.
      if (newX != prevX1stTouch || newY != prevY1stTouch)
      {
        moveGraphic(newX - prevX1stTouch, newY - prevY1stTouch);
        prevX1stTouch = newX;
        prevY1stTouch = newY;
        ShowInformation(-1, -1);
      }      
    }
    else if (touches.length == 2 && !zoomDone)
    {      // Pinch (zoom in and zoom out) gesture.
      touch2 = evt.touches[1];
           // Get square of distance between fingers.
      diffX = prevX2ndTouch - prevX1stTouch;
      diffY = prevY2ndTouch - prevY1stTouch;
      oldDist = diffX * diffX + diffY * diffY;
      diffX = Math.round(touch2.pageX) - newX;
      diffY = Math.round(touch2.pageY) - newY;
      newDist = diffX * diffX + diffY * diffY;
      if (newDist > oldDist)
      {    // Zoom in.
        zoomIn();
        zoomDone = 1;
      }
      if (newDist < oldDist)
      {    // Zoom out.
        zoomOut();
        zoomDone = 1;
      }
    }
  }, false);
  centerX.onkeydown = function(e)
  {
    keydown(e);
  };
  centerX.oninput = function()
  {
	if (beforeMinus != "" && centerX.value != "-" + beforeMinus)
	{     // Minus sign is not in the first character. Move it to first character.
	  centerX.value = "-" + beforeMinus;
	}
	beforeMinus = "";
    get("info").innerHTML = "";	
    updateGraphic(1);
  };
  centerY.onkeydown = function(e)
  {
    keydown(e);
  };
  centerY.oninput = function()
  {
	if (beforeMinus != "" && centerY.value != "-" + beforeMinus)
	{     // Minus sign is not in the first character. Move it to first character.
	  centerY.value = "-" + beforeMinus;
	}
	beforeMinus = "";
    get("info").innerHTML = "";
    updateGraphic(1);
  };
  window.onresize = function()
  {
    var domRect = canvas.getBoundingClientRect();
    canvas.width = domRect.width;
    canvas.height = domRect.height;
    updateGraphic(1);
  };
  if (asmjs)
  {
    updateGraphic(1);
  }
}

addEventListener("load", startUp);
})(this);
