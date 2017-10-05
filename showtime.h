#ifndef _SHOWTIME_H
#define _SHOWTIME_H
#ifdef __EMSCRIPTEN__
double tenths(void);
extern double originalTenthSecond;
extern int oldTimeElapsed;
void GetDHMS(char **pptrText, int seconds);
void GetDHMSt(char **pptrText, int tenths);
#endif
void showElapsedTime(char **pptrOutput);
#endif