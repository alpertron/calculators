#ifndef _SHOWTIME_H
#define _SHOWTIME_H
#ifdef __EMSCRIPTEN__
extern double originalTenthSecond;
void GetDHMS(char **pptrText, int seconds);
void GetDHMSt(char **pptrText, int tenths);
void showElapsedTime(char **pptrOutput, int lang);
double tenths(void);

extern double originalTenthSecond;
extern int oldTimeElapsed;
#endif
#endif