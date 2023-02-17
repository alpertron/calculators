perl csp.pl
perl csp2.pl
del *.gz
del *.br
echo ***** gzip *****
for %%a in (*.HTM) do gzip -9 -c %%a > %%a.gz
for %%a in (*00*js) do gzip -9 -c %%a > %%a.gz
for %%a in (*.webmanifest) do gzip -9 -c %%a > %%a.gz
echo **** brotli ****
for %%a in (*.HTM) do brotli -f -Z %%a
for %%a in (*00*js) do brotli -f -Z %%a
for %%a in (*.webmanifest) do brotli -f -Z %%a
del /q toweb\*.*
copy *.HTM toweb
copy *.webmanifest toweb
copy *%1* toweb
copy .htaccess toweb
copy *.br toweb
copy *.gz toweb
