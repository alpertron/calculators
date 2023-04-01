perl csp.pl
perl csp2.pl
del *.gz
del *.br
echo ***** gzip *****
cd toweb
for %%a in (*.HTM) do gzip -9 -c %%a > %%a.gz
cd ..
for %%a in (*00*js) do gzip -9 -c %%a > %%a.gz
for %%a in (*.webmanifest) do gzip -9 -c %%a > %%a.gz
echo **** brotli ****
cd toweb
for %%a in (*.HTM) do brotli -f -Z %%a
cd ..
for %%a in (*00*js) do brotli -f -Z %%a
for %%a in (*.webmanifest) do brotli -f -Z %%a
copy *.webmanifest toweb
copy *%1* toweb
copy .htaccess toweb
copy *.br toweb
copy *.gz toweb
