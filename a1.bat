copy .htaccess_orig toweb\.htaccess
copy *00*js toweb
del *00*js
copy *.webmanifest toweb
cd toweb
perl ..\csp.pl
perl ..\csp2.pl
echo ***** gzip *****
for %%a in (*.HTM) do gzip -9 -c %%a > %%a.gz
for %%a in (*00*js) do gzip -9 -c %%a > %%a.gz
for %%a in (*.webmanifest) do gzip -9 -c %%a > %%a.gz
echo **** brotli ****
for %%a in (*.HTM) do brotli -f -Z %%a
for %%a in (*00*js) do brotli -f -Z %%a
for %%a in (*.webmanifest) do brotli -f -Z %%a
