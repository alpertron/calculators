copy .htaccess_orig toweb\.htaccess
copy *00*js toweb
del *00*js
copy *.webmanifest toweb
cd toweb
perl ..\csp.pl
perl ..\csp2.pl
