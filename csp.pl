use Digest::SHA qw(sha256_base64);
use File::Copy;

my $start = "Header set Content-Security-Policy \"default-src 'none'; base-uri 'none'; form-action 'self'; child-src 'self' blob:;worker-src 'self' blob:; img-src 'self' www.google-analytics.com; connect-src 'self' www.google-analytics.com; frame-ancestors 'self'; manifest-src 'self';";
copy(".htaccess_orig", ".htaccess");
open(my $htaccess, '>>', '.htaccess');
opendir($dir, ".");
while (readdir $dir)
{
  my $dirEntry = $_;
  my $firstIndex = index($dirEntry, ".HTM");
  if ($firstIndex == -1)
  {
    next;
  }
  open(my $filehandle, '<', $dirEntry);
  my $data = do { local $/; <$filehandle> };
  my $extra = " 'self' blob:";
  if (index($data, "WebAssembly") >= 0)
  {
    if (index($data, "instantiate") >= 0)
    {
      $extra = " 'self' blob: 'unsafe-eval'";
    }
  }
  print $htaccess  "<Files ${dirEntry}>\n${start}";
  if (index($dirEntry, "ECM.HTM") != -1 || index($dirEntry, "ECMC.HTM") != -1)
  {
    print $htaccess " media-src 'self';"
  }
  getHashes($data, "style", $hash, "");
  print $htaccess  $hash;
  if ($hash ne "")
  {
    print $htaccess ";";
  }
  getHashes($data, "script", $hash, $extra);
  print $htaccess  $hash;
  if (index($dirEntry, "ECM.HTM") != -1)
  {
    getFileHash("blockly.js", $hash);
    print $htaccess  $hash;
    getFileHash("en.js", $hash);
    print $htaccess  $hash;
  }
  if (index($dirEntry, "ECMC.HTM") != -1)
  {
    getFileHash("blockly.js", $hash);
    print $htaccess  $hash;
    getFileHash("es.js", $hash);
    print $htaccess  $hash;
  }
  if ($hash ne "")
  {
    print $htaccess ";";
  }
  print $htaccess  " report-uri https://alpertron23.report-uri.com/r/d/csp/enforce\"\n";
  print $htaccess  "</Files>\n\n";
  close($filehandle);
}
closedir $dir;

sub getHashes
{
  my $data = $_[0];
  my $tagname = $_[1];
  my $extra = $_[3];
  $_[2] = "";
  for (;;)
  {
    my $firstIndex = index($data, "<${tagname}>");
    if ($firstIndex == -1)
    {
      last;
    }
    my $lastIndex = index($data, "</${tagname}>", $firstIndex);
    $firstIndex = $firstIndex + length($tagname) + 2;
    my $substr = substr($data, $firstIndex, $lastIndex - $firstIndex);
    my $hash = sha256_base64($substr);
    while (length($hash) % 4)
    {
      $hash .= '=';
    }
    if ($_[2] eq "")
    {
      $_[2] = " ${tagname}-src 'unsafe-inline'${extra}";
    }
    $_[2] .= " 'sha256-${hash}'";
    $data = substr($data, $lastIndex);
  }
  close $fh;
}

sub getFileHash
{
  my $filename = $_[0];
  open my $filehandle, '<', $filename or die $!; 
  my $string = do { local $/; <$filehandle> };
  close $filehandle;
  my $hash = sha256_base64($string);
  while (length($hash) % 4)
  {
    $hash .= '=';
  }
  $_[1] = " 'sha256-${hash}'";
}