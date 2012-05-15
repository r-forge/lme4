
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<! --- R-Forge Logo --- >
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<h2>Resources</h2>
<ul>
<li>A <a href="misc/lme4_conversion.html">guide</a> to changes in the upcoming release of <code>lme4</code> (<a href="misc/lme4_conversion.md">Markdown source</a>, <a href="misc/lme4_compat_report.html">downstream package compatibility report</a>)</li>
<li>The <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/">project summary page</a>, including links to the source code browser, bug and feature request trackers, etc.</li>
<li><a href="./doxygen/"><em>Doxygen</em> documentation</a> of the underlying C functions</li>
<li><a href="bib/lme4bib.html">References to articles and other research</a> using <code>nlme</code> or <code>lme4</code>, or the corresponding
<a href="bib/lme4bib.bib">BibTeX file</a>. (If you would like to add your work to this database, please email <code>vasishth.shravan at gmail dot com</code>.)</li>
<li><a href="./slides/">Slides</a> from short courses on <code>lme4</code></li>
<li><a href="./lMMwR/">Chapter drafts</a> of the book <em>lme4:
Mixed-effects Modeling with R</em></li>

<h2>Installation</h2>
<p>Binary and source packages of recent versions
of the development version  (<code>lme4</code>, previously <code>lme4Eigen</code>) are available
<a href="./repos"><strong>here</strong></a>,
as well as MacOS binaries for <code>RcppEigen</code>.
The current status for the various versions is as follows:
<ul>
<li>The stable release of <code>lme4</code> on CRAN (version 0.999375-42) is the last version of this branch.</li>
<li>R-forge has <code>lme4.0</code> (version 0.9999-1), which is the continuation of the old branch.  It should be perfectly compatible with the current CRAN version, <emph>except</emph> for its package name.</li>
<li>R-forge also has a new version of <code>lme4</code> (version 0.999902344-0), which is as far as possible <emph>externally</emph> compatible with the old branch, but its internal structures have changed, so packages dependent on the old structures will need to be updated to work with the new release.  The new version also produces objects of class <code>merMod</code> rather than class <code>mer</code>.</li>
<li>An older development branch, <code>lme4a</code>, is also available on R-forge, but we believe that it is now dominated by the new development branch (i.e. the r-forge version of <code>lme4</code>), and we strongly request that you work with the new version. If you encounter a situation where <code>lme4</code> can't do something that <code>lme4a</code> can, please inform the maintainers.
</li>
</ul>
<p>
<ul>
<li>Although installing packages into older versions of R is sometimes possible by downloading and manually installing the package files, installation will generally work <emph>much</emph> more smoothly if you are running the latest version of R.  If you want help installing <code>lme4</code> to an older version of R you will have to provide a very good reason ...
</li>
<li>You may have to install dependencies 
(e.g. <code>RcppEigen</code> and <code>minqa</code>
for <code>lme4</code>) from other repositories first.</li>
</ul>
<ol>
<li>If possible, you should install the most recent version of the development version built on r-forge, e.g.
<pre>
install.packages("lme4",repos="http://r-forge.r-project.org")
</pre>.
</li>
<li>If the appropriate version for your OS is not available 
(e.g. binaries for MacOS), the next best option is to install
the most recent version available at the <code>lme4</code>
repository, e.g.
<pre>
install.packages("lme4",repos="http://lme4.r-forge.r-project.org/repos")
</pre>
These versions will generally be slightly less up-to-date than the ones built by R-forge, because they must be updated by hand.
</li>
<li>
Alternatively (e.g. if you have network problems),
you should be able to install these versions
by descending through the directory hierarchy 
from <a href="http://lme4.r-forge.r-project.org/repos">here</a>,
downloading a copy and using
<pre>
install.packages("[name_of_file_here]",repos=NULL)
</pre>
or a menu item such as "install from local zip file"
</li>
<li>The most reliable way to install the package 
is to check out 
<a href="http://r-forge.r-project.org/scm/?group_id=60">the most recent copy</a> and build it from source on your machine; however, this may require extra tools (Subversion client, compilers, etc.).  Alternatively, if you have the compiler tools but don't want to mess with Subversion, you will often get farther with <code>install.packages(...,type="source")</code> than trying to find a binary package compiled for your OS.
</li>
<li>
If all else fails, please contact the maintainers.
</li>
</p>
</body>
</html>
