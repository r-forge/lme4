
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

<code>lme4</code> development has moved to <a href="https://github.com/lme4/lme4/">Github</a>; please go there for up-to-date development versions and information.  

<h3>Resources that remain on R-forge</h3>
<ul>
<li><a href="bib/lme4bib.html">References to articles and other research</a> using <code>nlme</code> or <code>lme4</code>, or the corresponding
<a href="bib/lme4bib.bib">BibTeX file</a>. (If you would like to add your work to this database, please email <code>vasishth.shravan at gmail dot com</code>.)</li>
<li><a href="./slides/">Slides</a> from short courses on <code>lme4</code></li>
<li><a href="./lMMwR/">Chapter drafts</a> of the book <em>lme4:
Mixed-effects Modeling with R</em></li>
<p><a href="./repos">Repository</a> containing (relatively) up-to-date binary and source packages of recent versions
of the development version (for most up-to-date code, which must be installed from source, see Github).
</ul>
</body>
</html>
