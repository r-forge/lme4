
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

<p> <strong>Doxygen</strong> documentation of the underlying C functions is <a href="./doxygen/"><strong>here</strong></a>. </p>

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

<p> <strong>References to articles and other research using nlme or lme4</strong> can be found <a href="bib/lme4bib.html"><strong>here</strong></a>. The LaTeX bibliography file can be accessed from <a href="bib/lme4bib.bib">here</a>. If you would like to add your work to this database, please email vasishth.shravan at gmail dot com</p>

<p>
<strong>Slides</strong> from short courses on lme4 are
<a href="./slides/"><strong>here</strong></a>.
</p>

<p> <strong>Chapter drafts</strong> of the book <em>lme4:
Mixed-effects Modeling with R</em> are available <a
href="./book/"><strong>here</strong></a>.
</p>
<p><strong>Binary and source packages</strong> of recent versions
of the development version <code>lme4a</code> are available
<href="./repos"><strong>here</strong></a> (<em>under development</em>).
</p>

</body>
</html>
