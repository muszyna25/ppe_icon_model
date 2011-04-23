<!--- PHP web page to show set of plots  -------------------------------------------->
<!--- Martin & Mashrab 201005, 201010, 201012    ------------------------------------>

<html>
<head>
<title>ICON plot catalog</title>
<link rel="stylesheet" href="style.css" type="text/css">
</head>
<body bgcolor="#002266">

<!--- Java Scripts ------------------------------------------------------------------>

<script type="text/javascript" src="plot-catalog.js">
</script>

<!--- PHP --------------------------------------------------------------------------->

<?php
  include "common.php";

// XML parser ----------------------------------------------------------------------->

  include("parser_php4.phps");
  show_header();
  echo "\n<table width=\"100%\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\" bgcolor=\"#FFFFFF\">".
       "<tr>\n";

// Read XML-catalog ----------------------------------------------------------------->

  $cat_filename = htmlspecialchars($_GET["cat"]);
  if ($cat_filename)
  {
    $xml_cat = file_get_contents($cat_filename);
    $catalog = new XMLParser($xml_cat);
    $catalog->Parse();
    $cat_author    = $catalog->document->author[0]->tagData;
    $cat_descr     = $catalog->document->description[0]->tagData;
    $cat_docu      = $catalog->document->docu[0]->tagData;
    $cat_docutext  = $catalog->document->docutext[0]->tagData;
    $cat_related   = $catalog->document->related[0]->tagData;
    $plot_dir      = $catalog->document->plotdir[0]->tagData;
    $plot_filename0= $catalog->document->filename[0]->tagData;
    $cat_ps        = $catalog->document->ps[0]->tagData;
    $cat_pdf       = $catalog->document->pdf[0]->tagData;
    echo "<td class=sidemenu id=\"td_sidemenu\">\n".
         "<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\" width=\"136\">\n".
         "<form  name=\"sidemenu\" method=\"get\" action=\"" . $PHP_SELF . "\">\n";
    echo "<tr><td colspan=2 class=imagebox> <div class=plotsum>\n".
         "<img id=\"menuhide\" class=menuhide src=\"leftarrow.png\" onclick=\"menuHide();\" title=\"Hide Menu\">".
         "</div>".
         "</td></tr>\n";

    $plot_filename = array($plot_filename0, $plot_filename0, $plot_filename0, $plot_filename0, $plot_filename0, $plot_filename0, $plot_filename0, $plot_filename0, $plot_filename0);
    $plot_filename = array_fill(0, 32, $plot_filename0) ;
    $all = array(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31) ;
    $index[0] = 0 ;
    $i = 0;          //number of plots

// Menu of dimensions left ---------------------------------------------------------->

    echo "<tr><td colspan=2 class=menuspace>".
       "</td></tr>\n";

    foreach($catalog->document->parameter as $c_entry)
    {
      $c_par_name = $c_entry->shortname[0]->tagData;
      eval("\$c_par_val = \$_GET[\"" . $c_par_name . "\"];");
      if (! $c_par_val)
      {
        $c_par_val = $c_entry->value[0]->short[0]->tagData;
      }

      // Create plot file names by replacing generic parameter names. --------------->

      if ( count($c_par_val) > 1)
      { 
        foreach ($c_par_val as $t) 
        {
          $plot_filename[$i] = str_replace("#" . $c_par_name . "#", $t, $plot_filename[$i]);
          $index[$i] = $i;
          $aaa[$i]=$t;
          $i++;
        }
      } else {
        foreach ($all as $ii) 
        {
          if ( is_array($c_par_val) )
          {
            $plot_filename[$ii] = str_replace("#" . $c_par_name . "#", $c_par_val[0], $plot_filename[$ii]) ; 
          } else {
            $plot_filename[$ii] = str_replace("#" . $c_par_name . "#", $c_par_val,    $plot_filename[$ii]) ; 
          }
        }
      }

      // Icons for next/previous plot and checkbox for enabling multi-plot feature -->
 
      echo "<tr><td class=menuheaderleft>" .
           $c_entry->fullname[0]->tagData .
           "</b>";
      echo "</td>\n";

      echo "<td class=menuheaderright align=right>" ;
      if ( $_GET["multiplot"] != $c_par_name )
      {
        echo "<img src=step_prev.gif width=16 onclick=\"onPrevious('". $c_par_name . "');\">".
	     "<img src=step_next.gif width=16 onclick=\"onNext('"    . $c_par_name . "');\">\n";
      } 
      echo "<input name=multiplot value=$c_par_name type=checkbox onChange=\"autoSubmit()\" title=\"Multiple plots for this parameter.\" ";
      if ($_GET["multiplot"] == $c_par_name)
      {
        echo "checked";
      } 
      echo ">";
      echo "</td></tr>\n";

      echo "<tr><td class=menubox colspan=2>\n";

      // Selection of parameter choice - including multi-plot option ---------------->
 
      if ($_GET["multiplot"] == $c_par_name)
      {
        echo "<select name=" . $c_par_name . "[] onChange=\"autoSubmit()\"   multiple size=5 >\n";
      //echo "<select name=" . $c_par_name . "[] onKeyup=\"autoSubmit()\"    multiple size=5 >\n";
      //echo "<select name=" . $c_par_name . "[] onMouseout=\"autoSubmit()\" multiple size=5 >\n";
      //echo "<select name=" . $c_par_name . "[] onClick=\"autoSubmit()\"    multiple size=5 >\n";
      }
      else
      {
        echo "<select id=" . $c_par_name . " name=" . $c_par_name . "[] onChange=\"autoSubmit()\"  >\n";
      }

      foreach($c_entry->value as $c_value)
      {
        $c_val = $c_value->short[0]->tagData;
        echo "<option value=\"" . $c_val . "\"";

        if ( $_GET["multiplot"] == $c_par_name )
	{
          foreach ($c_par_val as $t)
          {
            if ($t == $c_val)
	    {
	      echo " selected";
	    }
	  }
        }
        else
	{
          if ($c_par_val[0] == $c_val)
	  {
	    echo " selected";
	  }
	}

        echo ">" . $c_value->full[0]->tagData . "</option>\n";
      }
      echo "</select>\n".
           "</td></tr>\n";
      echo "<tr><td colspan=2 class=menuspace>".
           "</td></tr>\n";
    }

    echo "<tr><td colspan=2 class=menuspace>".
         "<br>".
         "</td></tr>\n";


// Download PS or PDF --------------------------------------------------------------->

    if ($catalog->document->ps or $catalog->document->pdf)
    {  
      echo "<tr><td colspan=2 class=menuheader>".
           "Download".
           "</td></tr>\n";
      echo "<tr><td colspan=2 class=menubox>\n";
      $find=array('.gif', '.png', '.jpg');

      if ($catalog->document->pdf and $cat_pdf == "1")
      {
        $plot_filename_pdf = str_replace($find , ".pdf" , $plot_filename[0]);
        echo "<a href=\"" . $plot_dir . "/" . $plot_filename_pdf  . "\">\n"  .
             "<img title=\"PDF\" alt=\"PDF\" border=\"0\" src=\"pdf.gif\"> PDF \n</a><br>";
      }

      if ($catalog->document->ps and $cat_ps == "1")
      {
        $plot_filename_ps = str_replace($find , ".ps" , $plot_filename[0]);
        echo "<a href=\"" . $plot_dir . "/" . $plot_filename_ps  . "\">\n" .
             "<img title=\"Postscript\" alt=\"Postscript\" border=\"0\" src=\"ps.gif\"> Postscript \n</a>";
      }
      echo "<tr><td colspan=2 class=menuspace>".
           "</td></tr>\n";
    }
     

// Related pages -------------------------------------------------------------------->

    if ($catalog->document->related)
    {  
      echo "<tr><td colspan=2 class=menuheader>".
           "Related Products".
           "</td></tr>\n";
      echo "<tr><td colspan=2 class=menubox>\n";
      foreach($catalog->document->related as $c_related)
      {
        echo "<li><a href=\"index.php?cat=" . $c_related->cat[0]->tagData . "\">".
             $c_related->text[0]->tagData . "<a><br>\n";
      }
      echo "<tr><td colspan=2 class=menuspace>".
           "</td></tr>\n";
    }

// Documentation -------------------------------------------------------------------->

    echo "<tr><td colspan=2 class=menuheader>".
         "Documentation".
         "</td></tr>\n";
    echo "<tr><td colspan=2 class=menubox>\n";
    if ($catalog->document->docu)
    {  
      foreach($catalog->document->docu as $c_docu)
      {
        echo "<li><a href=\"" . $c_docu->link[0]->tagData . "\">".
           $c_docu->text[0]->tagData . "<a>\n";
      }
      echo " <br>";
    }

// Text ----------------------------------------------------------------------------->

    echo "$cat_docutext \n";

    echo "</td></tr>\n".
         "</table>\n".
         "<input type=\"hidden\" name=\"cat\" value=\"" . $cat_filename . "\">".
         "</form>".
         "</td>\n";

// Image ---------------------------------------------------------------------------->

    echo "<td class=imagebox>".
         "<div class=plotsum>\n".
         "<img id=\"menushow\" class=menushow src=\"rightarrow.png\" onclick=\"menuShow();\" title=\"Show Menu\">".
         $cat_descr . " (" . $cat_author . ") </div>\n";
    foreach ($index as $ii) 
    {
      echo "<img alt=\"Missing file: " . $plot_dir . "/" . $plot_filename[$ii] . "\" .
                 style=\"border-left:1px solid grey; border-right:1px solid grey; border-bottom:1px solid grey\" . 
                 src=\""               . $plot_dir . "/" . $plot_filename[$ii] . "\">"; 
    } 

    echo "</td>\n";

  }
  else

// General page --------------------------------------------------------------------->

  {
    include("cat-links.html"); 
  }
  echo "</tr>".
       "</table>\n";
  show_footer();
?>

</body>
</html>
