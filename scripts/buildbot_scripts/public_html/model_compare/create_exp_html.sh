#!/bin/ksh

Part1="test_hat_jww-moist_cld-cnv-vdf"
Lists="iconR2B04-grid_0002_precip_APRC"

cat > ${Part1}.html <<EOF
<html>
<head>
<link rel="stylesheet" type="text/css" href="ScreenStyle.css">
<style type="text/css"></style>
</head>

<body>

<h3>exp.${Part1}</h3>
<a>To see the Postprocessing results please select the plot</a><br><br> 

<ol style="list-style-type: none;">
EOF

for List in ${Lists}
do
cat >> ${Part1}.html <<EOF
  <li>
    ${List}<br>

         <table id="UserCanvas" border="1" cellspacing="0" cellpadding="0" align="left">
           <tr>
             <td id="ZeileTL_2_icon">
                <a href="compare_${List}.html"> 
	        <img src="eps/reference_${List}/${Part1}_${List}.jpg" 
                     width="350" 
                     height="84" /></a> 
             </td>
           </tr>
         </table>
    <br><br><br><br><br>
  </li>
<br>
EOF
done

cat >> ${Part1}.html <<EOF
</ol>

<!-- ==================================================================== -->


</body>
</html>
EOF
