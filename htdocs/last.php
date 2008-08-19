<?php
    include_once 'config.php';
    
    $id = $_GET['id'];     
    
    $link= mysql_connect($db_server,$db_username,$db_password);
    mysql_select_db($db_database) or die( "Unable to select database");
    $query = "SELECT * FROM tasks WHERE id=$id";
    $result = mysql_query($query) or die(mysql_error($link));        
    mysql_close();
    
    if(mysql_numrows($result) == 0)
        $result = null;
    
    if($result == null)
    {
        echo "Invalid Task id $id";    
	return(0);
    }
    else
    {
        $email = mysql_result($result, 0, 'email');
        if( $email == '' ){
                $email="ref-$id.txt";
        }
        $maxisteps = mysql_result($result, 0, 'maxisteps');
        $data = mysql_result($result, 0, 'data');
        $progress = mysql_result($result, 0, 'progress');
        $date =  mysql_result($result, 0, 'date');
        $done =  mysql_result($result, 0, 'done');
        $started =  mysql_result($result, 0, 'started');
        $running =  mysql_result($result, 0, 'running');
        $error =  mysql_result($result, 0, 'error');
        $sigma0 =  mysql_result($result, 0, 'sigma0');
        $slab =  mysql_result($result, 0, 'slab');
        $nslab =  mysql_result($result, 0, 'nslab');
    }
?>





<html>
    <head>
        <TITLE></TITLE>
    
    
        <script type="text/javascript" language="JavaScript">
            function UpdateTimer()
            {
                window.location.reload(true);
            }
        </script>
        
    </head>
<BODY <?php if($started >=0 || $running ==1 ) echo "onload=\"setTimeout('UpdateTimer()', 20000);\"";  ?>>
        
<?php
		echo "Task $id($running) <b>";
                if ($started ==1 ) echo "<span style=\"color: rgb(51, 204, 0);\">Running</span>";
                if ($started ==0 ) echo "<span style=\"color: rgb(255, 204, 0);\">Waiting to be started</span>(if it waits forever, please ask Dongxu to start the dispatcher)<br>";
                if ($started ==-1 && $running==1 ){
		echo "<span style=\"color: rgb(255,204, 0);\">Stopping</span>";
		} 
                if ($started ==-1 && $running==0 ){
		echo "<span style=\"color: rgb(255,0,0);\">Stopped</span>";
		if($error == 1) echo "<p><span style=\"color: rgb(255,128,32);\">Couldn't start the task. Incorrect data format?</span></p>";
		}
		echo "</b><br>Progress: $progress/$maxisteps<br><BR>";
		echo "<center>";
		echo "<form method=\"get\" action=\"stop.php?id=$id\" name=\"Stop\"><INPUT type=\"hidden\" name=\"id\" value=\"$id\">";
		    if ($running == 1) {
		    echo "<button value=\"-1\" name=\"act\">Stop</button>";
		    } else {
		    if($progress>=$maxisteps) echo "<INPUT class=\"input\" TYPE=\"text\" NAME=\"isteps\" value=\"100\" size=\"3\" maxlength=\"3\"> more steps &nbsp";
		    if($started<0) echo "<button value=\"0\" name=\"act\">Start</button>";
		    }
		    echo "</form>";
		echo "<table style=\"text-align: center;\" border=\"0\"
 cellpadding=\"2\" cellspacing=\"2\"><tbody><tr><td>";
 		                    if($progress >= 1){ 
		    echo "<table cellpadding='2' cellspacing='10'><tbody><tr><td>";
		    echo "<img src='downloads/image-$id-$progress.png'><br></td>";
		    echo "<td style='vertical-align:top;text-align:right'><table border='2'><tbody><tr><center><h3>Results</h3> (all units in &Aring)</center></tr>";
		    echo "<tr><td><h4>Roughness</h4></td><td>$sigma0 &Aring</td></tr>";
		    echo "<tr><td><h4>R(q<sub>z</sub>)/R<sub>F</sub></h4> (Experimental)</td><td><a href='downloads/ref-$id.txt'>$email</a></td></tr>";
		    echo "<tr><td><h4>R(q<sub>z</sub>)/R<sub>F</sub></h4> (best-fit)</td><td><a href='downloads/rf-$id-$progress.txt'>rf-$id-$progress.txt</a></td></tr>";
                    echo "<tr><td><h4><a href='http://links.jstor.org/sici?sici=0027-8424(19870715)84%3A14%3C4709%3AROTLIO%3E2.0.CO%3B2-3'>Longitudinal Density Profile</a></h4></td><td><a href='downloads/rho-$id-$progress.txt'>rho-$id-$progress.txt</a></td></tr>";
                    echo "<tr><td>Interface thickness</td><td> $slab &Aring; </td></td>";
                    echo "<tr><td>Number of slabs</td><td>$nslab</td></td>";
                    echo "</tbody></table></td></tr></tbody></table>";
echo "<a href='view.php?id=$id'>View Progress</a><p>";
                    } else {
                            echo "Waiting for results (about 40 seconds each step) <br>";
                    }
		    echo "</td></tr></tbody></table>";
		?>
</BODY>
</html>
