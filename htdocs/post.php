<?php

    include_once 'config.php';
  
    //get posted email and maxisteps
    $email = $_POST['email'];
    $maxIsteps = $_POST['maxisteps'];
    $formula0 = $_POST['formula0'];
    $formula1 = $_POST['formula1'];
    $formula12 = $_POST['formula12'];
    $rhoin0 = $_POST['rhoin0'];
    $rhoin1 = $_POST['rhoin1'];
    $rhoin12 = $_POST['rhoin12'];
    if($_POST['top'] !="") {
	$formula0=$_POST['top'];
    	if($_POST['top'] == "H2O" ) {
	$rhoin0=0.998;
	}
    	if($_POST['top'] == "Vacuum/Air") {
	$formula0="";
	$rhoin0=0;
	}
    }
    if($_POST['maximum'] !="") {
	$formula12=$_POST['maximum'];
    	if($_POST['maximum'] == "H2O" ){
	$formula12="H2O";
	$rhoin12=0.998;
	}
    	if($_POST['maximum'] == "Lipids" ){
	$formula12="H2O";
	$rhoin12=2.;
	}
    }
    if($_POST['substrate'] !="") {
	$formula1=$_POST['substrate'];
    	if($_POST['substrate'] == "H2O" ){
	$rhoin1=0.998;
	}
    	if($_POST['substrate'] == "lipids") {
	$formula1="H2O";
	$rhoin1=2.;
	}
    }

    $slab = $_POST['slab'];
    $nslab = $_POST['nslab'];
    if(!($nslab >= 1 && $nslab <= 128)) $nslab=64; //use valid $nslab
    $energy = $_POST['energy'];
    $qmin = $_POST['qmin'];
    $qmax = $_POST['qmax'];
    if($formula1 =="") {
    	$formula1="H2O";
	$rhoin1=0.998;
    }
    if(! ($energy >= 1 && $energy <= 100)) {
    	$energy=8.;
    }
    if(! ($slab >= 10 && $slab<= 500)) {
    	$slab=45.;
    }
    if($qmin > $qmax) {
    	$qtmp=$qmin;
	$qmin=$qmax;
	$qmax=$qtmp;
    }
    //end get email and maxisteps;
    
    //get data from file or pasted content
	$_FILES['uploaded']['tmp_name'];
	$oldName=$_FILES['uploaded']['name'];
    if( $oldName != "" )
    {
        $fileSize = $_FILES['uploaded']['size'];
		//echo "Uploaded file size: $fileSize<br>";
	if($fileSize > 65535 ) {
		echo "Uploaded file too big.<br>Please Back and try again.<br>";
		exit;
	}
	$tmpName=$_FILES['uploaded']['tmp_name'];
        $fp = fopen($tmpName, 'r');
        $data = fread($fp, $fileSize);
        fclose($fp);
//	echo "uploaded<br>$data";
    }
    else
    {
	$oldName='ref.txt';
        $data = $_POST['pasted'];
//	echo "pasted<br>$data";
    }
    //end get data

    //make data insertable
    $email = addslashes($email);
//    echo "email=$email";
    $data = addslashes($data);
//    $formula0 = addslashes($formula0);
//    $formula1 = addslashes($formula1);
//    echo "formula1=$formula1";
//    $formula12 = addslashes($formula12);
//    $rhoin0 = addslashes($rhoin0);
//    $rhoin1 = addslashes($rhoin1);
//    $rhoin12 = addslashes($rhoin12);
    //end make data insertable
    //$fitbulk = $_POST['fitbulk'];
    if ( $fitbulk != '1' ) {
            $fitbulk = '0';
    }
    //insert form data into db
    $link= mysql_connect($db_server,$db_username,$db_password);
    mysql_select_db($db_database) or die( "Unable to select database");
    //$query = "INSERT INTO tasks(maxIsteps, email, data) VALUES($maxIsteps, '$email', '$data')";
    $date0 = date('U');
    $query = "INSERT INTO tasks(started,date,progress,done,maxIsteps, email, data,formula0,formula1,formula12,rhoin0,rhoin1,rhoin12,slab,energy,qmin,qmax,nslab,fitbulk) VALUES(0,'$date0','0','0','$maxIsteps','$oldName', '$data','$formula0','$formula1','$formula12','$rhoin0','$rhoin1','$rhoin12','$slab','$energy','$qmin','$qmax','$nslab','$fitbulk')";
    //die($query);
    mysql_query($query) or die( mysql_error($link) );
    mysql_close();
    $id = mysql_insert_id($link) ;
    //end insert

    //redirect to another page
    header("Location: last.php?id=$id");
?>
