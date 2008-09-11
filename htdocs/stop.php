<?php
    include_once 'config.php';
    
    $id = $_GET['id'];     
    $act = $_GET['act'];
    $isteps = $_GET['isteps'];
    
    if(!($isteps >=1 && $isteps <=500)) $isteps= 100;
    
    $link= mysql_connect($db_server,$db_username,$db_password);
    mysql_select_db($db_database) or die( "Unable to select database");
  
    if ($act == 0) {
    $query = "UPDATE tasks set maxisteps=maxisteps+$isteps WHERE id=$id and progress>=maxisteps";
    $result = mysql_query($query) or die(mysql_error($link));
    
    $qmin = $_GET['qmin'];
    $qmax = $_GET['qmax'];
    $slab = $_GET['slab']; 
    $nslab = $_GET['nslab']; 
    $forward = 0;
    $forward = $_GET['forward'];
    if( $forward == '1' ) {      
    $nq = $_GET['nq'];  
    $query = "UPDATE tasks set started=0,running=0,error=0,qmin=$qmin,qmax=$qmax,nq=$nq,forward=1,slab=$slab,nslab=$nslab WHERE id=$id";
    }else
    {        
    $fitbulk = $_GET['fitbulk']; 
    if ( $fitbulk != '1' ) {
            $fitbulk='0';
    }
    $query = "UPDATE tasks set started=0,running=0,error=0,qmin=$qmin,qmax=$qmax,slab=$slab,nslab=$nslab,fitbulk=$fitbulk WHERE id=$id";
    }
    $result = mysql_query($query) or die(mysql_error($link));
    }else {
      $query = "UPDATE tasks set started=-1,running=-1 WHERE id=$id";
    $result = mysql_query($query) or die(mysql_error($link));        
    }
            
    mysql_close();
    if($result == null)
    {
        echo 'invalid id';    
    }
     header("Location: last.php?id=$id");
?>
