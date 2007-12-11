<?php
    include_once 'config.php';
    
    $id = $_GET['id'];     
    $act = $_GET['act'];
    $isteps = $_GET['isteps'];
    if(!($isteps >=1 && $isteps <=500)) $isteps=20;
    
    $link= mysql_connect($db_server,$db_username,$db_password);
    mysql_select_db($db_database) or die( "Unable to select database");
    $query = "UPDATE tasks set started='$act' WHERE id=$id";
    $result = mysql_query($query) or die(mysql_error($link));        
    if ($act == 0) {
    $query = "UPDATE tasks set maxisteps=maxisteps+$isteps WHERE id=$id and progress>=maxisteps";
    $result = mysql_query($query) or die(mysql_error($link));        
    $query = "UPDATE tasks set running=0,error=0 WHERE id=$id";
    $result = mysql_query($query) or die(mysql_error($link));        
    mysql_close();
    }
    
    if($result == null)
    {
        echo 'invalid id';    
    }
    header("Location: last.php?id=$id");
?>
