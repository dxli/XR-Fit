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
        echo 'invalid id';    
    }
    else
    {
        $email = mysql_result($result, 0, 'email');
        $maxisteps = mysql_result($result, 0, 'maxisteps');
        $data = mysql_result($result, 0, 'data');
        $progress = mysql_result($result, 0, 'progress');
        $date =  mysql_result($result, 0, 'date');
        $done =  mysql_result($result, 0, 'done');
        $started =  mysql_result($result, 0, 'started');
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
<BODY <?php if($done != 1) echo "onload=\"setTimeout('UpdateTimer()', 20000);\"";  ?>>
<center><?php 
    echo "<a href=\"last.php?id=$id\">Back</a>"; ?>
</center>
       <hr style="width: 100%; height: 2px;"> 
<?php
        if($result != null)
        {
            if($started == 1)
            {
                echo "Task $id Started<br>Progress: $progress/$maxisteps<br><BR>";
                for($i = 1; $i<=$progress; $i++)
                {
                    echo "<img src='downloads/image-$id-$i.png'><br>";
                }
            }
            else
            {
                echo "Task $id is waiting in the queue<br><br>";
            }
        }
?>
<center><?php 
    echo "<a href=\"last.php?id=$id\">Back</a>"; ?>
</center>
       <hr style="width: 100%; height: 2px;"> 
</BODY>
</html>
