<?php
session_start();
$path = "uploads/";

if(isset($_POST) and $_SERVER['REQUEST_METHOD'] == "POST"){ 
    $sessionID=$_POST['sessionID'];
    $key=array_keys($_FILES)[0];
    if ($key=="inputFileSelector"){
        $valid_file_formats = array("csv", "CSV", "Csv");
    }else{
        $valid_file_formats = array("dat", "DAT", "Dat");
    }

    $name = $_FILES[$key]['name'];
    $size = $_FILES[$key]['size'];
    if(strlen($name)) {
        list($txt, $ext) = explode(".", $name);
        if(in_array($ext,$valid_file_formats)) {
            if($size<(1024*1024)) {
                $file_name = $sessionID."-".$txt.".".$ext;
                $tmp = $_FILES[$key]['tmp_name'];
                if(move_uploaded_file($tmp, $path.$file_name)){
                    //echo "Upload successful";
                    $response=$file_name;
                    $outcome=true;
                    $return_arr[] = array("outcome" => $outcome,
                    "response" => $response, "sessionID"=>$sessionID,
                    );
                    echo json_encode($return_arr);
                }else{
                    $response="Upload failed, not sure why. Try again.";
                    $outcome=false;
                    $return_arr[] = array("outcome" => $outcome,
                    "response" => $response, "sessionID"=>$sessionID,
                    );
                    echo json_encode($return_arr);
                }
            }else{
                $response="File size maximum 1 MB";
                $outcome=false;
                $return_arr[] = array("outcome" => $outcome,
                    "response" => $response, "sessionID"=>$sessionID,
                );
                echo json_encode($return_arr);
            }
        }else{
            $response="Invalid file format";
            $outcome=false;
            $return_arr[] = array("outcome" => $outcome,
                    "response" => $response, "sessionID"=>$sessionID,
            );
            echo json_encode($return_arr);
        }
    }else{
        $response="Please select file with required extension";
        $outcome=false;
        $return_arr[] = array("outcome" => $outcome,
                    "response" => $response, "sessionID"=>$sessionID,
        );
        echo json_encode($return_arr);
    }
}
?>
