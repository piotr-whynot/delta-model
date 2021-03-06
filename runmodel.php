<?php

if(isset($_POST) and $_SERVER['REQUEST_METHOD'] == "POST"){

    $inputFile = $_POST['runinputFile'];
    $paramFile = $_POST['runparamFile'];
    $initcondFile = $_POST['runinitcondFile'];
    if($paramFile=="default"){
        $paramFile="model/config/modpar.dat";
    }
    if($initcondFile=="default"){
        $initcondFile="model/config/init.dat";
    }
    $spinup = $_POST['spinup'];
    $outputs= $_POST['outputs'];

    $runID = $_POST['runID'];
    $sessionID = $_POST['sessionID'];
    
    $outfiles="";
    foreach($outputs as $output){
        $outfiles=$outfiles." ../results/".$sessionID."-".$runID."_".$output;
    }

    #$outfile="results/".$sessionID."-".$runID."_alloutflows.csv";

    chdir("./model");
    $command = escapeshellcmd("python ./hydro_model.py ../".$paramFile." ../".$initcondFile." ../".$inputFile." ".$spinup.$outfiles)." 2>&1";
#    $command= escapeshellcmd("python ./hydro_model.py ../model/config/modpar.dat ../model/config/init.dat ../uploads/31279404-input_1967-2014.csv 0 ../results/31279404-yy_alloutflows.csv");
    #$command= escapeshellcmd("./run.sh");

    $outcome = exec($command, $response);

    $return_arr[] = array("command"=> $command, "outcome" => $outcome,
                    "response" => $response,
    );
    echo json_encode($return_arr);
}
?>
