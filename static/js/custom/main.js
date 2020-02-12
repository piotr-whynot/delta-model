var inputQC=false;
var inputColNames={"PET":["Date", "Inflow-Mohembo", "Rainfall-Maun", "Rainfall-Shakawe", "PET-Maun"], "Temp":["Date", "Inflow", "Rainfall-Maun", "Rainfall-Shakawe", "MinTemp-Maun", "MaxTemp-Maun"]};
var paramFileDefault="default";
var initcondFileDefault="default";
var sessionID=Math.floor(Math.random()*100000000);
var resultsDir="results/";
var uploadsDir="uploads/";
var modelRun=false;

var outputsList=JSON.parse('{"alloutflows": ["Outflows from all model units", "checked", "csv", "Mm3/month"], "allinundation":["Inundated area of all model units", "checked", "csv", "km2"], "totalinundation":["Total inundated area", "checked", "csv", "km2"], "totalecoregions":["Eco-hydrological regions for the entire system", "", "csv", "km2"], "animatedinundation": ["Animation of inundated area", "", "gif", "km2"]}');


function initialize () {

    $("#intro_container").show();
    $('#showhide-intro').html("[close]");
    $("#input_container").show();
    $('#showhide-input').html("[close]");
    $("#parameters_container").hide();
    $("#initcond_container").hide();
    $("#run_container").hide();
    $("#results_container").hide();
//"<span class='tooltip'>?<span class='tooltiptext'>placeholder</span></span>"
    txt="<form id='input_form' method='post' enctype='multipart/form-data' action='upload.php' autocomplete='off'>";
    txt+="<input type=text id=inputFile name=inputFile hidden>";
    txt+="<input type=text id=firstDate name=firstDate hidden>";
    txt+="<input type=text id=lastDate name=lastDate hidden>";
    txt+="<input type=text id=sessionID name=sessionID hidden value="+sessionID+">";
    txt+="<input type=text id=nts name=nts hidden>";
    txt+="CSV file with input data:";
    txt+="<input type='file' name='inputFileSelector' id='inputFileSelector' class=activeInput>";
    txt+="</form>";
    txt+="<div id=input_outcome class=outcome></div>";
    $("#input_container").html(txt);

    txt="<form id='param_form' method='post' enctype='multipart/form-data' action='upload.php' autocomplete='off'>";
    txt+="Use default: <input type=text id=paramFile name=paramFile value="+paramFileDefault+" disabled class=activeInput><br>";
    txt+="or select customized parameters file: <input type=file id=paramFileSelector name=paramFileSelector >";
    txt+="</form>";
    txt+="<div id=param_outcome class=outcome><div>";
    $("#parameters_container").html(txt);
    
    txt="<form id='initcond_form' method='post' enctype='multipart/form-data' action='upload.php' autocomplete='off'>";
    txt+="Use default: <input type=text id=initcondFile name=initcondFile value="+initcondFileDefault+" disabled class=activeInput><br>";
    txt+="or select customized initial conditions file: <input type=file id=initcondFileSelector name=initcondFileSelector >";
    txt+="</form>";
    txt+="<div id=initcond_outcome class=outcome><div>";
    $("#initcond_container").html(txt);
   
    txt="<p>Nothing to see here yet. Load input data first.";
    txt+="<div id=run_outcome class=outcome><div>";
    $("#run_container").html(txt);

    txt="<p>Nothing to see here yet. Run model first.";
    $("#results_container").html(txt); 
}


$(document).on('click','.showhideDiv', function(event){
    var targetID=$(event.target).attr('id');
    var targetDiv=targetID.split("-")[1];
    targetDiv=targetDiv+"_container";
    var animate="slow";
    if($('#'+targetDiv).css('display')=='none'){
        // show element
        $('#'+targetDiv).show(animate);
        $(event.target).toggleClass('visible');
        $('#'+targetID).html("[close]");
    }else{
        $('#'+targetDiv).hide(animate);
        $(event.target).toggleClass('visible');
        $('#'+targetID).html("[open]");
    }
});




$(document).on('change','#inputFileSelector', function(event){
    $("#input_outcome").html('<img src="./static/img/ajax-loader.gif" alt="Uploading...." class=centered>');
    form=$("#input_form")[0];
    formdata=new FormData(form);
    $.ajax({
        url: "upload.php",
        method: "POST",
        data: formdata,
        processData: false,
        contentType: false,
        success: function (data) {
        // success callback
            $("#initcond_outcome").html(""); 
            data=JSON.parse(data);
            outcome=data[0]['outcome'];
            if (outcome!=true){
                $("#input_outcome").html("<span class=alert>"+data[0]['response']+"</span>");
            }else{
                uploadedFileName=data[0]['response'];
                uploadedFile=uploadsDir+uploadedFileName;
                readFile(uploadedFile, function(rawdata){
                    outcome=false;
                    arraydata=parseCSV(rawdata);
                    //QC of csv data;
                    ncol=arraydata[0].length;
                    if (ncol==5){
                        colnames=inputColNames['PET'];
                        outcome=true;
                    }else if(ncol==6){
                        colnames=inputColNames['Temp'];
                        outcome=true;
                    }else{
                        outcome=false;
                        response="expected 5 or 6 columns in csv file, got "+ncol;
                    }
                    if (outcome){
                        for(var c=0, len=ncol; c < len; c++){
                            if(colnames[c]!=arraydata[0][c]){
                                outcome=false;
                                col=c+1;
                                response="in column "+col+" expected "+colnames[c]+", got "+arraydata[0][c];
                            }
                        }
                    }

                    if (outcome!=true){
                       $("#input_outcome").html("<span class=alert>"+response+"</span>"); 
                    }else{
                        seriesData=[];
                        allData=[];
                        seriesNames=[];
                        ncol=arraydata[0].length;
                        nrow=arraydata.length;
                        $("#nts").val(nrow);
                        for(var c=1, len=ncol; c < len; c++){
                          //get var names
                            seriesNames.push(arraydata[0][c]);
                            allData.push(new Array());
                        }
                        for(var r=1, lenr=nrow; r < lenr; r++){
                            dte=new Date(arraydata[r][0]);
                            if (r==1){
                                var firstDate=dte.getFullYear()+"-"+(dte.getMonth()+1)+"-"+dte.getDate()
                                $("#firstDate").val(firstDate);
                            }
                            if (r==nrow-1){
                                lastDate=dte.getFullYear()+"-"+(dte.getMonth()+1)+"-"+dte.getDate()
                                $("#lastDate").val(lastDate);
                            }
                            dte=dte.valueOf();
                            for(var c=1, lenc=ncol; c < lenc; c++){
                                allData[c-1].push([dte, parseFloat(arraydata[r][c])])
                            }
                        }
                        for(var c=1, len=ncol; c < len; c++){
                            seriesData.push({id: c, name: seriesNames[c-1], data: allData[c-1], showInLegend:true});
                        }
                        plotTimeSeries("input_outcome", seriesData, "Input data", "Inflow: [Mm3/month] <br>others: [mm/month]");

                        $("#inputFile").val(uploadedFile);
                        populateRunData();
                    }
                });
            }
         }
    });
});


$(document).on('change','#paramFileSelector', function(event){
    $("#param_outcome").html('<img src="./static/img/ajax-loader.gif" alt="Uploading...." class=centered>');
    form=$("#param_form")[0];
    formdata=new FormData(form);
    $.ajax({
        url: "upload.php",
        method: "POST",
        data: formdata,
        processData: false,
        contentType: false,
        success: function (data) {
        // success callback
            $("#param_outcome").html(""); 
            data=JSON.parse(data);
            outcome=data[0]['outcome'];
            if (outcome!=true){
                $("#param_outcome").html("<span class=alert>"+data[0]['response']+"</span>");
            }else{
                uploadedFileName=data[0]['response'];
                uploadedFile=uploadsDir+uploadedFileName;
              //  readFile(uploadedFile, function(rawdata){
                    outcome=true;
                    //check on contents if required
                    if (outcome!=true){
                       $("#param_outcome").html("<span class=alert>"+response+"</span>"); 
                    }else{
                        $("#paramFile").val(uploadedFileName);
                        populateRunData();
                    }
              //  });
            }
         }
    });
});



$(document).on('change','#initcondFileSelector', function(event){
    $("#initcond_outcome").html('<img src="./static/img/ajax-loader.gif" alt="Processing...."/>');
    form=$("#initcond_form")[0];
    formdata=new FormData(form);
    $.ajax({
        url: "upload.php",
        method: "POST",
        data: formdata,
        processData: false,
        contentType: false,
        success: function (data) {
            $("#initcond_outcome").html(""); 
        // success callback
            data=JSON.parse(data);
            outcome=data[0]['outcome'];
            if (outcome!=true){
                $("#initcond_outcome").html("<span class=alert>"+data[0]['response']+"</span>");
            }else{
                uploadedFileName=data[0]['response'];
                uploadedFile=uploadsDir+uploadedFileName;
                //readFile(uploadedFile, function(rawdata){
                    outcome=true;
                    //check on contents if required

                    if (outcome!=true){
                       $("#initcond_outcome").html("<span class=alert>"+response+"</span>"); 
                    }else{
                        $("#initcondFile").val(uploadedFileName);
                        populateRunData();
                    }
                //});
            }
         }
    });
});


function populateRunData(){
    paramFile=$("#paramFile").val();
    initcondFile=$("#initcondFile").val();
    inputFile=$("#inputFile").val();
    firstDate=$("#firstDate").val();
    lastDate=$("#lastDate").val();
    nts=$("#nts").val();
    if (inputFile!="" && initcondFile!="" && paramFile!=""){
        txt="<ul>";
        txt+="<form id='run_form' method='POST' action='runmodel.php'>";
        txt+="<h3>Model parameters</h3>";
        txt+="<li>Parameters file:"+paramFile+"<input type=text name=runparamFile id=runparamFile value="+paramFile+" hidden>";
        txt+="<h3>Initial conditions</h3>";
        txt+="<li>Initial conditions file:"+initcondFile+"<input type=text name=runinitcondFile id=runinitcontFile value="+initcondFile+" hidden>";
        txt+="<h3>Input data</h3>";
        txt+="<li>input data file:"+inputFile+"<input type=text name=runinputFile id=runinputFile value="+inputFile+" hidden>";
        txt+="<li>for the period from: "+firstDate+" to "+lastDate+", covering "+nts+" monthly time steps";
        txt+="<h3>Model run properties</h3>";
        txt+="<li>User run ID: <input type=text name=runID id=runID class=activeInput>";
        txt+="<li>Model run Code: "+sessionID+"<input type=text name=sessionID id=sessionID value="+sessionID+" hidden>";
        txt+="<li>Spin-up months (will be skipped from output, 0 for no spin-up): <input type=text name=spinup id=spinup value=0 class=activeInput>";
        txt+="<h3>Model output</h3>";

        for (var key in outputsList){
            txt+="<li><label>"+outputsList[key][0]+"<input type=checkbox class=activeInput name=outputs[] value="+key+"."+outputsList[key][2]+" id="+key+" "+outputsList[key][1]+">";
        }

        txt+="</ul>";
        txt+="<input type=button value=Run id=run_model>";
        txt+="</form>";
        txt+="<div id=run_outcome class=outcome><div>";
     
        $("#run_container").html(txt);
        $("#showhide-run").html("[close]");
        $("#run_container").show();
    }else{
        console.log("some data are still missing");
    }
}


$(document).on('click','#run_model', function(event){
    runID=$("runID").val();
    $("#run_outcome").html("");
    outcome=checkRunData();
    if (outcome==true){
    form=$("#run_form")[0];
    formdata=new FormData(form);
    console.log(formdata);
    $("#run_outcome").html('<img src="./img/ajax-loader.gif" alt="Calculating...." class=centered>');
    $.ajax({
        url: "runmodel.php",
        method: "POST",
        data: formdata,
        processData: false,
        contentType: false,
        success: function (data) {
        // success callback
            data=JSON.parse(data);
            outcome=data[0]['outcome'];
            response=data[0]['response'];
            if (outcome!="success"){
                $("#run_outcome").html(outcome);
            }else{
                txt="";
                for (line in response){
                txt+=response[line]+"<br>";
                }
                $("#run_outcome").html(txt);
                populateResults(runID);
                modelRun=true;
            }
        }
    });
    }else{
        $("#run_outcome").html("<span class=alert>Some run configuration data are missing. RunID is obligatory. At least one output item must be selected</span>");
    }
});


function checkRunData(){
    runID=$("#runID").val();
    nboxes=$("input[name='outputs[]']:checked").length;
    isInt=Number.isInteger(parseInt($("#spinup").val()));
    nts=$("#nts").val();
    if(isInt){
        isLow=parseInt($("#spinup").val())<parseInt(nts);
    }else{
        isLow=false;
    }
    if(runID!="" && nboxes>0 && isInt && isLow){
        outcome=true;
    }else{
        outcome=false;
    }
    return outcome;
}




function populateResults(runID){
    txt="<h3>Available outputs</h3>";
    files=$("input[name='outputs[]']:checked");
    nfiles=files.length;
    for(var c=0, len=nfiles; c < len; c++){
        fileCode=files[c].value;
        key=fileCode.split(".")[0];
        ext=outputsList[key][2];
        resultsFile=resultsDir+"/"+sessionID+"-"+runID+"_"+key+"."+ext;
        txt+="<b>"+outputsList[key][0]+": </b><label> plot <input type=checkbox class=previewResults name=res-"+key+" id=res-"+key+" value="+resultsFile+"></label>";
        txt+=" or download output file: <a href="+resultsFile+">"+sessionID+"-"+runID+"_"+key+"."+ext+"</a><br>";
        txt+="<div id=results_outcome-"+key+" class=outcome></div>";
    }
 
    $("#results_container").html(txt);
    $("#showhide-results").html("[close]");
    $("#results_container").show();
}




$(document).on('change','.previewResults', function(event){
    resID=$(this).attr("id");
    resID=resID.split("-")[1];
    resultFile=$(this).attr("value");
    if ($(this).is(':checked')==false){
        $("#results_outcome-"+resID).html("");
    }else{
        ext=resultFile.substring(resultFile.length - 3, resultFile.length);
        console.log(ext);
        if(ext=="csv"){
            readFile(resultFile, function(rawdata){
                arraydata=parseCSV(rawdata);
                seriesData=[];
                allData=[];
                seriesNames=[];
                ncol=arraydata[0].length;
                nrow=arraydata.length;
                for(var c=1, len=ncol; c < len; c++){
                  //get var names
                    seriesNames.push(arraydata[0][c]);
                    allData.push(new Array());
                }
                for(var r=1, lenr=nrow; r < lenr; r++){
                    dte=new Date(arraydata[r][0]);
                    if (r==1){firstDate=dte;}
                    if (r==nrow-1){lastDate=dte;}
                    dte=dte.valueOf();
                    for(var c=1, lenc=ncol; c < lenc; c++){
                        allData[c-1].push([dte, parseFloat(arraydata[r][c])])
                    }
                }
                for(var c=1, len=ncol; c < len; c++){
                    seriesData.push({id: c, name: seriesNames[c-1], data: allData[c-1], showInLegend:true});
                }
                console.log(resID);

                label=outputsList[resID][0]+"<br>"+outputsList[resID][3];
                outcome=plotTimeSeries("results_outcome-"+resID, seriesData, "Model results", label);
                if (outcome != true){
                   $("#results_outcome-"+resID).html("<span class=alert>"+response+"</span>");
                };
            });
        }else{
                   $("#results_outcome-"+resID).html("<img src="+resultFile+">");
        }
    }
});



$(document).on('click','.activeInput', function(event){
    if (modelRun){
        var retVal = confirm("Changing this will reset model output. Do you want to continue ?");
        if( retVal == true ) {
            $("#results_container").html("");
            $("#run_outcome").html("");
            modelRun=false;
            return true;
        }else{
            return false;
        }
    }

});


function readFile(_file, callback){
    $.get(_file,function(data){
        callback(data);
    });
}


function parseCSV(str) {
    //from https://stackoverflow.com/questions/1293147/javascript-code-to-parse-csv-data
    var arr = [];
    var quote = false;  // true means we're inside a quoted field

    // iterate over each character, keep track of current row and column (of the returned array)
    for (var row = 0, col = 0, c = 0; c < str.length; c++) {
        var cc = str[c], nc = str[c+1];        // current character, next character
        arr[row] = arr[row] || [];             // create a new row if necessary
        arr[row][col] = arr[row][col] || '';   // create a new column (start with empty string) if necessary

        // If the current character is a quotation mark, and we're inside a
        // quoted field, and the next character is also a quotation mark,
        // add a quotation mark to the current column and skip the next character
        if (cc == '"' && quote && nc == '"') { arr[row][col] += cc; ++c; continue; }  

        // If it's just one quotation mark, begin/end quoted field
        if (cc == '"') { quote = !quote; continue; }

        // If it's a comma and we're not in a quoted field, move on to the next column
        if (cc == ',' && !quote) { ++col; continue; }

        // If it's a newline (CRLF) and we're not in a quoted field, skip the next character
        // and move on to the next row and move to column 0 of that new row
        if (cc == '\r' && nc == '\n' && !quote) { ++row; col = 0; ++c; continue; }

        // If it's a newline (LF or CR) and we're not in a quoted field,
        // move on to the next row and move to column 0 of that new row
        if (cc == '\n' && !quote) { ++row; col = 0; continue; }
        if (cc == '\r' && !quote) { ++row; col = 0; continue; }

        // Otherwise, append the current character to the current column
        arr[row][col] += cc;
    }
    return arr;
}



