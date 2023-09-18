//==============================================================
// The following code processes TIFF files to measure fluorophore intensity & count nuclei of cells
// stained for and expressing specific protein markers.
// User provides channel/marker information (channel ch_number/marker ch_name) & a threshold range
// for measuring flurophore intensity.
//
// Algorithm:
//******** (I) Measurement of protein expression: ***********************************
//1) Preprocessing:
//- If valid image file, open it with Bioformats importer
//- Split channels
//2) Locating expression region:
//- Select protein channel
//- Apply threshold
//- Create Protein Expressor Region Mask
//- Save Mask
//3) Measuring expression:
//- Set Measurements(area mean min max integrated limit redirect=None decimal=0)
//- Measure
//4) Saving result:
//- Collect Results in a custom results table
//
//******** (II) Measurement of background *************************
// !!!NOTE: Used Subtract background function, radius=50 instead of the following steps!!!
//5) Create background region mask
//- Duplicate Protein Expressor Region Mask
//- Rech_name as background masks
//- Invert LUT
//- Measure
//- Collect Results in a custom results table
//
//******** (III) Measurement of protein expressor nuclei count *************************
//
//6) Locating the expressor cell nuclei:
//- Select DAPI channel
//- Duplicate it
//- Make 8-bit
//- AutoThreshold
//- Binarize
//- Watershed segmentation to split connected regions
//- Create DAPI Mask
//- Create overlap (AND) of DAPI Mask with Protein Expressor Region to identify expressor nuclei region
//7) Measuring the expressor cell count:
//- Analyze Particles to count nuclei in this overlap (count appears in Summary window)
//8) Saving result
//- Select Summary window & save count in excel file
//- Save fluorophore measurements in custom results table to excel file
//- Save background fluorophore measurements in custom results table to excel file
//9) Postprocessing:
//- Close all windows
//
//==============================================================


// **************** ENABLE/DISABLE  TEST MODE ******************
TEST_MODE = 1;
DISABLE_BATCH_MODE = 1;
//*****************************************************************

// GLOBAL ACCESS
dir = getDirectory("Choose a Directory");
resultsDir = dir + "Results/";
File.makeDirectory(resultsDir);

masksDir = resultsDir + "Masks/";
File.makeDirectory(masksDir);

files = getFileList(dir);
fileCount = files.length;

// START
main();

// FUNCTION DEFINITIONS
function main() {
   if(TEST_MODE) {
     //********** Use this to skip dialog box and test the code **********
     test_process();
   }
   else {
      presentDialog();
      readInput();
   }
}

function presentDialog() {
    // variables
    experimentID = "";
    cellLine = "";

    ch_name = "";
    ch_num = 0;
    ch_min_threshold = 0;
    ch_max_threshold = 0;
    dapi_ch_num = 0;
    cd144_ch_num = 0;
    cd144_min_threshold = 0;
    cd144_max_threshold = 0;

    Dialog.create("EC Image Data Extractor");
    Dialog.addString("Enter Experiment ID", experimentID);
    Dialog.addString("Enter Cell Line", cellLine);

    // P1 Channel
    Dialog.addString("Enter Channel Name", ch_name);
    Dialog.addNumber("Enter Channel Number", ch_num);
    Dialog.addNumber("Enter Min Threshold", ch_min_threshold);
    Dialog.addNumber("Enter Max Threshold", ch_max_threshold);

    Dialog.addNumber("Enter DAPI Channel Number", dapi_ch_num);

    // CD144 Channel
    Dialog.addNumber("Enter CD144 Channel Number", cd144_ch_num);
    Dialog.addNumber("Enter Min Threshold", cd144_min_threshold);
    Dialog.addNumber("Enter Max Threshold", cd144_max_threshold);

    Dialog.show();
}

function readInput(){
    // Read input
    experimentID = Dialog.getString();
    cellLine = Dialog.getString();

    ch_name = Dialog.getString();
    ch_num = Dialog.getNumber();
    ch_min_threshold = Dialog.getNumber();
    ch_max_threshold = Dialog.getNumber();

    dapi_ch_num = Dialog.getNumber();

    cd144_ch_num = Dialog.getNumber();
    cd144_min_threshold = Dialog.getNumber();
    cd144_max_threshold = Dialog.getNumber();

    process(experimentID,
    cellLine,
    ch_name,
    ch_num,
    ch_min_threshold,
    ch_max_threshold,
    dapi_ch_num,
    cd144_ch_num,
    cd144_min_threshold,
    cd144_max_threshold);
}

function process(experimentID,
                  cellLine,
                  ch_name,
                  ch_num,
                  ch_min_threshold,
                  ch_max_threshold,
                  dapi_ch_num,
                  cd144_ch_num,
                  cd144_min_threshold,
                  cd144_max_threshold) {
    // Table to track fluorophore intensity & area measurements
    resultsTableName = "Results_" + ch_name + "_measurements";
    Table.create(resultsTableName);

    // Table to track nuclei counts
    dapiTableName = "DAPI_Counts";
    Table.create(dapiTableName);

    // Table to track protein expressor cells' nuclei counts
    expressorDapiTableName = ch_name + "_DAPI_Counts";
    Table.create(expressorDapiTableName);

    // Table to track cd144 nuclei counts
    cd144DapiTableName = "CD144_DAPI_Counts";
    Table.create(cd144DapiTableName);

    // Table to track image name and condition
    imageInfoTableName = "Image Metadata";
    Table.create(imageInfoTableName);

    setOption("ExpandableArrays", true);
    areas = newArray();
    means = newArray();
    intden = newArray();
    rawintden = newArray();
    minThr = newArray();
    maxThr = newArray();

    dapiCounts = newArray();
    expressorDapiCounts = newArray();
    cd144DapiCounts = newArray();

    fileNames = newArray();
    conditions = newArray();

    updateImageInfoTable(fileNames, conditions, imageInfoTableName);
    updateResultsTable(areas, means, intden, rawintden, minThr, maxThr, resultsTableName);
    updateDapiTable(dapiCounts, dapiTableName);
    updateCD144DapiTable(cd144DapiCounts, cd144DapiTableName);
    updateExpressorDapiTable(expressorDapiCounts, expressorDapiTableName, ch_name);

    // **************************** (I) Measurement of protein expression  ****************************
    if(ch_num > 0){

        if(DISABLE_BATCH_MODE) {
           // *********** Seting number of files to process ***************
           fileCount = 1;
        }

        for(i=0; i<fileCount; i++) {
            fileName = files[i];
            if (endsWith(fileName, ".tif")){
                fileNames[i] = fileName;

                run("Bio-Formats", "open=[" + dir+fileName +"] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

                imgName = getTitle();
                if (imgName != NaN) {
                    conditions[i] = getConditionInfo(imgName);
                }

               updateImageInfoTable(fileNames, conditions, imageInfoTableName);

                // Set Measurements
                run("Set Measurements...", "area mean integrated limit redirect=None decimal=0");

                // Split Channels
                run("Split Channels");

                if(ch_num > 0)
                {
                    // Select Channel
                    pWin = "C"+ ch_num +"-"+ imgName;
                    pWin2 = pWin + "_dup";
                    selectWindow(pWin);

                    run("Duplicate..."," ");
                    rename(pWin2);


                    // **************************** (II) Measurement of background *************************
                    selectWindow(pWin);
                    run("Subtract Background...", "rolling=50");

                    //**************************** Measure fluorophore intensity in expressors ****************
                    // Set Threshold
                    setThreshold(ch_min_threshold, ch_max_threshold);

                    // Measure
                    run("Measure");

                    areas[i] = getResult("Area",0);
                    means[i] = getResult("Mean",0);
                    intden[i] = getResult("IntDen",0);
                    rawintden[i] = getResult("RawIntDen",0);
                    minThr[i] = getResult("MinThr",0);
                    maxThr[i] = getResult("MaxThr",0);
                    updateResultsTable(areas, means, intden, rawintden, minThr, maxThr, resultsTableName);

                    selectWindow("Results");
                    run("Clear Results");
                }

                //**** Create Expressors region mask ***
                //Create Mask
                run("Convert to Mask");

                // Invert LUT
                run("Invert LUT");

                pmask = ch_name +"_Mask_"+i + ".tif";
                pmaskFilePath = masksDir+pmask;
                saveAs("tif", pmaskFilePath);

                // **************************** (III) Measurement of protein expressor nuclei count  ****************************
                if(dapi_ch_num > 0)
                {
                    // Select Channel
                    dWin = "C"+ dapi_ch_num +"-"+ imgName;
                    pMaskWin = "C"+ dapi_ch_num +"-" + ch_name + imgName;
                    selectWindow(dWin);
                    run("Duplicate..."," ");
                    rename(pMaskWin);

                    // 8-bit
                    run("8-bit");

                    // Set Threshold
                    setAutoThreshold();

                    //Make Binary
                    run("Make Binary");

                    // Segmentation
                    run("Watershed");

                    //Create Mask
                    run("Convert to Mask");

                    // Mask for counting all nuclei  => dapi_mask
                    dapi_mask = "dapi_" + i + ".tif";
                    dapiMaskFilePath = masksDir + dapi_mask;
                    saveAs("tif", dapiMaskFilePath);

                    // Mask for counting expressor nuclei  => pdMask
                    pdMask = ch_name + dapi_mask;
                    run("Duplicate..."," ");
                    rename(pdMask);

                    //************************ Total nuclei count  ****************************
                    selectWindow(dapi_mask);
                    run("Analyze Particles...","size=10-Infinity circularity=0.00-1.00 show=Nothing display clear summarize");
                    selectWindow("Summary");
                    dapiCounts[i] = Table.get("Count",0,"Summary");
                    updateDapiTable(dapiCounts, dapiTableName);
                    removeLastRowOfSummary();

                    //************************ CD144 nuclei count ****************************
                    if(cd144_ch_num>0){
                        // select cd144 channel
                        cd144Win = "C"+ cd144_ch_num +"-"+ imgName;
                        selectWindow(cd144Win);

                        // Set Threshold
                        setThreshold(cd144_min_threshold, cd144_max_threshold);

                        //Create Mask
                        run("Convert to Mask");

                        // Mask for counting cd144 nuclei  => cd144_region_mask
                        cd144_region_mask = "cd144_" + i + ".tif";
                        cd144RegionMaskFilePath = masksDir + cd144_region_mask;
                        saveAs("tif", cd144RegionMaskFilePath);

                        selectWindow(dapi_mask);
                        // Get overlap region of total nuclei and protein expressors
                        imageCalculator("AND create", dapi_mask, cd144_region_mask);
                        cd144_nuclei_mask = "cd144_nuclei_"+ i + ".tif";
                        cd144NucleiMaskFilePath = masksDir + cd144_nuclei_mask;
                        saveAs("tif", cd144NucleiMaskFilePath);

                        // count nuclei in the overlap region
                        selectWindow(cd144_nuclei_mask);
                        run("Analyze Particles...","size=10-Infinity circularity=0.00-1.00 show=Nothing display clear summarize");

                        selectWindow("Summary");
                        cd144DapiCounts[i] = Table.get("Count",0,"Summary");
                        updateCD144DapiTable(cd144DapiCounts, cd144DapiTableName);
                        removeLastRowOfSummary();
                        selectWindow("Results");
                        run("Clear Results");
                    }

                    // **************************** Expressor nuclei count  ****************************
                    selectWindow(pdMask);
                    // Get overlap region of total nuclei and protein expressors
                    imageCalculator("AND create", pdMask, pmask);

                    pn_area_mask = ch_name+"_nuclei_region_"+i + ".tif";
                    pnMaskFilePath = masksDir + pn_area_mask;
                    saveAs("tif", pnMaskFilePath);

                    // count nuclei in the overlap region
                    selectWindow(pn_area_mask);
                    run("Analyze Particles...","size=10-Infinity circularity=0.00-1.00 show=Nothing display clear summarize");

                    selectWindow("Summary");
                    expressorDapiCounts[i] = Table.get("Count",0,"Summary");
                    updateExpressorDapiTable(expressorDapiCounts, expressorDapiTableName, ch_name);
                    removeLastRowOfSummary();

                    pnAreaMaskFilePath = masksDir+ pn_area_mask;
                    saveAs("tif", pnAreaMaskFilePath);

                    selectWindow(pn_area_mask);
                    close();
                    selectWindow(dapi_mask);
                    close();
                    selectWindow(pmask);
                    close();
                    selectWindow("Results");
                    run("Clear Results");
                }
             } //end if
       }//end for
    }// end if

   saveData(resultsTableName, dapiTableName, expressorDapiTableName);
}

function saveData(resultsTableName, dapiTableName, expressorDapiTableName) {
    selectWindow(resultsTableName);
    resultsFileName = experimentID + "_" + cellLine + "_" + resultsTableName;
    pathToOutputFile = resultsDir + resultsFileName + ".csv";
    saveAs("Results", pathToOutputFile);
    close(resultsFileName+ ".csv");

    selectWindow(expressorDapiTableName);
    expressorNucleiFileName = experimentID + "_" + cellLine + "_" + expressorDapiTableName;
    pathToExpressorNucleiOutputFile = resultsDir + expressorNucleiFileName;
    saveAs("Results", pathToExpressorNucleiOutputFile+".csv");
    close(expressorNucleiFileName+".csv");

    selectWindow(cd144DapiTableName);
    cd144NucleiFileName = experimentID + "_" + cellLine + "_" + cd144DapiTableName;
    pathToCD144NucleiOutputFile = resultsDir + cd144NucleiFileName;
    saveAs("Results", pathToCD144NucleiOutputFile+".csv");
    close(cd144NucleiFileName+".csv");

    selectWindow(dapiTableName);
    allNucleiFileName = experimentID + "_" + cellLine + "_" + dapiTableName;
    pathToAllNucleiOutputFile = resultsDir + allNucleiFileName;
    saveAs("Results", pathToAllNucleiOutputFile +".csv");
    close(allNucleiFileName+".csv");

    selectWindow(imageInfoTableName);
    allFileInfoName = experimentID + "_" + cellLine + "_" + imageInfoTableName;
    pathToFileInfoOutputFile = resultsDir + allFileInfoName;
    saveAs("Results", pathToFileInfoOutputFile +".csv");
    close(allFileInfoName+".csv");

    close("Results");
    close("Summary");
    close("*");
}

function updateResultsTable(areas,means,intden,rawintden,minThr,maxThr, resultsTableName) {
    Table.setColumn("Area", areas, resultsTableName);
    Table.setColumn("Mean", means, resultsTableName);
    Table.setColumn("IntDen",intden, resultsTableName);
    Table.setColumn("RawIntDen",rawintden, resultsTableName);
    Table.setColumn("MinThr", minThr, resultsTableName);
    Table.setColumn("MaxThr", maxThr, resultsTableName);
    Table.update(resultsTableName);
}

function updateDapiTable(dapiCounts, dapiTableName) {
    Table.setColumn("Total Nuclei", dapiCounts, dapiTableName);
    Table.update(dapiTableName);
}

function updateExpressorDapiTable(expressorDapiCounts, expressorDapiTableName, markerName) {
    columnName = markerName +" Expressor Nuclei";
    Table.setColumn(columnName, expressorDapiCounts, expressorDapiTableName);
    Table.update(expressorDapiTableName);
}

function updateCD144DapiTable(cd144DapiCounts, cd144DapiTableName){
    Table.setColumn("CD144 nuclei count", cd144DapiCounts, cd144DapiTableName);
    Table.update(cd144DapiTableName);
}

function updateImageInfoTable(imageNames, conditions, imageInfoTableName) {
    Table.setColumn("Filename", imageNames, imageInfoTableName);
    Table.setColumn("Condition", conditions, imageInfoTableName);
    Table.update(imageInfoTableName);
}

function removeLastRowOfSummary(){
  nRows = Table.size("Summary");
  Table.deleteRows(0, nRows-1, "Summary");
}


function getConditionInfo(imgName){
   condition = "Untreated";

   if(imgName.contains("DAPT") || imgName.contains("dapt")) {
      return "DAPT";
   }
   else if(imgName.contains("CpE") || imgName.contains("cpe")) {
      return "CpE";
   }
   else if(imgName.contains("Plasma") || imgName.contains("plasma")) {
      return "10% vWF Deficient Plasma";
   }
   else if(imgName.contains("Control") || imgName.contains("UT") || imgName.contains("cntrl")) {
      return "Untreated";
   }
   return condition;
}

//==================== TESTS ========================

function test_process() {
    process("Hx69",
    "hPSEC",
    "VWF",
    1,
    0,
    3236,
    4,
   2,
   0,
   1867);
}
