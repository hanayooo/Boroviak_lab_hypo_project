// FIJI Macro for processing multi-channel z-stack lif files (2D version, no Fill Holes, with labels)
// This macro will:
// 1. Import all series from a lif file
// 2. Perform z-max projection for each channel
// 3. Apply Gaussian blur (sigma=1.5)
// 4. Apply manual thresholds
// 5. Create masks and apply watershed (no fill holes)
// 6. Count and label objects using Analyze Particles
// 7. Save results, screenshots and labeled images
// 8. Output summary CSV

#@ File (label="Select your lif file", style="file") lifFile
#@ Integer (label="Threshold for Channel 1", value=680) threshold1
#@ Integer (label="Threshold for Channel 2", value=20) threshold2
#@ Integer (label="Threshold for Channel 3", value=280) threshold3
#@ Integer (label="Threshold for Channel 4", value=21) threshold4

function processLifFile(lifFilePath, threshold1, threshold2, threshold3, threshold4) {
    thresholds = newArray(threshold1, threshold2, threshold3, threshold4);
    
    lifFileName = File.getName(lifFilePath);
    lifFileDir = File.getParent(lifFilePath);
    baseFileName = substring(lifFileName, 0, lastIndexOf(lifFileName, "."));
    
    resultsDir = lifFileDir + File.separator + baseFileName + "_results";
    File.makeDirectory(resultsDir);
    
    summaryFile = resultsDir + File.separator + baseFileName + "_summary.csv";
    File.append("Series,Channel,Count", summaryFile);
    
    run("Bio-Formats Macro Extensions");
    Ext.setId(lifFilePath);
    Ext.getSeriesCount(seriesCount);
    print("Found " + seriesCount + " series in this lif file");
    
    // 关闭可能存在的所有结果窗口，确保开始时干净
    run("Close All");
    if (isOpen("Results")) {selectWindow("Results"); run("Close");}
    if (isOpen("Summary")) {selectWindow("Summary"); run("Close");}
    if (isOpen("ROI Manager")) {selectWindow("ROI Manager"); run("Close");}
    
    for (s = 0; s < seriesCount; s++) {
        Ext.setSeries(s);
        Ext.getSeriesName(seriesName);
        print("Processing Series #" + (s+1) + ": " + seriesName);
        
        Ext.getSizeC(channelCount);
        Ext.getSizeZ(zSlices);
        
        run("Bio-Formats Importer", "open=[" + lifFilePath + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_" + (s+1));
        originalSeriesID = getImageID();
        
        for (c = 1; c <= channelCount; c++) {
            // 清理ROI Manager以防止累积
            if (isOpen("ROI Manager")) {
                roiManager("reset");
            }
            
            if (!isOpen(originalSeriesID)) {
                run("Bio-Formats Importer", "open=[" + lifFilePath + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_" + (s+1));
                originalSeriesID = getImageID();
            }
            
            selectImage(originalSeriesID);
            run("Duplicate...", "duplicate channels=" + c);
            channelID = getImageID();
            
            if (zSlices > 1) {
                run("Z Project...", "projection=[Max Intensity]");
                projectionID = getImageID();
                selectImage(channelID); close();
                channelID = projectionID;
            }
            
            run("Gaussian Blur...", "sigma=1.5");
            
            channelThreshold = thresholds[Math.min(c-1, thresholds.length-1)];
            setThreshold(channelThreshold, 65535);
            setOption("BlackBackground", true);
            run("Convert to Mask");
            run("Watershed");
            
            processed_filename = resultsDir + File.separator + baseFileName + "_series" + (s+1) + "_channel" + c + "_processed.tif";
            saveAs("Tiff", processed_filename);
            currentImageID = getImageID();

            // 确保彻底清除任何叠加层
            run("Remove Overlay");
            
            // 清理现有结果以防止累积
            run("Clear Results");
            
            // 清理现有的ROI
            if (isOpen("ROI Manager")) {
                roiManager("reset");
            }
        
            // 计数并标记对象
            run("Set Measurements...", "area centroid redirect=None decimal=3");
            run("Analyze Particles...", "size=60-Infinity show=Overlay display summarize add label");

            // 保存带有标签的截图
            labeled_image_filename = resultsDir + File.separator + baseFileName + "_series" + (s+1) + "_channel" + c + "_labeled.png";
            saveAs("PNG", labeled_image_filename);
            
            // 通过重新打开处理过的图像并再次分析来创建一个单独的副本用于保存标记图像（完全避免叠加问题）
            /*
            processingID = getImageID();
            open(processed_filename);
            cleanImageID = getImageID();
            run("Analyze Particles...", "size=60-Infinity show=Overlay display exclude add label");
            labeled_image_filename = resultsDir + File.separator + baseFileName + "_series" + (s+1) + "_channel" + c + "_labeled.png";
            saveAs("PNG", labeled_image_filename);
            selectImage(cleanImageID);
            close();
            selectImage(processingID);
            */
            
            countResult = 0;
            if (nResults() > 0) {
                countResult = nResults();
                run("Clear Results");
            }
            
            File.append(seriesName + "," + "Channel" + c + "," + countResult, summaryFile);
            print("Series #" + (s+1) + ", Channel #" + c + ": Detected " + countResult + " objects");
            
            // 关闭当前图像
            selectImage(currentImageID);
            close();
            
            // 确保除了原始系列图像外的所有图像都被关闭
            while (nImages > 1) {
                closeFlag = false;
                for (i = 1; i <= nImages; i++) {
                    selectImage(i);
                    tempID = getImageID();
                    if (tempID != originalSeriesID) {
                        close();
                        closeFlag = true;
                        break;
                    }
                }
                if (!closeFlag) break;
            }
            
            // 清理ROI Manager
            if (isOpen("ROI Manager")) {
                roiManager("reset");
            }
        }
        
        if (isOpen(originalSeriesID)) {
            selectImage(originalSeriesID);
            close();
        }
    }
    
    Ext.close();
    
    if (isOpen("Log")) {
        selectWindow("Log");
        saveAs("Text", resultsDir + File.separator + baseFileName + "_processing_log.txt");
    }
    
    print("Processing complete. Results saved to: " + resultsDir);
    
    // 确保所有窗口在宏结束时关闭
    run("Close All");
    if (isOpen("Results")) {selectWindow("Results"); run("Close");}
    if (isOpen("Summary")) {selectWindow("Summary"); run("Close");}
    if (isOpen("ROI Manager")) {selectWindow("ROI Manager"); run("Close");}
}

processLifFile(lifFile, threshold1, threshold2, threshold3, threshold4);