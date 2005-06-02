Readme for tightClust

Compile:
Put the following the files in the same directory. Link and compile using 
	gcc -o tightClust tightClust_ANSI.c ranlib.c linpack.c com.c cluster.c
An executable file "tightClust" will be generated.

Run:
Run the program using commands like:
	tightClust c:\temp\ testData.txt 15 25 35 0.1 0.6 7 2 10 0.7 3
	tightClust c:\temp\ testData.txt 15 35 45 0.1 0.6 7 2 10 0.7 3
For details, check the manual of in the folder "./windows". The file "tightClust.exe" in "./windows" is compiled in Windows using Visual C++.

                   Manual for TightClust (ANSI C version)
Data format:
        See the attached file, "testData.txt". It comes with the number of genes and then
number of samples in the first row. Then sample names at the second rows and then data
matrix. The first column is reserved for geneID and second for annotation.
Note there should be no space after the last number at each row.
Command:
        Go to the directory when "tightClust.exe" resides. Type commands like the
following.
tightClust c:\temp\ testData.txt 15 25 35 0.1 0.6 7 2 10 0.7 3
tightClust c:\temp\ testData.txt 15 35 45 0.1 0.6 7 2 10 0.7 3
        The tight clustering result will be saved in a file like
"tightClust_15_25_35_0.1_0.6_7_2_10_0.7_3.txt". The output information on the screen
will be saved to a log file like "temp_logFile_15_25_35_0.1_0.6_7_2_10_0.7_3.txt".
Clusters found are immediately saved to a temp cluster file named as
"temp_clustFile.txt". This is to let users to explore the first few tight clusters before the
program finishes. The following 12 parameters are needed for tightClust.
    1. workingDir: the working directory where data resides. All output files will also
        be saved in this directory.
    2. dataFileName: the data matrix. See data format for details. A file named
        "temp_confirmData.txt" is saved for users to confirm if their data are read in
        correctly.
    3. targetClustNum: the total number of clusters that the user aims to find
    4. min_k: the starting point of k_0. It should be larger than targetClustNum+10. The
        larger min_k will find smaller and tighter clusters.
    5. max_k: the ending point of k_0. It should be larger than min_k.
    6. alpha: the threshold of comembership index.
    7. beta: the threshold of clusters stably found in consecutive k_0.
    8. topNum: the number of top (size) candidate clusters for a specific k_0.
    9. seqNum: the number of subsequent k_0 that finds the tight cluster.
    10. resampNum: total number of resampling to obtain comembership matrix.
    11. subSampPercnt: percentage of subsamples selected
    12. npass: number of different random inital for
