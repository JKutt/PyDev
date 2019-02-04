H. LARNIER - 03 February 2019 - hugo.larnier@diasgeo.com for suggestions, inquiries..
Will be subject to revisions often!
------

How to use the field data viewer:

To open the software: Double click on "DIAS Field Data Viewer"

MISSING DATA:
To check for missing data, the software will compare a file database and the node messages:
	1. Click on "File", "Load DATA". A folder selection window will open. Select the folder containing the node data. (for the Las Bambas project, the folder is located at: C:\Users\FieldUser\Desktop\DATA_LAS_BAMBAS\DATA)
	2. In the main software window, click on "Tools", "Check missing data". A file selection window will open. Select the node messages file corresponding to the line you want to check (or the grid) (for the Las Bambas project, the node messages files are located at: C:\Users\FieldUser\Desktop\DATA_LAS_BAMBAS\logs, select for example "NodeMessages_73.txt" to check for missing files for injection line 73)
	3. A folder selection window will open. Select the folder where you want the list to be written.
	4. Check the list and see if node files are missing for nodes you just dumped.
	5. If node files are missing, do better.

GPX FILES:
To automatically make GPX files, you have two options: nodes files clustering, or nodes messages clustering.
1. For nodes files clustering:
	1. Click on "File", "Load DATA". A folder selection window will open. Select the folder containing the node data. (for the Las Bambas project, the folder is located at: C:\Users\FieldUser\Desktop\DATA_LAS_BAMBAS\DATA).
	2. In the main software window, click on "Tools", "GPS locations". A new window will open.
	3. Click on "Clustering", "Cluster nodes files".
	4. Wait a bit.
	5. A window will open and ask you where you want to save the GPX file.
	6. Check for locations.

2. For nodes messages clustering:
	1. In the main software window, click on "Tools", "GPS locations". A new window will open.
	2. Click on "Clustering", "Cluster nodes messages".
	3. Wait a bit.
	4. A window will open and ask you where you want to save the GPX file.
	5. Check for locations.

GPX files can be viewed in the "GPS locations" window ("File", "Load GPX files"), or in MapSource. The nodes files clustering will output as a comment the distribution of node ID at a same location.

TO VIEW CURRENT DATA:
	1. Click on "File", "Load DATA". A folder selection window will open. Select the folder containing the node data. (for the Las Bambas project, the folder is located at: C:\Users\FieldUser\Desktop\DATA_LAS_BAMBAS\DATA)
	2. Click on "File", "Load recording file", and select the recording file you want to check.
	3. Click on "File", "Read folder". The software will then read the DATA folder and assign node files to each injection in the recording file.
	4. Click "View", "Single injection map". A new window will open.
	5. On the left side panel, select the injection you wish to check. A map will be drawn, as well as current files available for this injection. WARNING: For now the injection number doesn't correspond to the MEM number of the injection. This will be solved in the future.

TO VIEW DATA:
	1. Click "View", "Data information". A new window will open.	
	2. Click "File", "Load single file".
	3. Data will be drawn.

----
To come:
 - Automatic check of simple stuff on the field after a day (current variations, current leaks between TX and CC)
 - Possibility to read node messages to check that everything was alright after the day (grid integrity, grid connection...)
 - Display of field Vp and Vs.
 - Get rid of MapSource for GPS processing.