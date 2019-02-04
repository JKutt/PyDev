import numpy as np
import matplotlib.pylab as plt
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
import io_lib as iolib
from io_lib import Record
from glob import glob
from mpl_toolkits import mplot3d
import gpxpy
import gpxpy.gpx
import utm


def gps_from_minutes_to_decimal(gps_value):
    # Get decimal coordinates from decimal minutes
    # Inputs:
    # gps_value: gps coordinates in decimal minutes
    # Output:
    # gps value in decimal value
    degree = np.floor(gps_value)
    tmp = (gps_value - degree) * 100
    minutes = np.floor(tmp)
    seconds = ((tmp - minutes) * 60.)

    return degree + minutes / 60. + seconds / 3600.

def cluster_nodes_messages(nodes_messages_file):
	#nodes_messages_file = 'C:/Users/FieldUser/Desktop/Tools/logs/NodeMessages_complete.txt'
	print('Starting clustering on file: ' + nodes_messages_file)
	fIn = open(nodes_messages_file, 'r')
	lines = fIn.readlines()
	fIn.close()

	list_files = []

	X = []
	Clines = []
	Cstations = []
	altitudes = []
	id = []
	mem = []
	name = []
	coordinates = []
	for line in lines:
		if 'i;' in line:
			#print(line)
			try:
				tmp = line.split(';G')
				lat = float(tmp[1].split(',')[0])
				lon = float(tmp[1].split(',')[2])
				lon = gps_from_minutes_to_decimal(lon / 100.)
				lat = gps_from_minutes_to_decimal(lat / 100.)
				if 	tmp[1].split(',')[1] == 'S':
					lat *= -1
				if 	tmp[1].split(',')[3] == 'W':
					lon *= -1
				#tmp = utm.from_latlon(lat, lon)
				alt = tmp[1].split(',')[4].split(';')[0]
				Line = line.split(';l')[1].split(';')[0]
				Station = line.split(';n')[1].split(';')[0]
				ID = line.split('ID')[1].split(':')[1][1:-1]
				
				#print(lat, lon, alt, Line, Station, ID)
				# loading in clustering variables
				tmp = utm.from_latlon(lat, lon)
				#print(tmp)			
				altitudes.append(float(alt))
				zone = [tmp[2], tmp[3]]
				X.append([tmp[0], tmp[1]])
				coordinates.append([tmp[0], tmp[1], float(alt)])
				id.append(ID)
				Clines.append(int(Line))
				Cstations.append(int(Station))
			except:
				pass
				
	# X = np.asarray(X)
	# Clines = np.asarray(Clines)
	# Cstations = np.asarray(Cstations)
	# id = np.asarray(id)
	# altitudes = np.asarray(altitudes)
	
	# print('Starting clustering')
	# db = DBSCAN(eps=50, min_samples=20).fit(X)
	# core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
	# core_samples_mask[db.core_sample_indices_] = True
	# labels = db.labels_

	# # Number of clusters in labels, ignoring noise if present.
	# n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

	# print('Estimated number of clusters: %d' % n_clusters_)
										
	size = 100
	std_clustering = 35

	X = np.asarray(X)
	altitudes = np.asarray(altitudes)
	lines = np.asarray(Clines)
	stations = np.asarray(Cstations)
	id = np.asarray(id)
	#mem = np.asarray(mem)
	db = DBSCAN(eps=std_clustering, min_samples=5).fit(X)
	core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
	core_samples_mask[db.core_sample_indices_] = True
	labels = db.labels_

	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

	print('Estimated number of clusters: %d' % n_clusters_)

	unique_labels = set(labels)
	colors = [plt.cm.Spectral(each)
			  for each in np.linspace(0, 1, len(unique_labels))]
	
	x = []
	y = []
	z = []
	line_cluster = []
	station_cluster = []
	ID = []
	MEM = []
	print(len(X), len(altitudes), len(lines), len(stations))
	for k in unique_labels:
		class_member_mask = (labels == k)
		xy = X[class_member_mask]
		alt = altitudes[class_member_mask& core_samples_mask]
		L = lines[class_member_mask& core_samples_mask]
		S = stations[class_member_mask& core_samples_mask]
		ID.append(id[class_member_mask& core_samples_mask])
		#MEM.append(mem[class_member_mask& core_samples_mask])
		x.append(np.median(xy[:, 0]))
		y.append(np.median(xy[:, 1]))
		z.append(np.median(alt))
		line_cluster.append(np.median(L))
		station_cluster.append(np.median(S))
	##ax.plot_trisurf(x, y, z, cmap='viridis')

	#ax.scatter(x, y, z, c=z, cmap='viridis')
	#plt.show()

	file = 'C:/Users/FieldUser/Desktop/Tools/clusteringTest/GPS_nodes_messages.csv'
	fOut = open(file, 'w')

	for i in range(len(x)):
		fOut.write(str(line_cluster[i]) + ',' + str(station_cluster[i]) + ',' + str(x[i]) + ',' + str(y[i]) + ',' + str(z[i]) + '\n')
	fOut.close()	



	gpx = gpxpy.gpx.GPX()

	for i in range(len(x)):
		coord = utm.to_latlon(x[i], y[i], zone[0], zone[1])		
		try:
			gpx.waypoints.append(gpxpy.gpx.GPXWaypoint(coord[0],\
																							  coord[1],\
																							  elevation=z[i],
																							  name=str(int(line_cluster[i])) + ' ' + str(int(station_cluster[i])) + ' ID'))
		except:
			pass
	return gpx
	#fOut = open('GPX_nodes_messages.gpx', 'w')
	#fOut.write(gpx.to_xml())
	#fOut.close()	

def cluster_data_folder(recording_file, data_path):

	records = []
	if data_path and recording_file:
	    print("Reading...")
	    logFile = open(data_path + 'checkDataFolder.log', 'w')

	    list_nodes = iolib.list_nodes(data_path)

	    for node in list_nodes:
	        print("Checking node: " + node)
	        logFile.write("Checking node: " + node + '\n')
	        datNode = iolib.get_list_files(data_path, node)
	        iolib.check_list_files(datNode, logFile)
	    logFile.close()
	    print("Done! You can now check the file " + data_path + "checkDataFolder.log")

	    if len(glob(data_path + 'nodeLibrary.dat')):
	        records_all = iolib.read_library_file(data_path + 'nodeLibrary.dat')
	        records_all, nodeList = iolib.add_new_nodes(records_all, list_nodes)
	    else:
	        records_all = [[] for i in range(len(list_nodes))]
	    # Reading recording file
	    injections = iolib.read_log_file(recording_file)

	    nodeIndex = 0
	    # Making library of records to compare with injections
	    for node in list_nodes:
	        print("Processing node: " + node)
	        datNode = iolib.get_list_files(data_path, node)
	        fileIndex = 0
	        for datFile in datNode:  # loop on DAT files in folder
	            test = iolib.is_file_in_library(records_all, data_path + node + '/' + datFile, nodeIndex)
	            if test == 1:
	                print("\tReading file: " + datFile)
	                fIn = open(data_path + node + '/' + datFile)
	                linesFIn = fIn.readlines()
	                fIn.close()
	                info = iolib.get_dat_info(linesFIn)               # Getting DAT file info from header
	                timeStamp = iolib.get_start_end_time(linesFIn)    # Getting time stamp
	                gpsLocations = iolib.get_average_gps(linesFIn)    # Getting average GPS location from DAT file
	                try:
	                    records_all[nodeIndex].append(Record(node_id=info[0],
	                                                      node_type=info[3],
	                                                      name=data_path + node + '/' + datFile,
	                                                      mem=info[1],
	                                                      start_date=timeStamp[0],
	                                                      end_date=timeStamp[1],
	                                                      relay_state=info[2],
	                                                      northing=gpsLocations[0],
	                                                      easting=gpsLocations[1],
	                                                      altitude=gpsLocations[2],
														  line=info[6],
														  station=info[7]))
	                except IndexError:
	                    pass
	            else:
	                print("\t File: " + datFile + " is already in the library")
	            fileIndex += 1
	        nodeIndex += 1
	    # Comparing recording file and list of records
	    records = [[] for i in range(len(records_all))]
	    for injection in injections:
	    		for nodeIndex in range(len(records_all)):
	        		for record in records_all[nodeIndex]:
	        		    test = iolib.is_node_active(injection, record)
	        		    if test:
	        		        records[nodeIndex].append(record)
				
	## Current stations
	# Clustering will issue only Tx sites, non clustered data are injections sites.
	X = []
	Clines = []
	Cstations = []
	altitudes = []
	id = []
	mem = []
	name = []
	for injection in injections:
	    if injection.valid == ' Good':
	        for nodeIndex in range(len(records_all)):
	            	for record in records_all[nodeIndex]:
	        	        if record.node_type == 'C':
	        	            test = iolib.is_node_active(injection, record)
	        	            if test:
	        	            		    try:
	        	            			    tmp = utm.from_latlon(record.northing, record.easting)
	        	            			    altitudes.append(record.altitude)
	        	            			    zone = [tmp[2], tmp[3]]
	        	            			    X.append([tmp[0], tmp[1]])
	        	            			    id.append(record.node_id)
	        	            			    Clines.append(record.line)
	        	            			    Cstations.append(record.station)
	        	            			    mem.append(injection.num)        	            										
	        	            			    name.append(record.name)        	            										
	        	            		    except:
	        	            			    pass


	X = np.asarray(X)
	Clines = np.asarray(Clines)
	Cstations = np.asarray(Cstations)
	id = np.asarray(id)
	mem = np.asarray(mem)
	name = np.asarray(name)
	altitudes = np.asarray(altitudes)
	db = DBSCAN(eps=50, min_samples=20).fit(X)
	core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
	core_samples_mask[db.core_sample_indices_] = True
	labels = db.labels_

	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

	print('Estimated number of clusters: %d' % n_clusters_)

	unique_labels = set(labels)
	colors = [plt.cm.Spectral(each)
	          for each in np.linspace(0, 1, len(unique_labels))]

	"""
	for k, col in zip(unique_labels, colors):
	    if k == -1:
	        # Black used for noise.
	        col = [0, 0, 0, 1]

	    class_member_mask = (labels == k)
	    xy = X[class_member_mask & core_samples_mask]
	    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
	             markeredgecolor='k', markersize=6)

	    xy = X[class_member_mask & ~core_samples_mask]
	    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
	             markeredgecolor='k', markersize=6)
	
	nonTx = X[class_member_mask & ~core_samples_mask]
	Clines = Clines[class_member_mask & ~core_samples_mask]
	Cstations = Cstations[class_member_mask & ~core_samples_mask]
	id = id[class_member_mask & ~core_samples_mask]
	name = name[class_member_mask & ~core_samples_mask]
	mem = mem[class_member_mask & ~core_samples_mask]
	altitudes = altitudes[class_member_mask & ~core_samples_mask] 

	gpx = gpxpy.gpx.GPX()
	for i in range(len(Clines)):
	    coord = utm.to_latlon(nonTx[i, 0], nonTx[i, 1], zone[0], zone[1])
	    gpx.waypoints.append(gpxpy.gpx.GPXWaypoint(coord[0], coord[1], elevation=altitudes[i], name=str(Clines[i]) + ' ' + str(Cstations[i]) + ' ' + str(id[i])))
		
	fOut = open('currentRecords.gpx', 'w')
	fOut.write(gpx.to_xml())
	fOut.close()
	file = 'C:/Users/FieldUser/Desktop/Tools/clusteringTest/GPS_currentRecorders.csv'
	fOut = open(file, 'w')
	for i in range(len(Clines)):
		fOut.write(str(Clines[i]) + ',' + str(Cstations[i]) + ',' + str(nonTx[i, 0]) + ',' + str(nonTx[i, 1]) + ',' + str(altitudes[i]) + ',' + str(id[i]) + '\n')
	fOut.close()	
	#plt.title('Estimated number of clusters: %d' % n_clusters_)
	#plt.show()
										
	"""
	
	# Using recording time
	for injection in injections:
	    id = []
	    line = []
	    station = []
	    coordinates = []
	    for nodeIndex in range(len(records_all)):
	        	for record in records_all[nodeIndex]:
	        	    if record.node_type == 'C':
	        	        test = iolib.is_node_active(injection, record)
	        	        if test:
	        	            if record.node_id in id:
	        	            		ind = id.index(id)
	        	            		line[ind].append(record.line)
	        	            		station[ind].append(record.station)
	        	            		try:
	        	            			tmp = utm.from_latlon(record.northing, record.easting)
	        	            			coordinates[ind].append([tmp[0], tmp[1], record.altitude])	
	        	            		except:
	        	            			pass
		

	print("Calculating distances")
	distances = []
	altitudes = []
	line = []
	stations = []
	coordinates = []
	current = ['EB', 'DM', 'DP']
	for i in range(len(records)):
	    for j in range(int(np.min([len(records[i]), 5]))):
	        #if not records[i][j].node_id in current:
	        if records[i][j].node_type == 'V':# or records[i][j].node_type == 'C':
	            try:
	                tmp = utm.from_latlon(records[i][j].northing, records[i][j].easting)
	                coordinates.append([tmp[0], tmp[1], records[i][j].altitude])
	            except:
	                pass

	for i in range(len(coordinates)):
	    for j in range(len(coordinates)):
	        if not i == j:
	            distances.append(np.sqrt((coordinates[i][0] - coordinates[j][0]) ** 2 + (coordinates[i][1] - coordinates[j][1]) ** 2 ))
	            altitudes.append(coordinates[i][2])

	distances = np.asarray(distances)
	distances = distances.reshape(-1, 1)
	print("Starting clustering")
	db = DBSCAN(eps=1, min_samples=10).fit(distances)
	print("End clustering")
	core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
	core_samples_mask[db.core_sample_indices_] = True
	labels = db.labels_

	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

	print('Estimated number of clusters: %d' % n_clusters_)
	unique_labels = set(labels)
	clusters = []
	var_clusters = []
	size_clusters = []
	for k in unique_labels:

	    class_member_mask = (labels == k)
	    xy = distances[class_member_mask & core_samples_mask]
	    clusters.append(np.median(xy))
	    var_clusters.append(np.std(xy))
	    size_clusters.append(len(xy))


	ind = np.argsort(clusters)

	for cluster, var_cluster, size_cluster in zip(clusters, var_clusters, size_clusters):
	    if cluster > 10:
	        size = size_cluster
	        std_clustering = var_cluster
	        break

	X = []
	altitudes = []
	lines = []
	stations = []
	id = []
	mem = []
	for i in range(len(records)):
	    for record in records[i]:
	        if record.node_type == 'V':
	            try:
	                tmp = utm.from_latlon(record.northing, record.easting)
	                X.append([tmp[0], tmp[1]])
	                altitudes.append(record.altitude)
	                lines.append(record.line)
	                stations.append(record.station)
	                id.append(record.node_id)
	                mem.append(record.mem)
	            except:
	                pass


	X = np.asarray(X)
	altitudes = np.asarray(altitudes)
	lines = np.asarray(lines)
	stations = np.asarray(stations)
	id = np.asarray(id)
	mem = np.asarray(mem)
	db = DBSCAN(eps=std_clustering, min_samples=5).fit(X)
	core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
	core_samples_mask[db.core_sample_indices_] = True
	labels = db.labels_

	# Number of clusters in labels, ignoring noise if present.
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

	print('Estimated number of clusters: %d' % n_clusters_)

	unique_labels = set(labels)
	colors = [plt.cm.Spectral(each)
	          for each in np.linspace(0, 1, len(unique_labels))]
	"""
	for k, col in zip(unique_labels, colors):
	    if k == -1:
	        # Black used for noise.
	        col = [0, 0, 0, 1]

	    class_member_mask = (labels == k)
	    xy = X[class_member_mask & core_samples_mask]
	    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
	             markeredgecolor='k', markersize=6)

	    xy = X[class_member_mask & ~core_samples_mask]
	    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
	             markeredgecolor='k', markersize=6)

	plt.title('Estimated number of clusters: %d' % n_clusters_)
	plt.show()
	"""
	"""
	plt.figure()
	for k, col in zip(unique_labels, colors):
	    if k == -1:
	        # Black used for noise.
	        col = [0, 0, 0, 1]

	    class_member_mask = (labels == k)
	    z = altitudes[class_member_mask & core_samples_mask]
	    xy = X[class_member_mask & core_samples_mask]
	    plt.plot(xy[:, 0], z, 'o', markerfacecolor=tuple(col),
	             markeredgecolor='k', markersize=6)
	plt.show()
	"""
	#fig = plt.figure()
	#ax = plt.axes(projection='3d')
	x = []
	y = []
	z = []
	line_cluster = []
	station_cluster = []
	ID = []
	MEM = []
	for k in unique_labels:
	    class_member_mask = (labels == k)
	    xy = X[class_member_mask]
	    alt = altitudes[class_member_mask& core_samples_mask]
	    L = lines[class_member_mask& core_samples_mask]
	    S = stations[class_member_mask& core_samples_mask]
	    ID.append(id[class_member_mask& core_samples_mask])
	    MEM.append(mem[class_member_mask& core_samples_mask])
	    x.append(np.median(xy[:, 0]))
	    y.append(np.median(xy[:, 1]))
	    z.append(np.median(alt))
	    line_cluster.append(np.median(L))
	    station_cluster.append(np.median(S))
	#ax.plot_trisurf(x, y, z, cmap='viridis')
	#ax.scatter(x, y, z, c=z, cmap='viridis')
	#plt.show()

	file = 'C:/Users/FieldUser/Desktop/Tools/clusteringTest/GPS_nodes.csv'
	fOut = open(file, 'w')

	for i in range(len(x)):
		fOut.write(str(line_cluster[i]) + ',' + str(station_cluster[i]) + ',' + str(x[i]) + ',' + str(y[i]) + ',' + str(z[i]) + '\n')
	fOut.close()	



	gpx = gpxpy.gpx.GPX()

	for i in range(len(x)):
	    MEM[i] = [MEM[i][j] for j in range(len(MEM[i])) if MEM[i][j] >= 0 ]
	    coord = utm.to_latlon(x[i], y[i], zone[0], zone[1])
	    name = ''
	    ID[i] = ID[i].tolist()
	    ID2 = np.unique(ID[i]).tolist()
	    for j in range(len(ID2)):
	        name += ID2[j] + '/'
	    description_gpx = ''
	    ind = np.argsort(MEM[i]).tolist()
	    if len(MEM[i]) > 0:
	        description_gpx = ID[i][0] + ':' + str(MEM[i][0]) + '/'
	        for j in range(len(MEM[i]) - 1):
	            if MEM[i][ind[j + 1]] > 0:
	                #print(ID[i][ind[j + 1]], str(MEM[i][ind[j + 1]]))	
	                description_gpx += ID[i][ind[j + 1]] + ':' + str(MEM[i][ind[j + 1]]) + '/'
	    try:
	        gpx.waypoints.append(gpxpy.gpx.GPXWaypoint(coord[0],\
	                                                                                          coord[1],\
	                                                                                          elevation=z[i],
	                                                                                          name=str(int(line_cluster[i])) + ' ' + str(int(station_cluster[i])) + ' ' + name[:-1],
	                                                                                          description=description_gpx[:-1]))
	    except:
	        pass
		
	return gpx
	#fOut = open('nodes.gpx', 'w')
	#fOut.write(gpx.to_xml())
	#fOut.close()

	