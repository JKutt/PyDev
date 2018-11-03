# Main code for friday
import multiprocessing
import sys
from classes import *
import stages_par as stages
import time

if __name__ == '__main__':


    print("Hello, this is Friday.")
    # Input
    project_path = sys.argv[1]
    print("Analyzing project: " + project_path)
    log_file = 'Recordings_36.txt'
    print("Loaded recording file: " + log_file)
    project = Project(path=project_path,
                      data_path=project_path + 'DATA\\',
                      log_file=[],
                      list_nodes=[],
                      injections=[],
                      records=[],
                      recording_file=project_path + 'LogFiles\\' + log_file)

    print("Seems like what you told me exists. Starting analysis.")
    # project_path + '/LogFiles/' + log_file
    # Organisation of project folders:
    # ZIP---|- DATA folder with node and current recorders data
    #       |- LOGS folder with nodeMessages/recordingFiles (nodeMessages.txt/Recordings_0.txt)
    #       |- DDN file

    # Getting content of data folder
    project.list_nodes = stages.get_content_folder(project)

    # First pass on checking data structure
    s_time = time.time()
    print("Folder and file analysis")
    project.records, project.injections = stages.get_records(project)
    print("\tElapsed time for library reading: " + str(time.time() - s_time) + ' seconds.')
    project.get_node_list()
    project.initialize_node_data()
    for ind in range(len(project.injections)):
        if 'Good' in project.injections[ind].valid:
            print("Analysing injection " + str(project.injections[ind].num) + "/" + str(project.injections[-1].num))

            print('Starting harmonic analysis of DATA folder')
            s_time = time.time()
            project.node_data[ind] = stages.harmonic_vp_analysis_injection(project, project.injections[ind])
            print("\tElapsed time for harmonic and Vp analysis: " + str(time.time() - s_time) + ' seconds.')

            print('Coherency analysis')
            s_time = time.time()
            project.node_data[ind] = stages.coherency_analysis(project, project.injections[ind], project.node_data[ind])
            print("\tElapsed time for coherency analysis: " + str(time.time() - s_time) + ' seconds.')

            print('Statistical analysis')
            project.node_data[ind] = stages.statistical_analysis(project, project.node_data[ind])

            print('Checking outliers on Vp and harmonics')
            project.node_data[ind] = stages.vp_harmonics_check(project, project.node_data[ind])

            print('Quarantining data')
            stages.quarantine(project, ind)
            # Potential room for DIASPro call to process the injection
            #exit()
# Second pass on

# Third pass on time series
# Using multiprocessing module to split and analyse injections separately

# Quarantining data
