% read data of an .h5 file
% replace '/group_name/dataset_name' with the corresponding group and
% dataset name you would like to read
data1 = h5read('h_mec_1000.h5', '/Matrix/pt');
data2 = h5read('h_mec_1000.h5', '/Vector/maxvxsmod');
data3 = h5read('h_mec_1000.h5', '/Value/values');
% data4 = h5read('h_mec_1000.h5', '/Antiplane/VZ0');


