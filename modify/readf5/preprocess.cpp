#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <set>
#include <iostream>
#include <fstream>
#include <sys/types.h>  
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "H5Cpp.h"


using namespace std;
using namespace H5;
#include "preprocess.h"
const int OUTPUT_INFO = 1;
const long long iinfo_int64 = 9223372036854775807;


extern "C"
int get_basecall(char* mfile_path, char* m_event_basecall, float* my_raw_signals) {

	// string mfile_path = "read116.fast5";
	int signal_len, events_len;
	H5File file = H5File(mfile_path, H5F_ACC_RDONLY);
    Group reads = file.openGroup("/Raw/Reads");
    H5std_string reads_idx = reads.getObjnameByIdx(0);
    string signal_path = "/Raw/Reads/" + reads_idx + "/Signal";
    DataSet signal = file.openDataSet(signal_path);
    signal_len = signal.getSpace().getSimpleExtentNpoints();

    DataSet events = file.openDataSet("/Analyses/Basecall_1D_000/BaseCalled_template/Events");
    events_len = events.getSpace().getSimpleExtentNpoints();
	// cout<<"signal_len: "<<signal_len<<" events_len:"<<events_len<<endl;
	short* raw_signals = new short[signal_len];
	// float* my_raw_signals = new float[signal_len];
	unsigned long long start_time = 0;
	int used_albacore_version = 2;
	H5std_string channel_number("");
	double digitisation = 0.0;
	double offset = 0.0;
	double range = 0.0;
	double sampling_rate = 0.0;
	H5std_string fastq_buf("");
	Events_t * events_data = new Events_t[events_len];
	readfast5(mfile_path, raw_signals, start_time, used_albacore_version, channel_number, digitisation,
	offset, range, sampling_rate, fastq_buf, events_data);
	vector<string> fastq_split;
	String_Split(fastq_buf, '\n', fastq_split);
	string fq_seq = fastq_split[1];
	string read_id = String_Replace(String_Replace(fastq_split[0], " ", ":::"), "\t", "|||");
	
	int moptions_outLevel = 2;
	string SignalGroup = "simple";
	string f5status = "";

	get_event(my_raw_signals, m_event_basecall, read_id, f5status, mfile_path, SignalGroup, moptions_outLevel, used_albacore_version, fq_seq, raw_signals, &sampling_rate, &start_time,  signal_len, events_len, events_data);

	fstream fout;
	fout.open("DataPreprocess.log", ios::out | ios::app);
	for (int j=0;j<signal_len;j++){
		fout<<my_raw_signals[j]<<" ";
	} 
	fout << endl;
	fout.close();

	delete []raw_signals;
	delete []events_data;
	
	return 0;
}

extern "C"
void get_f5len(char* f5name, int& signal_len, int& events_len){
    H5File file = H5File(f5name, H5F_ACC_RDONLY);
    Group reads = file.openGroup("/Raw/Reads");
    H5std_string reads_idx = reads.getObjnameByIdx(0);
    string signal_path = "/Raw/Reads/" + reads_idx + "/Signal";
    DataSet signal = file.openDataSet(signal_path);
    signal_len = signal.getSpace().getSimpleExtentNpoints();
    DataSet events = file.openDataSet("/Analyses/Basecall_1D_000/BaseCalled_template/Events");
    events_len = events.getSpace().getSimpleExtentNpoints();
	cout<<"signal_len: "<<signal_len<<" events_len:"<<events_len<<endl;
}
extern "C"
void get_AlbacoreVersion(char* f5name, int& used_albacore_version){
	H5File file = H5File(f5name, H5F_ACC_RDONLY);
	Group basecall = file.openGroup("/Analyses/Basecall_1D_000");
    Attribute albacorev = basecall.openAttribute("version");
    DataType albacore_type = albacorev.getDataType();
    H5std_string albacorev_buf("");
    albacorev.read(albacore_type, albacorev_buf);
    // cout<<"Albacore version= "<<albacorev_buf<<endl;
    if (albacorev_buf[0] == '2'){
        used_albacore_version = 2;
    }
    else {
        used_albacore_version = 1;
    }
}


extern "C"
int get_read_id(char* f5name, char* read_id_buf){
	H5File file = H5File(f5name, H5F_ACC_RDONLY);
	DataSet fastq = file.openDataSet("/Analyses/Basecall_1D_000/BaseCalled_template/Fastq");
    DataType fastq_type = fastq.getDataType();
	H5std_string fastq_buf("");
    fastq.read(fastq_buf, fastq_type);
	vector<string> fastq_split;
	String_Split(fastq_buf, '\n', fastq_split);
	string read_id_str = String_Replace(String_Replace(fastq_split[0], " ", ":::"), "\t", "|||");
	cout<<read_id_str.length()<<endl;
	const char* read_id = read_id_str.data();
	strcpy(read_id_buf, read_id);
	cout<<read_id_buf<<endl;
	return 0;
}

void read_f5len(string f5name, int& signal_len, int& events_len){
    H5File file = H5File(f5name, H5F_ACC_RDONLY);
    Group reads = file.openGroup("/Raw/Reads");
    H5std_string reads_idx = reads.getObjnameByIdx(0);
    string signal_path = "/Raw/Reads/" + reads_idx + "/Signal";
    DataSet signal = file.openDataSet(signal_path);
    signal_len = signal.getSpace().getSimpleExtentNpoints();

    DataSet events = file.openDataSet("/Analyses/Basecall_1D_000/BaseCalled_template/Events");
    events_len = events.getSpace().getSimpleExtentNpoints();
}
void readfast5(string f5name, short* signal_buf, unsigned long long& start_time, int& used_albacore_version, string& channel_number_buf, double& digitisation_buf,
double& offset_buf, double& range_buf, double& sampling_rate_buf, string& fastq_buf, Events_t *events_data){
    H5File file = H5File(f5name, H5F_ACC_RDONLY);
    /* read Signal */
    Group reads = file.openGroup("/Raw/Reads");
    H5std_string reads_idx = reads.getObjnameByIdx(0);
    string signal_path = "/Raw/Reads/" + reads_idx + "/Signal";
    DataSet signal = file.openDataSet(signal_path);
    DataType signal_type = signal.getDataType();
    signal.read(signal_buf, signal_type);
    
    /*read start time*/
    string start_time_path = "/Raw/Reads/" + reads_idx;
    Group start_time_group = file.openGroup(start_time_path);
    Attribute start_time_buf = start_time_group.openAttribute("start_time");
    DataType start_time_type = start_time_buf.getDataType();
    start_time_buf.read(start_time_type, &start_time);
	// cout<<"start_time: "<< start_time<<endl;
    
    /* read Albacore version */
    Group basecall = file.openGroup("/Analyses/Basecall_1D_000");
    Attribute albacorev = basecall.openAttribute("version");
    DataType albacore_type = albacorev.getDataType();
    H5std_string albacorev_buf("");
    albacorev.read(albacore_type, albacorev_buf);
    // cout<<"Albacore version= "<<albacorev_buf<<endl;
    if (albacorev_buf[0] == '2'){
        used_albacore_version = 2;
    }
    else {
        used_albacore_version = 1;
    }
   
    /*read channel_id attribute*/ 
    Group group1 = file.openGroup("UniqueGlobalKey");
    Group group2 = group1.openGroup("channel_id");
    /*channel_number*/
    Attribute channel_number = group2.openAttribute("channel_number");
    DataType channel_number_type = channel_number.getDataType();
    channel_number.read(channel_number_type, channel_number_buf);
    // cout << "channel_number: " << channel_number_buf << endl;
    
    /*digitisation*/
    Attribute digitisation = group2.openAttribute("digitisation");
    DataType digitisation_type = digitisation.getDataType();
    digitisation.read(digitisation_type, &digitisation_buf);
    // cout << "digitisation: " << digitisation_buf << endl;
    
    /*offset*/
    Attribute offset = group2.openAttribute("offset");
    DataType offset_type = offset.getDataType();
    offset.read(offset_type, &offset_buf);
    // cout << "offset: " << offset_buf << endl;
    /*range*/
    Attribute range = group2.openAttribute("range");
    DataType range_type = range.getDataType();
    range.read(range_type, &range_buf);
    // cout << "range: " << range_buf << endl;
    
    /*sampling_rate*/
    Attribute sampling_rate = group2.openAttribute("sampling_rate");
    DataType sampling_rate_type = sampling_rate.getDataType();
    sampling_rate.read(sampling_rate_type, &sampling_rate_buf);
    // cout << "sampling_rate: " << sampling_rate_buf << endl;

    /*read fastq*/
    DataSet fastq = file.openDataSet("/Analyses/Basecall_1D_000/BaseCalled_template/Fastq");
    DataType fastq_type = fastq.getDataType();
    fastq.read(fastq_buf, fastq_type);
    // cout<<"fastq:"<<fastq_buf<<endl;
    
    /*read events*/
    DataSet events = file.openDataSet("/Analyses/Basecall_1D_000/BaseCalled_template/Events");
    CompType ctype(sizeof(Events_t));
    hid_t tid_s = H5Tcopy(H5T_C_S1);
    H5Tset_size(tid_s, 6);

    
    ctype.insertMember("mean", HOFFSET(Events_t, mean), PredType::NATIVE_FLOAT);
    ctype.insertMember("start", HOFFSET(Events_t, start), PredType::NATIVE_ULONG);
    ctype.insertMember("stdv", HOFFSET(Events_t, stdv), PredType::NATIVE_FLOAT);
    ctype.insertMember("length", HOFFSET(Events_t, length), PredType::NATIVE_ULONG);
    ctype.insertMember("model_state", HOFFSET(Events_t, model_state), tid_s);
    ctype.insertMember("move", HOFFSET(Events_t, move), PredType::NATIVE_INT);
    ctype.insertMember("p_model_state", HOFFSET(Events_t, p_model_state), PredType::NATIVE_FLOAT);
    ctype.insertMember("weights", HOFFSET(Events_t, weights), PredType::NATIVE_FLOAT);
    
    events.read(events_data, ctype);
    
}

void String_Split(string s, char delim, vector<string>& ans)
{
	string::size_type pos_1, pos_2 = 0;
	while (pos_2 != s.npos) {
		pos_1 = s.find_first_not_of(delim, pos_2);
		if (pos_1 == s.npos) break;
		pos_2 = s.find_first_of(delim, pos_1);
		ans.push_back(s.substr(pos_1, pos_2 - pos_1));
	}
}

string String_Replace(string s, const char* oldstr, const char* newstr) {
	while (s.find(oldstr) != string::npos) {
		s = s.replace(s.find(oldstr), strlen(oldstr), newstr);
	}
	return s;
}

void showAllFiles(const char* dir_name, vector<string>& filename_vec)
{
	// check the parameter !
	if (NULL == dir_name)
	{
		cout << " dir_name is null ! " << endl;
		return;
	}

	// check if dir_name is a valid dir
	struct stat s;
	lstat(dir_name, &s);
	if (!S_ISDIR(s.st_mode))
	{
		cout << "dir_name is not a valid directory !" << endl;
		return;
	}

	struct dirent* filename;    // return value for readdir()
	DIR* dir;                   // return value for opendir()
	dir = opendir(dir_name);
	if (NULL == dir)
	{
		cout << "Can not open dir " << dir_name << endl;
		return;
	}
	cout << "Successfully opened the dir !" << endl;

	/* read all the files in the dir ~ */
	while ((filename = readdir(dir)) != NULL)
	{
		// get rid of "." and ".."
		if (strcmp(filename->d_name, ".") == 0 ||
			strcmp(filename->d_name, "..") == 0)
			continue;
		// cout<<filename ->d_name <<endl;
		filename_vec.push_back(filename->d_name);
	}
}

char* Upper(char* str, int length) {
	/* 该函数用于将字符串中所有字符变为大写 */
	int i;
	for (i = 0; i < length; i++) {
		if (str[i] >= 'a' && str[i] <= 'z') {
			str[i] -= 32;
		}
	}
	return str;
}

float myround(float num, int precision) {
	/* 该函数用于四舍五入精确到小数 */
	if (num > 0) {
		num = (float)(int)(num * pow(10, precision) + 0.5);
		num = num / pow(10, precision);
	}
	else {
		num = (float)(int)(num * pow(10, precision) - 0.5);
		num = num / pow(10, precision);
	}
	return num;
}

void get_event(float* my_raw_signals, char* m_event_basecall, string read_id, string f5status, string mfile_path, string SignalGroup, int moptions_outLevel, int used_albacore_version,const string& fq_seq, short* signals, double* sampling_rate, unsigned long long* start_time, int signal_len, int events_len, Events_t* events_data) {
	/* get events from a fast5 file */
	clock_t start, end;
	double  duration;
	start = clock();

	bool convertError = false;
	int i;
	string m_event_basecall_str = "";
	int left_right_skip_left = 0;
	int left_right_skip_right = 0;

	vector<M_event> m_event;
	if (used_albacore_version == 1) {
		int move0_left = 0, move0_right = events_len - 1;
		while (move0_left < move0_right) {
			if (events_data[move0_left].move == 0) {
				move0_left += 1;
			}
			else {
				break;
			}
		}
		if (move0_left > move0_right - 20) {
			fprintf(stderr, "Too many move0 at 3'(l%d, r%d)\n", move0_left, move0_right);
			exit(EXIT_FAILURE);
		}
		while (move0_left > move0_right) {
			/* get the last non-stay event at the right tail */
			if (events_data[move0_right].move == 0) {
				move0_right -= 1;
			}
			else {
				break;
			}
		}
		if (move0_right < move0_left + 20) {
			fprintf(stderr, "Too many move0 at 5'(l%d, r%d)\n", move0_left, move0_right);
			exit(EXIT_FAILURE);
		}

		// get the start time
		float based_ind = events_data[move0_left].start * *sampling_rate - *start_time;
		unsigned long long first_base_index_in_raw_signal = (unsigned long long)round(based_ind);

		// get the potential error of the starting time
		if (first_base_index_in_raw_signal < -2) {
			fprintf(stderr, "The index of the first base is less than -2(%llu=%llu*%.0f-%llu)\n", first_base_index_in_raw_signal, events_data[move0_left].start, *sampling_rate, *start_time);
			exit(EXIT_FAILURE);
		}
		else if (first_base_index_in_raw_signal < 0) {
			first_base_index_in_raw_signal = 0;
			if (moptions_outLevel <= OUTPUT_INFO) {
				printf("Warning!!! first_base_index_in_raw_signal less than 0 %s\n", mfile_path.c_str());
			}
		}

		int pre_i = move0_left;
		unsigned long long cur_length = events_data[pre_i].length * (int)*sampling_rate;
		double cal_st;
	
		for (i = move0_left + 1; i < move0_left + 1; i++) {
			if (events_data[i].move > 0) {
				// for non-stay event
				if (pre_i == move0_left) {
					M_event m = { myround(events_data[pre_i].mean, 3), myround(events_data[pre_i].stdv, 3), \
					first_base_index_in_raw_signal, \
					cur_length, };
					strcpy(m.model_state, Upper(events_data[pre_i].model_state, events_data[pre_i].length));
					m_event.push_back(m);
				}
				else {
					// calculate starting index in raw signal
					// calculated position
					cal_st = (events_data[pre_i].start - events_data[move0_left].start) * (double)*sampling_rate + based_ind;
					if (cal_st < 0) {
						printf("Warning less than 0\n");
					}
					if (cal_st > 0 && cal_st - (m_event.back().start + m_event.back().length) > 0\
						&& (unsigned long long)(cal_st - (m_event.back().start + m_event.back().length)) > 0) {
						if ((unsigned long long)(cal_st - (m_event.back().start + m_event.back().length)) > 2) {

							M_event m = { myround(events_data[pre_i].mean, 3), myround(events_data[pre_i].stdv, 3), \
							m_event.back().start + m_event.back().length, \
							(unsigned long long)(cal_st - (m_event.back().start + m_event.back().length)),};
							strcpy(m.model_state, Upper(events_data[pre_i].model_state, events_data[pre_i].length));
							m_event.push_back(m);

							m = { myround(events_data[pre_i].mean, 3), myround(events_data[pre_i].stdv, 3), \
							(unsigned long long)cal_st, \
							cur_length, };
							strcpy(m.model_state, Upper(events_data[pre_i].model_state, events_data[pre_i].length));
							m_event.push_back(m);
						}
						else {
							// for a normal event
							M_event m = { myround(events_data[pre_i].mean, 3), myround(events_data[pre_i].stdv, 3), \
							m_event.back().start + m_event.back().length, \
							(unsigned long long)(cal_st - (m_event.back().start + m_event.back().length)), };
							strcpy(m.model_state, Upper(events_data[pre_i].model_state, events_data[pre_i].length));
							m_event.push_back(m);
						}
					}
					else {
						M_event m = { myround(events_data[pre_i].mean, 3), myround(events_data[pre_i].stdv, 3), \
						m_event.back().start + m_event.back().length, \
						cur_length,  };
						strcpy(m.model_state, Upper(events_data[pre_i].model_state, events_data[pre_i].length));
						m_event.push_back(m);
					}
					if (m_event.back().start > iinfo_int64 - 2 || m_event.back().start < 0) {
						if (convertError == false) {
							printf("ex: %llu*%.0f=%llu\n%llu\n%s\n%llu\n%llu\n", events_data[move0_left].start, *sampling_rate, events_data[move0_left].start * (int)*sampling_rate, *start_time, mfile_path.c_str(), m_event.back().start, m_event.back().length);
						}
						convertError = true;
					}
				}
				pre_i = i;
				cur_length = (unsigned long long)(events_data[i].length * (double)*sampling_rate);
			}
			else {
				// for stay event
				cur_length += (unsigned long long)(events_data[i].length * (double)*sampling_rate);
			}
		}
		// for the last event
		cal_st = (events_data[pre_i].start - events_data[move0_left].start) * (double)*sampling_rate + based_ind;
		if (cal_st < 0) {
			printf("Warning less than 0\n");
		}
		if (cal_st > 0 && cal_st - (m_event.back().start + m_event.back().length) > 0\
			&& (unsigned long long)(cal_st - (m_event.back().start + m_event.back().length)) > 0) {
			if ((unsigned long long)(cal_st - (m_event.back().start + m_event.back().length)) > 2) {

				M_event m = { myround(events_data[pre_i].mean, 3), myround(events_data[pre_i].stdv, 3), \
				m_event.back().start + m_event.back().length, \
				(unsigned long long)(cal_st - (m_event.back().start + m_event.back().length)), };
				strcpy(m.model_state, Upper(events_data[pre_i].model_state, events_data[pre_i].length));
				m_event.push_back(m);

				m = { myround(events_data[pre_i].mean, 3), myround(events_data[pre_i].stdv, 3), \
				(unsigned long long)cal_st, \
				cur_length,};
				strcpy(m.model_state, Upper(events_data[pre_i].model_state, events_data[pre_i].length));
				m_event.push_back(m);
			}
			else {
				// for a normal event
				M_event m = { myround(events_data[pre_i].mean, 3), myround(events_data[pre_i].stdv, 3), \
				m_event.back().start + m_event.back().length, \
				(unsigned long long)(cal_st - (m_event.back().start + m_event.back().length)), };
				strcpy(m.model_state, Upper(events_data[pre_i].model_state, events_data[pre_i].length));
				m_event.push_back(m);
			}
		}
		else {
			M_event m = { myround(events_data[pre_i].mean, 3), myround(events_data[pre_i].stdv, 3), \
			m_event.back().start + m_event.back().length, \
			cur_length,};
			strcpy(m.model_state, Upper(events_data[pre_i].model_state, events_data[pre_i].length));
			m_event.push_back(m);
		}

		/* decode */

		// sp_param['m_event_basecall']
		// string m_event_basecall = "";
		for (i = 0; i < m_event.size(); i++) {
			m_event_basecall_str += string(1, m_event[i].model_state[2]);
		}
		const char* m_event_basecall_data = m_event_basecall_str.data();
		strcpy(m_event_basecall, m_event_basecall_data);

		// sp_param['left_right_skip']
		left_right_skip_left = move0_left;
		left_right_skip_right = events_len - move0_right - 1;
	}
	else if (used_albacore_version == 2) {
	
		if (SignalGroup == "simple") {
			
			int pre_i = 0;
			unsigned long long pre_length = events_data[pre_i].length;
			for (i = 1; i < events_len; i++) {
				//777
				if (events_data[i].move > 0) {
					M_event m = {
						myround(events_data[pre_i].mean, 3), myround(events_data[pre_i].stdv, 3), \
						events_data[pre_i].start, pre_length,};
					strcpy(m.model_state, events_data[pre_i].model_state);
					m_event.push_back(m);
					
					pre_i = i;
					pre_length = events_data[pre_i].length;
				}
				else {
					// for stay events
					pre_length += events_data[i].length;
				}
			}
			M_event m = {
				myround(events_data[pre_i].mean, 3), myround(events_data[pre_i].stdv, 3), \
				events_data[pre_i].start, pre_length,};
			strcpy(m.model_state, events_data[pre_i].model_state);
			m_event.push_back(m);
			

		}
		else {
			m_event = get_EventInfo(signal_len, events_len, signals, events_data, fq_seq);			
		}
		//sp_param['m_event'] = m_event

		// sp_param['m_event_basecall']
		for (i = 0; i < m_event.size(); i++) {
			m_event_basecall_str += string(1, m_event[i].model_state[2]);
		}
		const char* m_event_basecall_data = m_event_basecall_str.data();
		strcpy(m_event_basecall, m_event_basecall_data);

		// sp_param['left_right_skip']
		left_right_skip_left = 0;
		left_right_skip_right = 0;
	}
	else {
		fprintf(stderr, "This version of Albacore is not supported. Please use the version of Albacore 1.x or 2.x");
		exit(EXIT_FAILURE);
	}

	// my_raw_signal = mnormalized(my_raw_signals, mfile_path, signal_len, signals, events_len, events_data);
	mnormalized(my_raw_signals, mfile_path, signal_len, signals, events_len, events_data);

	for (int i = 0; i < m_event.size(); i++) {
		vector<float> Signal_range(my_raw_signals + m_event[i].start, my_raw_signals + m_event[i].start + m_event[i].length);
		if (Signal_range.size() == 0) {
			cout << "Signal out of range " << i << ":"
				<< events_data[i].start << "-" << events_data[i].length << ";" << events_len << "for" << mfile_path << endl;
			if (i > 500)
			{
				m_event.assign(m_event.begin(), m_event.begin() + i - 1);
			}
			else f5status = "Less event";
			break;
		}
		float cur_mean, cur_stdv;
		cal_mean_stdv(Signal_range, cur_mean, cur_stdv);
		m_event[i].mean = myround(cur_mean, 3);
		m_event[i].stdv = myround(cur_stdv, 3);
	}

	// end = clock();
	// duration = (double)(end - start) / CLOCKS_PER_SEC;
	// fstream fout;
	// fout.open("DataPreprocess.log", ios::out | ios::app);
	// for (int j=0;j<signal_len;j++){
	// 	fout<<my_raw_signals[j]<<" ";
	// } 
	// fout << endl;
	// fout.close();
	// printf("%f seconds\n", duration);

}

vector<M_event> get_EventInfo(int signal_len, int lines_event, short* signals, Events_t* events_data ,const string& fq_seq) {
	int i;
	int min_signal_num = 4;
	int mysignals_len = signal_len + 1;
	float* mysignals = new float[mysignals_len];
	mysignals[0] = 0.0;
	vector<M_event> event_info;
	for (i = 1; i < mysignals_len + 1; i++) {
		mysignals[i] = myround((float)(signals[(size_t)i - 1]) / 50.0, 5);
	}
	float* signal_sum = cumsum(mysignals_len + 1, mysignals);

	int m_signal_dif_len = mysignals_len - 2 * min_signal_num;
	float* m_signal_dif = new float[m_signal_dif_len];
	for (i = min_signal_num; i < mysignals_len - min_signal_num; i++) {
		m_signal_dif[i - min_signal_num] = signal_sum[i] * 2;
	}
	for (i = 0; i < m_signal_dif_len; i++) {
		m_signal_dif[i] -= signals[i];
	}
	for (i = 2 * min_signal_num; i < mysignals_len; i++) {
		m_signal_dif[i - 2 * min_signal_num] -= signal_sum[i];
	}
	for (i = 0; i < m_signal_dif_len; i++) {
		m_signal_dif[i] = fabsf(m_signal_dif[i]);
	}
	unsigned long long last_signal_i = events_data[0].start;
	int fq_seq_i = 2;
	int c_move_num = 1;
	// incrrt_event_list = []
	vector<int> incrrt_event_list;
	for (i = 1; i < lines_event; i++) {
		if (events_data[i].move != 0) {
			c_move_num += events_data[i].move;
			vector<int> split_points = get_extreme_N(m_signal_dif_len, m_signal_dif, c_move_num - 1, last_signal_i, events_data[i].start + events_data[i].length, min_signal_num);
			int c_m_i;
			for (c_m_i = 0; c_m_i < c_move_num - 1; c_m_i++) {
				int h_m_i = 0;
				unsigned long long c_e_p;
				if (c_m_i < split_points.size()) {
					h_m_i = c_m_i;
					c_e_p = split_points[h_m_i];
				}
				else {
					h_m_i = split_points.size() - 1;
					c_e_p = last_signal_i + min_signal_num;
					incrrt_event_list.push_back(event_info.size());
				}
				float c_mnn = mean(signal_len, signals, last_signal_i, c_e_p);
				float c_std = stdvariancec(signal_len, signals, last_signal_i, c_e_p, c_mnn);
				unsigned long long c_start = last_signal_i;
				unsigned long long c_length = c_e_p - last_signal_i;
				string c_mode(&fq_seq[(size_t)fq_seq_i-2], &fq_seq[(size_t)fq_seq_i+3]);
				M_event m = { c_mnn, c_std, c_start, c_length };
				strcpy(m.model_state, c_mode.c_str());
				event_info.push_back(m);

				last_signal_i = split_points[h_m_i];
				fq_seq_i += 1;
			}
			int c_move_num = 1;
			int ev_i = lines_event - 1;
			unsigned long long c_e_p = events_data[ev_i].start + events_data[ev_i].length;
			float c_mnn = mean(signal_len, signals, last_signal_i, c_e_p);
			float c_std = stdvariancec(signal_len, signals, last_signal_i, c_e_p, c_mnn);
			unsigned long long c_start = last_signal_i;
			unsigned long long c_length = c_e_p - last_signal_i;
			string c_mode(&fq_seq[(unsigned long long)fq_seq_i - 2], &fq_seq[fq_seq_i]);
			M_event m = { c_mnn, c_std, c_start, c_length };
			strcpy(m.model_state, c_mode.c_str());
			event_info.push_back(m);

			vector<int>::iterator pi;
			for (pi = incrrt_event_list.begin(); pi != incrrt_event_list.end(); pi++) {
				unsigned long long h_2 = (unsigned long long)((event_info[(size_t)*pi + 1].length + event_info[(size_t)*pi + 1].start - event_info[*pi].start) / 2.0 + 0.2);
				event_info[*pi].length = h_2;
				event_info[(size_t)*pi + 1].start = event_info[*pi].start + event_info[*pi].length;
				event_info[(size_t)*pi + 1].length = event_info[(size_t)*pi + 1].length - h_2;
			}
		}
	}

	delete[] mysignals;
	delete[] signal_sum;
	delete[] m_signal_dif;
	return event_info;
}

float* cumsum(int mysignals_len, float* mysignals) {
	float* cumsumsignals = new float[mysignals_len];
	int i, j;
	for (i = 0; i < mysignals_len; i++) {
		cumsumsignals[i] = 0.0;
		for (j = 0; j <= i; j++) {
			cumsumsignals[i] += mysignals[j];
		}
	}
	return cumsumsignals;
}

vector<int> get_extreme_N(int m_signal_dif_len, float* m_signal_dif, int n_splits, unsigned long long p_signal_start, unsigned long long p_signal_end, int min_signal_num) {
	int i;
	int cu_region_sort_pos_start = (int)(p_signal_start - min_signal_num + 0.5);
	int cu_region_sort_pos_end = (int)(p_signal_end - min_signal_num + 0.5);
	vector<float> m_signal_dif_v(m_signal_dif + cu_region_sort_pos_start, m_signal_dif + m_signal_dif_len - cu_region_sort_pos_end);
	vector<int> cu_region_sort_pos = reargsort(m_signal_dif_v);
	vector<int>::iterator pr;
	for (pr = cu_region_sort_pos.begin(); pr != cu_region_sort_pos.end(); pr++) {
		*pr += p_signal_start;
	}
	set<int> m_nb_pos;
	int m_nb_pos_end = (int)(p_signal_start + min_signal_num - 0.5);
	for (i = p_signal_start; i < m_nb_pos_end; i++) {
		m_nb_pos.insert(i);
	}
	int m_nb_pos_begin = (int)(p_signal_end - min_signal_num + 1.5);
	for (i = m_nb_pos_end; i < p_signal_end; i++) {
		m_nb_pos.insert(i);
	}
	vector<int> split_points_list;
	for (i = 0; i < cu_region_sort_pos.size(); i++) {
		if (m_nb_pos.count(i) == 0) {
			split_points_list.push_back(i);
			if (split_points_list.size() == n_splits) {
				break;
			}
			int c_pos_begin = i - min_signal_num + 1;
			int c_pos_end = i + min_signal_num + 1;
			int j;
			for (j = c_pos_begin; j < c_pos_end; j++) {
				m_nb_pos.insert(j);
			}
		}
	}
	sort(split_points_list.begin(), split_points_list.end());
	return split_points_list;
}

template <typename T>
vector<int> reargsort(const vector<T>& v) {
	/* for returning the index of reverse sorted array */
	// construct index array
	vector<int> idx(v.size());
	iota(idx.begin(), idx.end(), 0);
	// reverse sort
	sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] > v[i2];});
	return idx;
}

float mean(int signal_len,const short* raw_signals,int begin, int end) {
	/* for calculating the mean of raw_signals */
	float sum = 0.0;
	int i;
	for (i = begin; i < end; i++) {
		sum += raw_signals[i];
	}
	float mean = sum / signal_len;
	return mean;
}

float variance(int signal_len, const short* raw_signals, int begin, int end, float mean) {
	/* for calculating the variance of raw_signals */
	float variance = 0.0;
	for (int i = begin; i < end; i++) {
		variance += pow(raw_signals[i] - mean, 2);
	}
	return variance;
}

float stdvariancec(int signal_len, const short* raw_signals, int begin, int end, float mean) {
	/* for calculating the standard variance of raw_signals */
	float stdvariance = variance(signal_len, raw_signals, begin, end, mean);
	return sqrt(stdvariance);
}

void mnormalized(float* my_raw_signals, string mfile_path, int signal_len, short* raw_signals, int events_len, Events_t *events_data) {
	/* for normalizing raw_signals */
	if (events_data[0].start >= events_data[events_len-1].start + events_data[events_len-1].length) {
		cout << "Fatal error signal start position is less than the end position.\n" << " The path is " << mfile_path
			<< "\nstart[0] is " << events_data[0].start << "\nstarts[-1] is " << events_data[events_len-1].start << "\nlengths[-1] is " << events_data[events_len-1].length << endl;
	}

	vector<short> vsshift(raw_signals + events_data[0].start, raw_signals + events_data[events_len-1].start + events_data[events_len-1].length);

	float mshift = median(vsshift);
	// cout << mshift << endl;
	vsshift.erase(vsshift.begin(), vsshift.end());
	vector<short> vsscale(raw_signals + events_data[0].start, raw_signals + events_data[events_len-1].start + events_data[events_len-1].length);
	vector<float> vfscale;
	vector<short>::iterator ps;
	for (ps = vsscale.begin(); ps != vsscale.end(); ps++) {
		vfscale.push_back(fabsf(*ps-mshift));
	}
	float mscale = median(vfscale);
	vsscale.erase(vsscale.begin(), vsscale.end());
	vfscale.erase(vfscale.begin(), vfscale.end());
	// standardize
	// float* my_raw_signals = new float[signal_len];
	for (int i = 0; i < signal_len; i++) {
		my_raw_signals[i] = (raw_signals[i] - mshift) / mscale;
	}
	// get meand
	vector<float> vfmed(my_raw_signals + events_data[0].start, my_raw_signals + events_data[events_len-1].start + events_data[events_len-1].length);
	float read_med = median(vfmed);
	vfmed.erase(vfmed.begin(), vfmed.end());
	vector<float> vfmad(my_raw_signals + events_data[0].start, my_raw_signals + events_data[events_len-1].start + events_data[events_len-1].length);
	vector<float>::iterator pf;
	for (pf = vfmad.begin(); pf != vfmad.end(); pf++) {
		*pf = fabsf(*pf - read_med);
	}
	float read_mad = median(vfmad);
	vfmad.erase(vfmad.begin(), vfmad.end());
	float lower_lim = read_med - (read_mad * 5);
	float upper_lim = read_med + (read_mad * 5);
	float round_upper_lim = myround(upper_lim, 3);
	float round_lower_lim = myround(lower_lim, 3);
	//cout << upper_lim << " " << lower_lim << endl;
	for (int i = 0; i < signal_len; i++) {
		if (my_raw_signals[i] > upper_lim) {
			my_raw_signals[i] = round_upper_lim;
		}
		else {
			if (my_raw_signals[i] < lower_lim) {
				my_raw_signals[i] = round_lower_lim;
			}
			else {
				my_raw_signals[i] = myround(my_raw_signals[i], 3);
			}
		}
		// cout << my_raw_signals[i] << " ";
	}
//	cout << endl;
	// return my_raw_signals;
}

template <typename T>
float median(vector<T>& vt){	
	/* for getting the median */
	sort(vt.begin(), vt.end());
	float median;
	size_t n = vt.size();
	if (n % 2 == 0) {
		median = (vt[n / 2 - 1] + vt[n / 2]) / 2.0;
	}
	else {
		median = (float)vt[n / 2];
	}
	return median;
}
//lx: calculate mean and stdv
template <typename T>
void cal_mean_stdv(vector<T>& vt, float& mean, float& stdv){	
	float sum=accumulate(vt.begin(),vt.end(),0.0);
	mean = sum/vt.size();
	float accum = 0.0;
	for_each (vt.begin(),vt.end(), [&](const float d){
		accum  += (d-mean)*(d-mean);
	});
	stdv = sqrt(accum/(vt.size()));
}

