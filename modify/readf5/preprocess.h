#ifndef PROCESS_H
#define PROCESS_H


struct Events_t {
	float mean;
	unsigned long long start;
	float stdv;
	unsigned long long length;
	char model_state[6];
	int move;
	float p_model_state;
	float weights;
};
struct M_event {
	float mean;
	float stdv;
	unsigned long long start;
	unsigned long long length;
	char model_state[6];
};


void read_f5len(string f5name, int& signal_len, int& events_len);

void readfast5(string f5name, short* signal_buf, unsigned long long &start_time, int& used_albacore_version, string& channel_number_buf, double& digitisation_buf, double& offset_buf, 
double& range_buf, double& sampling_rate_buf, string& fastq_buf, Events_t *events_data);

// // void write_data(int fd, string f5status, vector<M_event>& m_event, short* raw_signals, int signal_len, string read_id, string mfile_path, string m_event_basecall, int left_right_skip_left, int left_right_skip_right);
// /* for writing data to python by pipe*/

void String_Split(string s, char delim, vector<string>& ans);

string String_Replace(string s, const char* oldstr, const char* newstr);

void showAllFiles(const char* dir_name, vector<string>& filename_vec);

string Upper(string str, int length);
/* for capitalizing every char */

float myround(float num, int precision);
/* for rounding up */

void get_event(float* my_raw_signals, char* m_event_basecall, string read_id, string f5status, string mfile_path, string SignalGroup, int moptions_outLevel, int used_albacore_version,const string& fq_seq, short* signals, double* sampling_rate, unsigned long long* start_time, int signal_len, int events_len, 
		Events_t* events_data, int& left_right_skip_left, int& left_right_skip_right);
/* get events from a fast5 file */

vector<M_event> get_EventInfo(int signal_len, int lines_event, short* signals, Events_t* events_data, const string& fq_seq);
/* for getting events information */

float* cumsum(int mysignals_len, float* mysignals);
/* for cumulative sum of the elements */

vector<int> get_extreme_N(int m_signal_dif_len, float* m_signal_dif, int n_splits, unsigned long long p_signal_start, unsigned long long p_signal_end, int min_signal_num);
/* for getting extreme N */

template <typename T>
vector<int> reargsort(const vector<T>& v);
/* for returning the index of reverse sorted array */

float mean(int signal_len, const short* raw_signals, int begin, int end);
/* for calculating the mean of raw_signals */

float variance(int signal_len, const short* raw_signals, int begin, int end, float mean);
/* for calculating the variance of raw_signals */

float stdvariancec(int signal_len, const short* raw_signals, int begin, int end, float mean);
/* for calculating the standard variance of raw_signals */

// float* mnormalized(string mfile_path, int signal_len, short* raw_signals, int events_len, Events_t* events_data);
/* for normalizing raw_signals */
void mnormalized(float* my_raw_signals, string mfile_path, int signal_len, short* raw_signals, int events_len, Events_t *events_data);

template <typename T>
float median(vector<T>& vt);
/* for getting the median */

template <typename T>
void cal_mean_stdv(vector<T>& vt, float& mean, float& stdv);
/* for getting the mean and stdv*/

#endif