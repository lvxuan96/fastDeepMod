preprocess:preprocess.h preprocess.cpp
	g++ -I /home/lvxuan/anaconda3/envs/mdeepmod/include -L /home/lvxuan/anaconda3/envs/mdeepmod/lib \
	 -fPIC -shared -o preprocess.so  preprocess.cpp -std=c++11 -lhdf5_cpp -lhdf5 -fopenmp
