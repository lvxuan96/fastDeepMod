import h5py
import numpy as np
import struct
import time
import glob
import os
import sys
import ctypes
from ctypes import *

import myCom
# from . import myCom   
from distutils.version import LooseVersion
from collections import defaultdict

import multiprocessing

fast5_channel_id = 'UniqueGlobalKey/channel_id'
fast5_analysis = ''.join(['/', myCom.analyses_base])
fast5_events = myCom.basecall_events_base
fast5_rawReads = ''.join(['/', myCom.raw_base, '/', myCom.reads_base])
fast5_basecall_fq = myCom.basecall_fastq_base
fast5_signal = myCom.signal_base

pre_base_str = 'rnn.pred.ind'


def mDetect_manager(moptions):
    # pmanager = multiprocessing.Manager();
    # get input folder
    while (not moptions['wrkBase'] == None) and len(moptions['wrkBase']) > 0 and moptions['wrkBase'][-1] in ['/', '\\']:
        moptions['wrkBase'] = moptions['wrkBase'][:-1]

    # need to make prediction of modification
    if moptions['predDet'] == 1:
        # get well-trained model
        if moptions['modfile'].rfind('/') == -1:
            moptions['modfile'] = [moptions['modfile'], './']
        else:
            moptions['modfile'] = [moptions['modfile'], moptions['modfile'][:moptions['modfile'].rfind('/') + 1]]
        start_time = time.time()

        # output folder
        if not os.path.isdir(moptions['outFolder'] + moptions['FileID']):
            os.system('mkdir -p ' + moptions['outFolder'] + moptions['FileID'])

        detect_handler(moptions)


def detect_handler(moptions):

    sp_options = defaultdict()
    sp_options['ctfolderid'] = 0
    sp_options['ctfolder'] = moptions['outFolder'] + moptions['FileID'] + '/' + str(sp_options['ctfolderid'])
    if not os.path.isdir(sp_options['ctfolder']):
        os.system('mkdir ' + sp_options['ctfolder'])
    sp_options['Mod'] = []

    mDetect1(moptions, sp_options)

def mDetect1(moptions, sp_options):
    f5data = get_Event_Signals(moptions, sp_options)

def get_Event_Signals(moptions, sp_options):

    f5data = {}
    sp_options["Error"] = defaultdict(list)
    sp_options["get_albacore_version"] = defaultdict(int)

    mylib = CDLL('preprocess.so')
    
    f5name = "read116.fast5"

    signal_len = c_int()
    events_len = c_int()
    mylib.get_f5len(c_char_p(bytes(f5name,'utf-8')), byref(signal_len), byref(events_len))
    print('signal_len=', signal_len.value, 'events_len=', events_len.value)
    

    
    class M_event(Structure):
        _fields_ = [("mean",c_float),
                ("stdv", c_float),
                ("start", c_ulonglong),
                ("length",c_ulonglong),
                ("model_state", c_char * 6)
                ]
    m_event_buf = (M_event * events_len.value)()
    m_event_basecall = create_string_buffer(events_len.value)
    read_id_buf = create_string_buffer(1000)
    raw_signals = (c_float * signal_len.value)()
    left_right_skip_left = c_int()
    left_right_skip_right = c_int()

    mylib.get_event_signals.restypes = c_int
    mylib.get_event_signals.argtypes = [c_char_p, c_char_p, c_char_p, POINTER(c_float), POINTER(c_int), POINTER(c_int), POINTER(M_event)]
    returnstatus = mylib.get_event_signals(c_char_p(bytes(f5name,'utf-8')), m_event_basecall, read_id_buf,
        raw_signals, byref(left_right_skip_left), byref(left_right_skip_right), m_event_buf)
    # print("read id=", read_id_buf.value.decode("utf-8"))
    # print("m_event_basecall=", m_event_basecall.value.decode("utf-8"))
    # print("raw_signals=", list(raw_signals))
    # print("left_right_skip_left=",left_right_skip_left.value, "left_right_skip_right=",left_right_skip_right.value)
    # m_event_buf = list(m_event_buf)
    print(m_event_buf[45836].mean, m_event_buf[45836].length, m_event_buf[45836].model_state.decode("utf-8"))
    
    # bufsize=memoryview(m_event_buf).itemsize
    # bufform = memoryview(m_event_buf)
    # print(bufsize)
    # print(bufform.format)
    # m_event = np.ctypeslib.as_array(m_event_buf, events_len.value)
   
    
    
    

    return f5data







moptions = {}
moptions['basecall_1d'] = 'Basecall_1D_000'
# moptions['basecall_1d'] = ['Basecall_1D_000']
moptions['basecall_2strand'] = 'BaseCalled_template'

moptions['outLevel'] = myCom.OUTPUT_WARNING
moptions['outLevel'] = myCom.OUTPUT_INFO

moptions['modfile'] = '../../mod_output/train1/2/mod_train'

moptions['fnum'] = 3;
moptions['hidden'] = 100;
moptions['windowsize'] = 21;

moptions['threads'] = 8
moptions['threads'] = 1
moptions['files_per_thread'] = 500

moptions['wrkBase'] = "../PRJEB31789-Chlamydomonas_0-reads-downloads-pass-0/0"
moptions["FileID"] = "pred"
moptions['predDet'] = 1
moptions['recursive'] = 1
moptions['outFolder'] = "../out"

# start_wall_time = time.time()
# start_cpu_time = time.clock()
mDetect_manager(moptions)
# end_wall_time = time.time()
# end_cpu_time = time.clock()
# print("wall time: ", end_wall_time - start_wall_time)
# print("cpu time: ", end_cpu_time - start_cpu_time)

