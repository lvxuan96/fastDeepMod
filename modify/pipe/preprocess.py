import h5py
import numpy as np
import struct
import time
import glob
import os
import sys

import myCom
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

    try:
        os.mkfifo("/tmp/numfile.pipe")

    except OSError as e:
        print("OS Error: {0}".format(e))

    rf = os.open("/tmp/numfile.pipe", os.O_RDONLY)

    numfile, = struct.unpack("i", os.read(rf, 4))

    os.close(rf)
    os.unlink("/tmp/numfile.pipe")

    read_path = "/tmp/data.pipe"
    try:
        os.mkfifo(read_path)

    except OSError as e:
        print("OS Error: {0}".format(e))

    rf = os.open(read_path, os.O_RDONLY)
    print("os.open successfully")

    for i in range(numfile):
        data2py(f5data, rf)
        
    os.close(rf)
    os.unlink("/tmp/data.pipe")
    return f5data

def data2py(f5data, rf):
    # tyh: need to change
    starttime = time.time()

    len_read_id, = struct.unpack("Q", os.read(rf, 8))
    read_id = bytes.decode(os.read(rf, len_read_id))

    len_mfile_path, = struct.unpack("Q", os.read(rf, 8))
    mfile_path = bytes.decode(os.read(rf, len_mfile_path))

    len_m_event_basecall, = struct.unpack("Q", os.read(rf, 8))
    all_m_event_basecall = ""
    # max size for a named pipe is 64kb
    while len_m_event_basecall:
        m_event_basecall = os.read(rf, len_m_event_basecall)

        all_m_event_basecall += bytes.decode(m_event_basecall)
        len_m_event_basecall -= len(m_event_basecall)
        # print("len ", len_m_event_basecall)

    left_right_skip_left, = struct.unpack("i", os.read(rf, 4))
    left_right_skip_right, = struct.unpack("i", os.read(rf, 4))

    len_m_event, = struct.unpack("Q", os.read(rf, 8))
    lines_signals, = struct.unpack("i", os.read(rf, 4))
    print("len_m_event ", len_m_event)
    print("lines_signals ", lines_signals)

    m_event = []
    for i in range(len_m_event):
        mean, = struct.unpack("f", os.read(rf, 4))
        stdv, = struct.unpack("f", os.read(rf, 4))
        start, = struct.unpack("Q", os.read(rf, 8))
        length, = struct.unpack("Q", os.read(rf, 8))
        model_state = os.read(rf, 5)
        # os.read(rf, 3)
        model_state = bytes.decode(model_state)
        m_event.append((mean, stdv, start, length, model_state))

    m_event = np.array(m_event, dtype=[('mean', '<f4'), ('stdv', '<f4'), ('start', np.uint64), ('length', np.uint64),
                                  ('model_state', 'U5')])
    raw_signals = []
    for i in range(lines_signals):
        raw_signal, =struct.unpack("h",os.read(rf, 2))
        raw_signals.append(raw_signal)

    print("Having got all data from CPP")

    f5data[read_id] = (m_event_basecall, m_event, raw_signals, mfile_path, (left_right_skip_left, left_right_skip_right))

    endtime = time.time()
    with open("DataTransmission.log", "at") as f:
        f.write("data2py " + str(endtime - starttime) + "\n")




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

starttime = time.time()
mDetect_manager(moptions)
endtime = time.time()
print(endtime - starttime)

