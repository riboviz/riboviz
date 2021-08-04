"""
Estimate the time of a task based on the sample size, gff size and
the number of process the task will run in parallel.
The estimation is based on experiments on some down sampled datasets,
and uses linear regression and interpolation to estimate.
The experiments shows that the estimation may be 50% underestimate,
and different CPU model has different IPC, so the result is for
reference only.

The input is:

* sample_size: The size of the sample file in Bytes.
* gff_lines: The number of lines of the gff file.
* num_process: The number of process that the task will run with
* task: The name of task, should be cutadapt, hisat2orf, hisat2rrna,
        samviewsort or bamtoh5. Case insensitive.
* experiment data file: `time_estimation.csv`, contains result of the
        down sampled files, uses for estimation.

The output should be the estimation time in seconds, by exit code:
"""
import math
import os

import numpy as np
from scipy import interpolate, stats


class record:
    def __init__(self, task, sample_size, core, gff_size):
        self.task = task
        self.sample_size = sample_size
        self.core = core
        self.gff_size = gff_size
        self.raw_time = []
        self.raw_memory = []

    @property
    def avg_time(self):
        assert len(self.raw_time) != 0
        return sum(self.raw_time) / len(self.raw_time)

    @property
    def avg_memory(self):
        if len(self.raw_memory) == 0:
            return 1
        if sum(self.raw_memory) / len(self.raw_memory) == 0:
            return 1
        return sum(self.raw_memory) / len(self.raw_memory)


def read_dataset():
    with open('riboviz/time_estimation.csv', 'r') as fin:
        all_record = []
        for i in fin.readlines()[1:]:
            tmp = i.split(',')
            all_record.append(record(tmp[0], float(tmp[1]), int(tmp[2]), int(tmp[3])))
            all_record[-1].raw_time.append(float(tmp[4]))
            all_record[-1].raw_memory.append(float(tmp[5]))
    return all_record


def search_result(records, task, core):
    ret = []
    for i in records:
        if i.task == task and i.core == core:
            ret.append(i)
    ret.sort(key=lambda x: x.sample_size * 10000 + x.gff_size)
    return ret


def estimate(task, sample_size, gff_size, core):
    gff_related = ['bamtoh5', 'samviewsort', 'hisat2orf']
    gff_not_related = ['cutadapt', 'hisat2rrna']
    if task in gff_not_related:
        time = _estimate_not_related(task, sample_size, core)
    else:
        time = _estimate_related(task, sample_size, gff_size, core)
    return time


def _estimate_related(task, sample_size, gff_size, core):
    # estimate time on gff related items
    all_record = read_dataset()
    estimate_list = []
    if core in [1, 2, 4, 8, 16, 32]:
        for j in search_result(all_record, task, core):
            estimate_list.append((j.sample_size, j.gff_size, j.avg_time, j.avg_memory))
    else:
        tmp = []
        for cores in [1, 2, 4, 8, 16, 32]:
            tt = []
            for j in search_result(all_record, task, cores):
                tt.append((j.sample_size, j.gff_size, j.avg_time, j.avg_memory))
            tmp.append(tt)
        # get possible time and memory value by interpolation
        for idx in range(len(tmp[0])):
            pp = []
            for idy, i in enumerate([1, 2, 4, 8, 16, 32]):
                assert tmp[0][idx][0] == tmp[idy][idx][0]
                assert tmp[0][idx][1] == tmp[idy][idx][1]
                pp.append(tmp[idy][idx][2])
            spl_t = interpolate.interp1d([1, 2, 4, 8, 16, 32], pp, bounds_error=False, fill_value="extrapolate")
            estimate_list.append((tmp[0][idx][0], tmp[0][idx][1], float(spl_t(core))))
    # print(estimate_list)
    min_sample = min([x[0] for x in estimate_list])
    max_sample = max([x[0] for x in estimate_list])
    min_gff = min([x[1] for x in estimate_list])
    max_gff = max([x[1] for x in estimate_list])
    if min_sample < sample_size < max_sample and min_gff < gff_size < max_gff:
        int1d_t = interpolate.interp2d([x[0] for x in estimate_list], [x[1] for x in estimate_list],
                                       [x[2] for x in estimate_list], bounds_error=False)
        time = float(int1d_t(sample_size, gff_size))
        print("Time: {0}".format(time))
    elif min_sample < sample_size < max_sample:
        # we first estimate the gff size
        ax = list(set([x[1] for x in estimate_list]))
        ay = []
        for gff in ax:
            tmpx = []
            tmpy = []
            for i in estimate_list:
                if i[1] == gff:
                    tmpx.append(i[0])
                    tmpy.append(i[2])
            int1d_tmp = interpolate.interp1d(tmpx, tmpy)
            ay.append(int1d_tmp(sample_size))

        int1d_tmp = interpolate.interp1d(ax, ay, bounds_error=False, fill_value="extrapolate")
        time = float(int1d_tmp(gff_size))
    elif min_gff < gff_size < max_gff:
        # we first estimate the sample size
        ax = list(set([x[0] for x in estimate_list]))
        ay = []
        for sample in ax:
            tmpx = []
            tmpy = []
            for i in estimate_list:
                if i[0] == sample:
                    tmpx.append(i[1])
                    tmpy.append(i[2])
            int1d_tmp = interpolate.interp1d(tmpx, tmpy)
            ay.append(int1d_tmp(gff_size))

        int1d_tmp = interpolate.interp1d(ax, ay, bounds_error=False, fill_value="extrapolate")
        time = float(int1d_tmp(sample_size))
    else:
        # we first estimate the gff
        ax = list(set([x[1] for x in estimate_list]))
        ay = []
        for gff in ax:
            tmpx = []
            tmpy = []
            for i in estimate_list:
                if i[1] == gff:
                    tmpx.append(i[0])
                    tmpy.append(i[2])
            int1d_tmp = interpolate.interp1d(tmpx, tmpy, bounds_error=False, fill_value="extrapolate")
            ay.append(int1d_tmp(sample_size))
        if task != 'hisat2orf':
            int1d_tmp = interpolate.interp1d(ax, ay, bounds_error=False, fill_value="extrapolate")
            time = float(int1d_tmp(gff_size))
        else:
            slope, intercept, r, p, stderr = stats.linregress(ax, ay)
            time = slope * gff_size + intercept
    if time < 120:
        return 120
    return time


def _estimate_not_related(task, sample_size, core):
    # estimate time on gff not related items
    all_record = read_dataset()
    estimate_list = []
    if core in [1, 2, 4, 8, 16, 32]:
        for j in search_result(all_record, task, core):
            estimate_list.append((j.sample_size, j.avg_time, j.avg_memory))
    else:
        tmp = []
        for cores in [1, 2, 4, 8, 16, 32]:
            tt = []
            for j in search_result(all_record, task, cores):
                tt.append((j.sample_size, j.avg_time))
            tmp.append(tt)
        # get possible time and memory value by interpolation
        for idx in range(len(tmp[0])):
            pp = []
            for idy, i in enumerate([1, 2, 4, 8, 16, 32]):
                assert tmp[0][idx][0] == tmp[idy][idx][0]
                pp.append(tmp[idy][idx][1])
            spl_t = interpolate.interp1d([1, 2, 4, 8, 16, 32], pp)
            estimate_list.append((tmp[0][idx][0], float(spl_t(core))))
    # print(estimate_list)
    # use 1d inerpolate
    print(estimate_list)
    int1d_t = interpolate.interp1d([x[0] for x in estimate_list],
                                   [x[1] for x in estimate_list], bounds_error=False, fill_value="extrapolate")
    time = float(int1d_t(sample_size))
    print("Time: {0}".format(time))
    return time


def estimate_time(sample_size, gff_lines, num_process, task):
    if task not in ['bamtoh5', 'samviewsort', 'hisat2orf', 'cutadapt', 'hisat2rrna']:
        return 7200

    return estimate(task, sample_size, gff_lines / 3 - 1, num_process)
