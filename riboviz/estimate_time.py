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
from scipy import interpolate, stats


class Record:
    """
    Class `Record` is a data structure to store the records in the csv dataset file.
    It contains the following properties:
    * task: The name of task, like `cutadapt`, `hisat2orf`.
    * sample_size: The size (file size in Bytes) of the sample file.
    * core: The number of cores (processes) that the task is running with.
    * gff_size: The number of items in the gff file.
    * raw_time: The time (list) when running the task under the given sample size, processes and gff file.
    * raw_memory: The memory (list) when running the task under the given sample size, processes and gff file.
    * avg_time: The average of raw_time. The length of raw_time must be bigger than zero when calling this property.
    * avg_memory: The average of raw memory.
    """
    def __init__(self, task, sample_size, core, gff_size):
        """
        Init function. Initialize the properties.

        :param task: The name of task
        :type task: str or unicode
        :param sample_size: The file size of sample in bytes
        :type sample_size: int
        :param core: The number of cores (processes)
        :type core: int
        :param gff_size: The number of items in the gff file
        :type gff_size: int
        """
        self.task = task
        self.sample_size = sample_size
        self.core = core
        self.gff_size = gff_size
        self.raw_time = []
        self.raw_memory = []

    @property
    def avg_time(self):
        """
        Function as a property.
        Return the average of self.raw_time.
        If there are no records in self.raw_time, then an error will be raised.

        :return: The average time.
        :rtype: float
        """
        assert len(self.raw_time) != 0
        return sum(self.raw_time) / len(self.raw_time)

    @property
    def avg_memory(self):
        """
        Function as a property.
        Return the average of self.raw_memory.
        If there are no records in self.raw_memory, then 1 will be returned.

        :return: The average memory.
        :rtype: float
        """
        if len(self.raw_memory) == 0:
            return 1
        if sum(self.raw_memory) / len(self.raw_memory) == 0:
            return 1
        return sum(self.raw_memory) / len(self.raw_memory)

    def add_item(self, time, memory):
        """
        Add a item in self.raw_time and self.raw_memory,
        so the average time can be visited from other functions.

        :param time: The execution time in seconds.
        :type time: int or float
        :param memory: The memory consumption in M Bytes.
        :type memory: int or float
        """
        self.raw_time.append(time)
        self.raw_memory.append(memory)


def read_dataset(file='riboviz/time_estimation.csv'):
    """
    Read the dataset and save the data into a list of Records.

    :param file: The path of dataset file, default is `riboviz/time_estimation.csv`
    :type file: str
    :return: The data in Record data structure
    :rtype: list(Record)
    """
    with open(file, 'r') as fin:
        all_record = []
        for line in fin.readlines()[1:]:
            # The first line of a CSV file is the header, so we omit it.
            tmp = line.split(',')
            all_record.append(Record(tmp[0], float(tmp[1]), int(tmp[2]), int(tmp[3])))
            all_record[-1].add_item(float(tmp[4]), float(tmp[5]))
    return all_record


def search_result(records, task, core):
    """
    Search for the records that have the corresponding task name and number of processes
    in the records list.

    :param file: The list that contains all records.
    :type file: list(Record)
    :param task: The name of task that need to filter.
    :type task: str
    :param core: The number of processes that need to filter.
    :type core: int
    :return: The target records
    :rtype: list(Record)
    """
    ret = []
    for i in records:
        if i.task == task and i.core == core:
            ret.append(i)
    # sort the return list, so that the list is in a correct order.
    ret.sort(key=lambda x: x.sample_size * 10000 + x.gff_size)
    return ret


def estimate(task, sample_size, gff_size, core):
    """
    Estimate the execution time based on the task name, sample size, gff items and number of cores.

    :param task: The name of task that need to filter.
    :type task: str
    :param sample_size: The file size of sample in Bytes.
    :type sample_size: int
    :param gff_size: The number of items in the gff file.
    :type gff_size: int
    :param core: The number of processes that run in parallel.
    :type core: int
    :return: The estimated time in seconds
    :rtype: int
    """
    # Three tasks are related to gff size, so they are called `_estimate_related`
    gff_related = ['bamtoh5', 'samviewsort', 'hisat2orf']
    gff_not_related = ['cutadapt', 'hisat2rrna']
    if task in gff_related:
        time = _estimate_related(task, sample_size, gff_size, core)
    if task in gff_not_related:
        time = _estimate_not_related(task, sample_size, core)
    return time


def _estimate_related(task, sample_size, gff_size, core):
    """
    Estimate the execution time of tasks that are related to the gff_size.
    We will first use a interpolation to fit the core, then use another interpolation to fit the sample size and gff size.

    :param task: The name of task that need to filter.
    :type task: str
    :param sample_size: The file size of sample in Bytes.
    :type sample_size: int
    :param gff_size: The number of items in the gff file.
    :type gff_size: int
    :param core: The number of processes that run in parallel.
    :type core: int
    :return: The estimated time in seconds
    :rtype: int
    """
    all_record = read_dataset()
    estimate_list = []
    # we first fit the core.
    # If core is already in [1, 2, 4, 8, 16, 32], we can then use the existed data.
    # otherwise, we need to estimate the value based on the data of core = [1, 2, 4, 8, 16, 32].
    if core in [1, 2, 4, 8, 16, 32]:
        for j in search_result(all_record, task, core):
            estimate_list.append((j.sample_size, j.gff_size, j.avg_time, j.avg_memory))
    else:
        # we first get all data into a 2d list tmp, then run a interpolation to estimate.
        tmp = []
        for cores in [1, 2, 4, 8, 16, 32]:
            # tmp_line is a temporial variable to store the second dimension of tmp
            tmp_line = []
            for j in search_result(all_record, task, cores):
                tmp_line.append((j.sample_size, j.gff_size, j.avg_time, j.avg_memory))
            tmp.append(tmp_line)
        # get possible time and memory value by interpolation
        for idx in range(len(tmp[0])):
            tmp_list = []
            for idy, i in enumerate([1, 2, 4, 8, 16, 32]):
                tmp_list.append(tmp[idy][idx][2])
            spl_t = interpolate.interp1d([1, 2, 4, 8, 16, 32], tmp_list, bounds_error=False, fill_value="extrapolate")
            estimate_list.append((tmp[0][idx][0], tmp[0][idx][1], float(spl_t(core))))

    # then, we need to do a 2d interpolation, to estimate the time based on the core and sample size.
    min_sample = min([x[0] for x in estimate_list])
    max_sample = max([x[0] for x in estimate_list])
    min_gff = min([x[1] for x in estimate_list])
    max_gff = max([x[1] for x in estimate_list])
    # if the target core and sample size is in the boundary of existing dataset, we can call the 2d interpolation
    # directly, otherwise, we have to manually do two 1d-interpolation.
    if min_sample < sample_size < max_sample and min_gff < gff_size < max_gff:
        int1d_t = interpolate.interp2d([x[0] for x in estimate_list], [x[1] for x in estimate_list],
                                       [x[2] for x in estimate_list], bounds_error=False)
        time = float(int1d_t(sample_size, gff_size))
        # print("Time: {0}".format(time))
    elif min_sample < sample_size < max_sample:
        # we first estimate the gff size, because gff size is not in the boundary
        ax = list(set([x[1] for x in estimate_list]))
        ay = []
        for gff in ax:
            # for different gff size, estimate the sample size vs execution time
            # tmpx is the x-axis of different gff sizes
            # tmpy is the y-axis (time)
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
        # similar as above, but we first estimate the sample size
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
        # if both sample size and gff size are out of boundary, then we first estimate the gff
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
        if True:
            int1d_tmp = interpolate.interp1d(ax, ay, bounds_error=False, fill_value="extrapolate")
            time = float(int1d_tmp(gff_size))
        else:
        # set to false to use the linear regression
            slope, intercept, r, p, stderr = stats.linregress(ax, ay)
            time = slope * gff_size + intercept

    if time < 120:
        return 120
    return time


def _estimate_not_related(task, sample_size, core):
    """
    Estimate the execution time of tasks that are not related to the gff_size.
    We will first use a interpolation to fit the core, then use another interpolation to fit the sample size.

    :param task: The name of task that need to filter.
    :type task: str
    :param sample_size: The file size of sample in Bytes.
    :type sample_size: int
    :param core: The number of processes that run in parallel.
    :type core: int
    :return: The estimated time in seconds
    :rtype: int
    """
    # we first fit the core.
    # If core is already in [1, 2, 4, 8, 16, 32], we can then use the existed data.
    # otherwise, we need to estimate the value based on the data of core = [1, 2, 4, 8, 16, 32].
    all_record = read_dataset()
    estimate_list = []
    if core in [1, 2, 4, 8, 16, 32]:
        for j in search_result(all_record, task, core):
            estimate_list.append((j.sample_size, j.avg_time, j.avg_memory))
    else:
        tmp = []
        for cores in [1, 2, 4, 8, 16, 32]:
            tmp_line = []
            for j in search_result(all_record, task, cores):
                tmp_line.append((j.sample_size, j.avg_time))
            tmp.append(tmp_line)
        # get possible time and memory value by interpolation
        for idx in range(len(tmp[0])):
            possible_time = []
            for idy, _ in enumerate([1, 2, 4, 8, 16, 32]):
                assert tmp[0][idx][0] == tmp[idy][idx][0]
                possible_time.append(tmp[idy][idx][1])
            spl_t = interpolate.interp1d([1, 2, 4, 8, 16, 32], possible_time)
            estimate_list.append((tmp[0][idx][0], float(spl_t(core))))
    # print(estimate_list)
    # use 1d inerpolate
    # print(estimate_list)
    int1d_t = interpolate.interp1d([x[0] for x in estimate_list],
                                   [x[1] for x in estimate_list], bounds_error=False, fill_value="extrapolate")
    time = float(int1d_t(sample_size))
    # print("Time: {0}".format(time))
    if time < 120:
        return 120
    return time


def estimate_time(sample_size, gff_lines, num_process, task):
    """
    Interface for time estimation, do some validity check.

    :param sample_size: The file size of sample in Bytes.
    :type sample_size: int
    :param gff_lines: The number of items in the gff file.
    :type gff_lines: int
    :param core: The number of processes that run in parallel.
    :type core: int
    :param task: The name of task that need to filter.
    :type task: str
    :return: The estimated time in seconds
    :rtype: int
    """
    # if task is not leagal, then return a default error time, 7200
    if task not in ['bamtoh5', 'samviewsort', 'hisat2orf', 'cutadapt', 'hisat2rrna']:
        return 7200
    # three lines in the gff file is a item, so the lines needs to be divided by three.
    return estimate(task, sample_size, gff_lines / 3 - 1, num_process)
