

import glob
import os
import datetime

class Benchmarks:

    def __init__(self, args):
        self.benchmarks_dir = args.benchmarks_dir
        self.outfn = args.outfn

        self.file_dic = dict()
        self.file_path_finder()

        self.od_step_s = self.step_ordering()
        self.write_table()

    def file_path_finder(self):
        for fp in glob.glob(os.path.join(self.benchmarks_dir, '*.txt')):
            fn = os.path.basename(fp)
            step = fn.split('.')[0]
            individual = '_'.join(fn.split('.')[1].split('_')[:-1])
            ctime = os.path.getctime(fp)
            rtime = datetime.datetime.fromtimestamp(ctime)
            self.file_dic.setdefault(step, {}).setdefault(individual, {}).setdefault('path', fp)
            self.file_dic.setdefault(step, {}).setdefault(individual, {}).setdefault('ctime', ctime)
            self.file_dic.setdefault(step, {}).setdefault(individual, {}).setdefault('rtime', rtime)

    def step_ordering(self):
        time_dic = dict()
        for step, individual_dic in self.file_dic.items():
            for individual, info_dic in individual_dic.items():
                time_dic.setdefault(step, float(info_dic['ctime']))
        #
        od_step_s = list()
        for step in sorted(time_dic, key=time_dic.get, reverse=False):
            od_step_s.append((step, time_dic[step]))
        return od_step_s


    def write_table(self):
        running_second = 0
        outfh = open(self.outfn, 'w')
        headers = ['timestamp']
        headers.append('start_time')
        headers.append('step_name')
        headers.append('running_second(avg)')
        headers.append('running_time(avg)')
        headers.append('use_memory(avg)')
        headers.append('io_in(avg)')
        headers.append('io_out(avg)')
        outfh.write('{0}\n'.format('\t'.join(headers)))
        for step, ctime in self.od_step_s:
            items = [ctime]
            items.append(datetime.datetime.fromtimestamp(ctime))
            items.append(step)
            items.append(round(self.cal_avg(step, 's'),2))
            items.append(datetime.timedelta(seconds=self.cal_avg(step, 's')))
            items.append(round(self.cal_avg(step, 'mem'),2))
            items.append(round(self.cal_avg(step, 'ioin'),2))
            items.append(round(self.cal_avg(step, 'ioout'),2))
            running_second += float(self.cal_avg(step, 's'))
            outfh.write('{0}\n'.format('\t'.join([str(x) for x in items])))

        running_time = datetime.timedelta(seconds=running_second)
        ppl_running_second = self.od_step_s[-1][1] - self.od_step_s[0][1]
        ppl_running_time = datetime.timedelta(seconds=ppl_running_second)
        items = ['.','.','Sum of working time for all steps in series.','.','.',running_second, running_time,'.','.','.']
        outfh.write('{0}\n'.format('\t'.join([str(x) for x in items])))
        items = ['.','.','Pipeline running time including reruns due to errors.','.','.',ppl_running_second, ppl_running_time,'.','.','.']
        outfh.write('{0}\n'.format('\t'.join([str(x) for x in items])))
        print('Pipeline Running Real Time : {0:,} {1}'.format(ppl_running_second, ppl_running_time))

        outfh.close()

    def cal_avg(self, step, key):
        values = list()
        for individual in self.file_dic[step]:
            fp = self.file_dic[step][individual]['path']
            data_dic = self.fp_parser(fp)
            values.append(float(data_dic[key]))
        return sum(values)/len(values)

    def cal_sum(self, step, key):
        values = list()
        for individual in self.file_dic[step]:
            fp = self.file_dic[step][individual]['path']
            data_dic = self.fp_parser(fp)
            values.append(float(data_dic[key]))
        return sum(values)

    def fp_parser(self, fp):
        data_dic = dict()
        data_dic.setdefault('s', 0.0)
        data_dic.setdefault('hms', 0.0)
        data_dic.setdefault('mem', 0.0)
        data_dic.setdefault('ioin', 0.0)
        data_dic.setdefault('ioout', 0.0)
        for line in open(fp):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['s']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            if items[idx_dic['s']] not in ['-']:
                data_dic['s'] = items[idx_dic['s']]
            if items[idx_dic['h:m:s']] not in ['-']:
                data_dic['hms'] = items[idx_dic['h:m:s']]
            if items[idx_dic['max_vms']] not in ['-']:
                data_dic['mem'] = items[idx_dic['max_vms']]
            if items[idx_dic['io_in']] not in ['-']:
                data_dic['ioin'] = items[idx_dic['io_in']]
            if items[idx_dic['io_out']] not in ['-']:
                data_dic['ioout'] = items[idx_dic['io_out']]
        return data_dic



def main(args):
    print(args)
    benchmarks = Benchmarks(args)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--benchmarks-dir')
    parser.add_argument('--outfn')
    args = parser.parse_args()
    main(args)

