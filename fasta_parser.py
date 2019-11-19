#!/usr/bin/python3
import os
import getopt, sys

if __name__ == '__main__':

    opts, args = getopt.getopt(sys.argv[1:], "hi:o:", ["help", "input=", "output="])
    for o, a in opts:
        if o in ("-h", "--help"):
            sys.exit()
        if o in ("-i", "--input"):
            fastA_path = a
            output_path = a + '.dump'
        if o in ("-o", "--output"):
            output_path = a

    amb_info = []
    with open(fastA_path + '.amb', mode='r') as f_amb:
        line = f_amb.readline()
        amb_line = line.split(' ')
        print('fa len:', amb_line[0], 'chr num:', amb_line[1])
        while True:
            line = f_amb.readline()
            if line:
                amb_line = line.split(' ')
                print('N pos:', amb_line[0], 'N len:', amb_line[1])
                amb_info.append(amb_line)
            else:
                break
    print(amb_info)

    ann_info = []
    ann_size = 0
    with open(fastA_path + '.ann', mode='r') as f_ann:
        line = f_ann.readline()
        ann_line = line.split(' ')
        print('fa len:', ann_line[0], 'chr num:', ann_line[1])
        while True:
            line = f_ann.readline()
            if line:
                ann_line1 = line.split(' ')
                line = f_ann.readline()
                ann_line2 = line.split(' ')
                print(ann_line1[1], ' pos:', ann_line2[0], ':', ann_line2[1])
                ann_info.append([ann_line1[1], int(ann_line2[0]), int(ann_line2[1])])
                ann_size += 1
            else:
                break
    print(ann_info)

    nst_dent4_table = ['A', 'C', 'G', 'T']

    pac_step = (64 * 1024 * 1024)
    with open(fastA_path + '.pac', mode='rb') as f_pac:
        with open(output_path + '.temp', mode='w') as f_tfa:
            while True:
                pac_srt = f_pac.read(pac_step)
                if pac_srt:
                    dump_str = ''
                    for pca_char in pac_srt:
                        dump_str += nst_dent4_table[int(pca_char) >> 6 & 3]
                        dump_str += nst_dent4_table[int(pca_char) >> 4 & 3]
                        dump_str += nst_dent4_table[int(pca_char) >> 2 & 3]
                        dump_str += nst_dent4_table[int(pca_char) >> 0 & 3]
                    f_tfa.write(dump_str)
                else:
                    break

            for amb_t in amb_info:
                f_tfa.seek(int(amb_t[0]))
                N_str = ''
                for i in range(int(amb_t[1])):
                    N_str += 'N'
                f_tfa.write(N_str)

    with open(output_path + '.temp', mode='r') as f_tfa:
        with open(output_path, mode='w') as f_dfa:
            for ann_t in ann_info:
                f_dfa.write('>')
                f_dfa.write(ann_t[0])
                f_dfa.write('\n')

                dump_str = f_tfa.read(ann_t[2])
                f_dfa.write(dump_str)
                f_dfa.write('\n')

    os.remove(output_path + '.temp')

    sys.exit(0)
