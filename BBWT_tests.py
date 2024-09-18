import argparse, string, csv, heapq, subprocess
from itertools import product
from pathlib import Path
from collections import Counter
from statistics import stdev

from count_runs import noRuns

import os, sys, shutil

cwd = Path(__file__).parent.absolute()

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("sigma", help="size of the alphabet", default=2, type=int)
    parser.add_argument("max_k", help="maximum strings length", default=10, type=int)
    parser.add_argument("-o", "--outdir", help="existing output directory for generated files",
                        default=cwd.joinpath("outdir"), type=Path)
    args = parser.parse_args()
    sigma = args.sigma
    max_k = args.max_k
    outdir = Path(args.outdir).joinpath("sigma_" + str(sigma) + "_" + str(args.max_k) )
    if not Path(outdir).exists():
        outdir.mkdir()
        print("A new output directory has been created.")
    tempdir = Path("temp"+str(sigma)+"_"+str(max_k))
    if not Path(tempdir).exists():
        tempdir.mkdir()
        print("A new temporary directory has been created.")
    allk_bbwt_bbwtrev(sigma, max_k, outdir, tempdir)
    shutil.rmtree(tempdir)
    exit()

# Lyndon factorization in linear time with Duval algorithm
def Duval(T):
    n,i,factors=len(T),0,[]
    while i < n:
        j = i+1
        k = i
        while (j<n and T[k]<=T[j]):
            if T[k] < T[j]: k = i-1
            k += 1
            j += 1
        while i <= k:
            factors.append(T[i:i+(j-k)])
            i += j-k
    return factors


def write_csv_w_a(filename, row, mode):
    with open(filename, mode) as f:
        wr = csv.writer(f, quoting=csv.QUOTE_ALL)
        row.append('\n')
        wr.writerow(row)
    return

def check_rev(seq):
    i = 0
    j = len(seq) - 1
    while seq[i] == seq[j]:
        i += 1
        j -= 1
        if i >= j:
            return False  # palindrome
    if seq[i] <= seq[j]:
        return True
    return False

def bbwt_bbwtrev(sigma,k, outfile, outfile2, outfile2_2, outfile3, outfile4, outfile4_2, outfile4_2_1, outfile4_2_2, outfile5, outfile6, outfile6_2, outfile6_2_1, outfile6_2_2, outfile7, max_rho, max_diff, min_diff, max_lyn, min_lyn, tempdir):
    rho_list = []
    no_rho_1 = 0
    max_rho_l = 1

    diff_list = []
    max_diff_l = 0
    min_diff_l = 0

    lyn_list = []
    max_lyn_l = 0
    min_lyn_l = 0
    no_lyn_0 = 0

    # create strings |k| with alphabet |sigma|
    it = product(string.ascii_uppercase[:sigma], repeat=k)
    for s in it:
        seq = ''.join(s)
        if check_rev(seq):
            # write forward to file
            f_file = str(tempdir)+"/forward_"+str(k)+".txt"
            r_file = str(tempdir)+"/reverse_"+str(k)+".txt"
            lr_file = str(tempdir)+"/lyn_reverse_"+str(k)+".txt"
            with open(f_file,"w") as file:
                file.write(seq)
            # compute bbwt forward
            bbwt_f = "{exe} -t {input}".format(exe = str(cwd)+"/cais/cais", input= f_file)
            subprocess.check_call(bbwt_f.split())
            # compute forward Lyndon factorization length
            lyn_f = Duval(seq)
            cf_f = len(set(lyn_f))
            # compute BBWT forward number of runs
            rf = noRuns(f_file+".bbwt")
            # compute reverse sequence
            rev = seq[::-1]
            # write reverse to file
            with open(r_file,"w") as file:
                file.write(rev)
            # compute bbwt reverse
            bbwt_r = "{exe} -t {input}".format(exe = str(cwd)+"/cais/cais", input= r_file)
            subprocess.check_call(bbwt_r.split())
            # compute reverse Lyndon factorization length
            lyn_r = Duval(rev)
            cf_r = len(set(lyn_r))
            # compute BBWT reverse number of runs
            rr = noRuns(r_file+".bbwt")

            # compute diff
            diff = rf-rr

            ## RUNS RATIO
            # compute rho
            r = round(max((rf / rr), (rr / rf)), 2)
            # chekf if rho == 1
            if r==1: no_rho_1 += 1
            # append rho to rho list
            rho_list.append(r)
            # check max_rho
            if max_rho_l < r: # if r not larger than the max_rho_l it cannot be larger than max_rho
                max_rho_l = r 
                # compute forward BWT to check standard words 
                bwt_f = "{exe} -b {input}".format(exe = str(cwd)+"/cais/cais", input= f_file)
                subprocess.check_call(bwt_f.split())
                nrf = noRuns(f_file+".bwt")
                std = (nrf == 2)
                if std == False:
                    for f in lyn_r:
                        if len(f)>2:
                            with open(lr_file,"w") as file:
                                file.write(f)
                                # compute forward BWT to check standard words 
                                bwt_lyn_r = "{exe} -b {input}".format(exe = str(cwd)+"/cais/cais", input= f_file)
                                subprocess.check_call(bwt_lyn_r.split())
                                if (noRuns(f_file+".bwt")!=2):
                                    print (f + " is not a standard word " + str(k) + " " + str(r))
                write_csv_w_a(outfile2, ["n","s","rho", "diff", "r_f", "r_r", "standard", "Lyn_r"], 'w')
                write_csv_w_a(outfile2, [str(k),seq, str(r), str(diff), str(rf), str(rr), str(std), lyn_r], 'a')
                
                if max_rho < r:
                    max_rho = r
                    write_csv_w_a(outfile2_2, ["n","s","rho", "diff", "r_f", "r_r", "standard", "Lyn_r"], 'w')
                    write_csv_w_a(outfile2_2, [str(k),seq, str(r), str(diff), str(rf), str(rr), str(std), lyn_r], 'a')
                
            elif max_rho_l == r:
                # compute forward BWT to check standard words 
                bwt_f = "{exe} -b {input}".format(exe = str(cwd)+"/cais/cais", input= f_file)
                subprocess.check_call(bwt_f.split())
                nrf = noRuns(f_file+".bwt")
                std = (nrf == 2)
                if std == False:
                    for f in lyn_r:
                        if len(f)>2:
                            with open(lr_file,"w") as file:
                                file.write(f)
                                # compute forward BWT to check standard words 
                                bwt_lyn_r = "{exe} -b {input}".format(exe = str(cwd)+"/cais/cais", input= f_file)
                                subprocess.check_call(bwt_lyn_r.split())
                                if (noRuns(f_file+".bwt")!=2):
                                    print (f + " is not a standard word " + str(k) + " " + str(r))
                                    
                write_csv_w_a(outfile2, [str(k),seq, str(r), str(diff), str(rf), str(rr), str(std), lyn_r], 'a')
                if max_rho == r:
                    write_csv_w_a(outfile2_2, [str(k),seq, str(r), str(diff), str(rf), str(rr), str(std), lyn_r], 'a')
            
            ## DIFFERENCE    
            # check if diff = 0 is equal to rho=1
            if diff == 0 and r != 1: print("BUG: rho and diff do not coincide")
            # append diff to diff_list
            diff_list.append(diff)
            # check max_diff
            if max_diff_l < diff:
                max_diff_l = diff
                # compute forward BWT to check standard words 
                command = "{exe} -b {input}".format(exe = str(cwd)+"/cais/cais", input= f_file)
                subprocess.check_call(command.split())
                nrf = noRuns(f_file+".bwt")
                std = (nrf == 2)
                write_csv_w_a(outfile4, ["n","s", "diff", "rho", "r_f", "r_r", "standard"], 'w')
                write_csv_w_a(outfile4, [str(k),seq, str(diff), str(r), str(rf), str(rr), str(std)], 'a')

                if max_diff < diff:
                    max_diff = diff
                    write_csv_w_a(outfile4_2_1, ["n","s", "diff", "rho", "r_f", "r_r", "standard"], 'w')
                    write_csv_w_a(outfile4_2_1, [str(k),seq, str(diff), str(r), str(rf), str(rr), str(std)], 'a')
            elif diff == max_diff_l:
                # compute forward BWT to check standard words 
                command = "{exe} -b {input}".format(exe = str(cwd)+"/cais/cais", input= f_file)
                subprocess.check_call(command.split())
                nrf = noRuns(f_file+".bwt")
                std = (nrf == 2)
                write_csv_w_a(outfile4, [str(k),seq, str(diff), str(r), str(rf), str(rr), str(std)], 'a')
                if max_diff == diff:
                    write_csv_w_a(outfile4_2_1, [str(k),seq, str(diff), str(r), str(rf), str(rr), str(std)], 'a')

            # check min_diff
            if min_diff_l > diff:
                min_diff = diff
                # compute forward BWT to check standard words 
                command = "{exe} -b {input}".format(exe = str(cwd)+"/cais/cais", input= f_file)
                subprocess.check_call(command.split())
                nrf = noRuns(f_file+".bwt")
                std = (nrf == 2)
                write_csv_w_a(outfile4_2, ["n","s", "diff", "rho", "r_f", "r_r", "standard"], 'w')
                write_csv_w_a(outfile4_2, [str(k),seq, str(diff), str(r), str(rf), str(rr), str(std)], 'a')
                if min_diff > diff:
                    min_diff = diff
                    write_csv_w_a(outfile4_2_2, ["n","s", "diff", "rho", "r_f", "r_r", "standard"], 'w')
                    write_csv_w_a(outfile4_2_2, [str(k),seq, str(diff), str(r), str(rf), str(rr), str(std)], 'a')
            elif diff == min_diff_l:
                # compute forward BWT to check standard words 
                command = "{exe} -b {input}".format(exe = str(cwd)+"/cais/cais", input= f_file)
                subprocess.check_call(command.split())
                nrf = noRuns(f_file+".bwt")
                std = (nrf == 2)
                write_csv_w_a(outfile4_2, [str(k),seq, str(diff), str(r), str(rf), str(rr), str(std)], 'a')
                if diff == min_diff:
                    write_csv_w_a(outfile4_2_2, [str(k),seq, str(diff), str(r), str(rf), str(rr), str(std)], 'a')

            ## LYNDON FACTORS
            # compute diff lyndon factors
            lyn = cf_f - cf_r
            if lyn == 0 : no_lyn_0 += 1
            lyn_list.append(lyn)
            # check max_lyn
            if max_lyn_l < lyn:
                max_lyn_l = lyn
                write_csv_w_a(outfile6, ["n","s", "lyn", "diff", "rho", "lyn_f", "lyn_r"], 'w')
                write_csv_w_a(outfile6, [str(k),seq, str(lyn), str(diff), str(r), str(cf_f), str(cf_r), str(cf_r)], 'a')
                if max_lyn < lyn:
                    max_lyn = lyn
                    write_csv_w_a(outfile6_2_1, ["n","s", "lyn", "diff", "rho", "lyn_f", "lyn_r"], 'w')
                    write_csv_w_a(outfile6_2_1, [str(k),seq, str(lyn), str(diff), str(r), str(cf_f), str(cf_r), str(cf_r)], 'a')
           
            elif max_lyn_l == lyn:
                write_csv_w_a(outfile6, [str(k),seq, str(lyn), str(diff), str(r), str(cf_f), str(cf_r), str(cf_r)], 'a')
                if max_lyn == lyn:
                    write_csv_w_a(outfile6_2_1, [str(k),seq, str(lyn), str(diff), str(r), str(cf_f), str(cf_r), str(cf_r)], 'a')
            
            # check min_lyn
            if min_lyn_l > lyn:
                min_lyn_l = lyn
                write_csv_w_a(outfile6_2, ["n","s", "lyn", "diff", "rho", "lyn_f", "lyn_r"], 'w')
                write_csv_w_a(outfile6_2, [str(k),seq, str(lyn), str(diff), str(r), str(cf_f), str(cf_r), str(cf_r)], 'a')
                if min_lyn > lyn:
                    min_lyn = lyn
                    write_csv_w_a(outfile6_2_2, ["n","s", "lyn", "diff", "rho", "lyn_f", "lyn_r"], 'w')
                    write_csv_w_a(outfile6_2_2, [str(k),seq, str(lyn), str(diff), str(r), str(cf_f), str(cf_r), str(cf_r)], 'a')

            elif min_lyn_l == lyn:
                write_csv_w_a(outfile6_2, [str(k),seq, str(lyn), str(diff), str(r), str(cf_f), str(cf_r), str(cf_r)], 'a')
                if min_lyn == lyn:
                    write_csv_w_a(outfile6_2_2, [str(k),seq, str(lyn), str(diff), str(r), str(cf_f), str(cf_r), str(cf_r)], 'a')

            

    # compute stats rho
    mean_rho_l = round(sum(rho_list)/len(rho_list),3)
    min_rho_l = min(rho_list)
    std_rho_l = round(stdev(rho_list),3)
    p_rho_1_l = round((no_rho_1 / len(rho_list) * 100), 2)

    # compute stats diff
    mean_diff_l = round(sum(diff_list)/len(diff_list),3)
    std_diff_l = round(stdev(diff_list),3)

    # compute stats lyndon factors
    mean_lyn_l = round(sum(lyn_list)/len(lyn_list),3)
    std_lyn_l = round(stdev(lyn_list),3)
    p_lyn_0_l = round((no_lyn_0 / len(rho_list) * 100), 2)


    # write stats for all sequences of length k
    write_csv_w_a(outfile, [str(k),str(max_rho_l),str(min_rho_l),str(mean_rho_l),str(std_rho_l),str(p_rho_1_l), str(len(rho_list))],'a')
    write_csv_w_a(outfile3, [str(k),str(max_diff_l),str(min_diff_l),str(mean_diff_l),str(std_diff_l),str(p_rho_1_l)],'a')
    write_csv_w_a(outfile7, [str(k),str(max_lyn_l),str(min_lyn_l),str(mean_lyn_l),str(std_lyn_l),str(p_lyn_0_l)],'a')
    return max_rho, max_diff, min_diff, max_lyn, min_lyn

def allk_bbwt_bbwtrev(sigma, k, outdir,tempdir):
    b_rho, b_diff, c_f = {}, {}, {}
    # open output files and write headings
    filename = str(outdir) + "/stats_rho" + "_" + str(sigma) + ".csv" 
    write_csv_w_a(filename, ["n","max","min","mean","std","equals1"], 'w')
    
    # max rho out of all strings
    filename2_2 = str(outdir) + "/max_rho" + "_" + str(sigma) + ".csv" 
    write_csv_w_a(filename2_2, ["n","s","rho", "diff", "r_f", "r_r", "standard", "Lyn_r"], 'w')

    filename3 = str(outdir) + "/stats_diff" + "_" + str(sigma) + ".csv" 
    write_csv_w_a(filename3, ["n","max","min","mean","std","equals0"], 'w')

    # max diff out of all strings
    filename4_2_1 = str(outdir) + "/max_diff" + "_" + str(sigma) + ".csv" 
    write_csv_w_a(filename4_2_1, ["n","s","diff_runs", "rho", "r_f", "r_r", "standard"], 'w')

    # min diff out of all strings
    filename4_2_2 = str(outdir) + "/min_diff" + "_" + str(sigma) + ".csv" 
    write_csv_w_a(filename4_2_2, ["n","s","diff_runs", "rho", "r_f", "r_r", "standard"], 'w')

    
    filename5 = str(outdir) + "/complete_lyn" + "_" + str(sigma) + ".csv" 
    write_csv_w_a(filename5, ["n","s","diff_lyn", "diff_runs", "rho", "lyn_f", "lyn_r"], 'w')

    filename6_2_1 = str(outdir) + "/max_lyn" + "_" + str(sigma) + ".csv" 
    write_csv_w_a(filename6_2_1, ["n","s","diff_lyn", "diff_runs", "rho", "r_f", "r_r", "standard"], 'w')

    filename6_2_2 = str(outdir) + "/min_lyn" + "_" + str(sigma) + ".csv" 
    write_csv_w_a(filename6_2_2, ["n","s","diff_lyn", "diff_runs", "rho", "lyn_f", "lyn_r", "standard"], 'w')

    filename7 = str(outdir) + "/stats_lyn" + "_" + str(sigma) + ".csv" 
    write_csv_w_a(filename7, ["n","max","min","mean","std","equals0"], 'w')

    max_rho = 1
    max_diff = 0
    min_diff = 0
    max_lyn = 0
    min_lyn = 0
    if k < 3:
        print("max_k must be >= 3")
        return
    for i in range(3, k+1):
        print(i)
        # max RHO of each string length
        filename2 = str(outdir) + "/max_rho" + "_" + str(sigma)+"_" + str(i)+ ".csv" 
        write_csv_w_a(filename2, ["n","s","rho", "diff", "r_f", "r_r", "standard", "Lyn_r"], 'w')
        
        filename4 = str(outdir) + "/max_diff" + "_" + str(sigma) +"_" + str(i)+ ".csv" 
        write_csv_w_a(filename4, ["n","s","diff_runs", "rho", "r_f", "r_r", "standard"], 'w')
        filename4_2 = str(outdir) + "/min_diff" + "_" + str(sigma) +"_" + str(i)+ ".csv" 
        write_csv_w_a(filename4_2, ["n","s","diff_runs", "rho", "r_f", "r_r", "standard"], 'w')
        
        filename6 = str(outdir) + "/max_lyn" + "_" + str(sigma) + "_" + str(i)+".csv" 
        write_csv_w_a(filename6, ["n","s","diff_lyn", "diff_runs", "rho", "r_f", "r_r", "standard"], 'w')
        filename6_2 = str(outdir) + "/min_lyn" + "_" + str(sigma) + "_" + str(i)+".csv" 
        write_csv_w_a(filename6_2, ["n","s","diff_lyn", "diff_runs", "rho", "lyn_f", "lyn_r", "standard"], 'w')

        max_rho, max_diff, min_diff, max_lyn, min_lyn = bbwt_bbwtrev(sigma, i, filename, filename2, filename2_2, filename3, filename4, filename4_2, filename4_2_1, filename4_2_2, filename5, filename6, filename6_2, filename6_2_1, filename6_2_2, filename7, max_rho, max_diff, min_diff, max_lyn, min_lyn, tempdir)
    return

if __name__ == '__main__':
    main()
