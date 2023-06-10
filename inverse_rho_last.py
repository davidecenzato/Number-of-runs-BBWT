import argparse, string, csv, heapq, subprocess
from itertools import product
from pathlib import Path
from collections import Counter
from statistics import stdev

from count_runs import noRuns

cwd = Path(__file__).parent.absolute()

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("sigma", help="size of the alphabet", default=2, type=int)
    parser.add_argument("max_k", help="maximum strings length", default=10, type=int)
    parser.add_argument("-o", "--outdir", help="existing output directory for generated files",
                        default=cwd.joinpath("outdir"), type=Path)
    parser.add_argument("-cr", "--countrho", choices=[0, 1], help="select if you want to print count of rho values",
                        default=0, type=int)
    parser.add_argument("-brho", "--bigrho",#, choices=[0, 1],
                        help="select if you want to output the strings with 10th higher rho values",
                        default=0, type=int)
    parser.add_argument("-bdiff", "--bigdiff",
                        help="select (-bdiff x) if you want to output the strings with xth higher and xth smaller diffrencee in the number of runs between bbwt(for) and bbwt(rev) ",
                        default=0, type=int)
    parser.add_argument("-c", "--cfactors", choices=[0, 1],
                       help="select if you want to output the number of lyndon factors for the forward and the reverse ",
                       default=0, type=int)
    parser.add_argument("--rep", choices=[0,1], help="disciminate repetitive and non repetitive sequences",
                        default=0, type=int)
    args = parser.parse_args()
    sigma = args.sigma
    outdir = Path(args.outdir).joinpath("sigma_" + str(sigma))
    if not Path(outdir).exists():
        outdir.mkdir()
        print("A new output directory has been created.")
    allk_bbwt_bbwtrev(sigma, args.max_k, outdir, args.countrho, args.bigrho, args.bigdiff, args.cfactors)
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

def write_csv(filename, rows):
    with open(filename, 'w', newline='') as f:
        wr = csv.writer(f, quoting=csv.QUOTE_ALL)
        for r in rows:
            if r: wr.writerow(r)
    return

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

def bbwt_bbwtrev(sigma,k, outfile, outfile2, count_rho=False, brho=False, bdiff=False, cfactors=False):
    d = []
    rho, diff , cf , repseq = {}, {}, {}, {}
    # create strings |k| with alphabet |sigma|
    it = product(string.ascii_uppercase[:sigma], repeat=k)
    no_seq = 0
    rho_list = []
    no_rho_1 = 0
    for s in it:
        seq = ''.join(s)
        if check_rev(seq):
            # write forward to file
            with open("forward.txt","w") as file:
                file.write(seq)
            # compute bbwt forward
            command = "{exe} -t {input}".format(exe = str(cwd)+"/cais/cais", input= "forward.txt")
            subprocess.check_call(command.split())
            # increase no. sequences
            no_seq += 1
            # compute forward Lyndon factorization length
            Lfact = Duval(seq)
            cf_f = len(Lfact)
            # compute BBWT forward number of runs
            rf = noRuns("forward.txt.bbwt")
            # compute reverse sequence
            rev = seq[::-1]
            # write reverse to file
            with open("reverse.txt","w") as file:
                file.write(rev)
            # compute bbwt reverse
            command = "{exe} -t {input}".format(exe = str(cwd)+"/cais/cais", input= "reverse.txt")
            subprocess.check_call(command.split())
            # compute reverse Lyndon factorization length
            Lfact = Duval(rev)
            # compute BBWT reverse number of runs
            rr = noRuns("reverse.txt.bbwt")
            # compute rho
            r = round(max((rf / rr), (rr / rf)), 2)
            # chekf if rho == 1
            if r==1: no_rho_1 += 1
            # append rho to rho list
            rho_list.append(r)
            # compute forward BWT
            command = "{exe} -b {input}".format(exe = str(cwd)+"/cais/cais", input= "forward.txt")
            subprocess.check_call(command.split())
            nrf = noRuns("forward.txt.bwt")
            '''
            bwt_for = BWT(seq)
            nrf = la.count_run(bwt_for)
            bwt_rev = BWT(rev)
            nrr = la.count_run(bwt_rev)
            '''
            '''
            dif = rf - rr
            cf_r = len(Lfact)
            cf_d = cf_f - cf_r
            rep = round(k/rf, 2)
            d.append([rf,
                      rr,
                      r,
                      dif,
                      cf_f,
                      lf_r,
                      cf_r,
                      cf_d,
                      rep
                      # cf_r := round(len(Duval(s))/len(Duval(rev)), 2),
                      ])
            print(seq)
            rf = la.count_run(la.bbwt(seq, False, False))
            print(rf)
            rev = seq[::-1]
            rr = la.count_run(la.bbwt(rev, False, False))
            r = round(max((rf / rr), (rr / rf)), 2)
            '''
            '''
            if r >= brho:
                
                cf_f = len(Duval(s))
                lf_r = Duval(rev)
                dif = rf - rr
                cf_d = cf_f - len(lf_r)
                
                rho[seq] = [k,r,
                      dif,
                      rf, rr, #la.count_run(la.bwt(seq, False, False))==2,
                      cf_d,
                      cf_f, cf_r, lf_r]
                
                rho[seq] = r
            '''

            # write single sequence stats
            write_csv_w_a(outfile2, [str(k),seq,str(r)], 'a')

    # compute stats
    mean_rho = sum(rho_list)/no_seq
    max_rho = max(rho_list)
    min_rho = min(rho_list)
    std_rho = stdev(rho_list)

    # write stats for all sequences of length k
    write_csv_w_a(outfile, [str(k),str(max_rho),str(min_rho),str(mean_rho),str(std_rho),str(no_rho_1)],'a')

                #if not count_rho: rho={k:v for k,v in heapq.nlargest(15, rho.items(), key=lambda i: i[1])}
            # if count_rho:
            #     rho[seq] = r
            # if bdiff:
            #     if dif >= bdiff or dif <= -(bdiff):
            #         diff[seq] = [dif, cf_d]
            #         #diff_max = {k: v for k, v in heapq.nlargest(bdiff, diff.items(), key=lambda i: i[1])}
                    #diff = {k: v for k, v in heapq.nsmallest(bdiff, diff.items(), key=lambda i: i[1])}
                    #diff.update(diff_max)
            # if True:
            #     repseq[seq]=rep
            #     repseq = {k: v for k, v in heapq.nlargest(10, repseq.items(), key=lambda i: i[1])}
            # if cfactors: cf[seq]=cf_d
    #print('k = ', k, " \par", sorted(Counter([s[3] for s in d]).items()))
    #if count_rho: print('k=',k, sorted(Counter([s[2] for s in d]).items()))
    #if cfactors: print('k = ', k, " \par", sorted(Counter([s[4] for s in d]).items()))
    #filename = str(outdir) + "/allruns" + str(k) + "_" + str(sigma) + ".csv"
    #write_csv(filename, d)
    return rho, diff, cf
'''
def min_max_diff(sigma,k, bdiff, b_diff, outdir):
    max_diff = sorted(b_diff.items(), key=lambda i: i[1], reverse=True)
    #max_diff = heapq.nlargest(bdiff, b_diff.items(), key=lambda i: i[1])
    filename = str(outdir) + "/maxMindf" + str(k) + "_" + str(sigma) + ".csv"
    #filename = str(outdir) + "/maxdiff_" + str(sigma) + ".csv"
    write_csv(filename, max_diff)
    #min_diff = sorted(b_diff.items(), key=lambda i: i[1], reverse=False)
    #min_diff = heapq.nsmallest(bdiff, b_diff.items(), key=lambda i: i[1])
    #filename = str(outdir) + "/mindf_" + str(sigma) + ".csv"
    #filename = str(outdir) + "/mindiff_" + str(sigma) + ".csv"
    #write_csv(filename, min_diff)
    return

def write_csv2(sigma, dict , fname, outdir):
    filename = str(outdir) + "/" + fname + str(sigma) + ".csv"
    write_csv(filename, dict)
    return

def big_rho(sigma,k, b_rho, outdir):
    max_rho = sorted(b_rho.items(), key=lambda i: i[1], reverse=True)
    filename = str(outdir) + "/maxrho" + str(k) + "_" + str(sigma) + ".csv"
    write_csv(filename, max_rho)
    return
'''
def allk_bbwt_bbwtrev(sigma, k, outdir, count_rho=False, brho=False, bdiff=False, cfactors=False):
    b_rho, b_diff, c_f = {}, {}, {}
    # open output files and write headings
    filename = str(outdir) + "/stats" + "_" + str(sigma) + ".csv"
    write_csv_w_a(filename, ["n","max","min","mean","std","equals1"], 'w')
    filename2 = str(outdir) + "/complete" + "_" + str(sigma) + ".csv"
    write_csv_w_a(filename2, ["n","s","rho"], 'w')

    if k < 3:
        print("max_k must be >= 3")
        return
    for i in range(3, k + 1):
        rho, diff, cf = bbwt_bbwtrev(sigma, i, filename, filename2, count_rho, brho, bdiff, cfactors)
        if brho or count_rho: b_rho.update(rho)
        if bdiff: b_diff.update(diff)
        if cfactors : c_f.update(cf)
    if brho: big_rho(sigma,k, b_rho, outdir)
    if bdiff: min_max_diff(sigma,k, bdiff, b_diff, outdir)
    if count_rho: print(sorted(Counter(b_rho.values()).items()))
    if cfactors:
        print(sorted(Counter(c_f.values()).items()))
        fname="dfactors_"
        min_max_diff(sigma, k, 50, c_f, outdir)
    return

if __name__ == '__main__':
    main()