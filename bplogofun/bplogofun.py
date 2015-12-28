from __future__ import print_function
from __future__ import division
from collections import defaultdict
from operator import itemgetter
from string import Template
from copy import deepcopy
import argparse
import re
import sys
import glob
import math as mt
import exact
import random
import statsmodels.api
import time
import pkgutil

def approx_expect2(H, k, N):
    return (H - ((k-1)/((mt.log(4)) * N)))

def permuted(items, pieces = 2):
    sublists = [[] for i in range(pieces)]
    for x in items:
        sublists[random.randint(0, pieces - 1)].append(x)
    permutedList = []
    for i in range(pieces):
        time.sleep(0.01)
        random.seed()
        random.shuffle(sublists[i])
        permutedList.extend(sublists[i])
    return permutedList

def weighted_dist(data):
    data_dict = defaultdict(int)
    for x in data:
        data_dict[x] += 1
    return data_dict

def rtp_binary(data, point, keys_sorted):
    low = 0
    high = len(keys_sorted)
    part = 0
    total = sum(data.values())
    while (low != high):
        mid = int((low + high) / 2)
        if (keys_sorted[mid] <= point):
            low = mid + 1
        else:
            high = mid
    for y in range(low, len(keys_sorted)):
        part += data[keys_sorted[y]]
    return part / total

def rtp(data, point):
    part = 0
    total = sum(data.values())
    results = filter(lambda x: x >= point, data.keys())
    for x in results:
        part += data[x]
    return part / total

def approx_expect(H, k, N):
    es = (k-1)/(2*N)
    es *= mt.log(mt.exp(1),2)
    return H - es

#logo_output(site_info, site_height_dict, adjusted_pvals[multipletesting[0]])
#    sprinzl = {'A': ["1:72", "2:71", "3:70", "4:69", "5:68", "6:67", "7:66"],
#               'D': ["10:25", "11:24", "12:23", "13:22"],
#               'C': ["27:43", "28:42", "29:41", "30:40", "31:39", "32:38"],
#               'T': ["49:65", "50:64", "51:63", "52:62", "53:61"]}
def logo_output(site_info, site_height_dict, pvals, p, P, pair_to_sprinzl, coord_to_pair, file_prefix):
    alpha = 0.05
    coord_length = 0 #used to determine eps height

    logo_outputDict = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
    for coord in range(1, len(site_info)+1):
        for base, info in sorted(site_info[coord].iteritems(), key = itemgetter(0)):
            for aainfo in sorted(site_height_dict[coord][base].iteritems(), key = itemgetter(1)):
                logo_outputDict[base][coord][aainfo[0]] = info * aainfo[1]
    
    for base in logo_outputDict:
        logodata = ""
        for coord in sorted(logo_outputDict[base].iterkeys()):
            if (len(str(coord)) > coord_length):
                coord_length = len(str(coord))
            logodata += "numbering {{({}) makenumber}} if\ngsave\n".format(coord)
            for aainfo in sorted(logo_outputDict[base][coord].iteritems(), key = itemgetter(1)):
                if (aainfo[1] == 0):
                    continue
                if P:
                    #pvalsP[bp][pairtype][aa_class[0]] = pv
                    if str(coord) in coord_to_pair: 
                        bp = coord_to_pair[str(coord)]
                        bp_list = bp.split(":")
                        if (bp_list.index(str(coord)) == 0):
                            sig = False
                            for pairs in filter(lambda k: k.startswith(base), pvals['P'][pair_to_sprinzl[bp]].keys()):
                                if (aainfo[0].upper() in pvals['P'][pair_to_sprinzl[bp]][pairs] and pvals['P'][pair_to_sprinzl[bp]][pairs][aainfo[0]] <= alpha):
                                    logodata += "/showingbox (s) def\n"
                                    logodata += "{:07.5f} ({}) numchar\n".format(aainfo[1], aainfo[0].upper())
                                    logodata += "/showingbox (n) def\n"
                                    sig = True
                                    break
                            if not sig:
                                logodata += "{:07.5f} ({}) numchar\n".format(aainfo[1], aainfo[0].upper())
                        else:
                            sig = False
                            for pairs in filter(lambda k: k.endswith(base), pvals['P'][pair_to_sprinzl[bp]].keys()):
                                if (aainfo[0].upper() in pvals['P'][pair_to_sprinzl[bp]][pairs] and pvals['P'][pair_to_sprinzl[bp]][pairs][aainfo[0]] <= alpha):
                                    logodata += "/showingbox (s) def\n"
                                    logodata += "{:07.5f} ({}) numchar\n".format(aainfo[1], aainfo[0].upper())
                                    logodata += "/showingbox (n) def\n"
                                    sig = True
                                    break
                            if not sig:
                                logodata += "{:07.5f} ({}) numchar\n".format(aainfo[1], aainfo[0].upper())

                    else:
                        logodata += "{:07.5f} ({}) numchar\n".format(aainfo[1], aainfo[0].upper())
                else:
                    logodata += "{:07.5f} ({}) numchar\n".format(aainfo[1], aainfo[0].upper())
            logodata += "grestore\nshift\n"
        #output logodata to template
        template_byte = pkgutil.get_data('bplogofun', 'eps/Template.eps')
        logo_template = template_byte.encode('utf-8')
        with open("{}_{}.eps".format(base, file_prefix), "w") as logo_output: 
            src = Template(logo_template)
            logodata_dict = {'logo_data': logodata, 'low': min(logo_outputDict[base].keys()), 'high': max(logo_outputDict[base].keys()), 'length': 15.28 * max(logo_outputDict[base].keys()), 'height': 735-(5*coord_length)}
            logo_output.write(src.substitute(logodata_dict))

def main():
    #Setup parser
    parser = argparse.ArgumentParser(description = "bpLogoFun")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-i", "--infernal", action="store_true")
    group.add_argument("-c", "--cove", action="store_true")
    parser.add_argument('--max', '-x', help='max exact entropy', type=int, default=5)
    parser.add_argument('--logo', help='Produce function logo ps files', action="store_true")
    parser.add_argument("-s", "--single", action="store_true")
    parser.add_argument("-p", help="Calculate permuation-based p-values for total information of CIFs",
                        action="store_true")
    parser.add_argument("-P",
                        help="calculate permutation-based p-values of individual CIF-class associations",
                        action="store_true")
    parser.add_argument("-B", help="Number of permutations", type=int, default=100)
    parser.add_argument("-o", "--stdout", action="store_true", help="Print results to STDOUT")
    parser.add_argument("-M", default="BY", help = "Specify method to correct p-values for multiple-comparisons. Current methods available: bonferroni, holm, hommel, BH, BY, and hochberg-simes. One or more can be specified separated by colons(:).")
    parser.add_argument("-d", action="store_true", help="Output the alignment coordinates that correspond to each base-pair")
    parser.add_argument("struct", help="Structure File")
    parser.add_argument("file_preifx", help="File prefix")
    args = parser.parse_args()
    permute = False
    multipletesting = []
    if (args.M):
        mp_testdict = {'BH': 'fdr_bh', 'BY': 'fdr_by', 'bonferroni': 'bonferroni',
                       'holm': 'holm', 'hommel': 'hommel',
                       'simes-hochberg': 'simes-hochberg'}
        multipletesting_list = args.M.split(":")
        for x in multipletesting_list:
            multipletesting.append(mp_testdict[x])

    #Dictionary of sprinzl coords
    sprinzl = {'A': ["1:72", "2:71", "3:70", "4:69", "5:68", "6:67", "7:66"],
               'D': ["10:25", "11:24", "12:23", "13:22"],
               'C': ["27:43", "28:42", "29:41", "30:40", "31:39", "32:38"],
               'T': ["49:65", "50:64", "51:63", "52:62", "53:61"]}
    print("# max_exact = {}".format(args.max), file=sys.stderr)
    
    if (args.p or args.P):
        permute = True
        pvalsP = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        pvalsp = defaultdict(lambda: defaultdict(float))
        if (len(multipletesting) == 0):
            multipletesting.append('fdr_bh')

        num_permutations = args.B
        print("# Number of permutations for p-values: {}".format(num_permutations), file=sys.stderr)


    print("Parsing base-pair coordinates", file=sys.stderr)
    #Parse Base-Pair Data
    ss = ""
    pairs = defaultdict(list)
    tarm = 0
    stack = []
    if args.infernal:
        with open(args.struct, "r") as infernal:
            for line in infernal:
                line = line.strip()
                ss += line.split()[2]

        state = "start"
        for count, i in enumerate(ss):
            if (i == "("):
                stack.append(count)
                if (state == "start"):
                    state = "A"
            elif (i == "<"):
                stack.append(count)
                if (state == "A"):
                    state = "D"
                elif (state == "cD"):
                    state = "C"
                elif (state == "cC"):
                    state = "T"
            elif (i == ">"):
                if (state == "D"):
                    state = "cD"
                elif (state == "C"):
                    state = "cC"
                elif (state == "T"):
                    state = "cT"

                arm = state.replace("c", "")
                pairs[arm].append([stack.pop(), count])
            elif (i == ")"):
                pairs['A'].append(stack.pop(), count)

    if args.cove:
        with open(args.struct, "r") as cove:
            for line in cove:
                line = line.strip()
                ss += line.split()[1]
        
        state = "start"
        for count, i in enumerate(ss):
            if (i == ">" and (state == "start" or state == "AD")):
                if (state == "start"):
                    state = "AD"
                stack.append(count)
    
            elif (i == "<" and (state == "AD" or state == "D")):
                if (state == "AD"):
                    state = "D"
                pairs[state].append([stack.pop(), count])
    
            elif (i == ">" and (state == "D" or state == "C")):
                if (state == "D"):
                    state = "C"
                stack.append(count)
    
            elif (i == "<" and (state == "C" or state == "cC")):
                if (state == "C"):
                    state = "cC"
                pairs["C"].append([stack.pop(), count])
    
            elif (i == ">" and (state == "cC" or state == "T")):
                if (state == "cC"):
                    state = "T"
                stack.append(count)
                tarm += 1
    
            elif (i == "<" and (state == "T" and tarm > 0)):
                pairs[state].append([stack.pop(), count])
                tarm -= 1
    
            elif (i == "<" and (state == "T" or state == "A") and tarm == 0):
                state = "A"
                pairs[state].append([stack.pop(), count])
    
    #reverse to start from 5' end
    for arm in pairs:
        pairs[arm].reverse()
    
    numclasses = 0
    size_dict = {} 
    size_int = 0
    summ = defaultdict(lambda : defaultdict(int))
    freq = defaultdict(lambda : defaultdict(lambda : defaultdict(int)))
    sitesum = defaultdict(lambda : defaultdict(int))
    sitefreq = defaultdict(lambda : defaultdict(lambda : defaultdict(int)))
    aa_classes = []
    seqs = []

    for fn in glob.glob("{}_?.aln".format(args.file_preifx)):
        match = re.search("_([A-Z])\.aln", fn)
        aa_class = match.group(1)
        numclasses += 1
        size_dict[aa_class] = 0
        with open(fn, "r") as ALN:
            good = False
            interleaved = False
            seq = {}
            for line in ALN:
                line = line.strip()
                match = re.search("^(\S+)\s+(\S+)", line)
                if (re.search("^CLUSTAL", line)):
                    good = True
                    continue
                elif (re.search("^[\s\*\:\.]+$", line) and not interleaved):
                    #print("Alignment data appears to be interleaved.", file = sys.stderr)
                    interleaved = True
                elif (re.search("^[\s\*\:\.]+$", line) and interleaved):
                    continue
                
                elif (match and not interleaved):
                    if (not good):
                        sys.exit("File {} appears not to be a clustal file".format(fn))
                    seq[match.group(1)] = match.group(2)
                    size_dict[aa_class] += 1
                    size_int += 1
                elif (match and interleaved):
                    seq[match.group(1)] += match.group(2)
        length = 0
        for sequence in seq.itervalues():
            if (permute):
                aa_classes.append(aa_class)
                seqs.append(sequence)

            if (not length > 0 or args.single):
                length = len(sequence)
            for arm in pairs:
                for i, pair in enumerate(pairs[arm]):
                    bp = sprinzl[arm][i]
                    pairtype = "{}{}".format(sequence[pair[0]].upper(),sequence[pair[1]].upper())
                    pairtype = pairtype.replace("T", "U")
                    summ[bp][pairtype] += 1
                    freq[bp][pairtype][aa_class] += 1
    
            if (args.single):
                for pos, char in enumerate(sequence):
                    char = char.upper()
                    state = char.replace("T", "U")
                    sitesum[pos+1][state] += 1
                    sitefreq[pos+1][state][aa_class] += 1
    
    print("{} alignments parsed".format(numclasses), file = sys.stderr)
    
    if (permute):
        indices = []
        pfreq = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda: defaultdict(int))))
        psitefreq = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda: defaultdict(int))))
        print("Generating permuted alignment data", file=sys.stderr)
        for p in range(num_permutations):
            indices.append(permuted(range(len(aa_classes))))
        for s, seq in enumerate(seqs):
            for arm in pairs:
                for j, coords in enumerate(pairs[arm]):
                    bp = sprinzl[arm][j]
                    pairtype = "{}{}".format(seq[coords[0]].upper(),seq[coords[1]].upper())
                    pairtype = pairtype.replace("T", "U")
                    for p in range(num_permutations):
                        pfreq[p][bp][pairtype][aa_classes[indices[p][s]]] += 1
            if (args.single):
                for pos, char in enumerate(seq):
                    char = char.upper()
                    state = char.replace("T", "U")
                    for p in range(num_permutations):
                        psitefreq[p][pos+1][state][aa_classes[indices[p][s]]] += 1         

        
    print("Computing exact expected entropies", file = sys.stderr)
    
    p = [x/sum(size_dict.values()) for x in size_dict.values()] # Holds Background
    exact_list = []
    start_sample_sz =1
    for n in range(start_sample_sz, args.max + 1):
        j = exact.calc_exact(n, p, numclasses)
        exact_list.append(j[1])
        print("{:2} {:07.5f}".format(n, exact_list[n-1]), file=sys.stderr)
    
    
    bg_entropy = 0
    info = defaultdict(lambda : defaultdict(float))
    height_dict = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
    for x in size_dict.itervalues():
        if (x != 0):
            bg_entropy -= (x/size_int) * mt.log(x/size_int, 2)
    
    if (permute):
        print("Computing permutation distributions ...", file = sys.stderr)
        bpinfodata = []
        bpheightdata = []
        for arm in ["A", "D", "C", "T"]:
            for i, bp in enumerate(sprinzl[arm]):
                if (not i < len(pairs[arm])):
                    continue
                for pairtype in sorted(summ[bp].keys()):
                    total = summ[bp][pairtype]
                    for p in range(num_permutations):
                        numpositives = len(pfreq[p][bp][pairtype])
                        fg_entropy = 0
                        for x in pfreq[p][bp][pairtype].itervalues():
                            fg_entropy -= (x/total)*mt.log(x/total,2)
                        if (total <= args.max):
                            expected_bg_entropy = exact_list[total - 1]
                        else:
                            expected_bg_entropy = approx_expect(bg_entropy, numclasses, total)

                        if ((expected_bg_entropy - fg_entropy) < 0):
                            info_int = 0.0
                        else:
                            info_int = expected_bg_entropy - fg_entropy
                        
                        if (args.p):
                            bpinfodata.append(info_int)

                        if (args.P):
                            pheightclass = {}
                            pheight = 0
                            for aa_class in pfreq[p][bp][pairtype].iterkeys():
                                pheightclass[aa_class] = ((pfreq[p][bp][pairtype][aa_class] / summ[bp][pairtype]) / (size_dict[aa_class] / size_int))
                                pheight += pheightclass[aa_class]
                            
                            for aa_class in pheightclass:
                                bpheightdata.append((pheightclass[aa_class] / pheight) * info_int)


        if (args.p):
            bpinfodist = weighted_dist(bpinfodata)
            bpinfodist_sortedKeys = sorted(bpinfodist.keys())
        if (args.P):
            bpheightdist = weighted_dist(bpheightdata)
            bpheightdist_sortedKeys = sorted(bpheightdist.keys())

        if (args.single):
            siteinfodata = []
            siteheightdata = []
            for i in sorted(sitesum.keys()):
                for state in sorted(sitesum[i].keys()):
                    total = sitesum[i][state]
                    for p in range(num_permutations):
                        numpositives = len(psitefreq[p][i][state])
                        fg_entropy = 0
                        for x in psitefreq[p][i][state].itervalues():
                            fg_entropy -= x/total * mt.log(x/total, 2)
                        if (total <= args.max):
                            expected_bg_entropy = exact_list[total - 1]
                        else:
                            expected_bg_entropy = approx_expect(bg_entropy, numclasses, total)
                        
                        if ((expected_bg_entropy - fg_entropy) < 0):
                            info_int = 0
                        else:
                            info_int = expected_bg_entropy - fg_entropy

                        if (args.p):
                            siteinfodata.append(info_int)

                        if (args.P):
                            pheightclass = {}
                            pheight = 0
                            for aa_class in psitefreq[p][i][state]:
                                pheightclass[aa_class] = ((psitefreq[p][i][state][aa_class] / sitesum[i][state]) / (size_dict[aa_class] / size_int))
                                pheight += pheightclass[aa_class]
                            for aa_class in pheightclass:
                                siteheightdata.append(((pheightclass[aa_class] / pheight) * info_int))
            if (args.p):
                siteinfodist = weighted_dist(siteinfodata)
                siteinfodist_sortedKeys = sorted(siteinfodist.keys())
            if (args.P):
                siteheightdist = weighted_dist(siteheightdata)
                siteheightdist_sortedKeys = sorted(siteheightdist.keys())
    
    #Compute information stats
    print("Computing information statistics", file = sys.stderr)
    pvals = []
    for arm in ["A", "D", "C", "T"]:
        for i, bp in enumerate(sprinzl[arm]):
            if (not i < len(pairs[arm])):
                continue
            for pairtype in sorted(summ[bp].keys()):
                total = summ[bp][pairtype]
                numpositives = len(freq[bp][pairtype])
                fg_entropy = 0
                for x in freq[bp][pairtype].itervalues():
                    fg_entropy -= (x/total) * mt.log(x/total, 2)
                if (total <= args.max):
                    expected_bg_entropy = exact_list[total - 1]
                else:
                    expected_bg_entropy = approx_expect(bg_entropy, numclasses, total)
                
                if ((expected_bg_entropy - fg_entropy) < 0):
                    info[bp][pairtype] = 0
                else:
                    info[bp][pairtype] = expected_bg_entropy - fg_entropy
                
                if (args.p):
                    #def rtp_binary(data, point, keys_sorted):
                    pv = rtp_binary(bpinfodist, info[bp][pairtype], bpinfodist_sortedKeys)
                    pvalsp[bp][pairtype] = pv

                height = 0
                for aa_class in freq[bp][pairtype].iterkeys():
                    height_dict[bp][pairtype][aa_class] = ((freq[bp][pairtype][aa_class] / summ[bp][pairtype]) / (size_dict[aa_class] / size_int))
    
                    height += height_dict[bp][pairtype][aa_class]
                for aa_class in sorted(height_dict[bp][pairtype].iteritems(), key = itemgetter(1), reverse = True):
                    height_dict[bp][pairtype][aa_class[0]] /= height
    
                    if (args.P):
                        pv = rtp_binary(bpheightdist, height_dict[bp][pairtype][aa_class[0]] * info[bp][pairtype], bpheightdist_sortedKeys)
                        pvalsP[bp][pairtype][aa_class[0]] = pv

    if (args.single):
        site_info = defaultdict(lambda : defaultdict(float))
        site_height_dict = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
                   # sitesum[pos+1][state] += 1
                    #sitefreq[pos+1][state][aa_class] += 1
        for i in sorted(sitesum.iterkeys()):
            for state in sorted(sitesum[i]):
                total = sitesum[i][state]
                numpositives = len(sitefreq[i][state])
                fg_entropy = 0
                fg_entropy = -(sum(map(lambda x: (x/total) * mt.log(x/total,2), sitefreq[i][state].itervalues())))
                if (total <= args.max):
                    expected_bg_entropy = exact_list[total-1]
                else:
                    expected_bg_entropy = approx_expect(bg_entropy, numclasses, total)
                if ((expected_bg_entropy - fg_entropy) < 0):
                    site_info[i][state] = 0.0
                else:
                    site_info[i][state] = expected_bg_entropy - fg_entropy

                if (args.p):
                    #def rtp_binary(data, point, keys_sorted):
                    pv = rtp_binary(siteinfodist, site_info[i][state], siteinfodist_sortedKeys)
                    pvalsp[i][state] = pv

                height = 0
                for aa_class in sitefreq[i][state].iterkeys():
                    site_height_dict[i][state][aa_class] = ((sitefreq[i][state][aa_class] / sitesum[i][state]) / (size_dict[aa_class] / size_int))
                    height += site_height_dict[i][state][aa_class]

                for aa_class in sorted(site_height_dict[i][state].iteritems(), key = itemgetter(1), reverse = True):
                    site_height_dict[i][state][aa_class[0]] /= height
                    if (args.P):
                        pv = rtp_binary(siteheightdist, site_height_dict[i][state][aa_class[0]] * site_info[i][state], siteheightdist_sortedKeys)
                        pvalsP[i][state][aa_class[0]] = pv

    if (permute):
        pvalsss = []
        adjusted_pvals = defaultdict(lambda: defaultdict(dict))
        print("Adjusting for multiple comparisons.", file=sys.stderr)
        if (args.p and not arg.P):
            for x in sorted(pvalsp.keys()):
                for y in sorted(pvalsp[x].items(), key=itemgetter(0)):
                    pvalsss.append(y[1])
        elif (args.P and args.p):
            #add args.p first
            for x in sorted(pvalsp.keys()):
                for y in sorted(pvalsp[x].items(), key=itemgetter(0)):
                    pvalsss.append(y[1])
            #add args.P
            for x in sorted(pvalsP.keys()):
                for y in sorted(pvalsP[x].keys()):
                    for z in sorted(pvalsP[x][y].items(), key=itemgetter(0)):
                        pvalsss.append(z[1])
        elif (args.P and not args.p):
            for x in sorted(pvalsP.keys()):
                for y in sorted(pvalsP[x].keys()):
                    for z in sorted(pvalsP[x][y].items(), key=itemgetter(0)):
                        pvalsss.append(z[1])
        
        for mtest in multipletesting:
            correctedp = statsmodels.api.stats.multipletests(pvalsss, method=mtest)
            correctedp = list(correctedp[1])
            if (args.p and not args.P):
                for x in sorted(pvalsp.keys()):
                    for y in sorted(pvalsp[x].keys()):
                        pvalsp[x][y] = correctedp.pop(0)
                adjusted_pvals[mtest]['p'] = deepcopy(pvalsp)
            elif (args.P and args.p):
                #process args.p first
                for x in sorted(pvalsp.keys()):
                    for y in sorted(pvalsp[x].keys()):
                        pvalsp[x][y] = correctedp.pop(0)
                adjusted_pvals[mtest]['p'] = deepcopy(pvalsp)
                #process args.P
                for x in sorted(pvalsP.keys()):
                    for y in sorted(pvalsP[x].keys()):
                        for z in sorted(pvalsP[x][y].keys()):
                            pvalsP[x][y][z] = correctedp.pop(0)
                adjusted_pvals[mtest]['P'] = deepcopy(pvalsP)
            elif (args.P and not args.p):
                for x in sorted(pvalsP.keys()):
                    for y in sorted(pvalsP[x].keys()):
                        for z in sorted(pvalsP[x][y].keys()):
                            pvalsP[x][y][z] = correctedp.pop(0)
                adjusted_pvals[mtest]['P'] = deepcopy(pvalsP)

    if (args.stdout):
        #restoring original p-values for output
        if (args.p and not args.P):
            for x in sorted(pvalsp.keys()):
                for y in sorted(pvalsp[x].keys()):
                    pvalsp[x][y] = pvalsss.pop(0)
        elif (args.P and args.p):
            #process args.p first
            for x in sorted(pvalsp.keys()):
                for y in sorted(pvalsp[x].keys()):
                    pvalsp[x][y] = pvalsss.pop(0)
            #process args.P
            for x in sorted(pvalsP.keys()):
                for y in sorted(pvalsP[x].keys()):
                    for z in sorted(pvalsP[x][y].keys()):
                        pvalsP[x][y][z] = pvalsss.pop(0)
        elif (args.P and not args.p):
            for x in sorted(pvalsP.keys()):
                for y in sorted(pvalsP[x].keys()):
                    for z in sorted(pvalsP[x][y].keys()):
                        pvalsP[x][y][z] = pvalsss.pop(0)
        
        #build output heading
        heading_dict = {}
        if (args.p):
            pstring = "\tp-value"
            for m in multipletesting:
                pstring += "\t{}".format(m)
            heading_dict['p'] = pstring
        else:
            heading_dict['p'] = ""
        if (args.P):
            Pstring = "\tclass:height:p-value"
            for m in multipletesting:
                Pstring += ":{}".format(m)
            heading_dict['P'] = Pstring
        else:
            heading_dict['P'] = "\tclass:height"
        if (args.d):
            heading_dict['d'] = "\tcoord"
        else:
            heading_dict['d'] = ""

        print("#bp\tarm\tsprinzl{d}\tbp\tN\tinfo{p}{P}".format(**heading_dict))
        for arm in ["A", "D", "C", "T"]:
            for i, coord in enumerate(sprinzl[arm]):
                if coord in info:
                    for pairtype in sorted(info[coord].iterkeys()):
                        output_string = "bp:\t{}\t{}".format(arm, coord)
                        if (args.d):
                            output_string += "\t{}".format(":".join([str(pairs[arm][i][0]+1),str(pairs[arm][i][1]+1)]))
                        output_string += "\t{}\t{}\t{:05.3f}\t".format(pairtype, summ[coord][pairtype],
                                                               info[coord][pairtype])
                        if (args.p):
                            output_string += "{:08.6f}".format(pvalsp[coord][pairtype])
                            for x in multipletesting:
                                output_string += "\t{:08.6f}".format(adjusted_pvals[x]['p'][coord][pairtype])
                                
                        output_string += "\t"
                        for aainfo in sorted(height_dict[coord][pairtype].iteritems(), key = itemgetter(1), reverse = True):
                            output_string += " {}:{:05.3f}".format(aainfo[0], aainfo[1])
                            if (args.P):
                                output_string += ":{:08.6f}".format(pvalsP[coord][pairtype][aainfo[0].upper()])
                                for x in multipletesting:
                                    output_string += ":{:08.6f}".format(adjusted_pvals[x]['P'][coord][pairtype][aainfo[0].upper()])

                        print(output_string)
        if (args.single):
            print("#ss\t\tcoord\tf\tN\tinfo{p}{P}".format(**heading_dict))
            for coord in range(1, len(site_info)+1):
                for base in sorted(site_info[coord].iterkeys()):
                    output_string = "ss:\t\t{}\t{}\t{}\t{:05.3f}".format(coord, base,
                                                                         sitesum[coord][base],
                                                                         site_info[coord][base])
                    if (args.p):
                        output_string += "{:08.6f}".format(pvalsp[coord][base])
                        for x in multipletesting:
                            output_string += "\t{}".format(adjusted_pvals[x]['p'][coord][base])
                    
                    output_string += "\t"
                    for aainfo in sorted(site_height_dict[coord][base].iteritems(), key = itemgetter(1), reverse = True):
                        output_string += " {}:{:05.3f}".format(aainfo[0], aainfo[1])
                        if (args.P):
                            output_string += ":{:08.6f}".format(pvalsP[coord][base][aainfo[0].upper()])
                            for x in multipletesting:
                                output_string += ":{:08.6f}".format(adjusted_pvals[x]['P'][coord][base][aainfo[0].upper()])

                    print(output_string)

    if (args.logo):
        print("Producing logo graphics")
        #def logo_output(site_info, site_height_dict, pvals = {}, p = False, P = False):
        pair_to_sprinzl = {}
        coord_to_pair = {}
        for arm in ["A", "D", "C", "T"]:
            for i, coord in enumerate(sprinzl[arm]):
                if (not i < len(pairs[arm])):
                    continue
                pair_to_sprinzl[":".join([str(pairs[arm][i][0] + 1), str(pairs[arm][i][1] + 1)])] = coord


        for key in pair_to_sprinzl:
            key_split = key.split(":")
            for x in key_split:
                coord_to_pair[x] = key

        if (permute):
            logo_output(site_info, site_height_dict, adjusted_pvals[multipletesting[0]], args.p, args.P,
                        pair_to_sprinzl, coord_to_pair, args.file_preifx)
        else:
            logo_output(site_info, site_height_dict, {}, args.p, args.P, pair_to_sprinzl, coord_to_pair, args.file_preifx)
