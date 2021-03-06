import exceptions
import logging as log

import sys
sys.path.append("../../")
from plancknull.reader import DXReader
from plancknull import utils

log.root.level = log.DEBUG

read_map = DXReader(sys.argv[1], debug=True)
full_survs = ["nominal", "full"]
survs = full_survs + list(range(1, 5+1))

print "FULL FREQUENCIES"
for freq in [30, 44, 70]:
    for surv in survs:
        for halfring in [0, 1, 2]:
            for bp_corr in [False, True]:
                if (bp_corr and halfring!=0):
                    pass
                else:
                    try:
                        read_map(freq, surv, chtag='', halfring=halfring, pol="IQU", bp_corr=bp_corr)
                    except exceptions.Exception as e:
                        print e


print "QUADRUPLETS"
freq = 70
for chtag in ["18_23", "19_22", "20_21"]:
    for surv in ["nominal", "full"] + list(range(1, 5+1)):
        for halfring in [0, 1, 2]:
            if isinstance(surv, int) and halfring != 0:
                pass
            else:
                try:
                    read_map(freq, surv, chtag=chtag, halfring=halfring, pol="I")
                except exceptions.Exception as e:
                    print e


print "CHANNELS"
for freq in [30, 44, 70]:
    for surv in range(1, 5+1):
        for chtag in utils.chlist(freq):
            try:
                read_map(freq, surv, chtag=chtag, halfring=0, pol="I")
            except exceptions.Exception as e:
                print e

print "horns"
for freq in [30, 44, 70]:
    for surv in range(1, 5+1):
        for horn in utils.HORNS[freq]:
            try:
                read_map(freq, surv, chtag="LFI%d" % horn, halfring=0, pol="I")
            except exceptions.Exception as e:
                print e

#hfi_folder = os.environ["DX9_HFI"]
#
#read_map = SingleFolderDXReader(hfi_folder)
#for freq in [100, 143, 217, 353, 545, 857]:
#    for surv in survs:
#        for chtag in ['', "detset_1", "detset_2"]:
#            for halfring in [0, 1, 2]:
#                if (isinstance(surv, int) and (halfring != 0 or chtag)) or (halfring!=0 and chtag):
#                    pass
#                else:
#                    try:
#                        read_map(freq, surv, chtag=chtag,  halfring=halfring, pol="I")
#                    except exceptions.IOError as e:
#                        pass
#
#for freq in [100, 143, 217, 353, 545, 857]:
#    for surv in full_survs:
#        for ch in PLANCK.f[freq].ch:
#            if not ch.tag[-1] in "ab":
#                try:
#                    read_map(freq, surv, chtag=ch.tag,  halfring=0, pol="I")
#                except exceptions.IOError as e:
#                    pass
