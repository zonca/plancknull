import exceptions
import logging as log
import os
from plancknull.reader import SingleFolderDXReader

from planck.Planck import Planck

PLANCK = Planck()
log.root.level = log.DEBUG

lfi_folder = os.environ["DX9_LFI"]

read_map = SingleFolderDXReader(lfi_folder)
full_survs = ["nominal", "full"]
survs = full_survs + list(range(1, 5+1))
for freq in [30, 44, 70]:
    for surv in survs:
        for halfring in [0, 1, 2]:
            for bp_corr in [False, True]:
                if (isinstance(surv, int) and halfring != 0) or (bp_corr and halfring!=0):
                    pass
                else:
                    try:
                        read_map(freq, surv, nside=1024, chtag='', halfring=halfring, pol="I", bp_corr=bp_corr)
                    except exceptions.IOError as e:
                        pass

freq = 70
for chtag in ["18_23", "19_22", "20_21"]:
    for surv in ["nominal", "full"] + list(range(1, 5+1)):
        for halfring in [0, 1, 2]:
            if isinstance(surv, int) and halfring != 0:
                pass
            else:
                try:
                    read_map(freq, surv, chtag=chtag, nside=1024, halfring=halfring, pol="I")
                except exceptions.IOError as e:
                    pass


for freq in [30, 44, 70]:
    for surv in range(1, 5+1):
        for ch in PLANCK.f[freq].ch:
            try:
                read_map(freq, surv, chtag=ch.tag,  nside=None,     halfring=0, pol="I")
            except exceptions.IOError as e:
                pass
hfi_folder = os.environ["DX9_HFI"]

read_map = SingleFolderDXReader(hfi_folder)
for freq in [100, 143, 217, 353, 545, 857]:
    for surv in survs:
        for chtag in ['', "detset_1", "detset_2"]:
            for halfring in [0, 1, 2]:
                if (isinstance(surv, int) and (halfring != 0 or chtag)) or (halfring!=0 and chtag):
                    pass
                else:
                    try:
                        read_map(freq, surv, chtag=chtag,  nside=None, halfring=halfring, pol="I")
                    except exceptions.IOError as e:
                        pass

for freq in [100, 143, 217, 353, 545, 857]:
    for surv in full_survs:
        for ch in PLANCK.f[freq].ch:
            if not ch.tag[-1] in "ab":
                try:
                    read_map(freq, surv, chtag=ch.tag,  nside=None, halfring=0, pol="I")
                except exceptions.IOError as e:
                    pass
