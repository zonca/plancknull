import os

import numpy as np
import healpy as hp
import json
from django.template import loader, Context
from django.conf import settings
import matplotlib.pyplot as plt
import exceptions

cm = plt.get_cmap('jet')
num_colors = 11
colors = [cm(float(i)/num_colors) for i in range(num_colors) ]

root_folder = "../ddx92/"
out_folder = "../dx9null/"

def write_html(filename, t, c):
    with open(os.path.join(out_folder, filename), 'w') as f:
        f.write(t.render(c))

HORNS = {30:[27,28], 44:[24,25,26], 70:list(range(18,23+1))}

def bin(y, samples):
    y_split = np.array_split(y, len(y)/samples)
    return map(np.mean, y_split), np.array(map(np.std, y_split))/np.sqrt(samples)


def chlist(freq):
    horns = HORNS[freq]
    chs = []
    for horn in horns:
        chs += ["LFI%dM" % horn, "LFI%dS" % horn]
    return chs

print "MENU"

freqs = [30,44,70]#, 100, 143, 217, 353, 545, 857]

menu = [{"title":"Halfrings", "file":"index.html",
"links" : [("%d" % freq, "%dGHz" % freq) for freq in freqs]}]

for comp in "IQU":
    menu_item = {"title":"Surv Diff %s" % comp, "file":"surveydiff_%s.html" % comp}
    menu_item["links"] = [("%d_%s" % (freq, comp), "%dGHz %s" % (freq, comp))  for freq in freqs]
    menu.append(menu_item)

menu_item = {"title":"Horn Diff", "file":"horndiff.html"}
menu_item["links"] = [("%d_SS1" % freq, "%dGHz" % freq) for freq in [30, 44, 70]]
menu.append(menu_item)

menu_item = {"title":"Spectra", "file":"spectra.html"}
menu_item["links"] = [("%d_%s" % (freq, comp), "%dGHz %s" % (freq, comp)) for comp in "TE" for freq in freqs]
menu.append(menu_item)

try:
    settings.configure(TEMPLATE_DEBUG=True, TEMPLATE_DIRS=("templates/",))
except exceptions.RuntimeError:
    pass

print "HALFRINGS"

survs = ["nominal", "full"]
table_list = []
summary_table = {"labels":[], "rows":[], "labels_done":False}
for freq in freqs:
    chtags = [""]
    if freq == 70:
        chtags += ["18_23", "19_22", "20_21"]
    for ch in chtags:
        if ch:
            chtag = ch
        else:
            chtag = str(freq)
        for surv in survs:
            table = {"title":chtag + " SS %s" % str(surv), "tag":chtag, "labels":"IQU"}
            table["rows"] = []
            summary_row = dict(tag=table["title"], numslinks=[])
            row = {"images":[]}
            f=os.path.join(root_folder, "halfrings", "%s_SS%s_map.json" % (chtag, str(surv)))
            print(f)
            metadata = json.load(open(f))
            for comp in "IQU":
                row["images"].append({ "title":metadata["title"] + " %s" % comp, 
                    "file_name" : metadata["base_file_name"] + "_map_%s" % comp,
                    "tag" : metadata["base_file_name"].replace("/","_") + "_%s" % comp,
                    })
                if not summary_table["labels_done"]:
                    summary_table["labels"].append(comp)
                try:
                    summary_row["numslinks"].append(("%.2f" % (metadata["map_chi2_%s" % comp]), row["images"][-1]["file_name"]))
                except exceptions.KeyError:
                    summary_row["numslinks"].append(("%.2f" % (metadata["map_chi2"]),row["images"][-1]["file_name"] ))
            table["rows"].append(row)
            table_list.append(table)
            summary_table["labels_done"] = True
            summary_table["rows"].append(summary_row)


t = loader.get_template(template_name="page.html")
c = Context({
            'page_title': 'Halfrings',
            'table_list': table_list,
            'summary_tables': [summary_table],
            'menu':menu
            })

write_html("index.html", t, c)

print "SurveyDiff"

def swap_surv(comb):
    if comb[1] % 2 != 0 and comb[0] % 2 == 0:
        return (comb[1], comb[0])
    else:
        return comb

bp_tag = {True:"_bpcorr", False:''}
survs = [1,2,3,4,5]
for comp in "IQU":
    table_list = []
    summary_table = {"labels":[], "rows":[], "labels_done":False}
    for freq in freqs:
        chtags = [""]
        if freq == 70:
            chtags += ["18_23", "19_22", "20_21"]
        if freq < 100:
            chtags += chlist(freq)
        for ch in chtags:
            if ch or freq>70:
                bp_corrs = [False]
            else:
                bp_corrs = [False, True]
            for bp_corr in bp_corrs:
                if ch and ch.find("_") < 0 and comp in "QU":
                    pass
                else:
                    if ch:
                        chtag = ch
                    else:
                        chtag = str(freq)
                    table = {"title":chtag + " %s" % comp, "tag":chtag  + "_%s" % comp, "labels":map(str, survs[1:])}
                    if bp_corr:
                        table["title"] += " BPCORR"
                    table["rows"] = []
                    summary_row = dict(tag=table["title"], numslinks=[])
                    for surv in survs[:-1]:
                        row = {"tag":str(surv), "images":[]}
                        for surv2 in survs[1:]:
                            if surv2 <= surv:
                                row["images"].append(None)
                            else:
                                comb = swap_surv((surv, surv2))
                                if not summary_table["labels_done"]:
                                    summary_table["labels"].append("SS%d-SS%d" % comb)
                                metadata_filename = os.path.join(root_folder, "surveydiff", "%s_SS%d-SS%d%s_map.json" % (chtag, comb[0], comb[1], bp_tag[bp_corr]))
                                metadata=json.load(open(metadata_filename))
                                row["images"].append({"file_name":metadata["base_file_name"] + "_map_%s" % comp, "title":metadata["title"] + " %s" % comp, 
                        "tag" : metadata["base_file_name"].replace("/","_")+ "_%s" % comp,
                                    })
                                try:
                                    summary_row["numslinks"].append(("%.2f" % (metadata["map_chi2_%s" % comp]), row["images"][-1]["file_name"]))
                                except exceptions.KeyError:
                                    summary_row["numslinks"].append(("%.2f" % (metadata["map_chi2"]),row["images"][-1]["file_name"] ))

                        table["rows"].append(row)
                    table_list.append(table)
                    summary_table["labels_done"] = True
                    summary_table["rows"].append(summary_row)

    t = loader.get_template(template_name="page.html")
    c = Context({
                'page_title': 'Survey differences %s' % comp,
                'table_list': table_list,
                'summary_tables': [summary_table],
                'menu':menu
                })

    write_html("surveydiff_%s.html" % comp, t, c)

bp_tag = {True:"_bpcorr", False:''}
survs = [1,2,3,4,5]

cl_comp = {"T":0, "E":1, "B":2}
table_list = []
for comp in "TE":
    tag = comp + comp
    summary_table = {"labels":[], "rows":[], "labels_done":False}
    for freq in freqs:
        chtags = [""]
        if freq == 70:
            chtags += ["18_23", "19_22", "20_21"]
        #if freq < 100:
        #    chtags += chlist(freq)
        for ch in chtags:
            if ch or freq>70:
                bp_corrs = [False]
            else:
                bp_corrs = [False, True]
            for bp_corr in bp_corrs:
                if ch and ch.find("_") < 0 and comp in "QU":
                    pass
                else:
                    if ch:
                        chtag = ch
                    else:
                        chtag = str(freq)

                    plt.figure()
                    i=0

                    table = {"title":chtag + " %s" % tag, "tag":chtag  + "_%s" % comp, "labels":["Spectra"]}
                    if bp_corr:
                        table["title"] += " BPCORR"
                    table["rows"] = []
                    summary_row = dict(tag=table["title"], numslinks=[])
                    for surv in survs[:-1]:
                        for surv2 in survs[1:]:
                            if surv2 <= surv:
                                pass
                            else:
                                comb = swap_surv((surv, surv2))
                                metadata_filename = os.path.join(root_folder, "surveydiff", "%s_SS%d-SS%d%s_cl.json" % (chtag, comb[0], comb[1], bp_tag[bp_corr]))
                                metadata=json.load(open(metadata_filename))
                                spec = hp.read_cl(root_folder + metadata["file_name"])[cl_comp[comp]]*1e12
                                binspec, binspecerr = bin(spec, 5)
                                ell, eller = bin(np.arange(len(spec)), 5)
                                ell = np.arange(len(spec))
                                binspec = spec
                                plt.loglog(ell, binspec, label="SS%d-SS%d" % comb, color=colors[i])
                                i += 1

                    f=os.path.join(root_folder, "halfrings", "%s_SS%s_cl.json" % (chtag, "nominal"))
                    metadata=json.load(open(f))
                    spec = hp.read_cl(root_folder + metadata["file_name"])[cl_comp[comp]]*1e12
                    binspec, binspecerr = bin(spec, 5)
                    ell, eller = bin(np.arange(len(spec)), 5)
                    ell = np.arange(len(spec))
                    binspec = spec
                    plt.loglog(ell, binspec, label="nom half", color=colors[i])
                    #wmap = np.loadtxt("wmap_binned_%s_spectrum_7yr_v4p1.txt" % tag.lower())
                    #wmap_ell = wmap[:,0]
                    #wmap_spec = wmap[:,1] / (wmap_ell * (wmap_ell + 1) / (2*np.pi))
                    ##plt.loglog(wmap_ell, wmap_spec, label="wmap " + tag, color='k')
                    whitenoise_cl_field = "whitenoise_cl"
                    if comp in "EB":
                        whitenoise_cl_field += '_P'

                    plt.loglog(ell, np.ones_like(ell) * metadata[whitenoise_cl_field]*1e12, 'k--', label="variance")
                    #print metadata[whitenoise_cl_field]*1e12
                    plt.ylim([1e-3, 1])
                    leg = plt.legend(loc=0, prop={"size":9})
                    leg.legendPatch.set_alpha(0.5)
                    plt.grid()
                    plt.title(table["title"])
                    plt.xlim([0, 2048])
                    plt.xlabel("ell"); plt.ylabel("C_ell [microK^2]")
                    out_filename = "spectra/cl_%s_%s" % (chtag, comp)
                    if bp_corr:
                        out_filename += "_bpcorr"
                    plt.savefig(out_folder + "/images/" + out_filename + "_thumb.jpg")
                    print out_filename
                    row = {"tag":"", "images":[]}
                    row["images"].append({"file_name":out_filename,
                        "tag" : metadata["base_file_name"].replace("/","_")+ "_%s" % comp,
                                    })
                    table["rows"] = [row]
                    table_list.append(table)

t = loader.get_template(template_name="page.html")
c = Context({
            'page_title': 'Spectra',
            'table_list': table_list,
            'menu':menu
            })

write_html("spectra.html", t, c)

print "CHDIFF"
survs = [1,2,3,4,5]
freqs = [30, 44, 70]
    
table_list = []
summary_tables = []
for freq in freqs:
    summary_table = {"tag":freq, "labels":[], "rows":[], "labels_done":False}
    for surv in survs:
        horns = ["LFI%d" % h for h in HORNS[freq]]
        table = {"title":"%dGHz SS%s" % (freq, str(surv)), "tag":"%d_SS%s" % (freq, str(surv)), "labels":map(str, horns[1:])}
        table["rows"] = []
        summary_row = dict(tag=table["title"], numslinks=[])
        for hornn, horn in enumerate(horns[:-1]):
            row = {"tag":str(horn), "images":[]}
            for horn2n, horn2 in enumerate(horns[1:]):
                horn2n += 1
                if horn2n <= hornn:
                    row["images"].append(None)
                else:
                    comb = (horn, horn2)
                    if not summary_table["labels_done"]:
                        summary_table["labels"].append("%s-%s" % tuple(map(str, comb)))
                    metadata=json.load(open(os.path.join(root_folder, "chdiff", "%s-%s_SS%s_map.json" % (comb[0], comb[1], str(surv)))))
                    row["images"].append({"file_name":metadata["base_file_name"] + "_map_I", "title":metadata["title"], 
                    "tag" : metadata["base_file_name"].replace("/","_"),
                        })
                    try:
                        summary_row["numslinks"].append(("%.2f" % (metadata["map_chi2_%s" % comp]), row["images"][-1]["file_name"]))
                    except exceptions.KeyError:
                        summary_row["numslinks"].append(("%.2f" % (metadata["map_chi2"]),row["images"][-1]["file_name"] ))
            table["rows"].append(row)
        table_list.append(table)
        summary_table["labels_done"] = True
        summary_table["rows"].append(summary_row)
    summary_tables.append(summary_table)

t = loader.get_template(template_name="page.html")
c = Context({
            'page_title': 'Horn differences',
            'table_list': table_list ,
            'summary_tables': summary_tables,
            'menu':menu

            })

write_html("horndiff.html", t, c)
