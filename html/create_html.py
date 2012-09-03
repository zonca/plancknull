import os
import operator
import json
from glob import glob
from django.template import loader, Context
from django.conf import settings
import exceptions

root_folder = "../dx9"
out_folder = "../dx9null/"

def write_html(filename, t, c):
    with open(os.path.join(out_folder, filename), 'w') as f:
        f.write(t.render(c))

HORNS = {30:[27,28], 44:[24,25,26], 70:list(range(18,23+1))}

def chlist(freq):
    horns = HORNS[freq]
    chs = []
    for horn in horns:
        chs += ["LFI%dM" % horn, "LFI%dS" % horn]
    return chs

print "MENU"

freqs = [30,44,70, 100, 143, 217, 353, 545, 857]

menu = [{"title":"Halfrings", "file":"index.html",
"links" : [("%d" % freq, "%dGHz" % freq) for freq in freqs]}]

for comp in "IQU":
    menu_item = {"title":"Survey Differences %s" % comp, "file":"surveydiff_%s.html" % comp}
    menu_item["links"] = [("%d_%s" % (freq, comp), "%dGHz %s" % (freq, comp))  for freq in freqs]
    menu.append(menu_item)

menu_item = {"title":"Horn Differences", "file":"horndiff.html"}
menu_item["links"] = [("%d_SS1" % freq, "%dGHz" % freq) for freq in [30, 44, 70]]
menu.append(menu_item)

try:
    settings.configure(TEMPLATE_DEBUG=True, TEMPLATE_DIRS=("templates/",))
except exceptions.RuntimeError:
    pass

print "HALFRINGS"

survs = ["nominal", "full"]
table_list = []
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
            row = {"images":[]}
            f=os.path.join(root_folder, "halfrings", "%s_SS%s_map.json" % (chtag, str(surv)))
            metadata = json.load(open(f))
            for comp in "IQU":
                row["images"].append({ "title":metadata["title"] + " %s" % comp, 
                    "file_name" : metadata["base_file_name"] + "_map_%s" % comp,
                    "tag" : metadata["base_file_name"].replace("/","_") + "_%s" % comp,
                    })
            table["rows"].append(row)
            table_list.append(table)


t = loader.get_template(template_name="page.html")
c = Context({
            'page_title': 'Halfrings',
            'table_list': table_list,
            'menu':menu
            })

write_html("index.html", t, c)

print "SurveyDiff"

def swap_surv(comb):
    if comb[1] % 2 != 0 and comb[0] % 2 == 0:
        return (comb[1], comb[0])
    else:
        return comb

survs = [1,2,3,4,5]
for comp in "IQU":
    table_list = []
    for freq in freqs:
        chtags = [""]
        if freq == 70:
            chtags += ["18_23", "19_22", "20_21"]
        if freq < 100:
            chtags += chlist(freq)
        for ch in chtags:
            if ch and ch.find("_") < 0 and comp in "QU":
                pass
            else:
                if ch:
                    chtag = ch
                else:
                    chtag = str(freq)
                table = {"title":chtag + " %s" % comp, "tag":chtag  + "_%s" % comp, "labels":map(str, survs[1:])}
                table["rows"] = []
                for surv in survs[:-1]:
                    row = {"tag":str(surv), "images":[]}
                    for surv2 in survs[1:]:
                        if surv2 <= surv:
                            row["images"].append(None)
                        else:
                            comb = swap_surv((surv, surv2))
                            metadata=json.load(open(os.path.join(root_folder, "surveydiff", "%s_SS%d-SS%d_map.json" % (chtag, comb[0], comb[1]))))
                            row["images"].append({"file_name":metadata["base_file_name"] + "_map_%s" % comp, "title":metadata["title"] + " %s" % comp, 
                    "tag" : metadata["base_file_name"].replace("/","_")+ "_%s" % comp,
                                })
                    table["rows"].append(row)
                table_list.append(table)

    t = loader.get_template(template_name="page.html")
    c = Context({
                'page_title': 'Survey differences %s' % comp,
                'table_list': table_list,
                'menu':menu
                })

    write_html("surveydiff_%s.html" % comp, t, c)

print "CHDIFF"
survs = [1,2,3,4,5]
freqs = [30, 44, 70]
    
table_list = []
for freq in freqs:
    for surv in survs:
        horns = ["LFI%d" % h for h in HORNS[freq]]
        table = {"title":"%dGHz SS%s" % (freq, str(surv)), "tag":"%d_SS%s" % (freq, str(surv)), "labels":map(str, horns[1:])}
        table["rows"] = []
        for hornn, horn in enumerate(horns[:-1]):
            row = {"tag":str(horn), "images":[]}
            for horn2n, horn2 in enumerate(horns[1:]):
                horn2n += 1
                if horn2n <= hornn:
                    row["images"].append(None)
                else:
                    comb = (horn, horn2)
                    metadata=json.load(open(os.path.join(root_folder, "chdiff", "%s-%s_SS%s_map.json" % (comb[0], comb[1], str(surv)))))
                    row["images"].append({"file_name":metadata["base_file_name"] + "_map_I", "title":metadata["title"], 
                    "tag" : metadata["base_file_name"].replace("/","_"),
                        })
            table["rows"].append(row)
        table_list.append(table)

t = loader.get_template(template_name="page.html")
c = Context({
            'page_title': 'Horn differences',
            'table_list': table_list ,
            'menu':menu

            })

write_html("horndiff.html", t, c)
