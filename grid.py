from PIL import Image, ImageFile
from PIL import ImageFont, ImageDraw
import pandas as pd
from os.path import join, isfile
from os import listdir
import numpy as np
import pysam
import sys

def highlight(img, coords, fill, width=5):
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle(coords, fill=fill)
    return img


def centered_pos(draw, x, y, txt, font, axis):
    tx, ty = draw.textsize(txt, font=font)
    if axis == "x":
        cx = x - tx / 2
        cy = y
    elif axis == "y":
        cx = x
        cy = y - ty / 2
    else:
        cx = x - tx / 2
        cy = y - ty / 2
    return (cx, cy)


def draw_labels(labels, positions, img, centering_axis=""):
    font_size = 40
    global_font = "/home/varelad/dash/LiberationSans-Regular.ttf"
    font = ImageFont.truetype(global_font, font_size, encoding="unic")
    for label, pos in zip(labels, positions):
        draw = ImageDraw.Draw(img)
        if centering_axis:
            draw.text(
                centered_pos(draw, pos[0], pos[1], label, font, centering_axis),
                label,
                fill="black",
                font=font,
            )
        else:
            draw.text((pos[0], pos[1]), label, fill="black", font=font)
    return img

def callers_color(sample_id, name_sv, df_out):
    list_svs = df_out["name"].tolist()
    index_sv = list_svs.index(name_sv)
    callers = df_out["callers"].tolist()
    callers_sv = callers[index_sv]
    called = sample_id in callers_sv
    return called


def create_imagen(outdir, rows_show, selected_svs, df_out, df_path, dict_annotation):
    columns_show = ["CIRCOS PLOT"] + selected_svs
    subset = df_path.loc[rows_show, columns_show]
    img_table = subset
    grid_size = 1000
    margin = 1 / 5
    margin = int(grid_size * margin)
    theight = grid_size * img_table.shape[0] + margin
    twidth = grid_size * img_table.shape[1] + margin
    mergedImg = Image.new("RGB", (twidth, theight), "white")
    for (i, j) in np.ndindex(img_table.shape):
        y, x = (i * grid_size + margin, j * grid_size + margin)
        mergedImg.paste(
            Image.open(img_table.iloc[i, j]).resize((grid_size, grid_size)), (x, y)
        )
        if j == 0:
            continue
        else:
            called = callers_color(img_table.index[i], img_table.columns[j], df_out)
        if called:
            mergedImg = highlight(
                mergedImg,
                [(x, y), (x + grid_size, y + grid_size)],
                (240, 215, 215, 100),
            )
        else:
            mergedImg = highlight(
                mergedImg, [(x, y), (x + grid_size, y + grid_size)], (180, 216, 230, 80)
            )
    annot = [dict_annotation.get(key) for key in selected_svs]
    labels = ["\n".join((label, annot)) for (label, annot) in zip(selected_svs, annot)]
    labels_circos = ["CIRCOS PLOT"] + labels
    positions = [
        (margin + i * grid_size + grid_size / 2, margin / 2)
        for i in range(len(labels_circos))
    ]
    mergedImg = draw_labels(labels_circos, positions, mergedImg, "xy")
    return mergedImg

def get_df_full(outdir, samples_ids):
    """Get a dataframe with ALL the information about the variants."""
    columns_df = (
        [
            "svs",
            "chr1",
            "start1",
            "end1",
            "chr2",
            "start2",
            "end2",
            "name",
            "score",
            "strand1",
            "strand2",
            "type",
            "annot",
        ]
        + ["COV {}".format(sample_id) for sample_id in samples_ids]
        + ["callers"]
    )
    svs_table = join(outdir, "merged_svs_annot.tsv")
    df_in = pd.read_csv(svs_table, sep="\t", header=None, names=columns_df)
    df_full = df_in
    return df_full

def add_names_files(outdir, samples_ids):
    path=join(outdir, samples_ids[0], "snapshots")
    allfiles = [f for f in listdir(path) if isfile(join(path, f))]    
    allfiles.remove('igvbatch_0.txt')
    #without .png
    files = [i.split(".")[0] for i in allfiles]
    d={}
    for i in files:
        vals=[float(x) if x != "X" else 24 for x in i.split("_")]
        vals= vals[0:4]
        d[i]=vals
# ORDER VALUES ACCORDING TO CHROMOSOME AND POSITION
    values = [ k for k in d.values() ]
    values.sort(key = lambda l: (l[0], l[1], l[2], l[3]))
    order_files=[]
    for value in values:
        for file in files:
            if value==d[file]:
                order_files.append(file)
    return order_files


def get_sv_path(outdir, sample_id, sv, df):
    list_svs = df["name"].tolist()
    files = df["names_files"].tolist()
    index_sv = list_svs.index(sv)
    base_name_sv = files[index_sv]
#     filename = "{}_{}.{}".format(sample_id, base_name_sv, "png")
    filename = "{}.png".format(base_name_sv)
    path_snap = join(outdir, sample_id, "snapshots", filename)
    return path_snap


def make_grids(outdir, samples_ids, df_out):
    list_svs = df_out["name"].tolist()
    rows = samples_ids + ["annotation"]
    df_svs = pd.DataFrame(index=rows, columns=list_svs)
    df_path = pd.DataFrame(index=samples_ids, columns=list_svs)
    circos_list = []
    for sample in range(len(samples_ids)):
        sample_id = samples_ids[sample]
        circos = join(outdir, sample_id, "%s_circos.png" % sample_id)
        circos_list.append(circos)
        for sv in range(len(list_svs)):
            name_sv = list_svs[sv]
            df_svs.iloc[sample, sv] = [sample_id, name_sv]
            df_path.iloc[sample, sv] = get_sv_path(outdir, sample_id, name_sv, df_out)
    df_svs.loc["annotation"] = df_out["annot"].tolist()
    df_svs.insert(0, "CIRCOS PLOT", circos_list + [""])
    df_path.insert(0, "CIRCOS PLOT", circos_list)
    df_svs.insert(0, "SAMPLE ID", samples_ids + [""])
    keys = list_svs
    values = df_out["annot"].tolist()
    dict_annotation = dict(zip(keys, values))
    return df_svs, df_path, dict_annotation

def filtered(outdir):
    table_filtered = join(outdir,  "filtered_table.tsv")
    df = pd.read_csv(table_filtered, sep="\t", header=None)
    new_header = df.iloc[0] #grab the first row for the header
    df = df[1:] #data less the header row
    df.columns = new_header
    indexes=df[df.columns[0]]
    table = df.set_index(indexes).iloc[:,1:]  #first columns as index
    cols=indexes.tolist()
    return cols, table




outdir="/work/isabl/home/mccartej/p230_isabl/analysis/I-H-136593/svs"
samples_ids=["I-H-136593-N1-1-D1-1", "I-H-136593-T1-1-D1-1","I-H-136593-T2-1-D1-1"]
df_out=get_df_full(outdir, samples_ids)
# print(df_out)
svs, table = filtered(outdir)

df_out_filt = df_out[df_out["name"].isin(svs)]
# print(df_out_filt)
order_files = add_names_files(outdir, samples_ids)
# print(order_files)
df_out_filt.insert(0, "names_files", order_files)
print(df_out_filt)
df_svs, df_path, dict_annotation=make_grids(outdir, samples_ids, df_out_filt)

grid=create_imagen(outdir, samples_ids, svs, df_out_filt, df_path, dict_annotation)
img_path="/home/varelad/136593.png"
grid.save(img_path)
print(img_path)