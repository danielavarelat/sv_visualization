from PIL import Image, ImageFile
from PIL import ImageFont, ImageDraw
import pandas as pd
from os.path import join
import numpy as np
import pysam


def get_df_full(outdir, samples_ids):
    """Get a dataframe with ALL the information about the variants."""
    with open(join(outdir, "names.txt")) as n:
        list_names = n.readlines()
    names = [word.strip() for word in list_names]
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
    df_in.insert(0, "names_files", names)
    callers = df_in["callers"].values.tolist()
    for sample_id in samples_ids:
        covs_sample = df_in["COV {}".format(sample_id)].values.tolist()
        max_covs_sample = []
        number_callers = []
        for sv, i in zip(covs_sample, callers):
            num_calls = len(i.split(","))
            number_callers.append(num_calls)
            cov1, cov2 = sv.replace("L", "").replace("(", "").split(",")
            cov2 = float(cov2.replace(")", ""))
            cov1 = float(cov1)
            max_cov = max(cov1, cov2)
            max_covs_sample.append(max_cov)
        df_in[sample_id] = max_covs_sample
    maxs = []
    for sv in range(df_in.shape[0]):
        max_cov = max(df_in.loc[sv, samples_ids].tolist())
        maxs.append(max_cov)
    df_full = df_in
    df_full["Max Cov"] = maxs
    df_full["Number calls"] = number_callers
    return df_full


def get_sv_path(outdir, sample_id, sv, df):
    list_svs = df["name"].tolist()
    files = df["names_files"].tolist()
    index_sv = list_svs.index(sv)
    base_name_sv = files[index_sv]
    filename = "{}_{}.{}".format(sample_id, base_name_sv, "png")
    path_snap = join(outdir, sample_id, "snapshots", filename)
    return path_snap


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
    font_size = 50
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


def get_df_filter(df_out, filters):
    df_to_filter = df_out[filters]
    return df_to_filter


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
            # df_path.iloc[sample, sv] = testing_img
    df_svs.loc["annotation"] = df_out["annot"].tolist()
    df_svs.insert(0, "CIRCOS PLOT", circos_list + [""])
    df_path.insert(0, "CIRCOS PLOT", circos_list)
    df_svs.insert(0, "SAMPLE ID", samples_ids + [""])
    keys = list_svs
    values = df_out["annot"].tolist()
    dict_annotation = dict(zip(keys, values))
    # df_trans = df_path.transpose()
    return df_svs, df_path, dict_annotation


def callers_color(sample_id, name_sv, df_out):
    list_svs = df_out["name"].tolist()
    index_sv = list_svs.index(name_sv)
    callers = df_out["callers"].tolist()
    callers_sv = callers[index_sv]
    called = sample_id in callers_sv
    return called


def create_imagen(
    outdir, rows_show, selected_svs, df_out, df_path, dict_annotation
):
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
#     mergedImg.save(join(outdir, "%s.png" % name_file))
#     status = "SAVED IMAGE: {}".format(join(outdir, "%s.png" % name_file))
    return mergedImg

