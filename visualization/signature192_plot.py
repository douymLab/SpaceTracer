import argparse
import itertools
import os
import pandas as pd
import io
import numpy as np
import matplotlib.pyplot as plt

class plot_profile:
    def __init__(self,input_matrix,title="Count",refine=False):
        self.input_matrix=input_matrix
        self.title=title
        self.refine=refine

    def _plot(self,savepath):
        plot_192(self.input_matrix,self.title,savepath,self.refine)


#### functions
def reindex_sbs192(data_f):
    first = ["A", "C", "G", "T"]
    inner_bracket = [[x] * 16 for x in ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G","A>C","A>G","A>T","G>A","G>C","G>T"]]
    inner_bracket = [item for sublist in inner_bracket for item in sublist]
    outter_bracket = [x for x in list(itertools.product(first, first))]
    result = [
        outter_bracket[f % 16][0]
        + "["
        + inner_bracket[f]
        + "]"
        + outter_bracket[f % 16][1]
        for f in range(0, 192)
    ]
    
    #data_f = data_f.reindex(result)
    data_f = data_f.fillna(0) 
    default_df=pd.DataFrame({"MutationType":result,"default":[0]*192})
    merge_df=default_df.merge(data_f,on="MutationType",how="outer")
    merge_df=merge_df.drop(columns=["default"])
    merge_df=merge_df.fillna(0) 
    merge_df = merge_df.set_index("MutationType", drop=True)
    merge_df = merge_df.reindex(result)
    
    return merge_df


def getylabels(ylabels):
    if max(ylabels) >= 10**9:
        ylabels = ["{:.2e}".format(x) for x in ylabels]
        ylabels[0] = "0.00"
    else:
        if max(ylabels) <= 1000:
            ylabels = ["{:,.0f}".format(x) for x in ylabels]
            ylabels[0] = "0"
        elif max(ylabels) < 10**5 and max(ylabels) > 1000:
            ylabels = ["{:,.0f}".format(x / 1000) + "k" for x in ylabels]
            ylabels[0] = "0"
        else:  # if max(ylabels)>= 10**5:
            ylabels = ["{:,.0f}".format(x / (10**6)) + "m" for x in ylabels]
            ylabels[0] = "0"
    return ylabels


def refine_counts(profile):
    # 创建新的 DataFrame，并初始化 Count 列为 0
    complement_pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}

    new_profile = pd.DataFrame({"Count": [0] * profile.shape[0]}, index=profile.index)
    
    for i in profile.index:
        char_up = i[0]
        char_before = i[2]
        char_after = i[4]
        char_down = i[6]
        
        # 获取互补碱基
        complements = [complement_pairs[base] for base in [char_down, char_before, char_after, char_up]]
        complement_row = f"{complements[0]}[{complements[1]}>{complements[2]}]{complements[3]}"
        
        if complement_row in profile.index:
            if profile.at[i, "Count"] > profile.at[complement_row, "Count"]:
                new_profile.at[i, "Count"] = profile.at[i, "Count"] - profile.at[complement_row, "Count"]
            else:
                new_profile.at[i, "Count"] = 0

    return new_profile


def plot_192(matrix_path,title,savepath,refine):

    # matrix_path = "matrix.txt"
    false_sig1=pd.read_csv(matrix_path,sep="\t")
    false_sig=reindex_sbs192(false_sig1)

    if refine:
        false_sig2=refine_counts(false_sig)
        false_sig=false_sig2

    seq192=false_sig.index.tolist()
    ctx = false_sig.index  # [seq[0]+seq[2]+seq[6] for seq in data.index]

    colors = ["darkgreen", "brown","#E41A1C", "grey", "orange","#377EB8","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"]
    colorsall = [
        [colors[j] for i in range(int(len(ctx) / 12))] for j in range(12)
    ]
    colors_flat_list = [item for sublist in colorsall for item in sublist]


    buf = io.BytesIO()
    figs = {}
    for sample in false_sig.columns:
        plot1 = plt.figure(figsize=(43.93, 9.92))
        panel1 = plt.axes([0.04, 0.09, 1.90, 0.77])
        total_count = np.sum(false_sig[sample].values)
        x = 0.4
        ymax = 0
        i = 0
        muts = false_sig[sample].values
        if total_count > 0:
            plt.bar(
                np.arange(len(ctx)) + x,
                muts ,
                width=0.4,
                color=colors_flat_list,
                align="center",
                zorder=1000,
            )
        xlabels = []
        x = 0.4
        ymax = 0

        colors=["darkgreen", "brown","#E41A1C", "grey", "orange","#377EB8","#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"]
        xlabels = [seq[0] + seq[2] + seq[6] for seq in seq192]
        i = 0
        x = 0.043
        y3 = 0.87
        y = int(ymax * 1.25)
        y2 = y + 2
        for i in range(0, 12, 1):
            panel1.add_patch(
                plt.Rectangle(
                    (x, y3),
                    0.15,
                    0.05,
                    facecolor=colors[i],
                    clip_on=False,
                    transform=plt.gcf().transFigure,
                )
            )
            x += 0.159
        yText = y3 + 0.06
        plt.text(
            0.1,
            yText,
            "C>A",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.255,
            yText,
            "C>G",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.415,
            yText,
            "C>T",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.575,
            yText,
            "T>A",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.735,
            yText,
            "T>C",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.885,
            yText,
            "T>G",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            1.035,
            yText,
            "A>C",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            1.185,
            yText,
            "A>G",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            1.335,
            yText,
            "A>T",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            1.485,
            yText,
            "G>A",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            1.635,
            yText,
            "G>C",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            1.79,
            yText,
            "G>T",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        if y <= 4:
            y += 4
        while y % 4 != 0:
            y += 1
        # ytick_offest = int(y/4)
        y = ymax / 1.025
        ytick_offest = float(y / 3)
        labs = np.arange(0.375, 192.375, 1)
        panel1.set_xlim([0, 192])
        # panel1.set_ylim([0, y])
        panel1.set_xticks(labs)
        # panel1.set_yticks(ylabs)
        count = 0
        m = 0
        for i in range(0, 192, 1):
            plt.text(
                i / 101 + 0.0415,
                0.02,
                xlabels[i][0],
                fontsize=30,
                color="gray",
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                transform=plt.gcf().transFigure,
            )
            plt.text(
                i / 101 + 0.0415,
                0.044,
                xlabels[i][1],
                fontsize=30,
                color=colors[m],
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                fontweight="bold",
                transform=plt.gcf().transFigure,
            )
            plt.text(
                i / 101 + 0.0415,
                0.071,
                xlabels[i][2],
                fontsize=30,
                color="gray",
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                transform=plt.gcf().transFigure,
            )
            count += 1
            if count == 16:
                count = 0
                m += 1
        # line_pos=seq96.index(aim_sub[0])
        # # print(aim_sub,line_pos+ 0.3)
        # plt.axvline(x=line_pos+ 0.35, color='r', linestyle='--')
        
        # plt.gca().yaxis.grid(True)
        # plt.gca().grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
        panel1.set_xlabel("")
        panel1.set_ylabel("")
        panel1.tick_params(
            axis="both",
            which="both",
            bottom=False,
            labelbottom=False,
            left=True,
            labelleft=True,
            right=True,
            labelright=False,
            top=False,
            labeltop=False,
            direction="in",
            length=25,
            colors="gray",
            width=2,
        )
        plt.text(
            0.045,
            0.75,
            sample+"+"+ title +": "+ "{:,}".format(int(total_count)) + " subs",
            fontsize=30,
            weight="bold",
            color="black",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )

    plt.savefig(savepath,bbox_inches='tight',dpi=300)
    plt.close()


def main():
    matrix_path=args.matrix
    savepath=args.outname
    plot_192(matrix_path=matrix_path,title=args.title,savepath=savepath,refine=args.refine)

    
## parameters
parser = argparse.ArgumentParser()
# parser.add_argument("--mutations","-m", required=False,default="",help="mutation identifier list")
parser.add_argument("--matrix", required=True,help="feature file")
parser.add_argument("--outname", required=True,help="output png file")
parser.add_argument("--refine", required=False,action='store_true',help="do you want to refine the signatures?")
parser.add_argument("--title", required=False,default="count",help="add a mark for out png")

args = parser.parse_args()
    
if __name__ == '__main__':
    main()  
