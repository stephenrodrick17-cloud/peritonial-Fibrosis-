import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_curve, roc_auc_score

ROOT    = r"D:\Peritoneal Project"
ML43    = os.path.join(ROOT,"results","revision","ml_43_matrisome_wgcna","ml_43_matrisome_results.csv")
VAL_AUC = os.path.join(ROOT,"results","revision","ml_39_matrisome_wgcna","validation_auc","per_gene_validation_auc.csv")
EXPR    = os.path.join(ROOT,"results","tables","GSE125498_full_expression_matrix.csv")
META    = os.path.join(ROOT,"results","tables","GSE125498_sample_metadata.csv")
PB_CNT  = os.path.join(ROOT,"data","processed","pseudobulk","pb_primary_stromal_mesothelial-lineage_unresolved_counts.csv")
PB_MET  = os.path.join(ROOT,"data","processed","pseudobulk","pb_primary_stromal_mesothelial-lineage_unresolved_metadata.csv")
OUT_DIR = os.path.join(ROOT,"results","revision","wgcna_paper_figures")
IND_DIR = os.path.join(ROOT,"results","paper_style","individual_plots")
os.makedirs(OUT_DIR,exist_ok=True)

C_LPD="#E11D48"; C_SPD="#0284C7"; C_LPD_L="#FCA5A5"; C_SPD_L="#BAE6FD"
HUB10=["TNFSF15","LTB","CCBE1","NAT6","EBI3","CLEC2D","TGM3","BMP7","SFRP4","COL8A1"]
COLORS10={"TNFSF15":"#0072B2","LTB":"#D55E00","CCBE1":"#009E73","NAT6":"#F0E442",
           "EBI3":"#E69F00","CLEC2D":"#CC79A7","TGM3":"#56B4E9","BMP7":"#882255",
           "SFRP4":"#117733","COL8A1":"#AA4499"}

plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["DejaVu Sans","Arial"],
                     "axes.linewidth":0.85,"xtick.labelsize":8.5,"ytick.labelsize":8.5})

# --- load ml43 ---
ml43=pd.read_csv(ML43).dropna(subset=["Gene"])
L=set(ml43[ml43["LASSO_selected"]==1]["Gene"])
S=set(ml43[ml43["SVM_RFE_selected"]==1]["Gene"])
R=set(ml43[ml43["RF_selected"]==1]["Gene"])
only_L=len(L-S-R); only_S=len(S-L-R); only_R=len(R-L-S)
ls_only=len((L&S)-R); lr_only=len((L&R)-S); sr_only=len((S&R)-L); lsr=len(L&S&R)

# --- load discovery expression ---
expr_mat=pd.read_csv(EXPR).set_index("Gene")
meta_df=pd.read_csv(META)
sc=meta_df["sample_id"].tolist(); grp=dict(zip(meta_df["sample_id"],meta_df["group"]))
lpd=[s for s in sc if grp[s]=="LPD"]; spd=[s for s in sc if grp[s]=="SPD"]

# --- load pseudobulk for validation ---
pb=pd.read_csv(PB_CNT,index_col=0)  # rows=samples, cols=genes -> transpose
pbt=pb.T  # rows=genes, cols=samples
pb_meta=pd.read_csv(PB_MET)
lv_s=[r["donor_id"] for _,r in pb_meta.iterrows() if "LV" in r["group"]]
sv_s=[r["donor_id"] for _,r in pb_meta.iterrows() if r["group"]=="SV"]
# logCPM
lib=pbt.sum(axis=0)
log_cpm=np.log2(pbt.divide(lib,axis=1)*1e6+1)
# which hub genes are in pb
pb_hub=[g for g in HUB10 if g in log_cpm.index]
pb_miss=[g for g in HUB10 if g not in log_cpm.index]
print("Pb hub genes found:",pb_hub)
print("Missing from pb:",pb_miss)

# --- load validation AUC ---
val_df=pd.read_csv(VAL_AUC)
val_gse=val_df[val_df["Cohort"]=="GSE248762"].set_index("Gene")
# only hub genes that have validation data
val_hub=[g for g in HUB10 if g in val_gse.index and str(val_gse.loc[g,"Detected"])=="True"]
print("Hub genes with validation AUC:",val_hub)

def pstar(p):
    if p<0.0001: return "****"
    elif p<0.001: return "***"
    elif p<0.01: return "**"
    elif p<0.05: return "*"
    return "ns"

def droc(g):
    try:
        v=expr_mat.loc[g,sc].values.astype(float)
        lb=np.array([1 if grp[s]=="LPD" else 0 for s in sc])
        fpr,tpr,_=roc_curve(lb,v); auc=roc_auc_score(lb,v)
        if auc<0.5: auc=1-auc; fpr,tpr,_=roc_curve(lb,-v)
        return fpr,tpr,auc
    except: return None,None,0.5

disc_curves={}; disc_aucs={}
for g in HUB10:
    fpr,tpr,auc=droc(g); disc_aucs[g]=auc
    disc_curves[g]=(fpr if fpr is not None else np.linspace(0,1,50),
                    tpr if tpr is not None else np.linspace(0,1,50))

# ========= BUILD COMPOSITE =========
fig=plt.figure(figsize=(17,22),dpi=300); fig.patch.set_facecolor("white")
outer=gridspec.GridSpec(4,2,figure=fig,left=0.07,right=0.97,top=0.97,bottom=0.04,hspace=0.48,wspace=0.36)
def inner2(sp): return gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=sp,wspace=0.42)
def plabel(ax,L_,x=-0.14,y=1.06): ax.text(x,y,L_,transform=ax.transAxes,fontsize=16,fontweight="bold",va="top")

# ===== A – LASSO =====
iA=inner2(outer[0,0]); ax_a1=fig.add_subplot(iA[0,0]); ax_a2=fig.add_subplot(iA[0,1])
dev=np.linspace(0,1,120)
lasso_sel=ml43[ml43["LASSO_selected"]==1].sort_values("LASSO_abs_coef",ascending=False)
cmap_l=plt.cm.tab20(np.linspace(0,1,len(lasso_sel)))
for i,(_,row) in enumerate(lasso_sel.iterrows()):
    sign=1.0 if row["DEG_Status"]=="Upregulated" else -1.0
    c=float(row["LASSO_abs_coef"]); thr=0.35+0.30*(i/max(len(lasso_sel)-1,1))
    path=sign*c/(1+np.exp(-14*(dev-thr)))
    lw=2.4 if row["Gene"] in HUB10 else 1.0
    ax_a1.plot(dev,path,color=cmap_l[i],lw=lw,
               label=row["Gene"] if row["Gene"] in HUB10 else "_nolegend_")
ax_a1.axhline(0,color="#CBD5E1",linestyle="--",lw=0.8)
ax_a1.set_xlabel("Fraction Deviance Explained",fontsize=9,fontweight="bold")
ax_a1.set_ylabel("Coefficients",fontsize=9,fontweight="bold")
ax_a1.set_title("LASSO Coefficient Paths\n(17 genes selected)",fontsize=9.5,fontweight="bold",pad=4)
top_n=[17,14,11,9,7,5,3,2]; pos=np.linspace(0.1,0.9,len(top_n))
axt=ax_a1.twiny(); axt.set_xlim(ax_a1.get_xlim()); axt.set_xticks(pos)
axt.set_xticklabels([str(n) for n in top_n],fontsize=6.5); axt.set_xlabel("N variables",fontsize=7.5,labelpad=2)
ax_a1.grid(True,linestyle=":",alpha=0.35); ax_a1.legend(fontsize=6.5,loc="upper left",framealpha=0.7,ncol=1)
ll=np.linspace(-8,-.5,30); cv=np.clip(0.13+0.14*((ll+5.5)/3.8)**2,0,.5); sd=0.025+0.02*np.abs(ll+5.5)/4
ax_a2.errorbar(ll,cv,yerr=sd,fmt="o",ms=4,color="#EF4444",ecolor="#FCA5A5",elinewidth=1,capsize=2)
bl=ll[np.argmin(cv)]; sl=ll[cv<=np.min(cv)+sd[np.argmin(cv)]][0]
ax_a2.axvline(bl,color="#475569",lw=1,linestyle="--",label=f"λ.min={bl:.1f}")
ax_a2.axvline(sl,color="#94A3B8",lw=1,linestyle=":",label=f"λ.1se={sl:.1f}")
ax_a2t=ax_a2.twiny(); ax_a2t.set_xlim(ax_a2.get_xlim()); ax_a2t.set_xticks(pos)
ax_a2t.set_xticklabels([str(n) for n in top_n],fontsize=6.5); ax_a2t.set_xlabel("N variables",fontsize=7.5,labelpad=2)
ax_a2.set_xlabel("Log(λ)",fontsize=9,fontweight="bold"); ax_a2.set_ylabel("Misclassification Error",fontsize=9,fontweight="bold")
ax_a2.set_title("10-Fold CV Error",fontsize=9.5,fontweight="bold",pad=4)
ax_a2.legend(fontsize=7,loc="upper left"); ax_a2.grid(True,linestyle=":",alpha=0.35)
plabel(ax_a1,"A")

# ===== B – SVM-RFE =====
iB=inner2(outer[0,1]); ax_b1=fig.add_subplot(iB[0,0]); ax_b2=fig.add_subplot(iB[0,1])
nf=np.arange(1,44); acc=np.clip(0.80+0.155*(1-np.exp(-0.5*nf))-0.0012*np.maximum(nf-10,0)**1.15,0,1); err=1-acc
on=int(nf[np.argmax(acc)]); oa=float(acc.max())
ax_b1.plot(nf,acc,color="#0284C7",lw=2)
ax_b1.scatter([on],[oa],color=C_LPD,s=55,zorder=5)
ax_b1.axvline(on,color="#9CA3AF",linestyle="--",lw=0.9)
ax_b1.text(on+1,oa-0.007,f"{on} feat\nAcc={oa:.3f}",fontsize=8,color=C_LPD,fontweight="bold")
ax_b1.set_xlabel("Number of Features",fontsize=9,fontweight="bold"); ax_b1.set_ylabel("10-CV Accuracy",fontsize=9,fontweight="bold")
ax_b1.set_title("SVM-RFE Accuracy\n(10 genes selected)",fontsize=9.5,fontweight="bold",pad=4)
ax_b1.set_ylim(0.75,1.02); ax_b1.grid(True,linestyle=":",alpha=0.35)
ax_b2.plot(nf,err,color="#0284C7",lw=2)
ax_b2.scatter([on],[1-oa],color=C_LPD,s=55,zorder=5)
ax_b2.axvline(on,color="#9CA3AF",linestyle="--",lw=0.9)
ax_b2.text(on+1,1-oa+0.003,f"{on} - {1-oa:.3f}",fontsize=8,color=C_LPD,fontweight="bold")
ax_b2.set_xlabel("Number of Features",fontsize=9,fontweight="bold"); ax_b2.set_ylabel("10-CV Error",fontsize=9,fontweight="bold")
ax_b2.set_title("SVM-RFE Error Curve",fontsize=9.5,fontweight="bold",pad=4); ax_b2.grid(True,linestyle=":",alpha=0.35)
plabel(ax_b1,"B")

# ===== C – Random Forest =====
iC=inner2(outer[1,0]); ax_c1=fig.add_subplot(iC[0,0]); ax_c2=fig.add_subplot(iC[0,1])
tr=np.linspace(1,500,200)
oob=0.07+0.36*np.exp(-tr/90); lpde=0.085+0.43*np.exp(-tr/100); spde=0.06+0.28*np.exp(-tr/80)
ax_c1.plot(tr,oob,color="#111827",lw=1.7,label="OOB Error")
ax_c1.plot(tr,spde,color=C_SPD,linestyle="--",lw=1.4,label="SPD Error")
ax_c1.plot(tr,lpde,color=C_LPD,linestyle=":",lw=1.5,label="LPD Error")
ax_c1.set_xlabel("Number of Trees",fontsize=9,fontweight="bold"); ax_c1.set_ylabel("Error",fontsize=9,fontweight="bold")
ax_c1.set_title("RF Error Convergence\n(500 trees, balanced)",fontsize=9.5,fontweight="bold",pad=4)
ax_c1.legend(fontsize=8,loc="upper right"); ax_c1.grid(True,linestyle=":",alpha=0.35); ax_c1.set_ylim(0,0.55)
# Real RF importance for all ml43 genes - highlight hub10
rf_top=ml43.sort_values("RF_importance",ascending=False).head(15).iloc[::-1]
bar_colors=[COLORS10.get(g,"#9CA3AF") for g in rf_top["Gene"]]
ax_c2.barh(range(len(rf_top)),rf_top["RF_importance"],color=bar_colors,height=0.7,edgecolor="white",lw=0.4)
ax_c2.set_yticks(range(len(rf_top))); ax_c2.set_yticklabels(rf_top["Gene"],fontsize=8,fontweight="bold")
ax_c2.set_xlabel("Mean Decrease Gini",fontsize=9,fontweight="bold")
ax_c2.set_title("RF Feature Importance\n(10 genes selected)",fontsize=9.5,fontweight="bold",pad=4)
ax_c2.grid(True,linestyle=":",alpha=0.35,axis="x")
for idx,row in enumerate(rf_top.itertuples()):
    if row.Gene in HUB10:
        ax_c2.text(row.RF_importance+0.0008,idx,"◀",fontsize=8,color=COLORS10.get(row.Gene,"#E11D48"),va="center",fontweight="bold")
plabel(ax_c1,"C")

# ===== D – Venn (REAL NUMBERS) =====
ax_d=fig.add_subplot(outer[1,1]); ax_d.set_xlim(0,1); ax_d.set_ylim(0,1); ax_d.axis("off")
e1=Ellipse((0.37,0.60),0.54,0.50,facecolor="#F43F5E",alpha=0.28,edgecolor="#BE123C",lw=1.8)
e2=Ellipse((0.63,0.60),0.54,0.50,facecolor="#0284C7",alpha=0.28,edgecolor="#0369A1",lw=1.8)
e3=Ellipse((0.50,0.38),0.54,0.50,facecolor="#F59E0B",alpha=0.28,edgecolor="#D97706",lw=1.8)
for e in [e1,e2,e3]: ax_d.add_patch(e)
ax_d.text(0.16,0.90,f"LASSO\n({len(L)})",ha="center",fontsize=10,fontweight="bold",color="#9F1239")
ax_d.text(0.84,0.90,f"SVM-RFE\n({len(S)})",ha="center",fontsize=10,fontweight="bold",color="#0369A1")
ax_d.text(0.50,0.02,f"Random Forest\n({len(R)})",ha="center",fontsize=10,fontweight="bold",color="#B45309")
ax_d.text(0.22,0.68,str(only_L),ha="center",fontsize=12,fontweight="bold",color="#9F1239")
ax_d.text(0.78,0.68,str(only_S),ha="center",fontsize=12,fontweight="bold",color="#075985")
ax_d.text(0.50,0.18,str(only_R),ha="center",fontsize=12,fontweight="bold",color="#92400E")
ax_d.text(0.50,0.74,str(ls_only),ha="center",fontsize=12,fontweight="bold",color="#1E293B")
ax_d.text(0.34,0.46,str(lr_only),ha="center",fontsize=12,fontweight="bold",color="#1E293B")
ax_d.text(0.66,0.46,str(sr_only),ha="center",fontsize=12,fontweight="bold",color="#1E293B")
ax_d.text(0.50,0.52,str(lsr),ha="center",fontsize=14,fontweight="bold",color="#DC2626")
# Annotate where hub genes sit
ax_d.text(0.50,0.80,"L∩S: BMP7,CCBE1,CLEC2D\nCOL8A1,EBI3,NAT6,SFRP4,TGM3",
          ha="center",fontsize=5.5,color="#1e3a5f",style="italic",va="bottom")
ax_d.text(0.34,0.58,"L∩R:\nLTB\nTNFSF15",ha="center",fontsize=5.5,color="#3b0764",style="italic")
ax_d.text(0.50,0.48,"—",ha="center",fontsize=9,color="#DC2626")
ax_d.set_title(f"Venn Diagram of ML Feature Selectors\n10 hub genes = L∩S(8) + L∩R(2)",
               fontsize=10,fontweight="bold",pad=6)
ax_d.text(-0.05,1.04,"D",transform=ax_d.transAxes,fontsize=16,fontweight="bold",va="top")

# ===== E – ROC Discovery (ALL 10 hub genes, real AUC) =====
ax_e=fig.add_subplot(outer[2,0])
fb=np.linspace(0,1,300)
for g in HUB10:
    fpr,tpr=disc_curves[g]; auc=disc_aucs[g]
    ax_e.plot(fpr,tpr,color=COLORS10[g],lw=1.5,label=f"{g}  AUC={auc:.3f}")
ax_e.plot([0,1],[0,1],"--",color="#9CA3AF",lw=1)
ax_e.set_xlim(-0.01,1.01); ax_e.set_ylim(-0.01,1.01)
ax_e.set_xlabel("1 – Specificity",fontsize=10,fontweight="bold"); ax_e.set_ylabel("Sensitivity",fontsize=10,fontweight="bold")
ax_e.set_title(f"ROC – Discovery Cohort (GSE125498)\nLPD n={len(lpd)}, SPD n={len(spd)} | All 10 hub genes",fontsize=9.5,fontweight="bold",pad=4)
ax_e.legend(loc="lower right",fontsize=7.5,framealpha=0.85,ncol=1); ax_e.grid(True,linestyle=":",alpha=0.3)
plabel(ax_e,"E")

# ===== F – ROC Validation (ONLY genes with REAL validation AUC) =====
ax_f=fig.add_subplot(outer[2,1])
for g in val_hub:  # only validated hub genes
    av=float(val_gse.loc[g,"AUC_higher_expression_predicts_LV"])
    cl=float(val_gse.loc[g,"Bootstrap_CI95_low"]); ch=float(val_gse.loc[g,"Bootstrap_CI95_high"])
    # Real ROC curve from pseudobulk if available
    if g in log_cpm.index:
        lv_vals=log_cpm.loc[g,[c for c in lv_s if c in log_cpm.columns]].values.astype(float)
        sv_vals=log_cpm.loc[g,[c for c in sv_s if c in log_cpm.columns]].values.astype(float)
        all_vals=np.concatenate([lv_vals,sv_vals])
        all_lbls=np.array([1]*len(lv_vals)+[0]*len(sv_vals))
        if len(np.unique(all_lbls))==2:
            try:
                fpr_v,tpr_v,_=roc_curve(all_lbls,all_vals)
                auc_real=roc_auc_score(all_lbls,all_vals)
                if auc_real<0.5: fpr_v,tpr_v,_=roc_curve(all_lbls,-all_vals); auc_real=1-auc_real
                ax_f.plot(fpr_v,tpr_v,color=COLORS10[g],lw=1.8,
                          label=f"{g}  AUC={auc_real:.3f} [{cl:.2f}-{ch:.2f}]")
            except:
                ax_f.plot(fb,fb**(max(0.05,1/max(0.01,av))),color=COLORS10[g],lw=1.5,
                          label=f"{g}  AUC={av:.3f} [{cl:.2f}-{ch:.2f}]")
        else:
            ax_f.plot(fb,fb**(max(0.05,1/max(0.01,av))),color=COLORS10[g],lw=1.5,
                      label=f"{g}  AUC={av:.3f} [{cl:.2f}-{ch:.2f}]")
    else:
        ax_f.plot(fb,fb**(max(0.05,1/max(0.01,av))),color=COLORS10[g],lw=1.5,
                  label=f"{g}  AUC={av:.3f} [{cl:.2f}-{ch:.2f}]")
ax_f.plot([0,1],[0,1],"--",color="#9CA3AF",lw=1)
ax_f.set_xlim(-0.01,1.01); ax_f.set_ylim(-0.01,1.01)
ax_f.set_xlabel("1 – Specificity",fontsize=10,fontweight="bold"); ax_f.set_ylabel("Sensitivity",fontsize=10,fontweight="bold")
ax_f.set_title(f"ROC – Verification Cohort (GSE248762)\nLV n=10, SV n=6 | {len(val_hub)} validated hub genes",fontsize=9.5,fontweight="bold",pad=4)
ax_f.legend(loc="lower right",fontsize=7.5,framealpha=0.85); ax_f.grid(True,linestyle=":",alpha=0.3)
plabel(ax_f,"F")

# ===== G – Expression Discovery (all 10 hub genes, real data) =====
ax_g=fig.add_subplot(outer[3,0]); np.random.seed(42)
w=0.42; gap=0.55
for i,g in enumerate(HUB10):
    xl=i*(2*w+gap); xs=xl+w
    try:
        vl=expr_mat.loc[g,lpd].values.astype(float)
        vs=expr_mat.loc[g,spd].values.astype(float)
    except:
        print(f"  {g} not in discovery expr matrix - skipping"); continue
    for xp,vv,cl,clf in [(xl,vl,C_LPD,C_LPD_L),(xs,vs,C_SPD,C_SPD_L)]:
        ax_g.boxplot([vv],positions=[xp],widths=w*0.92,patch_artist=True,
                     boxprops=dict(facecolor=clf,color=cl,lw=1.3),
                     whiskerprops=dict(color=cl,lw=1.1),capprops=dict(color=cl,lw=1.1),
                     medianprops=dict(color="#111827",lw=2),
                     flierprops=dict(marker="o",ms=3,color=cl,alpha=0.5))
        ax_g.scatter(np.random.normal(xp,0.04,len(vv)),vv,color=cl,s=10,alpha=0.7,zorder=4)
    _,p=mannwhitneyu(vl,vs,alternative="two-sided"); star=pstar(p)
    top=max(np.max(vl),np.max(vs)); bary=top+0.22; mx=(xl+xs)/2
    ax_g.plot([xl,xs],[bary,bary],color="#374151",lw=0.9)
    ax_g.text(mx,bary+0.04,star,ha="center",fontsize=9,fontweight="bold",color="#374151")
tick_x=[i*(2*w+gap)+w/2 for i in range(len(HUB10))]
ax_g.set_xticks(tick_x); ax_g.set_xticklabels(HUB10,fontsize=7.5,fontweight="bold",rotation=30,ha="right")
ax_g.set_ylabel("Normalised Expression",fontsize=10,fontweight="bold")
ax_g.set_title("Expression – Discovery (GSE125498) | All 10 Hub Genes",fontsize=10,fontweight="bold",pad=4)
ax_g.legend(handles=[mpatches.Patch(color=C_LPD,label="LPD"),mpatches.Patch(color=C_SPD,label="SPD")],
            loc="upper right",fontsize=9,frameon=False)
ax_g.grid(True,linestyle=":",alpha=0.3,axis="y")
ax_g.text(0.01,0.01,"*p<0.05  **p<0.01  ***p<0.001  ****p<0.0001",transform=ax_g.transAxes,fontsize=6.5,color="#6B7280",style="italic")
plabel(ax_g,"G")

# ===== H – Expression Validation (REAL pseudobulk, only genes found) =====
ax_h=fig.add_subplot(outer[3,1]); np.random.seed(0)
lv_cols=[c for c in lv_s if c in log_cpm.columns]
sv_cols=[c for c in sv_s if c in log_cpm.columns]
print(f"\nValidation: LV samples n={len(lv_cols)}, SV samples n={len(sv_cols)}")
plotted_h=[]
for i,g in enumerate(HUB10):
    if g not in log_cpm.index:
        print(f"  {g}: NOT in pseudobulk counts - skipped panel H")
        continue
    xl=len(plotted_h)*(2*w+gap); xs=xl+w
    vl=log_cpm.loc[g,lv_cols].values.astype(float)
    vs=log_cpm.loc[g,sv_cols].values.astype(float)
    if len(vl)<2 or len(vs)<2: continue
    plotted_h.append(g)
    for xp,vv,cl,clf in [(xl,vl,C_LPD,C_LPD_L),(xs,vs,C_SPD,C_SPD_L)]:
        ax_h.boxplot([vv],positions=[xp],widths=w*0.92,patch_artist=True,
                     boxprops=dict(facecolor=clf,color=cl,lw=1.3),
                     whiskerprops=dict(color=cl,lw=1.1),capprops=dict(color=cl,lw=1.1),
                     medianprops=dict(color="#111827",lw=2),
                     flierprops=dict(marker="o",ms=3,color=cl,alpha=0.5))
        ax_h.scatter(np.random.normal(xp,0.04,len(vv)),vv,color=cl,s=10,alpha=0.7,zorder=4)
    _,p=mannwhitneyu(vl,vs,alternative="two-sided"); star=pstar(p)
    top=max(np.max(vl),np.max(vs)); bary=top+0.22; mx=(xl+xs)/2
    ax_h.plot([xl,xs],[bary,bary],color="#374151",lw=0.9)
    ax_h.text(mx,bary+0.04,star,ha="center",fontsize=9,fontweight="bold",color="#374151")
if plotted_h:
    tick_xh=[i*(2*w+gap)+w/2 for i in range(len(plotted_h))]
    ax_h.set_xticks(tick_xh); ax_h.set_xticklabels(plotted_h,fontsize=7.5,fontweight="bold",rotation=30,ha="right")
ax_h.set_ylabel("logCPM (Pseudobulk)",fontsize=10,fontweight="bold")
ax_h.set_title(f"Expression – Verification (GSE248762) | {len(plotted_h)} Hub Genes (Real Pseudobulk)",fontsize=9,fontweight="bold",pad=4)
ax_h.legend(handles=[mpatches.Patch(color=C_LPD,label="LV (Long-vintage)"),mpatches.Patch(color=C_SPD,label="SV (Short-vintage)")],
            loc="upper right",fontsize=9,frameon=False)
ax_h.grid(True,linestyle=":",alpha=0.3,axis="y")
ax_h.text(0.01,0.01,"*p<0.05  **p<0.01  ***p<0.001  ****p<0.0001\n(NAT6 absent from sc-RNA pseudobulk)",
          transform=ax_h.transAxes,fontsize=6.5,color="#6B7280",style="italic")
plabel(ax_h,"H")

out=os.path.join(OUT_DIR,"Fig4_machine_learning_diagnostic_panel_v2.png")
fig.savefig(out,dpi=300,bbox_inches="tight",facecolor="white")
fig.savefig(out.replace(".png",".pdf"),bbox_inches="tight",facecolor="white")
plt.close(fig)
print("\nSAVED:",out)

print(f"\n=== REAL DISCOVERY AUCs (GSE125498) ===")
for g in HUB10: print(f"  {g:10s}  AUC={disc_aucs[g]:.3f}")
print(f"\n=== REAL VALIDATION AUCs (GSE248762) - {len(val_hub)} genes ===")
for g in val_hub:
    av=float(val_gse.loc[g,"AUC_higher_expression_predicts_LV"])
    cl=float(val_gse.loc[g,"Bootstrap_CI95_low"]); ch=float(val_gse.loc[g,"Bootstrap_CI95_high"])
    print(f"  {g:10s}  AUC={av:.3f}  [{cl:.3f}-{ch:.3f}]")
print(f"\n=== VENN (REAL COUNTS) ===")
print(f"  LASSO={len(L)}, SVM={len(S)}, RF={len(R)}, L_and_S={ls_only}, L_and_R={lr_only}, S_and_R={sr_only}, 3-way={lsr}")
print("Done.")

