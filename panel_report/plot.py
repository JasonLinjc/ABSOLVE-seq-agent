"""Render a two-page Dedud panel comparison with pointwise intervals and manuscript estimates."""
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('panel',type=Path);ap.add_argument('output',type=Path)
    args=ap.parse_args()
    with args.panel.open() as f: records=[r for r in csv.DictReader(f,delimiter='\t') if r['strategy']=='Dedud']
    lookup={(r['target'],r['editor']):r for r in records}
    colors=['#20725E','#943C68']
    with PdfPages(args.output) as pdf:
        for guide in ('1450','1617'):
            targets=sorted({r['target'] for r in records if r['target'].startswith(guide+'_')})
            fig,axs=plt.subplots(1,2,sharey=True,figsize=(11.3,max(8,len(targets)*.21+2.4)))
            fig.subplots_adjust(left=.185,right=.97,top=.89,bottom=.13,wspace=.12)
            fig.suptitle(f'sg{guide}   |   Background-corrected editing across the panel',x=.05,y=.975,ha='left',fontsize=16,fontweight='bold')
            fig.text(.05,.932,'Manuscript Dedud  ·  8 matched replicates  ·  3 donors  ·  model estimates and 95% confidence intervals',fontsize=10)
            for ax,arm,color in zip(axs,('WT','HiFi'),colors):
                for i,t in enumerate(targets):
                    r=lookup[t,arm]; v,lo,hi,q=map(float,[r[k] for k in ('delta_pp','ci95_low_pp','ci95_high_pp','q_bh')])
                    if i%2==0: ax.axhspan(i-.5,i+.5,color='#f2f6f5',zorder=-2)
                    ax.plot([lo,hi],[i,i],c=color,lw=1.3)
                    ax.plot(v,i,'o',mfc=color if q<.05 else 'white',mec=color,ms=4)
                    if 'manuscript_delta_pp' in r:ax.plot(float(r['manuscript_delta_pp']),i-.19,'|',color='#3c4144',ms=5)
                ax.axvline(0,color='#849899',lw=.85)
                for x in (.01,.1):ax.axvline(x,color='#A77254',lw=.6,ls=':',zorder=-1)
                ax.set_xscale('symlog',linthresh=.01,linscale=1)
                ax.set_xlim(-2,120); ax.set_xticks([-.1,-.01,0,.01,.1,1,10,100],
                    ['−0.1','−0.01','0','0.01','0.1','1','10','100'])
                ax.set_title(('WT SpCas9' if arm=='WT' else 'HiFi Cas9')+' − NoEP',fontsize=11,color=color)
                ax.set_yticks(np.arange(len(targets)),[t.replace(guide+'_OT_','') for t in targets])
                ax.tick_params(labelsize=8,length=2); ax.set_ylim(len(targets)-.4,-.6)
                ax.set_xlabel('Δindel rate (percentage points, symmetric log scale)',fontsize=9)
                for spine in ('top','right','left'):ax.spines[spine].set_visible(False)
            for lab,t in zip(axs[0].get_yticklabels(),targets):
                if t in ('1450_OT_0071_REF','1617_OT_0040_ALT') or '_0000_' in t:lab.set_fontweight('bold')
            handles=[Line2D([0],[0],ls='',marker='o',color='#275a5a',label='BH q < 0.05'),
                Line2D([0],[0],ls='',marker='o',mfc='white',mec='#275a5a',label='BH q ≥ 0.05')]
            if 'manuscript_delta_pp' in records[0]:handles.append(Line2D([0],[0],ls='',marker='|',color='#3c4144',label='Submitted manuscript estimate'))
            fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.12,.016),frameon=False,ncol=3,fontsize=9)
            pdf.savefig(fig);plt.close(fig)
if __name__=='__main__':main()
