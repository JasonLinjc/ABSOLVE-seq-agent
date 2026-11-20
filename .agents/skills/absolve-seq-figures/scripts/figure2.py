"""Redraw Figure 2 panels with verified, panel-specific source definitions."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from manuscript_figures import plt, digest, tsv, unique, ORDER, STRATEGIES, CONTROLS

ARMS = ['NoEP','WTCas9','HiFi']
COLORS = ['#808084','#369046','#ce50a0']
METHODS = dict(zip(['raw','raw_dedup','dedudS_filtered_v3','dedudS_filtered_v3_dedup'],STRATEGIES))
PUBLISHED = '2a0eb6d1a58652f234e49e47cf9f72807875d3287c40b5c8fc70a36c3b4effad'


def swarm(ax, values, x, color, outlined, log=False):
    """Separate nearby target markers deterministically in display units."""
    points = []
    for target,value in sorted(values,key=lambda pair:pair[1]):
        y = np.log10(value) if log else value
        if not np.isfinite(y):
            raise ValueError('Invalid target measurement: '+target)
        step = .043 if log else .016
        offsets = [0]
        for n in range(1, len(points)+2):
            offsets.extend([n,-n])
        offset = next(i for i in offsets if all(i!=old or abs(y-prev)>step for old,prev in points))
        points.append((offset,y))
        control = target in CONTROLS
        edge = '#d92935' if control else color
        ax.scatter(x+offset*.045,value,s=14,marker='^' if control else 'o',
                   facecolors='none' if outlined else edge,edgecolors=edge,lw=.7,zorder=3)


def load_source(panel, source, targets):
    hashes = {}
    if panel=='2d':
        if (source/'source.sha256').read_text().strip()!=PUBLISHED:
            raise ValueError('Figure 2d requires the original published Maxi index from setup_references.py')
        rows=[]
        for t in targets:
            path=source/(t+'.tsv')
            data=path.read_bytes()
            hashes[path.name]=hashlib.sha256(data).hexdigest()
            flags=[line.rsplit(b'\t',1)[-1] for line in data.splitlines()]
            if set(flags)-{b'0',b'1'}:
                raise ValueError('Invalid Maxi quality flag')
            rows.extend([dict(target=t,strategy='Unfiltered',value=len(flags)),
                         dict(target=t,strategy='High quality',value=flags.count(b'1'))])
        return rows, hashes
    hashes[source.name]=digest(source)
    if panel=='2a':
        rows=tsv(source)
        expected={(t,arm,s) for t in targets for arm in ARMS for s in STRATEGIES}
        if unique(rows,['target','arm','strategy']).keys()!=expected:
            raise ValueError('Expected 948 original Figure 2a coverage cells')
        return [dict(target=r['target'],arm=r['arm'],strategy=r['strategy'],value=int(r['coverage'])) for r in rows], hashes
    with source.open() as f:
        rows=[dict(target=r['Target'],strategy=METHODS[r['variable'].removeprefix('cInd_').removesuffix('/cTot')],
                   value=float(r['value'])) for r in csv.DictReader(f)]
    if unique(rows,['target','strategy']).keys()!={(t,s) for t in targets for s in STRATEGIES}:
        raise ValueError('Expected 316 original Figure 2b observations')
    return rows,hashes


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('panel',choices=['2a','2b','2d'])
    p.add_argument('source',type=Path)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    targets=[r['target'] for r in tsv(ORDER)]
    rows,hashes=load_source(args.panel,args.source,targets)
    if args.out.exists():
        raise FileExistsError(args.out)
    if args.panel=='2a':
        fig,axes=plt.subplots(1,4,figsize=(10,4),sharey=True,layout='constrained')
        for ax,s in zip(axes,STRATEGIES):
            groups=[[r for r in rows if r['strategy']==s and r['arm']==arm] for arm in ARMS]
            for i,group in enumerate(groups):
                swarm(ax,[(r['target'],r['value']) for r in group],i,COLORS[i],'Dedup' in s,True)
            ax.set_title(s)
            ax.set_yscale('log')
            ax.set_ylim(1e2,1.25e7)
            ax.set_ylabel('pBC observations' if 'Dedup' in s else 'Reads')
            ax.set_xticks(range(3),['Unedited','SpCas9','HiFi'],rotation=35)
    else:
        fig,ax=plt.subplots(figsize=(5.8,4.2),layout='constrained')
        axes=[ax]
        names=STRATEGIES if args.panel=='2b' else ['Unfiltered','High quality']
        for i,s in enumerate(names):
            values=[(r['target'],r['value']) for r in rows if r['strategy']==s]
            if args.panel=='2b':
                violin=ax.violinplot([v for _,v in values],[i],widths=.86,showextrema=False)
                for body in violin['bodies']:
                    body.set(facecolor='none',edgecolor='.52',linewidth=.75,alpha=.7)
            swarm(ax,values,i,COLORS[0],args.panel=='2d' or 'Dedup' in s,args.panel=='2d')
        if args.panel=='2d':
            ax.set_yscale('log')
            ax.set_ylim(2e3,1.25e6)
        else:
            ax.set_ylim(-.08,1.18)
        ax.set_ylabel('Unedited sample indels (%)' if args.panel=='2b' else 'Plasmid pBCs per target')
        ax.set_xticks(range(len(names)),names,rotation=25)
    for ax in axes:
        ax.grid(axis='y',lw=.45,color='.86')
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Figure '+args.panel+'   |   '+{'2a':'Effective data per candidate','2b':'Background indels','2d':'Plasmid barcode yield'}[args.panel])
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=args.out.parent) as tmp:
        dest=Path(tmp)/'result'
        dest.mkdir()
        for ext in ['pdf','png']:
            fig.savefig(dest/('Figure'+args.panel+'.'+ext),dpi=180)
        plt.close(fig)
        with (dest/'plotted_data.tsv').open('w') as f:
            w=csv.DictWriter(f,rows[0],delimiter='\t',lineterminator='\n')
            w.writeheader()
            w.writerows(rows)
        (dest/'provenance.json').write_text(json.dumps(dict(figure=args.panel,rows=len(rows),source_sha256=hashes),indent=2)+'\n')
        dest.rename(args.out)


if __name__=='__main__':
    main()
