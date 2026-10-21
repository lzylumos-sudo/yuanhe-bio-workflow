"""Prespecified production-window statistics; sample SD is not convergence evidence."""
from __future__ import annotations
import csv,json,math
from pathlib import Path
import numpy as np

def statistics(values):
    x=np.asarray(values,dtype=float)
    if len(x)<2 or not np.isfinite(x).all():raise ValueError('At least two finite observations required')
    n=len(x); centered=x-x.mean(); var=np.dot(centered,centered)/n
    tau=1.0
    if var>0:
        size=1 << (2*n-1).bit_length()
        f=np.fft.rfft(centered,n=size);ac=np.fft.irfft(f*f.conjugate(),n=size)[:n]
        ac=ac/(np.arange(n,0,-1)*var)
        # Initial-positive paired autocorrelation sequence.
        for i in range(1,n-1,2):
            pair=float(ac[i]+ac[i+1])
            if pair<=0:break
            tau+=2*pair
    neff=max(1.,min(float(n),n/tau));sd=float(np.std(x,ddof=1))
    return {'n':n,'mean':float(x.mean()),'sample_sd':sd,'naive_sem':sd/math.sqrt(n),
            'autocorrelation_time_frames':tau,'effective_samples':neff,'correlation_adjusted_sem':sd/math.sqrt(neff)}

def window_summary(rows, field='ligand_RMSD_A', window=(400.,500.), block_ns=10., production_ns=None):
    lo,hi=map(float,window)
    if not (0<=lo<hi) or block_ns<=0:raise ValueError('Invalid analysis window/block')
    times=np.array([r['time_ns'] for r in rows],dtype=float)
    if not len(times) or not np.isfinite(times).all() or np.any(np.diff(times)<=0):raise ValueError('Production time must be finite and strictly increasing')
    result={'window_ns':[lo,hi],'observable':'ligand heavy-atom RMSD after protein alignment','unit':'angstrom',
            'ddof':1,'observed_end_ns':float(times[-1]),'status':'insufficient_trajectory',
            'interpretation':'No binding-strength, activity or convergence verdict; may be unstable / require extension'}
    if production_ns is None:production_ns=times[-1]
    if production_ns<500-1e-6 or times[-1]<hi-1e-6:return result
    subset=[r for r in rows if lo-1e-9<=r['time_ns']<=hi+1e-9 and field in r]
    spacing=float(np.median(np.diff(times))) if len(times)>1 else math.inf
    if len(subset)<2 or subset[0]['time_ns']>lo+spacing+1e-6:return result
    result.update(status='computed_requires_review',**statistics([r[field] for r in subset]))
    blocks=[]
    for start in np.arange(lo,hi,block_ns):
        end=min(start+block_ns,hi)
        data=[r[field] for r in subset if start<=r['time_ns']<end or end==hi and r['time_ns']==hi]
        if len(data)>1:blocks.append({'start_ns':float(start),'end_ns':float(end),**statistics(data)})
    result['blocks']=blocks
    if len(blocks)>1:result['block_mean_sem']=float(np.std([b['mean'] for b in blocks],ddof=1)/np.sqrt(len(blocks)))
    return result

def compare_replicates(reports):
    complete=[r for r in reports if r.get('status')=='computed_requires_review']
    result={'replicates':reports,'required_replicates':len(reports),'complete_replicates':len(complete),'status':'incomplete' if len(complete)<len(reports) else 'requires_review',
            'interpretation':'Independent repeat means; never pooled into one 500 ns trajectory'}
    if len(complete)>1:result['between_repeat_means']=statistics([r['mean'] for r in complete])
    return result

def analyze(topology, trajectories, output, *, alignment='protein and backbone', ligand='resname LIG and not name H*',
            protein='protein', times_ns=None, window=(400,500), block_ns=10, reference=None):
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import minimize_vectors,distance_array
    from .trajectory import fit_to_reference,write_csv,plot_lines
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    u=mda.Universe(str(topology),list(map(str,trajectories)))
    align=u.select_atoms(alignment);prot=u.select_atoms(protein);lig=u.select_atoms(ligand)
    if len(align)<3 or not len(lig):raise ValueError('Empty ligand or insufficient stable protein selection')
    # Requires real topology bonds. Fragment reconstruction must not guess chemistry.
    if not len(prot.bonds) or not len(lig.bonds):raise ValueError('Bonded topology required for PBC repair')
    ref_a=ref_l=None;rows=[];contacts={};rmsf_mean=None;rmsf_m2=None
    for i,ts in enumerate(u.trajectory):
        if ts.dimensions is None or np.any(ts.dimensions[:3]<=0):raise ValueError('Periodic box missing')
        prot.unwrap(compound='fragments',reference=None,inplace=True)
        whole_lig=u.select_atoms('same residue as ('+ligand+')');whole_lig.unwrap(compound='fragments',reference=None,inplace=True)
        anchor=prot.fragments[0].center_of_geometry()
        for fragment in prot.fragments[1:]:
            delta=fragment.center_of_geometry()-anchor;fragment.positions+=minimize_vectors(delta[None,:],ts.dimensions)[0]-delta
        delta=whole_lig.center_of_geometry()-prot.center_of_geometry();whole_lig.positions+=minimize_vectors(delta[None,:],ts.dimensions)[0]-delta
        if ref_a is None:
            if reference:
                r=mda.Universe(str(reference));ref_a=r.select_atoms(alignment).positions.copy();ref_l=r.select_atoms(ligand).positions.copy()
            else:ref_a=align.positions.copy();ref_l=lig.positions.copy()
        rotation,center,target=fit_to_reference(align.positions,ref_a)
        a=(align.positions-center)@rotation+target;l=(lig.positions-center)@rotation+target
        if rmsf_mean is None:rmsf_mean=np.zeros_like(a);rmsf_m2=np.zeros_like(a)
        diff=a-rmsf_mean;rmsf_mean+=diff/(i+1);rmsf_m2+=diff*(a-rmsf_mean)
        d=distance_array(prot.positions,lig.positions)
        for residue in set(int(atom.resindex) for atom,dist in zip(prot,d.min(axis=1)) if dist<=4.5):contacts[residue]=contacts.get(residue,0)+1
        rows.append({'time_ns':float(times_ns[i] if times_ns is not None else ts.time/1000),
                     'protein_RMSD_A':float(np.sqrt(np.mean(np.sum((a-ref_a)**2,axis=1)))),
                     'ligand_RMSD_A':float(np.sqrt(np.mean(np.sum((l-ref_l)**2,axis=1)))),
                     'pocket_min_distance_A':float(d.min()),'protein_ligand_center_distance_A':float(np.linalg.norm(lig.center_of_geometry()-prot.center_of_geometry()))})
    write_csv(out/'timeseries.csv',rows)
    write_csv(out/'RMSF.csv',[{'atom_index':int(atom.index),'resid':int(atom.resid),'rmsf_A':float(value)} for atom,value in zip(align,np.sqrt(np.sum(rmsf_m2/len(rows),axis=1)))])
    result=window_summary(rows,window=window,block_ns=block_ns)
    result.update(alignment_selection=alignment,ligand_selection=ligand,reference=str(reference) if reference else 'first production frame',
                  contact_occupancy_4_5A={str(k):v/len(rows) for k,v in contacts.items()},pbc='bond-fragment unwrap, assemble protein chains, translate whole ligand to nearest protein image')
    # Hydrogen bonding uses topology charges/types rather than element-only guesses.
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
    hb=HydrogenBondAnalysis(u,between=[protein,'resname LIG'],d_a_cutoff=3.0,d_h_a_angle_cutoff=150)
    hb.run();events=hb.results.hbonds;counts={}
    for frame,donor,hydrogen,acceptor,dist,angle in events:
        key=(int(donor),int(acceptor));counts.setdefault(key,set()).add(int(frame))
    write_csv(out/'hydrogen_bonds.csv',[{'donor_atom':k[0],'acceptor_atom':k[1],'occupancy':len(v)/len(rows)} for k,v in counts.items()],['donor_atom','acceptor_atom','occupancy'])
    result['hydrogen_bonds']={'donor_acceptor_A':3.0,'angle_degrees':150,'chemical_definition':'MDAnalysis topology-charge donor/acceptor selection; review atom typing','events':len(events)}
    (out/'last100ns.json').write_text(json.dumps(result,indent=2)+'\n')
    for field in ['protein_RMSD_A','ligand_RMSD_A','pocket_min_distance_A']:plot_lines(out/(field+'.png'),[r['time_ns'] for r in rows],[r[field] for r in rows],field)
    u.trajectory.close();return result
