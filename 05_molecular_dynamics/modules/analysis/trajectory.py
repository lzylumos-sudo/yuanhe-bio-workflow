"""Streaming, explicit-role analysis of committed OpenMM production frames."""
from __future__ import annotations
import csv
import itertools
import shutil
from collections import Counter
from pathlib import Path
from ..utils.core import WorkflowError, read_json, atomic_json, duration_ps, digest


def fit_to_reference(x, reference):
    import numpy as np
    center = x.mean(axis=0)
    ref_center = reference.mean(axis=0)
    u, _, vt = np.linalg.svd((x-center).T @ (reference-ref_center))
    correction = np.eye(3)
    correction[-1,-1] = np.linalg.det(u @ vt)
    return u @ correction @ vt, center, ref_center


def write_csv(path, rows, fields=None):
    rows = list(rows)
    if not rows and not fields:
        raise ValueError("An empty table requires explicit columns.")
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_lines(path, x, y, ylabel, xlabel="Time (ns)"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with plt.rc_context({"path.simplify": False}):
        fig, ax = plt.subplots(figsize=(7,4), constrained_layout=True)
        ax.plot(x, y, lw=1)
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.grid(alpha=.2)
        fig.savefig(path, dpi=180)
        plt.close(fig)


def resolve_groups(config, meta):
    """Explicit roles when provided; otherwise topology chain labels only."""
    import numpy as np
    protein=set(meta['protein_indices'])
    chains=meta['chain_ids']
    definitions=config['analysis'].get('groups',{})
    groups={}
    if definitions:
        for role,definition in definitions.items():
            if role in ('ligand','glue'): continue
            selected=definition.get('chains',[])
            if not selected: raise WorkflowError('Analysis role needs a nonempty chains list.')
            ids=sorted(i for i in protein if chains[i] in selected)
            if not ids: raise WorkflowError('Analysis role does not match a protein chain.')
            if any(set(ids)&set(v) for v in groups.values()): raise WorkflowError('Analysis roles overlap.')
            groups[role]=ids
        if set(i for ids in groups.values() for i in ids)!=protein:
            raise WorkflowError('Assign every protein chain to an analysis role.')
    elif config.get('system_type')=='protein_ligand':
        groups['protein']=sorted(protein)
    else:
        for chain in sorted({chains[i] for i in protein}):
            groups['protein_'+(chain or 'blank')]=sorted(i for i in protein if chains[i]==chain)
    if meta.get('ligand_indices'): groups['ligand']=list(meta['ligand_indices'])
    if len(groups)<2: raise WorkflowError('Complex interface analysis requires at least two components.')
    return {role:np.asarray(ids,dtype=int) for role,ids in groups.items()}


def sasa_area(positions, elements, probe=1.4, points=100):
    """Heavy-atom Shrake-Rupley SASA, Angstrom squared; solvent is omitted."""
    import numpy as np
    from Bio.PDB import Structure, Model, Chain, Residue, Atom
    from Bio.PDB.SASA import ShrakeRupley
    structure = Structure.Structure("sasa")
    model, chain, residue = Model.Model(0), Chain.Chain("A"), Residue.Residue((" ",1," "), "MOL", " ")
    structure.add(model)
    model.add(chain)
    chain.add(residue)
    for i,(position,element) in enumerate(zip(positions,elements)):
        name = str(i)
        residue.add(Atom.Atom(name, np.asarray(position,dtype=float), 0., 1., " ", name, i, element=str(element).upper()))
    ShrakeRupley(probe_radius=probe,n_points=points).compute(structure,level="S")
    return float(structure.sasa)


def pair_sasa(positions, elements, group_a, group_b, probe=1.4, points=100):
    import numpy as np
    group_a,group_b = np.asarray(group_a,dtype=int),np.asarray(group_b,dtype=int)
    elements = np.asarray(elements)
    area_a = sasa_area(positions[group_a],elements[group_a],probe,points)
    area_b = sasa_area(positions[group_b],elements[group_b],probe,points)
    together = np.concatenate((group_a,group_b))
    area_ab = sasa_area(positions[together],elements[together],probe,points)
    buried = area_a+area_b-area_ab
    return {"SASA_A_A2":area_a,"SASA_B_A2":area_b,"SASA_AB_A2":area_ab,
            "buried_SASA_A2":buried,"interface_area_A2":buried/2}


def orientation_degrees(vector, reference):
    import numpy as np
    denominator = np.linalg.norm(vector)*np.linalg.norm(reference)
    if denominator < 1e-8:
        raise WorkflowError("Orientation atoms must define a nonzero vector.")
    return float(np.rad2deg(np.arccos(np.clip(np.dot(vector,reference)/denominator,-1,1))))


def plot_contact_map(path, counts, frames):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    interfaces = sorted({key[0] for key in counts})
    if not interfaces:
        fig,ax = plt.subplots(figsize=(7,3))
        ax.text(.5,.5,"No inter-group contacts at the configured cutoff",ha="center",va="center")
        ax.set_axis_off()
    else:
        fig,axes = plt.subplots(len(interfaces),1,figsize=(9,max(4,4*len(interfaces))),squeeze=False,constrained_layout=True)
        for ax,interface in zip(axes[:,0],interfaces):
            selected = {key:value for key,value in counts.items() if key[0]==interface}
            row_labels = sorted({key[1] for key in selected})
            column_labels = sorted({key[2] for key in selected})
            matrix = np.zeros((len(row_labels),len(column_labels)))
            rows = {label:i for i,label in enumerate(row_labels)}
            columns = {label:i for i,label in enumerate(column_labels)}
            for (_,a,b),value in selected.items():
                matrix[rows[a],columns[b]] = value/frames
            display = ax.imshow(matrix,vmin=0,vmax=1,aspect="auto",origin="lower",cmap="viridis")
            # All values stay in CSV; limit labels in large interfaces for readability.
            ri = np.arange(0,len(row_labels),max(1,len(row_labels)//20))
            ci = np.arange(0,len(column_labels),max(1,len(column_labels)//20))
            ax.set_yticks(ri,[row_labels[i] for i in ri],fontsize=7)
            ax.set_xticks(ci,[column_labels[i] for i in ci],rotation=90,fontsize=7)
            ax.set_title(interface+" residue-pair contact occupancy")
            fig.colorbar(display,ax=ax,label="Fraction of analyzed frames")
    fig.savefig(path,dpi=180,bbox_inches="tight")
    plt.close(fig)


def charged_atoms(meta, atom_names):
    """Conservative standard-sidechain ions; no histidine/ligand charge guessing."""
    from collections import defaultdict
    hydrogens = defaultdict(list)
    for a,b in meta["bonds"]:
        if meta["elements"][a] == "H":
            hydrogens[b].append(a)
        if meta["elements"][b] == "H":
            hydrogens[a].append(b)
    residues = defaultdict(list)
    for i in meta["protein_heavy"]:
        residues[meta["residues"][i]].append(i)
    charge = {}
    for label,indices in residues.items():
        name = label.split(":",1)[-1][:3]
        if name in ("ASP","GLU"):
            terminal = [i for i in indices if atom_names[i] in ("OD1","OD2","OE1","OE2")]
            if len(terminal)==2 and not any(hydrogens[i] for i in terminal):
                charge.update({i:-1 for i in terminal})
        elif name == "LYS":
            terminal = [i for i in indices if atom_names[i]=="NZ" and len(hydrogens[i])==3]
            charge.update({i:1 for i in terminal})
        elif name == "ARG":
            terminal = [i for i in indices if atom_names[i] in ("NE","NH1","NH2")]
            if len(terminal)==3 and sum(len(hydrogens[i]) for i in terminal)==5:
                charge.update({i:1 for i in terminal})
    return charge


def analyze_replica(config, paths, replica):
    import numpy as np
    import pandas as pd
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import distance_array, minimize_vectors, calc_angles, capped_distance
    settings = config["analysis"]
    source = paths["md"]/f"replica_{replica:02d}"
    out = paths["analysis"]/f"replica_{replica:02d}"
    out.mkdir(parents=True,exist_ok=True)
    progress = read_json(source/"progress.json")
    if progress["phase"] != "complete":
        raise WorkflowError("Analysis requires a completed replica.")
    segments = progress["segments"]
    for segment in segments:
        for key in ("dcd","csv"):
            if digest(source/segment[key]) != segment[key+"_sha256"]:
                raise WorkflowError("Trajectory/state checksum mismatch.")
    base = paths["md"]/"system"
    meta = read_json(base/"topology.json")
    u = mda.Universe(str(base/"topology.pdb"),[str(source/s["dcd"]) for s in segments])
    if u.atoms.n_atoms != meta["atoms"]:
        raise WorkflowError("Topology/trajectory atom count mismatch.")
    if hasattr(u.atoms,"bonds"):
        u.delete_bonds(u.bonds)
    u.add_TopologyAttr("bonds",meta["bonds"])
    protein = u.atoms[meta["protein_indices"]]
    ligand = u.atoms[meta.get("ligand_indices",[])]
    ca = np.asarray(meta["ca_indices"],dtype=int)
    fit_ids = u.select_atoms(settings.get("alignment_selection","protein and backbone")).indices
    if len(fit_ids)<3: raise WorkflowError("Stable protein alignment selection has fewer than three atoms")
    reference_fit = None
    heavy_l = np.asarray(meta.get("ligand_heavy",[]),dtype=int)
    residues = meta["residues"]
    elements = np.asarray(meta["elements"])
    groups = resolve_groups(config,meta)
    heavy_groups = {name:ids[elements[ids]!="H"] for name,ids in groups.items()}
    interfaces = list(itertools.combinations(groups,2))
    membership = {int(i):name for name,ids in groups.items() for i in ids}
    donor_h = {}
    for donor,hydrogen in meta["donor_hydrogen_pairs"]:
        donor_h.setdefault(donor,[]).append(hydrogen)
    donors = np.asarray(list(donor_h),dtype=int)
    acceptors = np.asarray(meta["acceptors"],dtype=int)
    charged = charged_atoms(meta,u.atoms.names)
    custom_distances = []
    for entry in settings.get("distance_pairs",[]):
        a = u.select_atoms(entry["selection1"]).indices
        b = u.select_atoms(entry["selection2"]).indices
        if not len(a) or not len(b):
            raise WorkflowError("A configured distance selection is empty.")
        custom_distances.append((entry["name"],a,b))
    orientation_ids = settings.get("glue_orientation_atom_indices")
    if orientation_ids is not None:
        if len(orientation_ids)!=2 or len(set(orientation_ids))!=2 or any(type(i) is not int or i<0 or i>=len(heavy_l) for i in orientation_ids):
            raise WorkflowError("glue_orientation_atom_indices requires two distinct ligand-local heavy-atom indices.")
        orientation_ids = heavy_l[np.asarray(orientation_ids)]
    reference_orientation = None
    interval_ns = duration_ps(config["MD"]["report_interval"])/1000
    expected_frames = round(progress["production_ns"]/interval_ns)
    if len(u.trajectory)!=expected_frames:
        raise WorkflowError("Missing or unexpected production frames.")
    rows,sasa_rows = [],[]
    contact_counts,pair_contact_counts,hbond_counts,salt_counts = Counter(),Counter(),Counter(),Counter()
    reference_ca = reference_lig = mean_ca = m2_ca = None
    count = 0
    sasa_stride = settings.get("sasa_stride",10)
    if type(sasa_stride) is not int or sasa_stride < 1:
        raise WorkflowError("analysis.sasa_stride must be a positive integer.")
    if type(settings.get("sasa_points",100)) is not int or settings.get("sasa_points",100) < 20:
        raise WorkflowError("analysis.sasa_points must be an integer of at least 20.")
    if not 0 < settings.get("sasa_probe_angstrom",1.4) < 10:
        raise WorkflowError("analysis.sasa_probe_angstrom must be positive and below 10 A.")
    for frame_index,ts in enumerate(u.trajectory):
        time_ns = (frame_index+1)*interval_ns
        if time_ns < settings["discard_ns"] or frame_index % settings["stride"]:
            continue
        if not np.isfinite(ts.positions).all() or ts.dimensions is None or not np.isfinite(ts.dimensions).all():
            raise WorkflowError("Non-finite coordinates or missing periodic box.")
        raw,box = ts.positions.copy(),ts.dimensions.copy()
        if min(box[:3])<=0:
            raise WorkflowError("Invalid periodic box.")
        protein.unwrap(compound="fragments",reference=None,inplace=True)
        if len(ligand):
            ligand.unwrap(compound="fragments",reference=None,inplace=True)
        fragments = protein.fragments
        anchor = fragments[0].center_of_geometry()
        for fragment in fragments[1:]:
            delta = fragment.center_of_geometry()-anchor
            fragment.positions += minimize_vectors(delta[None,:],box)[0]-delta
        if len(ligand):
            delta = ligand.center_of_geometry()-protein.center_of_geometry()
            ligand.positions += minimize_vectors(delta[None,:],box)[0]-delta
        positions = u.atoms.positions.copy()
        if reference_ca is None:
            reference_ca = positions[ca].copy()
        if reference_fit is None: reference_fit=positions[fit_ids].copy()
        rotation,center,ref_center = fit_to_reference(positions[fit_ids],reference_fit)
        aligned_ca = (positions[ca]-center) @ rotation + ref_center
        aligned_ligand = (positions[heavy_l]-center) @ rotation + ref_center
        if reference_lig is None and len(heavy_l):
            reference_lig = aligned_ligand.copy()
        count += 1
        if mean_ca is None:
            mean_ca,m2_ca = aligned_ca.copy(),np.zeros_like(aligned_ca)
        else:
            delta_ca = aligned_ca-mean_ca
            mean_ca += delta_ca/count
            m2_ca += delta_ca*(aligned_ca-mean_ca)
        row = {"time_ns":time_ns}
        if settings["rmsd"]:
            row["protein_CA_RMSD_A"] = float(np.sqrt(np.mean(np.sum((aligned_ca-reference_ca)**2,axis=1))))
            if len(heavy_l):
                row["ligand_RMSD_A"] = float(np.sqrt(np.mean(np.sum((aligned_ligand-reference_lig)**2,axis=1))))
        frame_residues,frame_pairs,frame_salt = set(),set(),set()
        for a_name,b_name in interfaces:
            a,b = heavy_groups[a_name],heavy_groups[b_name]
            interface = a_name+"_"+b_name
            distances = distance_array(raw[a],raw[b],box=box)
            if settings["distance"]:
                row["minimum_"+interface+"_distance_A"] = float(distances.min())
                delta = positions[a].mean(axis=0)-positions[b].mean(axis=0)
                row["COG_"+interface+"_distance_A"] = float(np.linalg.norm(minimize_vectors(delta[None,:],box)[0]))
            if settings["contacts"]:
                for i,j in zip(*np.where(distances<=settings["contact_cutoff_angstrom"])):
                    frame_pairs.add((interface,residues[a[i]],residues[b[j]]))
                    frame_residues.update((residues[a[i]],residues[b[j]]))
            for i,j in zip(*np.where(distances<=settings.get("salt_bridge_cutoff_angstrom",4.0))):
                if charged.get(int(a[i]),0)*charged.get(int(b[j]),0)==-1:
                    frame_salt.add((interface,residues[a[i]],residues[b[j]]))
            if settings.get("sasa",True) and (count-1)%sasa_stride==0:
                values = pair_sasa(positions,elements,a,b,settings.get("sasa_probe_angstrom",1.4),settings.get("sasa_points",100))
                sasa_rows.append({"time_ns":time_ns,"interface":interface,**values})
        if settings["distance"]:
            for name,a,b in custom_distances:
                row["distance_"+name+"_A"] = float(distance_array(raw[a],raw[b],box=box).min())
        if settings["contacts"]:
            contact_counts.update(frame_residues)
            pair_contact_counts.update(frame_pairs)
            row["contact_residues"] = len(frame_residues)
        salt_counts.update(frame_salt)
        row["salt_bridge_residue_pairs"] = len(frame_salt)
        if orientation_ids is not None:
            vector = (positions[orientation_ids[1]]-positions[orientation_ids[0]]) @ rotation
            if reference_orientation is None:
                reference_orientation = vector.copy()
            row["glue_orientation_change_degrees"] = orientation_degrees(vector,reference_orientation)
        if settings["hbond"]:
            occupied = set()
            if len(donors) and len(acceptors):
                near = capped_distance(raw[donors],raw[acceptors],settings["hbond_distance_angstrom"],box=box,return_distances=False)
                triples = []
                for di,ai in near:
                    donor,acceptor = int(donors[di]),int(acceptors[ai])
                    if membership.get(donor) is not None and membership.get(acceptor) is not None and membership[donor]!=membership[acceptor]:
                        triples.extend((donor,h,acceptor) for h in donor_h[donor])
                if triples:
                    candidates = np.asarray(triples,dtype=int)
                    angles = np.rad2deg(calc_angles(raw[candidates[:,0]],raw[candidates[:,1]],raw[candidates[:,2]],box=box))
                    for donor,hydrogen,acceptor in candidates[angles>=settings["hbond_angle_degrees"]]:
                        occupied.add((int(donor),int(acceptor)))
            hbond_counts.update(occupied)
            row["hydrogen_bonds"] = len(occupied)
        rows.append(row)
    if count<2:
        raise WorkflowError("At least two analyzed frames are required.")
    write_csv(out/"timeseries.csv",rows)
    distance_fields = ["time_ns"]+[key for key in rows[0] if "distance_A" in key]
    write_csv(out/"distance.csv",[{key:row[key] for key in distance_fields} for row in rows],distance_fields)
    if settings["rmsf"]:
        rmsf = np.sqrt(np.maximum(m2_ca/count,0).sum(axis=1))
        write_csv(out/"rmsf.csv",[{"residue":residues[i],"CA_RMSF_A":float(v)} for i,v in zip(ca,rmsf)])
        plot_lines(out/"RMSF.png",np.arange(1,len(ca)+1),rmsf,"C-alpha RMSF (A)","C-alpha index")
    fields = {
        "RMSD.png":("protein_CA_RMSD_A","Protein C-alpha RMSD (A)"),
        "Ligand_RMSD.png":("ligand_RMSD_A","Ligand heavy-atom RMSD (A)"),
        "Hbond.png":("hydrogen_bonds","Inter-group hydrogen bonds"),
        "contacts.png":("contact_residues","Inter-group contacted residues"),
        "glue_orientation.png":("glue_orientation_change_degrees","Defined glue axis change (degrees)")}
    distance_column = next((key for key in rows[0] if key.startswith("minimum_") and key.endswith("_distance_A")),None)
    if distance_column:
        fields["distance.png"] = (distance_column,"Minimum inter-group heavy-atom distance (A)")
    for name,(column,label) in fields.items():
        if column in rows[0]:
            plot_lines(out/name,[r["time_ns"] for r in rows],[r[column] for r in rows],label)
    if (out/"RMSD.png").exists():
        shutil.copyfile(out/"RMSD.png",out/"protein_RMSD.png")
    if not len(heavy_l):
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(7,4)); ax.axis('off')
        ax.text(.5,.5,'NOT APPLICABLE: no small-molecule ligand',ha='center',va='center')
        fig.savefig(out/'Ligand_RMSD.png',dpi=180); plt.close(fig)
    contact_rows = [{"residue":residue,"frames":n,"occupancy":n/count} for residue,n in contact_counts.most_common()]
    write_csv(out/"contact_occupancy.csv",contact_rows,["residue","frames","occupancy"])
    contact_pairs = [{"interface":interface,"residue_A":a,"residue_B":b,"frames":n,"occupancy":n/count}
                     for (interface,a,b),n in pair_contact_counts.most_common()]
    write_csv(out/"contact_map.csv",contact_pairs,["interface","residue_A","residue_B","frames","occupancy"])
    plot_contact_map(out/"Contact_map.png",pair_contact_counts,count)
    hb_rows = [{"donor_index":d,"acceptor_index":a,"donor_group":membership[d],"acceptor_group":membership[a],
                "donor_residue":residues[d],"acceptor_residue":residues[a],"frames":n,"occupancy":n/count}
               for (d,a),n in hbond_counts.most_common()]
    write_csv(out/"hbond_occupancy.csv",hb_rows,["donor_index","acceptor_index","donor_group","acceptor_group","donor_residue","acceptor_residue","frames","occupancy"])
    shutil.copyfile(out/"hbond_occupancy.csv",out/"Hbond.csv")
    write_csv(out/"salt_bridge_occupancy.csv",[{"interface":interface,"residue_A":a,"residue_B":b,"frames":n,"occupancy":n/count}
              for (interface,a,b),n in salt_counts.most_common()],["interface","residue_A","residue_B","frames","occupancy"])
    write_csv(out/"SASA.csv",sasa_rows,["time_ns","interface","SASA_A_A2","SASA_B_A2","SASA_AB_A2","buried_SASA_A2","interface_area_A2"])
    atomic_json(out/"optional_analysis.json",{
        "MMGBSA":{"status":"PENDING","reason":"AmberTools backend runs after geometric analysis."},
        "glue_orientation":{"status":"COMPUTED" if orientation_ids is not None else "NOT_COMPUTED",
           "reason":"Researcher-defined directed atom pair after whole-protein C-alpha alignment." if orientation_ids is not None else "No human-defined glue_orientation_atom_indices; no automatic orientation choice."},
        "SASA":{"status":"COMPUTED" if sasa_rows else "NOT_COMPUTED","frames_per_interface":len(sasa_rows)//len(interfaces)}})
    state_df = pd.concat([pd.read_csv(source/s["csv"]) for s in segments],ignore_index=True)
    numeric = state_df.select_dtypes(include="number")
    finite = bool(np.isfinite(numeric.to_numpy()).all())
    if not finite:
        raise WorkflowError("Non-finite thermodynamic observables.")
    summary = {"replica":replica,"frames":count,"production_ns":progress["production_ns"],
               "analyzed_start_ns":rows[0]["time_ns"],"analyzed_end_ns":rows[-1]["time_ns"],"finite_observables":finite}
    for key in rows[0]:
        if key!="time_ns":
            summary[key+"_mean"] = float(np.mean([r[key] for r in rows]))
            summary[key+"_std"] = float(np.std([r[key] for r in rows],ddof=1))
    for column in state_df:
        if "Temperature" in column:
            summary["temperature_mean_K"] = float(state_df[column].mean())
        if "Density" in column:
            summary["density_mean_g_mL"] = float(state_df[column].mean())
    key_residues = [r for r in contact_rows if r["occupancy"]>=settings["key_residue_occupancy"]]
    from .research import window_summary
    final_window=window_summary(rows,window=settings.get("formal_window_ns",[400,500]),block_ns=settings.get("block_ns",10),production_ns=progress.get("production_ns"))
    final_window.update(alignment_selection=settings.get("alignment_selection","protein and backbone"), reference="first analyzed production frame", ligand_selection="topology ligand heavy indices", atom_indices=heavy_l.tolist())
    atomic_json(out/"last100ns.json",final_window)
    atomic_json(out/"summary.json",summary)
    atomic_json(out/"quality.json",{
        "frames_expected":expected_frames,"frames_observed":len(u.trajectory),"frames_analyzed":count,
        "finite_observables":finite,"key_residues":key_residues,"system_type":config.get("system_type","protein_ligand"),
        "groups":{key:{"atoms":len(ids)} for key,ids in groups.items()},
        "reference":"First analyzed production frame per replica; all protein C-alpha least-squares fit.",
        "ligand_rmsd":"Heavy atoms after whole-protein fit; no ligand self-fit or symmetry correction; nearest periodic image." if len(heavy_l) else "NOT_APPLICABLE: no ligand.",
        "rmsf":"C-alpha fluctuation around its time-average after fitting to the reference.",
        "distance":"Minimum heavy-atom and geometric-center distances, minimum image; centers are not mass weighted.",
        "hbond":"Standard protein atom types and RDKit ligand Donor/Acceptor features; inter-role D-A distance and D-H-A angle.",
        "salt_bridge":"Inter-role standard deprotonated ASP/GLU oxygen vs fully protonated LYS/ARG nitrogen, 4 A by default; histidine/ligand ions excluded.",
        "occupancy":"Fraction of all analyzed frames; residue-pair contacts/salt bridges and donor-acceptor hydrogen bonds deduplicated per frame.",
        "SASA":{"method":"Biopython Shrake-Rupley; heavy atoms only; no solvent; reconstructed nearest-image groups.",
            "probe_angstrom":settings.get("sasa_probe_angstrom",1.4),"sphere_points":settings.get("sasa_points",100),
            "stride_analyzed_frames":sasa_stride,"buried_SASA":"SASA(A)+SASA(B)-SASA(A+B); interface area is half.",
            "ternary":"Each pair is computed in isolation from the third component; pairwise buried areas are not additive."},
        "limitations":"No automatic convergence, affinity, correct-pose, recruitment or degradation claim. Review protonation, PBC assembly and replica disagreement."})
    u.trajectory.close()
    return summary


def analyze(config, paths):
    paths["analysis"].mkdir(parents=True,exist_ok=True)
    summaries = [analyze_replica(config,paths,i) for i in range(1,config["replicates"]+1)]
    from .research import compare_replicates
    atomic_json(paths["analysis"]/"repeat_comparison.json",compare_replicates([read_json(paths["analysis"]/f"replica_{i:02d}/last100ns.json") for i in range(1,config["replicates"]+1)]))
    write_csv(paths["analysis"]/"summary.csv",summaries)
    atomic_json(paths["analysis"]/"analysis.json",{"replicates":len(summaries),"status":"complete",
                "MMGBSA":"PENDING","system_type":config.get("system_type","protein_ligand")})
