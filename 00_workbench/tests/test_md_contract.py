from pathlib import Path
import sys,json
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scientific.md_cli import export_bundle,validate

@pytest.mark.parametrize('problem,status', [
 (RuntimeError('parameterization failed'), 'failed'),
 (KeyboardInterrupt(), 'cancelled'),
 (SystemExit(143), 'cancelled'),
])
def test_membrane_construction_records_failure_and_preserves_evidence(tmp_path,problem,status):
 from scientific.membrane import construction_state
 original={'id':'MD_example','status':'building','source_hashes':{'ligand':'original'},'parameters':{'production_ns':.01}}
 path=tmp_path/'system.json';path.write_text(json.dumps(original))
 with pytest.raises(type(problem)) as captured:
  with construction_state(tmp_path):raise problem
 record=json.loads(path.read_text())
 assert captured.value is problem
 assert record['status']==status and record['error_type']==type(problem).__name__
 assert record['source_hashes']==original['source_hashes']
 assert record['parameters']==original['parameters']
 assert record['ended']>0

def test_membrane_construction_preserves_successful_record(tmp_path):
 from scientific.membrane import construction_state
 path=tmp_path/'system.json';path.write_text(json.dumps({'status':'building'}))
 expected={'status':'built','atom_count':100,'files':{'system.top':'verified'}}
 with construction_state(tmp_path):path.write_text(json.dumps(expected))
 assert json.loads(path.read_text())==expected

def test_handoff_preserves_identity_and_rejects_modified_input(tmp_path):
 r=tmp_path/'receptor.pdb';r.write_text('ATOM source')
 l=tmp_path/'ligand.sdf';l.write_text('source sdf')
 candidate=dict(compound_id='CMP1',pose_id='POS1',smiles='CCO',qualification='preliminary',molecule_id='generated1')
 bundle=export_bundle(tmp_path/'bundle',candidate,r,l,{'exploratory_ack':True},{'seed':17},{'run_id':'source'})
 assert validate(bundle)['candidate']==candidate
 (bundle/'inputs/ligand.sdf').write_text('changed')
 with pytest.raises(ValueError,match='checksum'):validate(bundle)

def test_handoff_refuses_path_escape_and_unacknowledged_candidate(tmp_path):
 r=tmp_path/'receptor.pdb';r.write_text('ATOM');l=tmp_path/'ligand.sdf';l.write_text('sdf')
 candidate=dict(molecule_id='generated1',compound_id='CMP1',pose_id='POS1',smiles='CCO',qualification='preliminary')
 with pytest.raises(ValueError,match='Exploratory'):export_bundle(tmp_path/'bundle',candidate,r,l,{}, {}, {})
 m=json.loads((tmp_path/'bundle/handoff.json').read_text());m['artifacts']['ligand']['path']='../ligand.sdf';(tmp_path/'bundle/handoff.json').write_text(json.dumps(m))
 with pytest.raises(ValueError,match='escapes'):validate(tmp_path/'bundle')

def test_rigid_water_export_preserves_settle_geometry(tmp_path):
 openmm=pytest.importorskip('openmm');parmed=pytest.importorskip('parmed')
 from openmm import app,unit
 import io
 pdb=app.PDBFile(io.StringIO('ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\nATOM      2  H1  HOH A   1       0.957   0.000   0.000  1.00  0.00           H\nATOM      3  H2  HOH A   1      -0.240   0.927   0.000  1.00  0.00           H\nEND\n'))
 ff=app.ForceField('amber14/tip3p.xml')
 system=ff.createSystem(pdb.topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
 structure=parmed.openmm.load_topology(pdb.topology,system,xyz=pdb.positions)
 structure.save(str(tmp_path/'water.top'))
 text=(tmp_path/'water.top').read_text()
 assert '[ settles ]' in text
 assert '0.09572000' in text and '0.15139' in text


def test_production_protocol_and_early_stop_detection(tmp_path):
 from scientific.membrane import mdp,check_completion
 p={'temperature_K':310,'equilibration_ns':1,'production_ns':1}
 recipe=mdp('production',p,17)
 assert 'nsteps = 500000' in recipe and 'pcoupltype = semiisotropic' in recipe
 assert 'define = -DPOSRES' not in recipe
 (tmp_path/'production.log').write_text('Step Time\n 250000 500.0\n')
 (tmp_path/'production.gro').touch();(tmp_path/'production.xtc').touch()
 with pytest.raises(RuntimeError,match='before the planned'):check_completion(tmp_path,'production',recipe)
 (tmp_path/'production.log').write_text('Step Time\n 500000 1000.0\n')
 check_completion(tmp_path,'production',recipe)

def test_lipid_rounding_is_bounded_and_recorded():
 mm=pytest.importorskip('openmm')
 from openmm import app,unit
 from scientific.membrane import normalize_lipid_rounding
 top=app.Topology();r=top.addResidue('POP',top.addChain());top.addAtom('P',app.element.phosphorus,r)
 system=mm.System();system.addParticle(31);force=mm.NonbondedForce();force.addParticle(-.000007,.3,.1);system.addForce(force)
 record=normalize_lipid_rounding(top,system,'POPE')
 assert record['total_delta_e']==pytest.approx(.000007)
 assert force.getParticleParameters(0)[0].value_in_unit(unit.elementary_charge)==pytest.approx(0)
 force.setParticleParameters(0,-.01,.3,.1)
 with pytest.raises(ValueError,match='materially'):normalize_lipid_rounding(top,system,'POPE')

def test_lipid_rounding_preserves_coulomb_exception_scaling():
 mm=pytest.importorskip('openmm')
 from openmm import app,unit
 from scientific.membrane import normalize_lipid_rounding
 top=app.Topology();chain=top.addChain();r=top.addResidue('POP',chain)
 top.addAtom('P',app.element.phosphorus,r);top.addAtom('O',app.element.oxygen,r)
 top.addAtom('H',app.element.hydrogen,top.addResidue('OTHER',chain))
 system=mm.System();force=mm.NonbondedForce()
 for q in [1.,-1.000007,.1]:system.addParticle(1);force.addParticle(q,.3,.1)
 force.addException(0,2,1*.1/1.2,.3,.05);system.addForce(force)
 normalize_lipid_rounding(top,system,'POPE')
 q1=force.getParticleParameters(0)[0].value_in_unit(unit.elementary_charge);q2=force.getParticleParameters(2)[0].value_in_unit(unit.elementary_charge)
 product=force.getExceptionParameters(0)[2].value_in_unit(unit.elementary_charge**2)
 assert product/(q1*q2)==pytest.approx(1/1.2)
