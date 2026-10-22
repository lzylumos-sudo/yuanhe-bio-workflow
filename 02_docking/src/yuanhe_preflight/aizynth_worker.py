"""Separate local AiZynthFinder process. No download or external model providers."""
import argparse, importlib.metadata, json
from pathlib import Path

def search(request):
    from aizynthfinder.aizynthfinder import AiZynthFinder
    cfg = request['retro']
    actual=importlib.metadata.version('aizynthfinder')
    if actual!='4.4.1':raise ValueError('Validated adapter requires AiZynthFinder 4.4.1; found '+actual)
    settings = {'search': {'max_transforms':cfg['max_depth'], 'time_limit':cfg['time_limit_seconds'],
                           'iteration_limit':cfg['iteration_limit'], 'exclude_target_from_stock':True},
                'expansion':{'local':[cfg['model_path'],cfg['templates_path']]},
                'stock':{'trusted':cfg['stock_path']}}
    finder = AiZynthFinder(configdict=settings)
    finder.expansion_policy.select('local'); finder.stock.select('trusted')
    finder.target_smiles = request['smiles']; finder.prepare_tree(); finder.tree_search(); finder.build_routes()
    return {'backend':'aizynthfinder', 'backend_version':importlib.metadata.version('aizynthfinder'),
            'target_smiles':request['smiles'], 'status':'searched', 'statistics':finder.extract_statistics(),
            'route_scores':finder.routes.scores,'all_route_scores':finder.routes.all_scores,
            'routes':[tree.to_dict(include_metadata=True) for tree in finder.routes.reaction_trees], 'stock_info':finder.stock_info(),
            'settings':settings, 'stereochemistry_capability':'not_certified; stereogenic routes require human review',
            'network':{'enabled':False,'sent_structures':[]}}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--request',required=True);ap.add_argument('--output',required=True);a=ap.parse_args()
    try: result=search(json.loads(Path(a.request).read_text()))
    except ModuleNotFoundError as e: result={'status':'dependency_missing','error':str(e)}
    except FileNotFoundError as e: result={'status':'model_or_stock_missing','error':str(e)}
    except Exception as e: result={'status':'backend_error','error':type(e).__name__+': '+str(e)}
    Path(a.output).write_text(json.dumps(result,ensure_ascii=False,default=str))
if __name__=='__main__':main()
