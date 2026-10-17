"""Find inputs only within a designated input directory, never a whole disk."""
from pathlib import Path

class SelectionRequired(ValueError):
    def __init__(self, message, candidates=()):
        super().__init__(message); self.candidates=list(map(str,candidates))

def candidates(directory):
    directory=Path(directory).expanduser().resolve()
    if not directory.is_dir(): return []
    return sorted(p for p in directory.rglob('*') if p.is_file() and p.suffix.lower()=='.pdb'
                  and not set(p.relative_to(directory).parts[:-1]) & {'results','runs','.backups','__pycache__'})

def resolve_input(value, directory, interactive=False, ask=input):
    if value:
        p=Path(value).expanduser().resolve()
        if not p.is_file() or p.suffix.lower()!='.pdb': raise ValueError('输入必须是存在的 PDB 文件：'+str(p))
        return p
    found=candidates(directory)
    if len(found)==1: return found[0]
    if not interactive:
        raise SelectionRequired('请选择一个输入；没有自动选取第一个文件。' if found else '输入目录没有 PDB，请使用 --input 指定文件。',found)
    for i,p in enumerate(found,1): print(f'{i}. {p}')
    answer=ask('输入文件序号或完整路径：').strip()
    if answer.isdigit() and 1<=int(answer)<=len(found): return found[int(answer)-1]
    return resolve_input(answer,directory,False)
